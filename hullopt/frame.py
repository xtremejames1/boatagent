"""Small 3D frame solver (Euler-Bernoulli beams, direct stiffness).

Nodes carry 6 DOF. Each element is a straight beam between two nodes with a
hollow-rectangle Section and a "depth direction" that fixes its orientation.
Supports clamp all 6 DOF. Results: nodal displacements and, per element, the
axial force, bending moments and torque at both ends, the combined stress and
an Euler buckling ratio.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .strength import Material, Section


@dataclass
class Element:
    name: str
    n0: int
    n1: int
    sec: Section
    depth_dir: np.ndarray      # roughly the section's depth direction (projected perpendicular)
    buckling_k: float = 1.0    # effective length factor


@dataclass
class ElementResult:
    name: str
    axial: float               # + tension (N)
    m_depth: float             # max |moment| bending in the depth plane (N m)
    m_thick: float             # max |moment| bending in the thickness plane (N m)
    torque: float
    sigma: float               # Pa
    tau: float
    utilisation: float
    buckling: float


@dataclass
class Frame:
    nodes: list = field(default_factory=list)       # positions (m)
    elements: list = field(default_factory=list)
    fixed: set = field(default_factory=set)

    def node(self, p) -> int:
        self.nodes.append(np.asarray(p, float))
        return len(self.nodes) - 1

    def add(self, name, n0, n1, sec, depth_dir, k=1.0) -> None:
        self.elements.append(Element(name, n0, n1, sec, np.asarray(depth_dir, float), k))

    # ------------------------------------------------------------------
    def _axes(self, e: Element):
        p0, p1 = self.nodes[e.n0], self.nodes[e.n1]
        L = float(np.linalg.norm(p1 - p0))
        e1 = (p1 - p0) / L
        d = e.depth_dir - e1 * (e.depth_dir @ e1)
        if np.linalg.norm(d) < 1e-9:     # depth along the axis: pick any perpendicular
            d = np.cross(e1, [0, 0, 1.0]) if abs(e1[2]) < 0.9 else np.cross(e1, [1.0, 0, 0])
        e2 = d / np.linalg.norm(d)
        return L, np.vstack([e1, e2, np.cross(e1, e2)])

    @staticmethod
    def _k_local(L, sec: Section, m: Material) -> np.ndarray:
        E, G = m.E, m.G
        A, Iz, Iy, J = sec.area, sec.I_h, sec.I_b, sec.J   # z: bending in the 1-2 (depth) plane
        k = np.zeros((12, 12))
        k[0, 0] = k[6, 6] = E * A / L
        k[0, 6] = k[6, 0] = -E * A / L
        k[3, 3] = k[9, 9] = G * J / L
        k[3, 9] = k[9, 3] = -G * J / L
        a = E * Iz / L**3
        for (i, j, v) in [(1, 1, 12), (1, 5, 6 * L), (1, 7, -12), (1, 11, 6 * L), (5, 5, 4 * L * L),
                          (5, 7, -6 * L), (5, 11, 2 * L * L), (7, 7, 12), (7, 11, -6 * L), (11, 11, 4 * L * L)]:
            k[i, j] = k[j, i] = a * v
        b = E * Iy / L**3
        for (i, j, v) in [(2, 2, 12), (2, 4, -6 * L), (2, 8, -12), (2, 10, -6 * L), (4, 4, 4 * L * L),
                          (4, 8, 6 * L), (4, 10, 2 * L * L), (8, 8, 12), (8, 10, 6 * L), (10, 10, 4 * L * L)]:
            k[i, j] = k[j, i] = b * v
        return k

    def solve(self, loads: dict, m: Material) -> tuple[np.ndarray, list[ElementResult]]:
        """loads: {node: (force(3), moment(3))}. Returns displacements and element results."""
        n = len(self.nodes)
        K = np.zeros((6 * n, 6 * n))
        cache = []
        for e in self.elements:
            L, R = self._axes(e)
            T = np.kron(np.eye(4), R)
            kl = self._k_local(L, e.sec, m)
            kg = T.T @ kl @ T
            idx = np.r_[6 * e.n0:6 * e.n0 + 6, 6 * e.n1:6 * e.n1 + 6]
            K[np.ix_(idx, idx)] += kg
            cache.append((L, T, kl, idx))
        F = np.zeros(6 * n)
        for i, (f, mo) in loads.items():
            F[6 * i:6 * i + 3] += f
            F[6 * i + 3:6 * i + 6] += mo
        free = np.array([d for i in range(n) if i not in self.fixed for d in range(6 * i, 6 * i + 6)])
        u = np.zeros(6 * n)
        u[free] = np.linalg.solve(K[np.ix_(free, free)], F[free])
        out = []
        for e, (L, T, kl, idx) in zip(self.elements, cache):
            f = kl @ (T @ u[idx])                    # local end forces on the element
            N = -f[0]                                # tension positive
            tq = max(abs(f[3]), abs(f[9]))
            m_thick = max(abs(f[4]), abs(f[10]))     # about local 2: bends in the 1-3 plane
            m_depth = max(abs(f[5]), abs(f[11]))     # about local 3: bends in the 1-2 plane
            s = e.sec
            sig = abs(N) / s.area + m_depth * (s.h / 2) / s.I_h + m_thick * (s.b / 2) / s.I_b
            tau = tq * s.tau_factor
            util = np.sqrt(sig**2 + 3 * tau**2) / m.sigma_allow
            pcr = np.pi**2 * m.E * min(s.I_h, s.I_b) / (e.buckling_k * L) ** 2
            out.append(ElementResult(e.name, N, m_depth, m_thick, tq, sig, tau, util, max(-N, 0.0) / pcr))
        return u, out

    def reactions(self, u: np.ndarray, m: Material) -> dict:
        """Support reaction (force, moment) acting on the frame at each fixed node."""
        n = len(self.nodes)
        K = np.zeros((6 * n, 6 * n))
        for e in self.elements:
            L, R = self._axes(e)
            T = np.kron(np.eye(4), R)
            kg = T.T @ self._k_local(L, e.sec, m) @ T
            idx = np.r_[6 * e.n0:6 * e.n0 + 6, 6 * e.n1:6 * e.n1 + 6]
            K[np.ix_(idx, idx)] += kg
        r = K @ u
        return {i: (r[6 * i:6 * i + 3], r[6 * i + 3:6 * i + 6]) for i in self.fixed}
