"""Strength check for the printed thruster mounts.

Each mount is a small 3D frame: straight members fixed at the hull plate and
joined rigidly at the thruster plate. Loads act at a point on that plate (with
any offset as a moment). Per member we report root axial force, bending and
torsion, the combined stress against an allowable for printed PLA, and Euler
buckling. The plate on the hull is checked as a 4-bolt group (M4 into the hull).

Sections are hollow rectangles: perimeter walls only (infill ignored, which is
conservative), `wall` thick.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .layout import Strut


@dataclass
class Material:
    name: str = "PLA (printed)"
    E: float = 3.0e9            # Pa
    G: float = 1.1e9
    sigma_allow: float = 25e6   # Pa, across-layer tensile (conservative for printed parts)


@dataclass
class Section:
    b: float      # thickness (m), across the member's thin direction
    h: float      # depth (m)
    t: float      # wall

    @property
    def area(self) -> float:
        return self.b * self.h - max(self.b - 2 * self.t, 0) * max(self.h - 2 * self.t, 0)

    @property
    def I_h(self) -> float:
        """Second moment for bending in the depth plane (strong axis)."""
        bi, hi = max(self.b - 2 * self.t, 0), max(self.h - 2 * self.t, 0)
        return (self.b * self.h**3 - bi * hi**3) / 12

    @property
    def I_b(self) -> float:
        """Second moment for bending in the thickness plane (weak axis)."""
        bi, hi = max(self.b - 2 * self.t, 0), max(self.h - 2 * self.t, 0)
        return (self.h * self.b**3 - hi * bi**3) / 12

    @property
    def J(self) -> float:
        """Thin-walled closed section (Bredt)."""
        am = (self.b - self.t) * (self.h - self.t)
        pm = 2 * ((self.b - self.t) + (self.h - self.t))
        return 4 * am**2 * self.t / pm

    @property
    def tau_factor(self) -> float:
        """Shear stress per unit torque: 1 / (2 A_m t)."""
        return 1.0 / (2 * (self.b - self.t) * (self.h - self.t) * self.t)


def member_frame(s: Strut) -> np.ndarray:
    """Rows: e1 along the member, e2 depth direction, e3 thickness direction."""
    u = (s.p1 - s.p0) / s.length
    if s.aligned:
        d = np.array([1.0, 0.0, 0.0]) - u * u[0]
    elif s.flat:
        d = np.array([0.0, 1.0, 0.0]) - u * u[1]
    else:
        d = np.cross([0.0, 1.0, 0.0], u)
    d /= np.linalg.norm(d)
    return np.vstack([u, d, np.cross(u, d)])


def tip_stiffness(L: float, sec: Section, m: Material) -> np.ndarray:
    """6x6 tip stiffness of a cantilever, local axes (1 axial, 2 depth, 3 thickness).

    Bending in the 1-2 plane (deflection along 2) uses I_h; in the 1-3 plane, I_b.
    """
    E, G = m.E, m.G
    Ih, Ib = sec.I_h, sec.I_b
    k = np.zeros((6, 6))
    k[0, 0] = E * sec.area / L
    k[1, 1] = 12 * E * Ih / L**3
    k[1, 5] = k[5, 1] = -6 * E * Ih / L**2
    k[5, 5] = 4 * E * Ih / L
    k[2, 2] = 12 * E * Ib / L**3
    k[2, 4] = k[4, 2] = 6 * E * Ib / L**2
    k[4, 4] = 4 * E * Ib / L
    k[3, 3] = G * sec.J / L
    return k


def _skew(r):
    return np.array([[0, -r[2], r[1]], [r[2], 0, -r[0]], [-r[1], r[0], 0]])


@dataclass
class MemberResult:
    name: str
    axial: float          # N, + tension
    m_depth: float        # N m, bending in the depth plane (strong)
    m_thick: float        # N m, bending in the thickness plane (weak)
    torque: float         # N m
    sigma: float          # Pa, |N|/A + bending (both axes, corner)
    tau: float            # Pa, torsion
    utilisation: float    # von Mises-ish / allowable
    buckling: float       # compression / Euler load (0 if in tension)


@dataclass
class CaseResult:
    case: str
    members: list = field(default_factory=list)
    bolt_tension_max: float = 0.0   # N per bolt
    bolt_shear_max: float = 0.0


def analyse(members: list[Strut], sections: list[Section], load_point: np.ndarray, force: np.ndarray,
            moment: np.ndarray, m: Material, bolt_plate: dict | None = None, case: str = "") -> CaseResult:
    """Solve the rigid-joint frame and evaluate member roots and the bolt group.

    bolt_plate: {"origin": point, "normal": unit vector out of the hull,
                 "bolts": list of in-plane offsets (u, v), "axes": (eu, ev)}.
    """
    K = np.zeros((6, 6))
    Ts, ks, Rs = [], [], []
    for s, sec in zip(members, sections):
        R = member_frame(s)
        r = s.p1 - load_point
        # rigid body (u_P, theta) -> tip (u, theta): u_tip = u_P + theta x r
        T = np.zeros((6, 6))
        T[:3, :3] = np.eye(3)
        T[:3, 3:] = -_skew(r)
        T[3:, 3:] = np.eye(3)
        Rb = np.zeros((6, 6))
        Rb[:3, :3] = R
        Rb[3:, 3:] = R
        k = tip_stiffness(s.length, sec, m)
        kg = Rb.T @ k @ Rb
        K += T.T @ kg @ T
        Ts.append(T), ks.append(kg), Rs.append(R)
    q = np.linalg.solve(K, np.concatenate([force, moment]))
    res = CaseResult(case)
    F_plate, M_plate = np.zeros(3), np.zeros(3)
    origin = bolt_plate["origin"] if bolt_plate else np.zeros(3)
    for s, sec, T, kg, R in zip(members, sections, Ts, ks, Rs):
        f_tip = kg @ (T @ q)                 # force/moment the member takes at its tip (global)
        F, Mt = f_tip[:3], f_tip[3:]
        # root reaction: equilibrium of the member
        L_vec = s.p1 - s.p0
        M_root = Mt + np.cross(L_vec, F)
        # internal forces at the root in local axes
        Fl, Ml = R @ F, R @ M_root
        N = Fl[0]
        m_depth, m_thick, tq = Ml[2], Ml[1], Ml[0]    # about e3 bends in 1-2 plane; about e2 in 1-3
        sig = abs(N) / sec.area + abs(m_depth) * (sec.h / 2) / sec.I_h + abs(m_thick) * (sec.b / 2) / sec.I_b
        tau = abs(tq) * sec.tau_factor
        util = np.sqrt(sig**2 + 3 * tau**2) / m.sigma_allow
        pcr = np.pi**2 * m.E * min(sec.I_h, sec.I_b) / (0.7 * s.length) ** 2
        buck = max(-N, 0.0) / pcr
        res.members.append(MemberResult(s.name, N, m_depth, m_thick, tq, sig, tau, util, buck))
        # load the member puts on the hull plate (equal and opposite to the root reaction on it)
        F_plate += F
        M_plate += M_root + np.cross(s.p0 - origin, F)
    if bolt_plate:
        n = bolt_plate["normal"]
        eu, ev = bolt_plate["axes"]
        bolts = np.asarray(bolt_plate["bolts"])
        fn = float(F_plate @ n)                  # + pulls the plate off the hull
        mu, mv = float(M_plate @ eu), float(M_plate @ ev)
        su, sv = (bolts[:, 0] ** 2).sum(), (bolts[:, 1] ** 2).sum()
        # moment about eu lifts the +v side when (eu x ev) == n
        tens = fn / len(bolts) + mu * bolts[:, 1] / sv - mv * bolts[:, 0] / su
        res.bolt_tension_max = float(max(tens.max(), 0.0))
        f_inplane = F_plate - fn * n
        mn = float(M_plate @ n)
        r2 = (bolts**2).sum()
        shear = np.linalg.norm(f_inplane) / len(bolts) + abs(mn) * np.linalg.norm(bolts, axis=1) / r2
        res.bolt_shear_max = float(shear.max())
    return res
