"""Hydrostatics of a closed triangle mesh, by surface integrals.

The mesh is placed in the earth frame (water surface at z=0) and clipped to
z <= 0. The divergence theorem with fields that vanish at z=0 gives the
submerged volume and its centroid without building the waterplane cap. The
waterplane itself comes from the clipped edges (Green's theorem).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


# ---------------------------------------------------------------- mesh integrals
def area_vectors(tri: np.ndarray) -> np.ndarray:
    return 0.5 * np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])


def volume_props(tri: np.ndarray) -> tuple[float, np.ndarray]:
    """Volume and centroid enclosed by outward triangles.

    Valid for a closed mesh, or for an open mesh whose missing part lies on
    z=0 (every field used here contains a factor z).
    """
    if len(tri) == 0:
        return 0.0, np.zeros(3)
    Az = area_vectors(tri)[:, 2]
    m = np.stack(
        [0.5 * (tri[:, 0] + tri[:, 1]), 0.5 * (tri[:, 1] + tri[:, 2]), 0.5 * (tri[:, 2] + tri[:, 0])],
        axis=1,
    )  # edge midpoints: exact for quadratic integrands
    x, y, z = m[..., 0], m[..., 1], m[..., 2]
    V = np.sum(Az * tri[:, :, 2].mean(1))
    Mx = np.sum(Az * (x * z).mean(1))
    My = np.sum(Az * (y * z).mean(1))
    Mz = np.sum(Az * (0.5 * z * z).mean(1))
    if abs(V) < 1e-15:
        return 0.0, np.zeros(3)
    return float(V), np.array([Mx, My, Mz]) / V


def area_props(tri: np.ndarray) -> tuple[float, np.ndarray]:
    """Surface area and area centroid."""
    a = np.linalg.norm(area_vectors(tri), axis=1)
    A = a.sum()
    return float(A), (tri.mean(1) * a[:, None]).sum(0) / max(A, 1e-15)


def _cut(p, q):
    t = p[:, 2] / (p[:, 2] - q[:, 2])
    r = p + t[:, None] * (q - p)
    r[:, 2] = 0.0
    return r


def clip_below(tri: np.ndarray):
    """Part of each triangle with z <= 0, plus waterline edges.

    Returns (triangles, seg_start, seg_end). Edges run along the boundary of
    the submerged surface in the same sense as its triangles.
    """
    z = tri[:, :, 2]
    below = z <= 0
    nb = below.sum(1)
    out = [tri[nb == 3]]
    s0, s1 = [], []

    sel = nb == 1
    if sel.any():
        t, bl = tri[sel], below[sel]
        k = np.argmax(bl, axis=1)
        idx = (k[:, None] + np.arange(3)[None]) % 3
        r = np.take_along_axis(t, idx[:, :, None], axis=1)
        p0, p1, p2 = r[:, 0], r[:, 1], r[:, 2]
        a, b = _cut(p0, p1), _cut(p0, p2)
        out.append(np.stack([p0, a, b], 1))
        s0.append(a)
        s1.append(b)

    sel = nb == 2
    if sel.any():
        t, ab = tri[sel], ~below[sel]
        k = np.argmax(ab, axis=1)
        idx = (k[:, None] + np.arange(3)[None]) % 3
        r = np.take_along_axis(t, idx[:, :, None], axis=1)
        p0, p1, p2 = r[:, 0], r[:, 1], r[:, 2]   # p0 is above water
        a, b = _cut(p0, p1), _cut(p0, p2)
        out.append(np.stack([a, p1, p2], 1))
        out.append(np.stack([a, p2, b], 1))
        s0.append(b)
        s1.append(a)

    tri_b = np.concatenate(out) if out else np.zeros((0, 3, 3))
    if s0:
        return tri_b, np.concatenate(s0), np.concatenate(s1)
    return tri_b, np.zeros((0, 3)), np.zeros((0, 3))


@dataclass
class Waterplane:
    area: float = 0.0
    centroid: np.ndarray = field(default_factory=lambda: np.zeros(2))
    i_t: float = 0.0   # second moment about the longitudinal axis through the centroid
    i_l: float = 0.0   # second moment about the transverse axis through the centroid
    x_min: float = 0.0
    x_max: float = 0.0
    y_min: float = 0.0
    y_max: float = 0.0


def waterplane(s0: np.ndarray, s1: np.ndarray) -> Waterplane:
    if len(s0) == 0:
        return Waterplane()
    # The cap's boundary runs opposite to the submerged surface's: use s1 -> s0.
    ax, ay, bx, by = s1[:, 0], s1[:, 1], s0[:, 0], s0[:, 1]
    cr = ax * by - bx * ay
    A = 0.5 * cr.sum()
    Sx = (cr * (ax + bx)).sum() / 6
    Sy = (cr * (ay + by)).sum() / 6
    Iyy = (cr * (ax * ax + ax * bx + bx * bx)).sum() / 12   # int x^2 dA
    Ixx = (cr * (ay * ay + ay * by + by * by)).sum() / 12   # int y^2 dA
    xc, yc = Sx / A, Sy / A
    pts = np.vstack([s0, s1])
    return Waterplane(
        area=float(A),
        centroid=np.array([xc, yc]),
        i_t=float(Ixx - A * yc * yc),
        i_l=float(Iyy - A * xc * xc),
        x_min=float(pts[:, 0].min()),
        x_max=float(pts[:, 0].max()),
        y_min=float(pts[:, 1].min()),
        y_max=float(pts[:, 1].max()),
    )


# ---------------------------------------------------------------- floating pose
@dataclass
class Pose:
    """Body -> earth: p_e = R (p_b - pivot) - [0, 0, heave], R = Rx(heel) Ry(trim).

    trim > 0 is bow down; heel > 0 lifts the port side. heave is the depth of
    the pivot below the water surface.
    """

    heave: float
    trim: float = 0.0
    heel: float = 0.0
    pivot: np.ndarray = field(default_factory=lambda: np.zeros(3))

    def R(self) -> np.ndarray:
        ct, st = np.cos(self.trim), np.sin(self.trim)
        ch, sh = np.cos(self.heel), np.sin(self.heel)
        Ry = np.array([[ct, 0, st], [0, 1, 0], [-st, 0, ct]])
        Rx = np.array([[1, 0, 0], [0, ch, -sh], [0, sh, ch]])
        return Rx @ Ry

    def to_earth(self, p: np.ndarray) -> np.ndarray:
        return (np.asarray(p) - self.pivot) @ self.R().T - np.array([0, 0, self.heave])

    def to_body(self, p: np.ndarray) -> np.ndarray:
        return (np.asarray(p) + np.array([0, 0, self.heave])) @ self.R() + self.pivot


@dataclass
class Appendages:
    """Small fully-enclosed parts treated as point volumes (skids, struts, thrusters)."""

    volumes: np.ndarray = field(default_factory=lambda: np.zeros(0))
    positions: np.ndarray = field(default_factory=lambda: np.zeros((0, 3)))

    def add(self, vol: float, pos) -> None:
        self.volumes = np.append(self.volumes, vol)
        self.positions = np.vstack([self.positions, np.asarray(pos, float)[None]])


@dataclass
class Hydro:
    volume: float            # mesh + submerged appendages (m^3)
    mesh_volume: float
    cob: np.ndarray          # centre of buoyancy, earth frame
    wetted_area: float       # mesh only
    wp: Waterplane
    lateral_area: float      # projected submerged area on the xz plane (sum over hulls)
    lateral_centroid: np.ndarray  # earth (x, z) of that area


def hydro(tri_body: np.ndarray, pose: Pose, app: Appendages | None = None) -> Hydro:
    tri_e = pose.to_earth(tri_body.reshape(-1, 3)).reshape(-1, 3, 3)
    # cheap reject of triangles fully above water
    tri_e = tri_e[tri_e[:, :, 2].min(1) <= 0]
    tb, s0, s1 = clip_below(tri_e)
    V, c = volume_props(tb)
    av = area_vectors(tb)
    wet = float(np.linalg.norm(av, axis=1).sum())
    ay = np.abs(av[:, 1])
    lat = 0.5 * float(ay.sum())
    if lat > 0:
        cen = tb.mean(1)
        latc = np.array([(ay * cen[:, 0]).sum(), (ay * cen[:, 2]).sum()]) / ay.sum()
    else:
        latc = np.zeros(2)
    Vt, M = V, V * c
    if app is not None and len(app.volumes):
        pe = pose.to_earth(app.positions)
        sub = pe[:, 2] < 0
        Vt += app.volumes[sub].sum()
        M = M + (app.volumes[sub, None] * pe[sub]).sum(0)
    cob = M / Vt if Vt > 0 else np.zeros(3)
    return Hydro(Vt, V, cob, wet, waterplane(s0, s1), lat, latc)


def solve_heave(tri, mass, rho, pose: Pose, app=None, h_lo=-1.0, h_hi=1.0, tol=1e-9):
    """Heave that floats `mass` with trim/heel held fixed. None if it sinks."""
    Vt = mass / rho

    def f(h):
        return hydro(tri, Pose(h, pose.trim, pose.heel, pose.pivot), app).volume - Vt

    if f(h_hi) < 0:
        return None
    lo, hi = h_lo, h_hi
    h = 0.5 * (lo + hi)
    # Newton on V(h) with dV/dh = waterplane area, bracketed by bisection.
    for _ in range(60):
        hy = hydro(tri, Pose(h, pose.trim, pose.heel, pose.pivot), app)
        r = hy.volume - Vt
        if abs(r) < tol * max(Vt, 1e-9):
            break
        if r > 0:
            hi = h
        else:
            lo = h
        awp = hy.wp.area * np.cos(pose.trim) * np.cos(pose.heel)
        hn = h - r / awp if awp > 1e-9 else 0.5 * (lo + hi)
        h = hn if lo < hn < hi else 0.5 * (lo + hi)
    return h


def equilibrium(tri, mass, cg_body, rho, pivot, app=None, h_hi=1.0, iters=30):
    """Free heave and trim, zero heel. Returns (pose, hydro) or (None, None)."""
    L = max(np.ptp(tri[:, :, 0]), 1e-6)
    Vt = mass / rho
    pose = Pose(0.0, 0.0, 0.0, np.asarray(pivot, float))
    h = solve_heave(tri, mass, rho, pose, app, h_lo=-h_hi, h_hi=h_hi)
    if h is None:
        return None, None
    pose.heave = h

    def resid(p):
        hy = hydro(tri, p, app)
        g = p.to_earth(np.asarray(cg_body)[None])[0]
        return np.array([(hy.volume - Vt) / Vt, (hy.cob[0] - g[0]) / L]), hy

    r, hy = resid(pose)
    for _ in range(iters):
        if np.max(np.abs(r)) < 1e-8:
            break
        J = np.empty((2, 2))
        for k, (dh, dt) in enumerate([(1e-5, 0.0), (0.0, 1e-5)]):
            rp, _ = resid(Pose(pose.heave + dh, pose.trim + dt, 0.0, pose.pivot))
            J[:, k] = (rp - r) / (dh + dt)
        try:
            step = np.linalg.solve(J, -r)
        except np.linalg.LinAlgError:
            return None, None
        step[1] = np.clip(step[1], -0.05, 0.05)
        pose = Pose(pose.heave + step[0], pose.trim + step[1], 0.0, pose.pivot)
        r, hy = resid(pose)
    if np.max(np.abs(r)) > 1e-5:
        return None, None
    return pose, hy


def gz_curve(tri, mass, cg_body, rho, pose0: Pose, heel_deg, app=None, h_hi=1.0):
    """Righting arm GZ (m) at fixed trim for each heel angle. Positive = righting."""
    out = []
    for a in np.atleast_1d(heel_deg):
        p = Pose(0.0, pose0.trim, np.radians(a), pose0.pivot)
        h = solve_heave(tri, mass, rho, p, app, h_lo=-h_hi, h_hi=h_hi)
        if h is None:
            out.append(np.nan)
            continue
        p.heave = h
        hy = hydro(tri, p, app)
        g = p.to_earth(np.asarray(cg_body)[None])[0]
        out.append(float(g[1] - hy.cob[1]))
    return np.array(out)
