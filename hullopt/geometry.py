"""Parametric demihull and catamaran meshes.

A demihull is a set of stations along x. Each station has the same number of
section points, so neighbouring stations join with quads. Section shape: a
deadrise bottom from the keel, a bilge arc, a vertical side and a flat deck.
Section shape is convex, so the end caps are simple fans.
"""
from __future__ import annotations

import struct
from dataclasses import dataclass

import numpy as np

from .params import Design

N_BOT, N_ARC, N_SIDE, N_DECK = 6, 10, 5, 4
# Index (in a half section) of the deck-edge point; points 0..N_EDGE run keel -> deck edge.
N_EDGE = N_BOT + (N_ARC - 1) + (N_SIDE - 1) - 1
B_MIN = 1e-3  # half-beam at the bow tip, keeps the mesh non-degenerate


def entry_exponent(d: Design) -> float:
    """Exponent of the forebody half-beam curve that gives the requested entry angle."""
    bh = d.beam / 2
    n = np.tan(np.radians(d.entry_half_angle_deg)) * d.loa * (1 - d.xmax_frac) / bh
    return float(np.clip(n, 1.0, 8.0))


def actual_entry_angle_deg(d: Design) -> float:
    """Entry half-angle after the exponent is clipped."""
    bh = d.beam / 2
    n = entry_exponent(d)
    return float(np.degrees(np.arctan(bh * n / (d.loa * (1 - d.xmax_frac)))))


def halfbeam(xi, d: Design) -> np.ndarray:
    """Deck/side half-beam at xi = x / LOA."""
    xi = np.atleast_1d(np.asarray(xi, float))
    bh, xm = d.beam / 2, d.xmax_frac
    b = np.empty_like(xi)
    f = xi >= xm
    s = np.clip((xi[f] - xm) / (1 - xm), 0, 1)
    b[f] = bh * (1 - s ** entry_exponent(d))
    s = np.clip((xm - xi[~f]) / xm, 0, 1)
    t = d.transom_frac
    b[~f] = bh * (t + (1 - t) * (1 - s ** 2))
    return np.maximum(b, B_MIN)


def keel_z(xi, d: Design) -> np.ndarray:
    """Keel height (rocker). Lowest (z=0) at xmax_frac, parabolic rise to the ends."""
    xi = np.atleast_1d(np.asarray(xi, float))
    xk = d.xmax_frac
    z = np.where(
        xi >= xk,
        d.rocker_fwd * (np.clip(xi - xk, 0, None) / (1 - xk)) ** 2,
        d.rocker_aft * (np.clip(xk - xi, 0, None) / xk) ** 2,
    )
    return z


def half_section(b, zk, zdeck: float, deadrise_deg: float, bilge_frac: float):
    """Port half section for each station.

    Returns Y, Z of shape (n_stations, n_half), ordered keel (y=0) -> bottom ->
    bilge -> side -> deck edge -> deck centreline (y=0).
    """
    b = np.atleast_1d(np.asarray(b, float))[:, None]
    zk = np.atleast_1d(np.asarray(zk, float))[:, None]
    beta = np.radians(deadrise_deg)
    sb, cb = np.sin(beta), np.cos(beta)
    r = bilge_frac * b
    zc = zk + ((b - r) * sb + r) / cb          # bilge circle centre height
    p1y = b - r + r * sb                       # bottom/bilge tangent point

    t = np.linspace(0, 1, N_BOT)[None, :]
    yb = t * p1y
    zb = zk + yb * np.tan(beta)

    th = np.linspace(beta - np.pi / 2, 0, N_ARC)[None, 1:]
    ya = b - r + r * np.cos(th)
    za = zc + r * np.sin(th)

    ts = np.linspace(0, 1, N_SIDE)[None, 1:]
    ys = np.broadcast_to(b, (b.shape[0], ts.shape[1]))
    zs = zc + ts * (zdeck - zc)

    td = np.linspace(0, 1, N_DECK)[None, 1:]
    yd = b * (1 - td)
    zd = np.full_like(yd, zdeck)

    Y = np.concatenate([yb, ya, ys, yd], axis=1)
    Z = np.concatenate([zb, za, zs, zd], axis=1)
    return Y, np.minimum(Z, zdeck)


def station_xi(n: int) -> np.ndarray:
    """Cosine spacing: more stations near the bow and transom."""
    t = np.linspace(0, 1, n)
    return 0.5 * (1 - np.cos(np.pi * t))


@dataclass
class Mesh:
    vertices: np.ndarray  # (nv, 3)
    faces: np.ndarray     # (nf, 3) int, outward normals

    @property
    def triangles(self) -> np.ndarray:
        return self.vertices[self.faces]

    def translated(self, dx=0.0, dy=0.0, dz=0.0) -> "Mesh":
        return Mesh(self.vertices + np.array([dx, dy, dz]), self.faces.copy())

    @staticmethod
    def merge(*meshes: "Mesh") -> "Mesh":
        verts, faces, off = [], [], 0
        for m in meshes:
            verts.append(m.vertices)
            faces.append(m.faces + off)
            off += len(m.vertices)
        return Mesh(np.vstack(verts), np.vstack(faces))


def loft(X, Y, Z) -> Mesh:
    """Mesh from closed section loops. X: (ns,), Y/Z: (ns, N), loops CCW about +x."""
    ns, N = Y.shape
    V = np.empty((ns * N + 2, 3))
    V[: ns * N, 0] = np.repeat(X, N)
    V[: ns * N, 1] = Y.ravel()
    V[: ns * N, 2] = Z.ravel()
    i = np.arange(ns - 1)[:, None]
    j = np.arange(N)[None, :]
    a = i * N + j
    b = i * N + (j + 1) % N
    c = (i + 1) * N + (j + 1) % N
    dd = (i + 1) * N + j
    side = np.concatenate(
        [np.stack([a, b, c], -1).reshape(-1, 3), np.stack([a, c, dd], -1).reshape(-1, 3)]
    )
    c0, c1 = ns * N, ns * N + 1
    V[c0] = [X[0], Y[0].mean(), Z[0].mean()]
    V[c1] = [X[-1], Y[-1].mean(), Z[-1].mean()]
    jj = np.arange(N)
    cap0 = np.stack([np.full(N, c0), (jj + 1) % N, jj], -1)
    last = (ns - 1) * N
    cap1 = np.stack([np.full(N, c1), last + jj, last + (jj + 1) % N], -1)
    return Mesh(V, np.concatenate([side, cap0, cap1]).astype(np.int64))


def demihull_mesh(d: Design, n_stations: int = 61) -> Mesh:
    """Demihull centred on y=0."""
    xi = station_xi(n_stations)
    Y, Z = half_section(halfbeam(xi, d), keel_z(xi, d), d.depth, d.deadrise_deg, d.bilge_frac)
    # Close the loop with the mirrored half (skip the shared keel/deck-centre points).
    Yf = np.concatenate([Y, -Y[:, -2:0:-1]], axis=1)
    Zf = np.concatenate([Z, Z[:, -2:0:-1]], axis=1)
    return loft(xi * d.loa, Yf, Zf)


def hull_centre_y(d: Design) -> float:
    """Port demihull centreline y; starboard is the mirror."""
    return d.spacing_inner / 2 + d.beam / 2


def catamaran_mesh(d: Design, n_stations: int = 61) -> Mesh:
    m = demihull_mesh(d, n_stations)
    yc = hull_centre_y(d)
    return Mesh.merge(m.translated(dy=yc), m.translated(dy=-yc))


def halfbreadth_at(x, z, d: Design) -> np.ndarray:
    """Demihull half-breadth at body (x, z). NaN where there is no hull."""
    x = np.atleast_1d(np.asarray(x, float))
    z = np.atleast_1d(np.asarray(z, float))
    out = np.full(x.shape, np.nan)
    xi = x / d.loa
    ok = (xi >= 0) & (xi <= 1)
    if not ok.any():
        return out
    idx = np.nonzero(ok)[0]
    Y, Z = half_section(halfbeam(xi[ok], d), keel_z(xi[ok], d), d.depth, d.deadrise_deg, d.bilge_frac)
    Y, Z = Y[:, : N_EDGE + 1], Z[:, : N_EDGE + 1]
    for k, i in enumerate(idx):
        zz = z[i]
        if zz < Z[k, 0] or zz > d.depth:
            continue
        out[i] = np.interp(zz, Z[k], Y[k])
    return out


def write_stl(path: str, tri: np.ndarray, name: str = "hull") -> None:
    """Binary STL."""
    tri = np.asarray(tri, dtype=np.float32)
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    ln = np.linalg.norm(n, axis=1, keepdims=True)
    n = np.where(ln > 0, n / np.where(ln > 0, ln, 1), 0).astype(np.float32)
    rec = np.zeros(len(tri), dtype=[("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")])
    rec["n"], rec["v"] = n, tri
    with open(path, "wb") as f:
        f.write(name.encode()[:80].ljust(80, b"\0"))
        f.write(struct.pack("<I", len(tri)))
        f.write(rec.tobytes())
