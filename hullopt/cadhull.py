"""Smooth CAD solid of a demihull (B-spline surfaces), for STEP export and editing.

Each station is a wire of three edges: a spline from the keel up the port side to
the deck edge, a straight deck line, and the mirrored starboard spline. A smooth
(non-ruled) loft through the stations gives continuous B-spline side and bottom
surfaces, a flat deck, a flat transom, and sharp creases only where the hull has
them (keel line when there is deadrise, deck edges). Units in CAD: mm.
"""
from __future__ import annotations

import numpy as np

from .geometry import N_EDGE, half_section, halfbeam, keel_z
from .params import Design


def _half_points(d: Design, x: float, n: int = 21) -> np.ndarray:
    """Keel -> deck edge (y >= 0) at station x, resampled evenly by arc length."""
    xi = np.clip(x / d.loa, 0, 1)
    b = halfbeam(np.array([xi]), d)
    k = keel_z(np.array([xi]), d)
    Y, Z = half_section(b, k, d.depth, d.deadrise_deg, d.bilge_frac)
    P = np.c_[Y[0, : N_EDGE + 1], Z[0, : N_EDGE + 1]]
    keep = np.r_[True, np.linalg.norm(np.diff(P, axis=0), axis=1) > 1e-7]
    P = P[keep]
    s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]
    t = np.linspace(0, s[-1], n)
    return np.c_[np.interp(t, s, P[:, 0]), np.interp(t, s, P[:, 1])]


def section_wire(d: Design, x: float, y0: float = 0.0, n: int = 21):
    import cadquery as cq
    P = _half_points(d, x, n) * 1000.0
    X = x * 1000.0
    port = [cq.Vector(X, y0 * 1000 + y, z) for y, z in P]
    stbd = [cq.Vector(X, y0 * 1000 - y, z) for y, z in P[::-1]]
    e1 = cq.Edge.makeSpline(port)
    e2 = cq.Edge.makeLine(port[-1], stbd[0])           # deck
    e3 = cq.Edge.makeSpline(stbd)
    return cq.Wire.assembleEdges([e1, e2, e3])


def smooth_demihull(d: Design, n_sections: int = 31, y0: float = 0.0, bow_stop: float = 0.996):
    """Smooth solid demihull from the transom (x=0) to just short of the bow tip
    (a knife edge cannot close a solid; the last section is ~1-2 mm wide)."""
    import cadquery as cq
    t = np.linspace(0, 1, n_sections)
    xi = 0.5 * (1 - np.cos(np.pi * t)) * bow_stop      # more sections near both ends
    wires = [section_wire(d, float(x) * d.loa, y0) for x in xi]
    return cq.Workplane("XY").add(cq.Solid.makeLoft(wires, False))
