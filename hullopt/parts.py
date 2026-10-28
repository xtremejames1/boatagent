"""Printable parametric parts (CadQuery), placed in the boat's body frame (mm in CAD).

Front, port side (starboard is the mirror image):
  front_frame    one part: a bow cap that slides over a spigot on the blunt-cut hull
                 (4x M4 sideways into heat-set inserts in the hull) and narrows to a
                 fine entry; a forefoot fin below it; an upper strut in two segments
                 (cap -> T200 plate -> top fairing lug) that veers outboard so the
                 4 M3 thruster screws are reachable from above; a level rail under the
                 thruster to the bottom fairing lug. Prints nose-up with no supports.
  front_fairing  separate nose fairing / bumper: starts at the thruster inlet, opening
                 where the inlet keep-out cone passes, closes at <= 35 deg to a small
                 flat tip. Bolts to the two lugs with nylon M3 screws (the break-off
                 point). Prints nose-up on its open aft rim.
  hull_bow_spigot  reference only: the last part of the hull with the spigot and
                 the insert holes the cap screws into.

Head-on hits go fairing -> lugs -> struts and rail -> cap -> hull; the thruster
hangs from the plate and is not in that path. Undo the 4 M3 screws from above,
lower the thruster a little and slide it out sideways under the upper strut.

Every build is checked: overhangs (print direction +x), keep-out cones at the
inlet and outlet, interference, the thruster's removal path, screw access, bed
fit. Parts are exported as STEP and STL, sliced with PrusaSlicer, and the frame
is checked for strength (hullopt.frame).

  python -m hullopt.parts --design results/stage2/designs/finalist_1.json --out results/parts
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import os
import re
import shutil
import subprocess
from dataclasses import dataclass, replace

import numpy as np

from .frame import Frame
from .geometry import Mesh, half_section, halfbeam, keel_z, write_stl
from .layout import cap_cut_x, cap_top, cap_z_range, mount_plate_centre
from .params import Config, Design
from .stage1 import evaluate
from .strength import Material, Section

MM = 1000.0


@dataclass
class PrintSettings:
    """Slicer settings (PrusaSlicer CLI)."""

    material: str = "PLA"
    density: float = 1.24          # g/cm^3
    perimeters: int = 4            # 4 x 0.4 mm = 1.6 mm walls
    infill: int = 10               # %
    pattern: str = "gyroid"
    layer: float = 0.2
    nozzle: float = 0.4
    top_layers: int = 5
    bottom_layers: int = 4
    bed: tuple = (250, 250, 250)   # mm

    @property
    def wall(self) -> float:
        return self.perimeters * self.nozzle / MM   # m


@dataclass
class FrontSpec:
    """Hardware and checks for the front frame and fairing (metres, newtons)."""

    hull_m4: float = 0.0045        # M4 clearance in the cap
    insert_m4: float = 0.0056      # M4 heat-set insert hole in the hull spigot
    insert_m4_depth: float = 0.008
    t200_m3: float = 0.0034        # M3 clearance in the T200 plate
    cable_hole: float = 0.012      # cable pass-through in the T200 plate
    insert_m3: float = 0.0040      # M3 heat-set insert in the lugs (fairing screws)
    insert_m3_depth: float = 0.006
    lug: float = 0.016             # lug block size (y and z)
    fairing_screws: int = 2        # nylon M3, the break-off point
    fuse_shear_each: float = 210.0  # N, nylon (PA66) M3 in single shear, approx.
    hold_factor: float = 3.0       # frame must hold this x full thrust
    impact_factor: float = 1.25    # and this x the fairing breakaway load head-on
    side_factor: float = 0.5       # side hits are secondary: report capacity at this x breakaway
    insert_shear: float = 1000.0   # N per M4 heat-set insert in shear, placeholder (test it)
    removal_clear: float = 0.002   # thruster drops this much before sliding out
    access_d: float = 0.008        # screwdriver/key clearance above each M3


# ------------------------------------------------------------------ CAD helpers
def _cq():
    import cadquery as cq
    return cq


def _vec(p):
    cq = _cq()
    return cq.Vector(*(float(c) * MM for c in p))


def _poly(points3d):
    cq = _cq()
    return cq.Wire.makePolygon([_vec(p) for p in points3d], close=True)


def _loft(wires, ruled=True):
    cq = _cq()
    return cq.Workplane("XY").add(cq.Solid.makeLoft(wires, ruled))


def bar3d(p0, p1, depth_dir, depth: float, thick: float):
    """Rectangular bar from p0 to p1 (m); `depth` along depth_dir (projected
    perpendicular to the bar), `thick` across."""
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    e1 = (p1 - p0) / np.linalg.norm(p1 - p0)
    d = np.asarray(depth_dir, float) - e1 * (np.asarray(depth_dir, float) @ e1)
    d /= np.linalg.norm(d)
    t = np.cross(e1, d)
    def rect(c):
        return _poly([c + sd * d * depth / 2 + st * t * thick / 2
                      for sd, st in ((1, 1), (-1, 1), (-1, -1), (1, -1))])
    return _loft([rect(p0), rect(p1)])


def hull_loop(d: Design, x: float, y0: float, inset: float = 0.0, scale_y: float = 1.0,
              zk: float | None = None, z_top: float | None = None, n_pts: int = 64):
    """Closed hull section at station x (m), centred on y0, optionally inset and
    cut off at z_top."""
    xi = np.clip(x / d.loa, 0, 1)
    b = float(halfbeam(np.array([xi]), d)[0])
    k = float(keel_z(np.array([xi]), d)[0]) if zk is None else zk
    Y, Z = half_section(np.array([b]), np.array([k]), d.depth, d.deadrise_deg, d.bilge_frac)
    Y, Z = Y[0] * scale_y, Z[0]
    top = d.depth if z_top is None else min(z_top, d.depth)
    Y = np.maximum(Y - inset, 0.0)
    Z = np.clip(Z, k + inset, top - inset)
    # mirror, then resample to a fixed number of evenly spaced points so that every
    # section has the same vertex count (robust ruled lofts, no repeated points)
    pts = list(zip(Y, Z))
    loop = np.array(pts + [(-y, z) for y, z in reversed(pts[1:-1])])
    closed = np.vstack([loop, loop[:1]])
    seg = np.linalg.norm(np.diff(closed, axis=0), axis=1)
    s_ = np.concatenate([[0], np.cumsum(seg)])
    t = np.linspace(0, s_[-1], n_pts, endpoint=False)
    yy, zz = np.interp(t, s_, closed[:, 0]), np.interp(t, s_, closed[:, 1])
    return _poly([(x, y0 + a, b) for a, b in zip(yy, zz)])


def ellipse_wire(x, y, z, a, b):
    cq = _cq()
    return cq.Wire.makeEllipse(a * MM, b * MM, _vec((x, y, z)), cq.Vector(1, 0, 0), cq.Vector(0, 1, 0))


def teardrop(centre, axis, d: float, length: float, up=(1.0, 0, 0), max_deg: float = 35.0):
    """Hole of diameter d along `axis` (unit, perpendicular to `up`), with a point
    toward `up` so the top of the hole prints without overhang."""
    cq = _cq()
    axis = np.asarray(axis, float) / np.linalg.norm(axis)
    up = np.asarray(up, float)
    r = d / 2
    a = math.radians(max_deg)
    tip = r / math.sin(a)                     # apex distance from the centre, along `up`
    # arc from one tangent point (angle 90-a from `up`) round the bottom to the other
    ang = np.linspace(math.pi / 2 - a, 2 * math.pi - (math.pi / 2 - a), 25)
    poly = [(r * math.cos(t), r * math.sin(t)) for t in ang] + [(tip, 0.0)]
    pl = cq.Plane(origin=tuple(np.asarray(centre) * MM), xDir=tuple(up), normal=tuple(axis))
    return cq.Workplane(pl).polyline([(u * MM, v * MM) for u, v in poly]).close().extrude(length * MM, both=True)


def frustum(base, axis, r0, r1, length):
    cq = _cq()
    return cq.Workplane("XY").add(cq.Solid.makeCone(r0 * MM, r1 * MM, length * MM, _vec(base),
                                                   cq.Vector(*axis)))


def thruster_solid(t, drop: float = 0.0, dy: float = 0.0):
    cq = _cq()
    c = t.centre + np.array([0.0, dy, -drop])
    return cq.Workplane("XY").add(cq.Solid.makeCylinder(t.radius * MM, t.length * MM,
                                                       _vec(c - t.axis * t.length / 2),
                                                       cq.Vector(*t.axis)))


def _vol(a, b) -> float:
    """Intersection volume (mm^3); 0 if none."""
    try:
        r = a.intersect(b)
        return float(sum(s.Volume() for s in r.solids().vals())) if r.solids().size() else 0.0
    except Exception:
        return 0.0


# ------------------------------------------------------------------ geometry
def front_geometry(d: Design, cfg: Config):
    row, det = evaluate(d, cfg, detail=True)
    if det is None:
        raise ValueError(f"design does not float: {row.get('fail_reason')}")
    fp = next(t for t in det.thrusters if t.name == "front_port")
    pod = det.pods[0]
    members = {s.name.replace("front_port_", ""): s for s in det.struts if s.name.startswith("front_port")}
    return row, det, fp, pod, members


def cap_screws(d: Design, cfg: Config) -> list[tuple[float, float]]:
    """(x, z) of the 4 cap screws (2 per side, staggered) through the sleeve."""
    xc = cap_cut_x(d, cfg)
    x0 = xc - cfg.cap_sleeve_len
    zk, zd = cap_z_range(d, cfg)
    zd = min(zd, zk + 0.12)          # screws in the lower, solid part of the spigot
    H = zd - zk
    return [(x0 + 0.30 * cfg.cap_sleeve_len, zk + 0.35 * H), (x0 + 0.70 * cfg.cap_sleeve_len, zk + 0.70 * H)]


def keep_out_cones(cfg: Config, t):
    """Inlet (forward end, +axis) and outlet (jet, -axis) keep-out frustums."""
    D = 2 * t.radius
    a_in, a_out = math.radians(cfg.inlet_cone_deg), math.radians(cfg.outlet_cone_deg)
    Li, Lo = cfg.inlet_cone_len_D * D, cfg.outlet_cone_len_D * D
    face_in = t.centre + t.axis * t.length / 2
    face_out = t.centre - t.axis * t.length / 2
    inlet = frustum(face_in, t.axis, t.radius, t.radius + Li * math.tan(a_in), Li)
    outlet = frustum(face_out, -t.axis, t.radius, t.radius + Lo * math.tan(a_out), Lo)
    return inlet, outlet


def front_frame(d: Design, cfg: Config, spec: FrontSpec, fp, pod, m: dict):
    cq = _cq()
    xc = cap_cut_x(d, cfg)
    x0 = xc - cfg.cap_sleeve_len
    zk, _ = cap_z_range(d, cfg)
    y0 = pod.y
    rail_ = m["rail"]
    zt = cap_top(d, cfg, rail_.p0[2] - rail_.thickness / 2)
    # --- bow cap: hull lines over the sleeve, then a wedge at cap_entry_deg to the tip
    b_cut = float(halfbeam(np.array([xc / d.loa]), d)[0])
    tan_e = math.tan(math.radians(cfg.cap_entry_deg))
    outer = [hull_loop(d, x0, y0, z_top=zt), hull_loop(d, xc, y0, z_top=zt)]
    xs = np.linspace(xc, d.loa, 6)[1:]
    for x in xs:
        s = max(b_cut - (x - xc) * tan_e, 0.0015) / b_cut
        outer.append(hull_loop(d, x, y0, scale_y=s, zk=zk, z_top=zt))
    cap = _loft(outer)
    spigot = spigot_solid(d, cfg, y0, zt + cfg.cap_wall + 0.002, x_start=x0 - 0.001)
    cap = cap.cut(spigot)      # struts start at the cut and the fin stays below, so none enter it
    # cap screw holes (through both sleeve walls, along y)
    for x, z in cap_screws(d, cfg):
        cap = cap.cut(teardrop((x, y0, z), (0, 1.0, 0), spec.hull_m4, 0.1, max_deg=cfg.max_overhang_deg))
    body = cap
    # --- forefoot fin: from the bed (sleeve's aft end) forward, keel down to the rail
    ff, rail = m["forefoot"], m["rail"]
    z_bot = rail.p0[2] - rail.thickness / 2
    x_ff1 = ff.p0[0] + ff.chord / 2
    fin = (cq.Workplane("XZ").workplane(offset=-(y0 + ff.thickness / 2) * MM)
           .polyline([(x0 * MM, (zk + 0.002) * MM), (x0 * MM, z_bot * MM), (x_ff1 * MM, z_bot * MM),
                      (x_ff1 * MM, (zk + 0.002) * MM)]).close().extrude(ff.thickness * MM))
    body = body.union(fin)
    # --- struts and rail (run a little into the cap / plate to fuse)
    ua, ub = m["upper_a"], m["upper_b"]
    up = np.array([0, 0, 1.0])
    u = (ua.p1 - ua.p0) / ua.length
    # upper_a runs a little past the joint so upper_b's start face is not left hanging
    body = body.union(bar3d(ua.p0, ua.p1 + u * 0.012, up, ua.chord, ua.thickness))
    body = body.union(bar3d(ub.p0, ub.p1, up, ub.chord, ub.thickness))
    # rail starts on the bed (the sleeve's aft end), alongside the forefoot fin
    body = body.union(bar3d(np.array([x0, rail.p0[1], rail.p0[2]]), rail.p1, [0, 1.0, 0], rail.chord,
                            rail.thickness))
    # --- T200 plate (square, rotated with the thruster) with 4x M3 and a cable hole
    pc = mount_plate_centre(cfg, fp)
    # kite: covers the screw pattern and the strut landing; its aft edges sit at the
    # overhang angle (a square turned 45 deg would overhang 45 deg nose-up)
    hs = abs(ua.p1[1] - pc[1]) + ua.thickness / 2 + 0.002
    ha = hs / math.tan(math.radians(cfg.max_overhang_deg))
    hf = max(hs, cfg.t200_bolt_square / math.sqrt(2) + 0.006)
    plate = (cq.Workplane("XY").workplane(offset=(pc[2] - cfg.mount_plate_t / 2) * MM)
             .polyline([((pc[0] + hf) * MM, pc[1] * MM), (pc[0] * MM, (pc[1] + hs) * MM),
                        ((pc[0] - ha) * MM, pc[1] * MM), (pc[0] * MM, (pc[1] - hs) * MM)]).close()
             .extrude(cfg.mount_plate_t * MM))
    body = body.union(plate)
    c, s_ = math.cos(fp.heading), math.sin(fp.heading)
    q = cfg.t200_bolt_square / 2
    m3 = [(pc[0] + c * uu - s_ * vv, pc[1] + s_ * uu + c * vv) for uu in (-q, q) for vv in (-q, q)]
    for (x, y), dia in [(p, spec.t200_m3) for p in m3] + [((pc[0], pc[1]), spec.cable_hole)]:
        body = body.cut(teardrop((x, y, pc[2]), (0, 0, 1.0), dia, cfg.mount_plate_t + 0.004,
                                 max_deg=cfg.max_overhang_deg))
    # --- lugs at the strut ends, with M3 inserts facing forward (the fairing screws on)
    lugs = [ub.p1, rail.p1]
    for p in lugs:
        blk = cq.Workplane("XY").box(cfg.hoop_len * MM, spec.lug * MM, spec.lug * MM).translate(tuple(p * MM))
        body = body.union(blk)
        ins = (cq.Workplane("YZ").workplane(offset=(p[0] + cfg.hoop_len / 2 - spec.insert_m3_depth) * MM)
               .center(p[1] * MM, p[2] * MM).circle(spec.insert_m3 * MM / 2).extrude(spec.insert_m3_depth * MM))
        body = body.cut(ins)
    return body, m3, lugs


def front_fairing(d: Design, cfg: Config, spec: FrontSpec, pod, lugs, wall: float, inlet_cone):
    cq = _cq()
    xb = pod.x_bumper                     # aft rim, just ahead of the lugs
    xn = xb + d.pod_nose_len              # start of the nose closure
    tan_o = math.tan(math.radians(cfg.max_overhang_deg))
    def semi(x):
        if x <= xn:
            return pod.half_w, pod.half_h
        w = max(pod.half_w - (x - xn) * tan_o, cfg.bumper_tip_r)
        return w, pod.half_h * w / pod.half_w
    xs = [xb, xn] + list(np.linspace(xn, pod.x1, 14)[1:])
    outer = _loft([ellipse_wire(x, pod.y, pod.z, *semi(x)) for x in xs])
    xi = [x for x in xs if min(semi(x)) - wall > 0.003]
    inner = _loft([ellipse_wire(x, pod.y, pod.z, semi(x)[0] - wall, semi(x)[1] - wall)
                   for x in [xb - 0.001] + xi[1:]])
    shell = outer.cut(inner)
    # tabs at the rim, in line with the lugs, with M3 clearance holes along x
    for p in lugs:
        tab = cq.Workplane("XY").box(0.006 * MM, spec.lug * MM, spec.lug * MM).translate(
            ((xb + 0.003) * MM, p[1] * MM, p[2] * MM))
        tab = tab.intersect(outer)
        shell = shell.union(tab)
        hole = (cq.Workplane("YZ").workplane(offset=(xb - 0.001) * MM).center(p[1] * MM, p[2] * MM)
                .circle(spec.t200_m3 * MM / 2).extrude(0.01 * MM))
        shell = shell.cut(hole)
    # inlet window: tall slots with pointed tops wherever the inlet keep-out cone
    # passes through the skin (a single big opening would need a long bridge)
    for sl in inlet_slots(cfg, pod, semi, xb, wall):
        shell = shell.cut(sl)
    return shell


def inlet_slots(cfg: Config, pod, semi, xb: float, wall: float, width: float = 0.008, rib: float = 0.004):
    """Slots through the fairing skin covering the inlet keep-out cone's footprint."""
    cq = _cq()
    t = pod.thr if hasattr(pod, "thr") else None
    if t is None:
        return []
    D = 2 * t.radius
    a_in = math.radians(cfg.inlet_cone_deg)
    Li = cfg.inlet_cone_len_D * D
    face = t.centre + t.axis * t.length / 2

    def in_cone(p):
        v = p - face
        s = v @ t.axis
        if s < 0 or s > Li:
            return False
        return np.linalg.norm(v - s * t.axis) <= t.radius + s * math.tan(a_in)

    xs = np.linspace(xb + 0.004, pod.x1, 120)
    out = []
    per = np.pi * (3 * (pod.half_w + pod.half_h) - np.sqrt((3 * pod.half_w + pod.half_h) * (pod.half_w + 3 * pod.half_h)))
    n_slots = int(per / (width + rib))
    tan_o = math.tan(math.radians(cfg.max_overhang_deg))
    for j in range(n_slots):
        th = 2 * np.pi * (j + 0.5) / n_slots
        hit = []
        for x in xs:
            a, b = semi(x)
            p = np.array([x, pod.y + a * math.cos(th), pod.z + b * math.sin(th)])
            if in_cone(p):
                hit.append(x)
        if len(hit) < 3:
            continue
        x0, x1 = min(hit), max(hit)
        if x1 - x0 < 0.012:
            continue
        xm = 0.5 * (x0 + x1)
        a, b = semi(xm)
        S = np.array([xm, pod.y + a * math.cos(th), pod.z + b * math.sin(th)])
        n = np.array([0.0, math.cos(th) / a, math.sin(th) / b])
        n /= np.linalg.norm(n)
        top = (width / 2) / tan_o          # pointed top at the overhang angle
        half = (x1 - x0) / 2
        pl = cq.Plane(origin=tuple(S * MM), xDir=(1, 0, 0), normal=tuple(n))
        pts = [(-half * MM, -width / 2 * MM), ((half - top) * MM, -width / 2 * MM), (half * MM, 0),
               ((half - top) * MM, width / 2 * MM), (-half * MM, width / 2 * MM)]
        out.append(cq.Workplane(pl).polyline(pts).close().extrude(0.015 * MM, both=True))
    return out


def spigot_solid(d: Design, cfg: Config, y0: float, z_top: float, x_start: float):
    """The hull spigot the cap slides over: hull lines inset by the cap wall, ending
    in a wedge (plan view) at the overhang angle, so the cap's cavity has no flat
    ceiling when it prints nose-up."""
    xc = cap_cut_x(d, cfg)
    # wedge starts where the inset half-width would close at the overhang angle by xc
    tan_o = math.tan(math.radians(cfg.max_overhang_deg))
    xs_ = np.linspace(xc - 0.06, xc, 241)
    bw_ = halfbeam(xs_ / d.loa, d) - cfg.cap_wall
    ok = xs_ + bw_ / tan_o >= xc
    xw = float(xs_[ok][0]) if ok.any() else xc - 0.02
    bw = float(halfbeam(np.array([xw / d.loa]), d)[0]) - cfg.cap_wall
    taper = xc - xw
    secs = [hull_loop(d, x_start, y0, inset=cfg.cap_wall, z_top=z_top),
            hull_loop(d, xw, y0, inset=cfg.cap_wall, z_top=z_top)]
    for x in np.linspace(xw, xc, 6)[1:]:
        half = max(1.0 - (x - xw) / taper, 0.01) * bw
        secs.append(_scaled_loop(d, x, y0, cfg.cap_wall, z_top, half))
    return _loft(secs)


def _scaled_loop(d: Design, x: float, y0: float, inset: float, z_top: float, half: float):
    """Inset hull section at x squeezed in y to the given half-width."""
    xi = np.clip(x / d.loa, 0, 1)
    b = float(halfbeam(np.array([xi]), d)[0]) - inset
    return hull_loop(d, x, y0, inset=inset, scale_y=1.0, z_top=z_top) if half >= b else \
        hull_loop(d, x, y0, inset=inset, scale_y=(half + inset) / (b + inset), z_top=z_top)


def hull_bow_spigot(d: Design, cfg: Config, spec: FrontSpec, y0: float, length: float = 0.06):
    """Reference: the end of the printed hull, cut blunt, with the spigot the cap slides
    over and the M4 heat-set insert holes (from both sides)."""
    cq = _cq()
    xc = cap_cut_x(d, cfg)
    x0 = xc - cfg.cap_sleeve_len
    hull = _loft([hull_loop(d, x0 - length, y0), hull_loop(d, x0, y0)])
    spig = spigot_solid(d, cfg, y0, d.depth, x_start=x0)
    part = hull.union(spig)
    for x, z in cap_screws(d, cfg):
        for side in (-1, 1):
            # insert hole from each side face inward
            xi = np.clip(x / d.loa, 0, 1)
            b = float(halfbeam(np.array([xi]), d)[0]) - cfg.cap_wall
            ys = y0 + side * b
            h = (cq.Workplane("XZ").workplane(offset=-(ys) * MM).center(x * MM, z * MM)
                 .circle(spec.insert_m4 * MM / 2).extrude(-side * spec.insert_m4_depth * MM))
            part = part.cut(h)
    return part


# ------------------------------------------------------------------ checks
def overhang_report(shape, max_deg: float, up=(1.0, 0, 0), tol=0.3) -> dict:
    """Area (mm^2) of downward-facing surface steeper than max_deg from vertical,
    printing with `up` as the build direction; faces on the bed are excluded."""
    v, f = shape.val().tessellate(0.3, 0.3)
    V = np.array([[p.x, p.y, p.z] for p in v])
    F = np.array(f)
    up = np.asarray(up, float)
    tri = V[F]
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    area = 0.5 * np.linalg.norm(n, axis=1)
    nz = (n @ up) / np.maximum(2 * area, 1e-12)
    h = tri @ up
    on_bed = h.max(1) < h.min() + tol
    bad = (nz < -math.sin(math.radians(max_deg + 0.5))) & ~on_bed
    return {"overhang_mm2": round(float(area[bad].sum()), 1),
            "worst_deg": round(float(np.degrees(np.arcsin(np.clip(-nz[~on_bed].min(), 0, 1)))), 1)}


def fits_bed(shape, ps: PrintSettings, up_axis: int = 0) -> tuple[bool, tuple]:
    bb = shape.val().BoundingBox()
    dims = [bb.xlen, bb.ylen, bb.zlen]
    h = dims[up_axis]
    foot = sorted(dims[i] for i in range(3) if i != up_axis)
    ok = h <= ps.bed[2] and foot[1] <= max(ps.bed[:2]) and foot[0] <= min(ps.bed[:2])
    return ok, tuple(round(x, 1) for x in dims)


def front_checks(cfg: Config, spec: FrontSpec, fp, frame_part, fairing, inlet, outlet, m3) -> dict:
    cq = _cq()
    thr = thruster_solid(fp)
    out = {
        "frame_vs_thruster_mm3": _vol(frame_part, thr),
        "fairing_vs_thruster_mm3": _vol(fairing, thr),
        "frame_vs_fairing_mm3": _vol(frame_part, fairing),
        "frame_in_inlet_cone_mm3": _vol(frame_part, inlet),
        "frame_in_outlet_cone_mm3": _vol(frame_part, outlet),
        "fairing_in_outlet_cone_mm3": _vol(fairing, outlet),
    }
    # removal: drop, then slide outboard (+y for port) clear of everything
    worst = 0.0
    side = 1.0 if fp.centre[1] >= 0 else -1.0
    for dy in np.linspace(0.0, 0.16, 9):
        t = thruster_solid(fp, drop=spec.removal_clear, dy=side * dy)
        worst = max(worst, _vol(frame_part, t), _vol(fairing, t))
    out["removal_path_max_mm3"] = worst
    # screw access straight down onto each M3 (above the plate top)
    pc = mount_plate_centre(cfg, fp)
    blocked = 0
    for x, y in m3:
        z0 = pc[2] + cfg.mount_plate_t / 2 + 0.001
        cyl = cq.Workplane("XY").circle(spec.access_d * MM / 2).extrude(0.15 * MM).translate(
            (x * MM, y * MM, z0 * MM))
        if _vol(frame_part, cyl) + _vol(fairing, cyl) > 1.0:
            blocked += 1
    out["m3_screws_blocked"] = blocked
    return {k: (round(v, 1) if isinstance(v, float) else v) for k, v in out.items()}


# ------------------------------------------------------------------ strength
def front_frame_model(cfg: Config, fp, pod, m: dict, wall: float):
    """Frame model of the printed front part: two fixed roots in the cap (upper strut,
    forefoot), the T200 hanging from the plate node, the fairing tying the two lugs."""
    f = Frame()
    ua, ub, ff, rl = m["upper_a"], m["upper_b"], m["forefoot"], m["rail"]
    n_root = f.node(ua.p0)
    n_plate = f.node(ua.p1)
    n_top = f.node(ub.p1)
    n_ff_top = f.node(ff.p1)
    n_rail0 = f.node(rl.p0)
    n_rail1 = f.node(rl.p1)
    n_thr = f.node(fp.centre)
    n_nose = f.node([pod.x1 - 0.01, pod.y, pod.z])
    f.fixed |= {n_root, n_ff_top}
    sec_u = Section(ua.thickness, ua.chord, wall)
    f.add("upper_a", n_root, n_plate, sec_u, [0, 0, 1])
    f.add("upper_b", n_plate, n_top, sec_u, [0, 0, 1])
    f.add("forefoot", n_ff_top, n_rail0, Section(ff.thickness, ff.chord, wall), [1, 0, 0])
    f.add("rail", n_rail0, n_rail1, Section(rl.thickness, rl.chord, wall), [0, 1, 0])
    stiff = Section(0.04, 0.04, 0.02)        # plate + T200 base: treated as rigid
    f.add("plate_link", n_plate, n_thr, stiff, [1, 0, 0])
    shell = Section(0.010, 0.060, wall)      # fairing skin between lug and nose (stiff-ish)
    f.add("fairing_top", n_top, n_nose, shell, [0, 0, 1])
    f.add("fairing_bottom", n_rail1, n_nose, shell, [0, 0, 1])
    return f, {"thruster": n_thr, "nose": n_nose, "root": n_root, "forefoot": n_ff_top}


def front_strength(d: Design, cfg: Config, spec: FrontSpec, fp, pod, m: dict, ps: PrintSettings):
    mat = Material()
    f, n = front_frame_model(cfg, fp, pod, m, ps.wall)
    g = cfg.g
    fwd = cfg.thr_force_kgf * g * spec.hold_factor
    rev = cfg.thr_rev_force_kgf * g * spec.hold_factor
    breakaway = spec.fairing_screws * spec.fuse_shear_each
    cases = [
        (f"forward thrust x{spec.hold_factor:g} ({fwd:.0f} N)", n["thruster"], fwd * fp.axis, True),
        (f"reverse thrust x{spec.hold_factor:g} ({rev:.0f} N)", n["thruster"], -rev * fp.axis, True),
        (f"head-on hit x{spec.impact_factor:g} breakaway ({breakaway * spec.impact_factor:.0f} N)",
         n["nose"], np.array([-breakaway * spec.impact_factor, 0, 0]), True),
        (f"side hit, secondary, x{spec.side_factor:g} breakaway ({breakaway * spec.side_factor:.0f} N)",
         n["nose"], np.array([0, -breakaway * spec.side_factor, 0]), False),
    ]
    # cap screws: (x, z) pairs on both sides -> 4 points in the cap's mid plane
    pts = [(x, z) for x, z in cap_screws(d, cfg)] * 2
    cx, cz = np.mean([p[0] for p in pts]), np.mean([p[1] for p in pts])
    r = np.array([[x - cx, z - cz] for x, z in pts])
    out = []
    for name, node, F, primary in cases:
        u, res = f.solve({node: (F, np.zeros(3))}, mat)
        reac = f.reactions(u, mat)
        Ft, Mt = np.zeros(3), np.zeros(3)
        for i, (fr, mo) in reac.items():
            p = f.nodes[i] - np.array([0, 0, 0])
            Ft += -fr                                   # load the frame puts on the cap
            Mt += -mo + np.cross(p - np.array([cx, pod.y, cz]), -fr)
        # screws take the in-plane (x, z) force and the pitch moment as shear
        f_xz = np.hypot(Ft[0], Ft[2]) / len(pts)
        my = abs(Mt[1]) * np.linalg.norm(r, axis=1).max() / (r**2).sum()
        out.append({"case": name, "primary": primary,
                    "members": [{"name": e.name, "axial_N": round(e.axial, 1), "sigma_MPa": round(e.sigma / 1e6, 2),
                                 "utilisation": round(e.utilisation, 3), "buckling": round(e.buckling, 3)}
                                for e in res if e.name in ("upper_a", "upper_b", "forefoot", "rail")],
                    "screw_shear_max_N": round(f_xz + my, 1),
                    "cap_moment_Nm": [round(x, 1) for x in Mt]})
    return out, breakaway, mat


def size_members(d: Design, cfg: Config, spec: FrontSpec, ps: PrintSettings) -> tuple[Config, list]:
    """Smallest upper strut, forefoot and rail (>= 12 x 24 mm) that pass the primary cases."""
    best = None
    tried = []
    for ut, uh, ft, lw in itertools.product([0.012, 0.016, 0.020, 0.024], [0.024, 0.032, 0.040, 0.048],
                                            [0.012, 0.016, 0.020], [0.024, 0.040]):
        lt = 0.012
        c = replace(cfg, front_upper_thick=ut, front_upper_depth=uh, front_lower_thick=lt,
                    front_lower_width=lw, forefoot_thick=ft)
        try:
            _, _, fp, pod, m = front_geometry(d, c)
        except ValueError:
            continue
        cases, _, _ = front_strength(d, c, spec, fp, pod, m, ps)
        worst = max(max(x["utilisation"] for x in k["members"]) for k in cases if k["primary"])
        buck = max(max(x["buckling"] for x in k["members"]) for k in cases if k["primary"])
        screws = max(k["screw_shear_max_N"] for k in cases if k["primary"])
        area = ut * uh + lt * lw + 0.5 * ft * c.forefoot_chord
        ok = worst <= 1.0 and buck <= 0.5 and screws <= spec.insert_shear
        tried.append((ut, uh, ft, lw, worst, buck, screws, ok))
        if ok and (best is None or area < best[0]):
            best = (area, c)
    return (best[1] if best else cfg), tried


# ------------------------------------------------------------------ slicing
def slice_mass(stl: str, ps: PrintSettings) -> dict:
    """Slice with PrusaSlicer; return grams, cm^3 and print time (or None if unavailable)."""
    exe = shutil.which("prusa-slicer") or shutil.which("PrusaSlicer")
    if not exe:
        return {"g": None, "cm3": None, "time": None, "note": "PrusaSlicer not installed"}
    gcode = stl[:-4] + ".gcode"
    cmd = [exe, "--export-gcode", "--perimeters", str(ps.perimeters),
           "--fill-density", f"{ps.infill}%", "--fill-pattern", ps.pattern,
           "--layer-height", str(ps.layer), "--nozzle-diameter", str(ps.nozzle),
           "--filament-diameter", "1.75", "--filament-density", str(ps.density),
           "--top-solid-layers", str(ps.top_layers), "--bottom-solid-layers", str(ps.bottom_layers),
           "--bed-shape", f"0x0,{ps.bed[0]}x0,{ps.bed[0]}x{ps.bed[1]},0x{ps.bed[1]}",
           "--max-print-height", str(ps.bed[2]),
           "--center", f"{ps.bed[0] / 2},{ps.bed[1] / 2}", "-o", gcode, stl]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    out = {"g": None, "cm3": None, "time": None, "note": ""}
    if r.returncode != 0 or not os.path.exists(gcode):
        lines = (r.stderr or r.stdout).strip().splitlines()
        out["note"] = lines[-1] if lines else "slicing failed"
        return out
    txt = open(gcode, errors="ignore").read()[-200_000:]
    for key, pat in (("g", r"filament used \[g\] = ([\d.]+)"), ("cm3", r"filament used \[cm3\] = ([\d.]+)")):
        mm = re.search(pat, txt)
        out[key] = float(mm.group(1)) if mm else None
    mm = re.search(r"estimated printing time \(normal mode\) = (.+)", txt)
    out["time"] = mm.group(1).strip() if mm else None
    os.remove(gcode)
    return out


def export_print_oriented(shape, path: str) -> None:
    """STL rotated so the print direction (+x body) is +z, as it sits on the bed."""
    cq = _cq()
    s = shape.rotate((0, 0, 0), (0, 1, 0), -90)   # +x -> +z
    cq.exporters.export(s, path, tolerance=0.05, angularTolerance=0.2)


def _to_mesh(shape, tol=0.2) -> Mesh:
    v, f = shape.val().tessellate(tol)
    V = np.array([[p.x, p.y, p.z] for p in v]) / MM
    return Mesh(V, np.array(f, dtype=np.int64))


# ------------------------------------------------------------------ run
def build_front(design: Design, cfg: Config, out: str, ps: PrintSettings | None = None,
                spec: FrontSpec | None = None, slice_parts: bool = True, size: bool = True) -> dict:
    cq = _cq()
    ps = ps or PrintSettings()
    spec = spec or FrontSpec()
    os.makedirs(out, exist_ok=True)
    tried = []
    if size:
        cfg, tried = size_members(design, cfg, spec, ps)
    row, det, fp, pod, m = front_geometry(design, cfg)
    wall = max(design.pod_wall, ps.wall)
    inlet, outlet = keep_out_cones(cfg, fp)
    frame_part, m3, lugs = front_frame(design, cfg, spec, fp, pod, m)
    pod.thr = fp
    fairing = front_fairing(design, cfg, spec, pod, lugs, wall, inlet)
    spigot = hull_bow_spigot(design, cfg, spec, pod.y)
    parts = {"front_frame": frame_part, "front_fairing": fairing}
    report = {"design": design.to_dict(),
              "sizes_mm": {"upper": [cfg.front_upper_thick * MM, cfg.front_upper_depth * MM],
                           "rail": [cfg.front_lower_width * MM, cfg.front_lower_thick * MM],
                           "forefoot": cfg.forefoot_thick * MM},
              "sized_ok": (not size) or any(t[-1] for t in tried),
              "parts": {}}
    for name, shp in parts.items():
        cq.exporters.export(shp, os.path.join(out, f"{name}.step"))
        stl = os.path.join(out, f"{name}.stl")
        cq.exporters.export(shp, stl, tolerance=0.05, angularTolerance=0.2)       # boat frame
        pstl = os.path.join(out, f"{name}_print.stl")
        export_print_oriented(shp, pstl)                                          # as printed
        ok, dims = fits_bed(shp, ps, up_axis=0)
        info = {"bbox_mm": list(dims), "fits_bed_nose_up": ok, "valid": shp.val().isValid(),
                "volume_cm3": round(shp.val().Volume() / 1000, 1),
                **overhang_report(shp, cfg.max_overhang_deg)}
        if slice_parts:
            info["slicer"] = slice_mass(pstl, ps)
        report["parts"][name] = info
    cq.exporters.export(spigot, os.path.join(out, "hull_bow_spigot.step"))
    cq.exporters.export(spigot, os.path.join(out, "hull_bow_spigot.stl"), tolerance=0.05, angularTolerance=0.2)
    report["checks"] = front_checks(cfg, spec, fp, frame_part, fairing, inlet, outlet, m3)
    # viewing assembly: frame, fairing, T200 stand-in, hull spigot, keep-out cones
    from .export import thruster_mesh
    asm = Mesh.merge(_to_mesh(frame_part), _to_mesh(fairing), thruster_mesh(fp), _to_mesh(spigot))
    write_stl(os.path.join(out, "front_assembly.stl"), asm.triangles, "front assembly")
    write_stl(os.path.join(out, "front_keepout_cones.stl"),
              Mesh.merge(_to_mesh(inlet), _to_mesh(outlet)).triangles, "keep-out cones")
    cases, breakaway, mat = front_strength(design, cfg, spec, fp, pod, m, ps)
    side = next(c for c in cases if not c["primary"])
    su = max(x["utilisation"] for x in side["members"])
    side_load = breakaway * spec.side_factor
    cap_side = min(side_load / max(su, 1e-9),
                   side_load * spec.insert_shear / max(side["screw_shear_max_N"], 1e-9))
    report["strength"] = {"material": mat.name, "sigma_allow_MPa": mat.sigma_allow / 1e6,
                          "breakaway_N": breakaway, "insert_shear_N": spec.insert_shear,
                          "side_capacity_N": round(cap_side), "cases": cases,
                          "sizing_tried": len(tried)}
    report["model_estimates_kg"] = {"pod": row["pod_mass_each"], "mounts_total_all_four": row["mass_mounts"]}
    with open(os.path.join(out, "front_report.json"), "w") as fh:
        json.dump(report, fh, indent=2, default=float)
    with open(os.path.join(out, "front_report.md"), "w") as fh:
        fh.write(report_md(report, ps, spec, cfg))
    return report


def report_md(r: dict, ps: PrintSettings, spec: FrontSpec, cfg: Config) -> str:
    sz = r["sizes_mm"]
    L = ["# Front frame and nose fairing", "",
         f"Print: {ps.material}, {ps.perimeters} walls ({ps.wall * MM:.1f} mm), {ps.infill} % {ps.pattern}, "
         f"{ps.layer} mm layers, {ps.bed[0]} x {ps.bed[1]} x {ps.bed[2]} mm bed, nose-up (boat +x up), "
         f"no supports, overhang limit {cfg.max_overhang_deg:g}° from vertical.", "",
         f"Sized members: upper strut {sz['upper'][0]:g} x {sz['upper'][1]:g} mm (thick x deep); "
         f"rail {sz['rail'][0]:g} x {sz['rail'][1]:g} mm (wide x thick); forefoot {sz['forefoot']:g} mm thick."
         + ("" if r.get("sized_ok", True) else " NO SIZE IN THE SEARCH PASSED: these are the defaults."), "",
         "## Parts", "",
         "| part | size x,y,z (mm) | fits bed nose-up | overhang > limit (mm²) | slicer mass (g) | print time |",
         "|---|---|---|---|---|---|"]
    for name, p in r["parts"].items():
        s = p.get("slicer", {})
        L.append(f"| {name} | {' x '.join(str(x) for x in p['bbox_mm'])} | {'yes' if p['fits_bed_nose_up'] else 'NO'} | "
                 f"{p['overhang_mm2']} | {s.get('g')} | {s.get('time')} |")
    est = r["model_estimates_kg"]
    L += ["", f"Optimizer's estimate: pod {est['pod'] * 1000:.0f} g; all four mounts "
              f"{est['mounts_total_all_four'] * 1000:.0f} g.", "", "## Checks", ""]
    names = {"frame_vs_thruster_mm3": "frame overlaps thruster", "fairing_vs_thruster_mm3": "fairing overlaps thruster",
             "frame_vs_fairing_mm3": "frame overlaps fairing", "frame_in_inlet_cone_mm3": "frame inside inlet keep-out",
             "frame_in_outlet_cone_mm3": "frame inside jet keep-out",
             "fairing_in_outlet_cone_mm3": "fairing inside jet keep-out",
             "removal_path_max_mm3": "thruster removal path blocked by", "m3_screws_blocked": "M3 screws without tool access"}
    L += [f"- {names[k]}: {v}" + (" mm³" if k.endswith("mm3") else "") for k, v in r["checks"].items()]
    st = r["strength"]
    L += ["", "## Strength", "",
          f"Allowable {st['sigma_allow_MPa']:.0f} MPa (printed PLA across layers, walls only). "
          f"Fairing breaks off at about {st['breakaway_N']:.0f} N ({spec.fairing_screws} nylon M3). "
          f"Cap screws: M4 into heat-set inserts, {st['insert_shear_N']:.0f} N each assumed in shear (test this). "
          f"Side hits are secondary: the frame takes about {st['side_capacity_N']} N sideways at the nose.", "",
          "| load case | member | axial (N) | stress (MPa) | utilisation | buckling | max cap screw shear (N) |",
          "|---|---|---|---|---|---|---|"]
    for c in st["cases"]:
        for i, mm in enumerate(c["members"]):
            L.append(f"| {c['case'] if i == 0 else ''} | {mm['name']} | {mm['axial_N']} | {mm['sigma_MPa']} | "
                     f"{mm['utilisation']} | {mm['buckling']} | {c['screw_shear_max_N'] if i == 0 else ''} |")
    return "\n".join(L) + "\n"


def main(argv=None):
    p = argparse.ArgumentParser(prog="hullopt.parts", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--design", help="design JSON (defaults if omitted)")
    p.add_argument("--config")
    p.add_argument("--out", default="results/parts")
    p.add_argument("--no-slice", action="store_true")
    p.add_argument("--no-size", action="store_true", help="keep the configured member sizes")
    a = p.parse_args(argv)
    d = Design()
    if a.design:
        with open(a.design) as f:
            d = Design.from_dict(json.load(f))
    cfg = Config.load(a.config)
    r = build_front(d, cfg, a.out, slice_parts=not a.no_slice, size=not a.no_size)
    print(report_md(r, PrintSettings(), FrontSpec(), cfg))
    print(f"wrote {a.out}/")


if __name__ == "__main__":
    main()
