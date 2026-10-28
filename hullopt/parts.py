"""Printable parametric parts (CadQuery), placed in the boat's body frame.

Front (port side; starboard is the mirror image):
  front_bracket  back plate on the flat bow face (4x M4, 40 x 60 mm), two members
                 (12 x 24 mm min) down to a plate the T200 bolts to (4x M3 on a
                 19 mm square + M6 centre), and two M3 heat-set holes for the pod.
  front_pod      faired shell over the thruster and bracket: solid nose that takes
                 head-on hits, side windows for inflow, open tail and a top slot for
                 the bracket members, two bosses that screw onto the bracket plate.
                 Fixed with nylon M3 screws, which shear first: the pod breaks off
                 before the bracket or the bow is damaged.

Every part is exported as STEP (edit in CAD) and STL, sliced with PrusaSlicer
for its real mass, and checked for strength (hullopt.strength).

  python -m hullopt.parts --design results/stage2/designs/finalist_1.json --out results/parts
"""
from __future__ import annotations

import argparse
import json
import math
import os
import re
import shutil
import subprocess
from dataclasses import dataclass

import numpy as np

from .export import pod_profile, thruster_mesh
from .geometry import Mesh, write_stl
from .layout import bow_face, mount_plate_centre
from .params import Config, Design
from .stage1 import evaluate
from .strength import Material, Section, analyse

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
    """Hardware and part dimensions for the front bracket and pod (metres)."""

    hull_bolt: float = 0.0045      # M4 clearance
    t200_m3: float = 0.0034        # M3 clearance
    t200_m6: float = 0.0065        # M6 clearance
    insert_m3: float = 0.0040      # M3 heat-set insert hole
    insert_depth: float = 0.0055
    fuse_offset: float = 0.022     # pod screws, either side of the plate centre (body y)
    fuse_boss_d: float = 0.009
    fuse_screws: int = 2           # nylon M3, shear pins
    fuse_shear_each: float = 210.0  # N, nylon (PA66) M3 in single shear, approx.
    slot_clear: float = 0.001      # each side of the members in the pod's top slot
    plate_corner_r: float = 0.003
    hold_factor: float = 3.0       # bracket must hold this x full thrust
    impact_factor: float = 1.25    # and this x the fuse breakaway load
    insert_pullout: float = 800.0  # N, M4 heat-set insert in the hull, placeholder (test it)


# ------------------------------------------------------------------ CAD helpers
def _cq():
    import cadquery as cq
    return cq


def _v(p) -> tuple:
    return tuple(float(c) * MM for c in p)


def bar(p0, p1, width_y: float, depth: float):
    """Rectangular bar from p0 to p1 (body frame, m), width across y, depth in the
    vertical plane through the bar; rounded long edges."""
    cq = _cq()
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    d = p1 - p0
    L = float(np.linalg.norm(d))
    ang = math.degrees(math.atan2(d[2], d[0]))
    b = (cq.Workplane("XY").box(L * MM, width_y * MM, depth * MM)
         .edges("|X").fillet(min(width_y, depth) * MM * 0.24))
    b = b.rotate((0, 0, 0), (0, 1, 0), -ang).translate(_v(0.5 * (p0 + p1)))
    return b


def ellipse_loft(stations, half_w, half_h, y, z, scale):
    """Solid lofted through ellipses at x stations (m) with section scale factors."""
    cq = _cq()
    wires = []
    for x, s in zip(stations, scale):
        wires.append(cq.Wire.makeEllipse(half_w * s * MM, half_h * s * MM,
                                         cq.Vector(x * MM, y * MM, z * MM),
                                         cq.Vector(1, 0, 0), cq.Vector(0, 1, 0)))
    # ruled (straight between stations): a smooth loft overshoots where the profile
    # changes quickly (short tail, nose)
    return cq.Workplane("XY").add(cq.Solid.makeLoft(wires, True))


# ------------------------------------------------------------------ front parts
def front_geometry(d: Design, cfg: Config):
    row, det = evaluate(d, cfg, detail=True)
    if det is None:
        raise ValueError(f"design does not float: {row.get('fail_reason')}")
    fp = next(t for t in det.thrusters if t.name == "front_port")
    pod = det.pods[0]
    members = [s for s in det.struts if s.name.startswith("front_port")]
    return row, det, fp, pod, members


def front_bracket(d: Design, cfg: Config, spec: FrontSpec, fp, pod, members):
    cq = _cq()
    zb, zt = bow_face(d)
    zc = 0.5 * (zb + zt)
    W, H, t = d.bow_face_width, (zt - zb) - 0.01, cfg.bow_plate_t
    bw, bh = cfg.hull_bolt_pattern
    # back plate on the bow face
    plate = (cq.Workplane("YZ").workplane(offset=d.loa * MM)
             .center(pod.y * MM, zc * MM).rect(W * MM, H * MM).extrude(t * MM)
             .edges("|X").fillet(spec.plate_corner_r * MM))
    holes = [(pod.y + sy * bw / 2, zc + sz * bh / 2) for sy in (-1, 1) for sz in (-1, 1)]
    plate = plate.faces(">X").workplane(centerOption="ProjectedOrigin", origin=(0, 0, 0)) \
        .pushPoints([(y * MM, z * MM) for y, z in holes]).hole(spec.hull_bolt * MM)
    body = plate
    for s in members:
        # bars run into the back plate and onto the top of the T200 plate
        u = (s.p1 - s.p0) / s.length
        p0 = s.p0 - u * cfg.bow_plate_t
        if s.flat:
            body = body.union(bar(p0, s.p1, s.chord, s.thickness))
        else:
            body = body.union(bar(p0, s.p1, s.thickness, s.chord))
    # T200 plate, square rotated with the thruster
    pc = mount_plate_centre(cfg, fp)
    a = cfg.mount_plate_size
    tp = (cq.Workplane("XY").box(a * MM, a * MM, cfg.mount_plate_t * MM)
          .edges("|Z").fillet(spec.plate_corner_r * MM)
          .rotate((0, 0, 0), (0, 0, 1), math.degrees(fp.heading)).translate(_v(pc)))
    body = body.union(tp)
    # T200 holes (in the thruster's frame), M6 centre, pod fuse inserts (body y)
    c, s_ = math.cos(fp.heading), math.sin(fp.heading)
    q = cfg.t200_bolt_square / 2
    m3 = [(pc[0] + c * u - s_ * v, pc[1] + s_ * u + c * v) for u in (-q, q) for v in (-q, q)]
    cut = cq.Workplane("XY")
    zbot = (pc[2] - cfg.mount_plate_t / 2 - 0.05) * MM
    hgt = (cfg.mount_plate_t + 0.10) * MM
    for (x, y), dia in [((x, y), spec.t200_m3) for x, y in m3] + [((pc[0], pc[1]), spec.t200_m6)]:
        cut = cut.union(cq.Workplane("XY").circle(dia * MM / 2).extrude(hgt)
                        .translate((x * MM, y * MM, zbot)))
    ztop = (pc[2] + cfg.mount_plate_t / 2) * MM
    for sy in (-1, 1):
        cut = cut.union(cq.Workplane("XY").circle(spec.insert_m3 * MM / 2).extrude(spec.insert_depth * MM)
                        .translate((pc[0] * MM, (pc[1] + sy * spec.fuse_offset) * MM,
                                    ztop - spec.insert_depth * MM)))
    return body.cut(cut)


def front_pod(d: Design, cfg: Config, spec: FrontSpec, pod, members, wall: float):
    cq = _cq()
    n = 48
    xs = pod.x0 + (pod.x1 - pod.x0) * 0.5 * (1 - np.cos(np.linspace(0, np.pi, n)))
    sc = pod_profile(pod, xs)
    outer = ellipse_loft(xs, pod.half_w, pod.half_h, pod.y, pod.z, sc)
    # cavity: from just behind the tail to the start of the solid nose (open tail)
    nose0 = pod.x1 - pod.nose_len
    xi = np.linspace(pod.x0 - 0.002, nose0, n)
    si = np.maximum(pod_profile(pod, np.clip(xi, pod.x0, None)), 0.3)
    inner = cq.Workplane("XY")
    wires = [cq.Wire.makeEllipse((pod.half_w * s - wall) * MM, (pod.half_h * s - wall) * MM,
                                 cq.Vector(x * MM, pod.y * MM, pod.z * MM), cq.Vector(1, 0, 0), cq.Vector(0, 1, 0))
             for x, s in zip(xi, si)]
    inner = inner.add(cq.Solid.makeLoft(wires, True))
    shell = outer.cut(inner)
    # side windows on the thruster bay
    bay_len = pod.bay_len
    xc = nose0 - pod.wall - bay_len / 2
    wl = 0.8 * bay_len
    wh = min(pod.window_area / wl, 1.6 * pod.half_h)
    for sy in (-1, 1):
        box = cq.Workplane("XY").box(wl * MM, pod.half_w * MM, wh * MM).translate(
            (xc * MM, (pod.y + sy * pod.half_w * 0.75) * MM, pod.z * MM))
        shell = shell.cut(box)
    # top slot for the bracket members, from the tail forward to where they enter
    x_end = max(float(s.p1[0]) for s in members) + max(s.chord for s in members) / 2
    th = max((s.chord if s.flat else s.thickness) for s in members) + 2 * spec.slot_clear
    slot = cq.Workplane("XY").box((x_end - pod.x0 + 0.004) * MM, th * MM, pod.half_h * MM).translate(
        ((0.5 * (pod.x0 + x_end) - 0.002) * MM, pod.y * MM, (pod.z + pod.half_h * 0.75) * MM))
    shell = shell.cut(slot)
    return shell


def front_pod_bosses(cfg: Config, spec: FrontSpec, fp, pod, wall: float):
    """Bosses from the pod's inner top down to the T200 plate, with M3 clearance."""
    cq = _cq()
    pc = mount_plate_centre(cfg, fp)
    z0 = pc[2] + cfg.mount_plate_t / 2
    out = None
    for sy in (-1, 1):
        y = pc[1] + sy * spec.fuse_offset
        # top of the ellipse at this y
        zt = pod.z + pod.half_h * math.sqrt(max(1 - ((y - pod.y) / pod.half_w) ** 2, 0))
        h = zt - z0
        b = (cq.Workplane("XY").circle(spec.fuse_boss_d * MM / 2).circle(spec.t200_m3 * MM / 2)
             .extrude(h * MM).translate((pc[0] * MM, y * MM, z0 * MM)))
        out = b if out is None else out.union(b)
    return out


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
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    out = {"g": None, "cm3": None, "time": None, "note": ""}
    if r.returncode != 0 or not os.path.exists(gcode):
        out["note"] = (r.stderr or r.stdout).strip().splitlines()[-1:] or ["slicing failed"]
        out["note"] = out["note"][0]
        return out
    txt = open(gcode, errors="ignore").read()[-200_000:]
    for key, pat in (("g", r"filament used \[g\] = ([\d.]+)"), ("cm3", r"filament used \[cm3\] = ([\d.]+)")):
        m = re.search(pat, txt)
        out[key] = float(m.group(1)) if m else None
    m = re.search(r"estimated printing time \(normal mode\) = (.+)", txt)
    out["time"] = m.group(1).strip() if m else None
    os.remove(gcode)
    return out


def fits_bed(shape, ps: PrintSettings) -> tuple[bool, tuple]:
    bb = shape.val().BoundingBox()
    dims = tuple(sorted([bb.xlen, bb.ylen, bb.zlen]))
    return all(dd <= b for dd, b in zip(dims, sorted(ps.bed))), (bb.xlen, bb.ylen, bb.zlen)


# ------------------------------------------------------------------ strength
def front_strength(d: Design, cfg: Config, spec: FrontSpec, fp, pod, members, ps: PrintSettings):
    mat = Material()
    secs = [Section(s.thickness, s.chord, ps.wall) for s in members]   # b thick, h depth
    zb, zt = bow_face(d)
    bw, bh = cfg.hull_bolt_pattern
    bolt = {"origin": np.array([d.loa, pod.y, 0.5 * (zb + zt)]), "normal": np.array([1.0, 0, 0]),
            "axes": (np.array([0, 1.0, 0]), np.array([0, 0, 1.0])),
            "bolts": [(sy * bw / 2, sz * bh / 2) for sy in (-1, 1) for sz in (-1, 1)]}
    pc = mount_plate_centre(cfg, fp)
    g = cfg.g
    fwd = cfg.thr_force_kgf * g * spec.hold_factor
    rev = cfg.thr_rev_force_kgf * g * spec.hold_factor
    breakaway = spec.fuse_screws * spec.fuse_shear_each
    imp = breakaway * spec.impact_factor
    cases = []
    for name, F, at in [
        (f"forward thrust x{spec.hold_factor:g} ({fwd:.0f} N)", fwd * fp.axis, fp.centre),
        (f"reverse thrust x{spec.hold_factor:g} ({rev:.0f} N)", -rev * fp.axis, fp.centre),
        (f"head-on impact at breakaway x{spec.impact_factor:g} ({imp:.0f} N)", np.array([-imp, 0, 0]),
         np.array([pod.x1, pod.y, pod.z])),
        (f"side hit at breakaway x{spec.impact_factor:g} ({imp:.0f} N)", np.array([0, -imp, 0]),
         np.array([pod.x1 - pod.nose_len / 2, pod.y + pod.half_w, pod.z])),
    ]:
        # move the load to the plate centre as force + moment
        M = np.cross(at - pc, F)
        cases.append(analyse(members, secs, pc, F, M, mat, bolt, name))
    return cases, breakaway, mat


# ------------------------------------------------------------------ run
def _to_mesh(shape, tol=0.2) -> Mesh:
    v, f = shape.val().tessellate(tol)
    V = np.array([[p.x, p.y, p.z] for p in v]) / MM
    return Mesh(V, np.array(f, dtype=np.int64))


def build_front(design: Design, cfg: Config, out: str, ps: PrintSettings | None = None,
                spec: FrontSpec | None = None, slice_parts: bool = True) -> dict:
    cq = _cq()
    ps = ps or PrintSettings()
    spec = spec or FrontSpec()
    os.makedirs(out, exist_ok=True)
    row, det, fp, pod, members = front_geometry(design, cfg)
    wall = max(design.pod_wall, ps.wall)
    bracket = front_bracket(design, cfg, spec, fp, pod, members)
    shell = front_pod(design, cfg, spec, pod, members, wall)
    pod_part = shell.union(front_pod_bosses(cfg, spec, fp, pod, wall))
    parts = {"front_bracket": bracket, "front_pod": pod_part}
    report = {"design": design.to_dict(), "parts": {}}
    for name, shp in parts.items():
        step = os.path.join(out, f"{name}.step")
        stl = os.path.join(out, f"{name}.stl")
        cq.exporters.export(shp, step)
        cq.exporters.export(shp, stl, tolerance=0.05, angularTolerance=0.2)
        ok, dims = fits_bed(shp, ps)
        info = {"bbox_mm": [round(x, 1) for x in dims], "fits_bed": ok,
                "volume_cm3": round(shp.val().Volume() / 1000, 1)}
        if slice_parts:
            info["slicer"] = slice_mass(stl, ps)
        report["parts"][name] = info
    # assembly for viewing: both parts + T200 stand-in
    asm = Mesh.merge(_to_mesh(bracket), _to_mesh(pod_part), thruster_mesh(fp))
    write_stl(os.path.join(out, "front_assembly.stl"), asm.triangles, "front assembly")
    cases, breakaway, mat = front_strength(design, cfg, spec, fp, pod, members, ps)
    report["strength"] = {
        "material": mat.name, "sigma_allow_MPa": mat.sigma_allow / 1e6,
        "breakaway_N": breakaway, "insert_pullout_N": spec.insert_pullout,
        "cases": [{"case": c.case, "bolt_tension_max_N": round(c.bolt_tension_max, 1),
                   "bolt_shear_max_N": round(c.bolt_shear_max, 1),
                   "members": [{"name": m.name, "axial_N": round(m.axial, 1),
                                "sigma_MPa": round(m.sigma / 1e6, 2), "tau_MPa": round(m.tau / 1e6, 2),
                                "utilisation": round(m.utilisation, 3), "buckling": round(m.buckling, 3)}
                               for m in c.members]} for c in cases],
    }
    report["model_estimates_kg"] = {"pod": row["pod_mass_each"],
                                    "mounts_total_all_four": row["mass_mounts"]}
    with open(os.path.join(out, "front_report.json"), "w") as f:
        json.dump(report, f, indent=2, default=float)
    with open(os.path.join(out, "front_report.md"), "w") as f:
        f.write(report_md(report, ps, spec))
    return report


def report_md(r: dict, ps: PrintSettings, spec: FrontSpec) -> str:
    L = ["# Front pod and bracket", "",
         f"Print settings: {ps.material}, {ps.perimeters} walls ({ps.wall * MM:.1f} mm), {ps.infill} % "
         f"{ps.pattern}, {ps.layer} mm layers, {ps.bed[0]} x {ps.bed[1]} x {ps.bed[2]} mm bed.", "",
         "## Parts", "", "| part | size (mm) | fits bed | slicer mass (g) | print time | CAD volume (cm³) |",
         "|---|---|---|---|---|---|"]
    for name, p in r["parts"].items():
        s = p.get("slicer", {})
        L.append(f"| {name} | {' x '.join(str(x) for x in p['bbox_mm'])} | {'yes' if p['fits_bed'] else 'NO'} | "
                 f"{s.get('g')} | {s.get('time')} | {p['volume_cm3']} |")
    est = r["model_estimates_kg"]
    L += ["", f"Optimizer's estimate: pod {est['pod'] * 1000:.0f} g; all four mounts {est['mounts_total_all_four'] * 1000:.0f} g.",
          "", "## Strength", "",
          f"Allowable stress {r['strength']['sigma_allow_MPa']:.0f} MPa (printed PLA, across layers). "
          f"Utilisation = combined stress / allowable; buckling = compression / Euler load. "
          f"Both must stay under 1.", "",
          f"Pod breakaway: {spec.fuse_screws} nylon M3 screws in shear ≈ {r['strength']['breakaway_N']:.0f} N. "
          f"Hull inserts: assumed {r['strength']['insert_pullout_N']:.0f} N pull-out each (test this).", "",
          "| load case | member | axial (N) | stress (MPa) | utilisation | buckling | max bolt tension (N) | max bolt shear (N) |",
          "|---|---|---|---|---|---|---|---|"]
    for c in r["strength"]["cases"]:
        for i, m in enumerate(c["members"]):
            L.append(f"| {c['case'] if i == 0 else ''} | {m['name']} | {m['axial_N']} | {m['sigma_MPa']} | "
                     f"{m['utilisation']} | {m['buckling']} | {c['bolt_tension_max_N'] if i == 0 else ''} | "
                     f"{c['bolt_shear_max_N'] if i == 0 else ''} |")
    return "\n".join(L) + "\n"


def main(argv=None):
    p = argparse.ArgumentParser(prog="hullopt.parts", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--design", help="design JSON (defaults if omitted)")
    p.add_argument("--config")
    p.add_argument("--out", default="results/parts")
    p.add_argument("--no-slice", action="store_true")
    a = p.parse_args(argv)
    d = Design()
    if a.design:
        with open(a.design) as f:
            d = Design.from_dict(json.load(f))
    r = build_front(d, Config.load(a.config), a.out, slice_parts=not a.no_slice)
    print(report_md(r, PrintSettings(), FrontSpec()))
    print(f"wrote {a.out}/")


if __name__ == "__main__":
    main()
