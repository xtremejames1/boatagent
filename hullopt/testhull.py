"""Scaled test-hull variant: the same demihull shape at a different length.

All hull lengths (loa, beam, depth, rocker) scale by one factor; angles and
fractions stay. The gap between the hulls' inner faces is set separately. No
thrusters, mounts or skids: just the hull shape, printed in pieces that fit the
bed, and a hydrostatics table (draft against load).

  python -m hullopt.testhull --design results/stage2/designs/finalist_1.json \\
      --loa 0.6096 --spacing 0.3048 --out results/test_hull_2ft
"""
from __future__ import annotations

import argparse
import json
import math
import os
from dataclasses import replace

import numpy as np

from .geometry import catamaran_mesh, demihull_mesh, hull_centre_y, write_stl
from .hydrostatics import Pose, equilibrium, gz_curve, hydro, solve_heave
from .params import Config, Design

INCH = 0.0254


def scaled(d: Design, loa: float, spacing_inner: float) -> tuple[Design, float]:
    k = loa / d.loa
    return replace(d, loa=loa, beam=d.beam * k, depth=d.depth * k, rocker_fwd=d.rocker_fwd * k,
                   rocker_aft=d.rocker_aft * k, spacing_inner=spacing_inner), k


def hull_solid(d: Design, n: int = 61):
    """Closed demihull solid (CadQuery), centred on y = 0, mm."""
    from .parts import _loft, hull_loop
    from .geometry import station_xi
    xs = station_xi(n) * d.loa
    xs[0], xs[-1] = 0.0, d.loa * 0.999      # the very tip is a knife edge
    return _loft([hull_loop(d, float(x), 0.0, n_pts=64) for x in xs])


def pieces(solid, d: Design, bed_len: float = 0.25, margin: float = 0.005):
    """Cut the hull along x into equal pieces no longer than the bed height."""
    import cadquery as cq
    n = math.ceil(d.loa / (bed_len - margin))
    edges = np.linspace(0.0, d.loa, n + 1)
    out = []
    big = 1.0
    for a, b in zip(edges[:-1], edges[1:]):
        box = cq.Workplane("XY").box((b - a) * 1000, big * 1000, big * 1000).translate(
            ((a + b) / 2 * 1000, 0, d.depth / 2 * 1000))
        out.append(((a, b), solid.intersect(box)))
    return out


def hydro_table(d: Design, cfg: Config, masses) -> list[dict]:
    tri = catamaran_mesh(d, cfg.n_stations).triangles
    pivot = np.array([d.loa / 2, 0, 0])
    rows = []
    for m in masses:
        h = solve_heave(tri, m, cfg.rho, Pose(0.0, pivot=pivot), h_lo=-0.02, h_hi=d.depth)
        if h is None:
            rows.append({"mass_kg": round(m, 3), "draft_mm": None})
            continue
        hy = hydro(tri, Pose(h, pivot=pivot))
        lwl = hy.wp.x_max - hy.wp.x_min
        bm_t = hy.wp.i_t / hy.volume
        kb = h + hy.cob[2]                    # keel to centre of buoyancy
        rows.append({"mass_kg": round(m, 3), "draft_mm": round(h * 1000, 1),
                     "freeboard_mm": round((d.depth - h) * 1000, 1), "lwl_mm": round(lwl * 1000),
                     "KB_mm": round(kb * 1000, 1), "BM_T_mm": round(bm_t * 1000),
                     "GM_T_mm_if_KG_at_deck": round((kb + bm_t - d.depth) * 1000),
                     "wetted_m2": round(hy.wetted_area, 4)})
    return rows


def build(d_full: Design, loa: float, spacing: float, out: str, cfg: Config | None = None,
          slice_parts: bool = True) -> dict:
    import cadquery as cq
    from .parts import PrintSettings, overhang_report, slice_mass
    cfg = cfg or Config()
    os.makedirs(out, exist_ok=True)
    d, k = scaled(d_full, loa, spacing)
    with open(os.path.join(out, "design.json"), "w") as f:
        json.dump(d.to_dict(), f, indent=2)
    # meshes for viewing / CFD
    write_stl(os.path.join(out, "demihull.stl"), demihull_mesh(d, cfg.n_stations).triangles, "demihull")
    write_stl(os.path.join(out, "catamaran.stl"), catamaran_mesh(d, cfg.n_stations).triangles, "catamaran")
    solid = hull_solid(d)
    cq.exporters.export(solid, os.path.join(out, "demihull.step"))
    ps = PrintSettings()
    rep = {"scale": k, "design": d.to_dict(), "pieces": []}
    for i, ((a, b), p) in enumerate(pieces(solid, d)):
        name = f"demihull_piece{i + 1}"
        cq.exporters.export(p, os.path.join(out, f"{name}.step"))
        # print standing on its aft cut face: +x up
        stl = os.path.join(out, f"{name}_print.stl")
        cq.exporters.export(p.rotate((0, 0, 0), (0, 1, 0), -90), stl, tolerance=0.05, angularTolerance=0.2)
        bb = p.val().BoundingBox()
        info = {"name": name, "x_mm": [round(a * 1000, 1), round(b * 1000, 1)],
                "size_mm": [round(bb.xlen, 1), round(bb.ylen, 1), round(bb.zlen, 1)],
                **overhang_report(p, cfg.max_overhang_deg)}
        if slice_parts:
            info["slicer"] = slice_mass(stl, ps)
            info["slicer_hollow"] = slice_mass(stl, replace(ps, infill=0))   # walls + end caps only
        rep["pieces"].append(info)
    # loads: around the Froude-scaled design mass (full-scale mass x k^3)
    full_mass = 25.7
    m_design = full_mass * k**3
    rep["froude_mass_kg"] = round(m_design, 2)
    rep["hydrostatics"] = hydro_table(d, cfg, [m_design * f for f in (0.5, 0.75, 1.0, 1.25, 1.5, 2.0)])
    with open(os.path.join(out, "report.json"), "w") as f:
        json.dump(rep, f, indent=2, default=float)
    with open(os.path.join(out, "README.md"), "w") as f:
        f.write(report_md(rep, d, d_full, ps))
    return rep


def report_md(r: dict, d: Design, d_full: Design, ps) -> str:
    k = r["scale"]
    yc = hull_centre_y(d)
    L = [f"# Test hull, {d.loa / INCH:.0f} in ({d.loa * 1000:.0f} mm)", "",
         f"The full-size hull (loa {d_full.loa * 1000:.0f} mm) scaled by {k:.3f}. Hull gap (inner "
         f"faces) {d.spacing_inner * 1000:.1f} mm = {d.spacing_inner / INCH:.1f} in; hull centrelines "
         f"{2 * yc * 1000:.1f} mm apart; overall width {(d.spacing_inner + 2 * d.beam) * 1000:.0f} mm.", "",
         "| | test hull | full size |", "|---|---|---|",
         f"| length | {d.loa * 1000:.0f} mm | {d_full.loa * 1000:.0f} mm |",
         f"| beam (each hull) | {d.beam * 1000:.1f} mm | {d_full.beam * 1000:.1f} mm |",
         f"| depth (keel to deck) | {d.depth * 1000:.1f} mm | {d_full.depth * 1000:.1f} mm |",
         f"| rocker bow / stern | {d.rocker_fwd * 1000:.1f} / {d.rocker_aft * 1000:.1f} mm | "
         f"{d_full.rocker_fwd * 1000:.1f} / {d_full.rocker_aft * 1000:.1f} mm |",
         f"| deadrise, bilge, entry | {d.deadrise_deg:.1f}°, {d.bilge_frac:.2f}, {d.entry_half_angle_deg:.1f}° | same |",
         "", "Angles and proportions are unchanged; the hull gap is not scaled (fixed at 12 in), so "
             "the two hulls sit relatively further apart than on the full-size boat.", "",
         f"Model-scale (Froude) testing: speed scales by √{k:.3f} = {math.sqrt(k):.3f}, so a full-size "
         f"1.0 m/s is {math.sqrt(k):.2f} m/s here. The matching total mass is about "
         f"{r['froude_mass_kg']} kg (full size 25.7 kg × {k:.3f}³).", "",
         "## Draft against load (both hulls, level)", "",
         "| total mass (kg) | draft (mm) | freeboard (mm) | waterline length (mm) | KB (mm) | BM_T (mm) | GM_T if KG at deck (mm) |",
         "|---|---|---|---|---|---|---|"]
    for h in r["hydrostatics"]:
        if h.get("draft_mm") is None:
            L.append(f"| {h['mass_kg']} | sinks | | | | | |")
            continue
        L.append(f"| {h['mass_kg']} | {h['draft_mm']} | {h['freeboard_mm']} | {h['lwl_mm']} | {h['KB_mm']} | "
                 f"{h['BM_T_mm']} | {h['GM_T_mm_if_KG_at_deck']} |")
    L += ["", f"## Print pieces ({ps.material}, {ps.perimeters} walls, {ps.infill} % {ps.pattern}, "
              f"{ps.bed[0]} mm bed)", "",
          "Each hull prints in pieces, standing on the aft cut face (no supports). Print two sets "
          "(port and starboard are identical). The cut faces are flat; glue them, or add pins.", "",
          f"| piece | x range (mm) | size (mm) | overhang > 35° (mm²) | mass, {ps.infill} % infill (g) | "
          "print time | mass, hollow (g) | print time |",
          "|---|---|---|---|---|---|---|---|"]
    tot, toth = 0.0, 0.0
    for p in r["pieces"]:
        s, h = p.get("slicer", {}), p.get("slicer_hollow", {})
        tot += s.get("g") or 0
        toth += h.get("g") or 0
        L.append(f"| {p['name']} | {p['x_mm'][0]}–{p['x_mm'][1]} | {' x '.join(str(v) for v in p['size_mm'])} | "
                 f"{p['overhang_mm2']} | {s.get('g')} | {s.get('time')} | {h.get('g')} | {h.get('time')} |")
    L += ["", f"One hull: {tot:.0f} g with {ps.infill} % infill, {toth:.0f} g hollow (4 walls, closed at the "
              f"cut faces). Two hulls: {2 * tot:.0f} g / {2 * toth:.0f} g, against the "
              f"{r['froude_mass_kg']} kg model-scale total."]
    L += ["", "Files: `demihull.step` (whole hull, solid), `demihull_pieceN.step` / `_print.stl` "
              "(print pieces, standing on the aft cut), `demihull.stl`, `catamaran.stl` (both hulls at "
              "the 12 in gap, for viewing or CFD), `design.json`."]
    return "\n".join(L) + "\n"


def main(argv=None):
    p = argparse.ArgumentParser(prog="hullopt.testhull", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--design", default="results/stage2/designs/finalist_1.json")
    p.add_argument("--loa", type=float, default=24 * INCH, help="test hull length (m)")
    p.add_argument("--spacing", type=float, default=12 * INCH, help="gap between hull inner faces (m)")
    p.add_argument("--out", default="results/test_hull_2ft")
    p.add_argument("--no-slice", action="store_true")
    a = p.parse_args(argv)
    with open(a.design) as f:
        d = Design.from_dict(json.load(f))
    r = build(d, a.loa, a.spacing, a.out, slice_parts=not a.no_slice)
    print(report_md(r, scaled(d, a.loa, a.spacing)[0], d, __import__("hullopt.parts", fromlist=["x"]).PrintSettings()))


if __name__ == "__main__":
    main()
