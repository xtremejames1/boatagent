import json
import math
from collections import Counter

import numpy as np
import pytest

from hullopt.geometry import (
    catamaran_mesh,
    demihull_mesh,
    halfbeam,
    halfbreadth_at,
    keel_z,
    write_stl,
)
from hullopt.hydrostatics import volume_props
from hullopt.mass import print_pieces
from hullopt.params import BOUNDS, Config, Design
from hullopt.stage1 import CONSTRAINTS, evaluate, lhs, main, nondominated_fronts, select_candidates


def test_demihull_is_closed_and_outward():
    m = demihull_mesh(Design())
    edges = Counter()
    for f in m.faces:
        for a, b in ((f[0], f[1]), (f[1], f[2]), (f[2], f[0])):
            edges[(a, b)] += 1
    # every directed edge appears once and its reverse once: closed, consistently oriented
    assert all(n == 1 for n in edges.values())
    assert all((b, a) in edges for (a, b) in edges)
    V, c = volume_props(m.triangles)
    assert V > 0
    d = Design()
    assert 0 < V < d.loa * d.beam * d.depth
    assert c[1] == pytest.approx(0, abs=1e-9)


def test_box_like_section_volume():
    # hard chine, no deadrise, no rocker, parallel sides: volume -> L*B*D as the ends fill out
    d = Design(deadrise_deg=0, bilge_frac=0, rocker_fwd=0, rocker_aft=0,
               transom_frac=1.0, entry_half_angle_deg=89, xmax_frac=0.5)
    V, _ = volume_props(demihull_mesh(d, 201).triangles)
    full = d.loa * d.beam * d.depth
    assert 0.9 * full < V < full


def test_planform_and_keel():
    d = Design()
    assert halfbeam([d.xmax_frac], d)[0] == pytest.approx(d.beam / 2)
    assert halfbeam([0.0], d)[0] == pytest.approx(d.transom_frac * d.beam / 2)
    assert halfbeam([1.0], d)[0] < 0.01
    assert keel_z([d.xmax_frac], d)[0] == 0
    assert keel_z([1.0], d)[0] == pytest.approx(d.rocker_fwd)
    assert keel_z([0.0], d)[0] == pytest.approx(d.rocker_aft)
    # side is vertical: full half-beam at mid-depth amidships
    x = d.xmax_frac * d.loa
    assert halfbreadth_at([x], [d.depth * 0.8], d)[0] == pytest.approx(d.beam / 2)
    assert np.isnan(halfbreadth_at([x], [-0.01], d)[0])
    assert np.isnan(halfbreadth_at([-0.1], [0.1], d)[0])


def test_catamaran_spacing():
    d = Design(spacing_inner=0.30, beam=0.2)
    v = catamaran_mesh(d).vertices
    port = v[v[:, 1] > 0]
    assert port[:, 1].min() == pytest.approx(0.15, abs=1e-9)
    assert v[:, 1].max() - v[:, 1].min() == pytest.approx(0.30 + 0.4)


def test_print_pieces():
    d = Design(loa=1.5, beam=0.2, depth=0.24)
    assert print_pieces(d, (0.256, 0.256, 0.256)) == 6
    assert print_pieces(d, (0.35, 0.35, 0.35)) == 5


def test_evaluate_default_is_consistent():
    cfg = Config()
    row = evaluate(Design(), cfg)
    assert set(CONSTRAINTS) <= set(row)
    assert row["mass_total"] == pytest.approx(sum(v for k, v in row.items()
                                                  if k.startswith("mass_") and k != "mass_total"))
    # floats: displacement equals mass, near-level trim from payload placement
    assert abs(row["trim_deg"]) < 0.5
    assert 0 < row["draft_hull"] < Design().depth
    assert row["gm_t"] > 0
    assert 20 < row["mass_total"] < 40
    # thruster cover hits the requested value at the design waterline
    assert row["thr_cover_min"] == pytest.approx(Design().front_thr_cover, abs=2e-3)
    # skid bottom is below every thruster by at least the margin
    assert row["skid_bottom_z"] <= min(row["front_thr_z"], row["rear_thr_z"]) - 0.05
    # slenderness definition
    disp = row["mass_total"] / (2 * cfg.rho)
    assert row["slenderness"] == pytest.approx(row["lwl"] / disp ** (1 / 3))


def test_heavier_payload_sinks_deeper():
    a = evaluate(Design(), Config())
    b = evaluate(Design(), Config(payload_kg=25.0))
    assert b["draft_hull"] > a["draft_hull"]
    assert b["freeboard"] < a["freeboard"]


def test_random_designs_do_not_crash():
    rows = [evaluate(d) for d in lhs(25, np.random.default_rng(1))]
    for r in rows:
        assert all(math.isfinite(r[k]) for k in CONSTRAINTS)
        assert r["feasible"] == all(r[k] <= 0 for k in CONSTRAINTS)


def test_nondominated_sort():
    F = np.array([[1, 1], [2, 2], [1, 2], [0, 3], [3, 0]], float)
    fronts = nondominated_fronts(F)
    assert sorted(fronts[0]) == [0, 3, 4]
    assert sorted(fronts[1]) == [2]
    assert sorted(fronts[2]) == [1]


def test_select_candidates_only_feasible():
    rows = [evaluate(d) for d in lhs(30, np.random.default_rng(2))]
    pareto, cands = select_candidates(rows, 10)
    assert all(r["feasible"] for r in cands)
    assert len(cands) <= 10


def test_bounds_match_design_fields():
    assert set(BOUNDS) == set(Design().to_dict())
    for k, (lo, hi) in BOUNDS.items():
        assert lo < hi
        assert lo <= getattr(Design(), k) <= hi, k


def test_cli_run_lhs(tmp_path):
    main(["run", "--sampler", "lhs", "--trials", "20", "--select", "5", "--out", str(tmp_path)])
    for f in ["all.csv", "feasible.csv", "pareto.csv", "stage2_candidates.csv", "summary.md", "config.json"]:
        assert (tmp_path / f).exists()
    cfg = json.loads((tmp_path / "config.json").read_text())
    assert Config.from_dict(cfg) == Config()


def test_stl_roundtrip_size(tmp_path):
    tri = demihull_mesh(Design()).triangles
    p = tmp_path / "h.stl"
    write_stl(str(p), tri)
    assert p.stat().st_size == 84 + 50 * len(tri)


def test_thruster_layout_front_ahead_rear_between():
    d = Design()
    cfg = Config()
    row, det = evaluate(d, cfg, detail=True)
    th = {t.name: t for t in det.thrusters}
    yc = d.spacing_inner / 2 + d.beam / 2
    # front: thruster inside a pod on the hull centreline, pod tail ahead of the bow by the gap
    pod = det.pods[0]
    assert pod.y == pytest.approx(yc) and th["front_port"].centre[1] == pytest.approx(yc)
    assert pod.x0 == pytest.approx(d.loa + d.front_thr_gap)
    # thruster envelope sits inside the pod (with the inner gap)
    pts = th["front_port"].envelope_points()
    assert pts[:, 0].min() >= pod.x0 + d.pod_wall - 1e-9
    assert pts[:, 0].max() <= pod.x1 - d.pod_nose_len - d.pod_wall + 1e-9
    assert np.abs(pts[:, 1] - pod.y).max() <= pod.half_w - d.pod_wall + 1e-9
    assert np.abs(pts[:, 2] - pod.z).max() <= pod.half_h - d.pod_wall + 1e-9
    # rear: whole envelope inside the gap between the hulls (local hull width, any depth)
    for name in ("rear_port", "rear_stbd"):
        pts = th[name].envelope_points()
        inner = yc - halfbeam(pts[:, 0] / d.loa, d)
        assert (np.abs(pts[:, 1]) <= inner - cfg.clearance_min).all()
        assert th[name].radius == cfg.thr_radius   # no shroud at the rear
    # overall length runs to the pod nose; pod prints in one piece on the 250 mm bed
    assert row["loa_overall"] == pytest.approx(pod.x1)
    assert row["pod_print_ratio"] <= 1.0
    assert row["impact_energy"] == pytest.approx(0.5 * row["mass_total"] * cfg.impact_speed ** 2)


def test_pod_too_long_for_bed_is_flagged():
    row = evaluate(Design(pod_nose_len=0.12, pod_tail_len=0.06))
    assert row["c_pod_print"] > 0


def test_export_meshes(tmp_path):
    from hullopt.export import build, export, foil_mesh, pod_mesh, thruster_mesh
    d = Design()
    row, det = evaluate(d, Config(), detail=True)
    pod = det.pods[0]
    # closed skin volume is close to the elliptic-cylinder bound, and below it
    V, _ = volume_props(pod_mesh(pod, windows=False).triangles)
    box = np.pi * pod.half_w * pod.half_h * pod.length
    assert 0.6 * box < V < box
    # windows remove faces only from the sides
    assert len(pod_mesh(pod).faces) < len(pod_mesh(pod, windows=False).faces)
    for t in det.thrusters:
        V, c = volume_props(thruster_mesh(t).triangles)
        assert V == pytest.approx(np.pi * t.radius**2 * t.length, rel=0.02)
        assert c == pytest.approx(t.centre, abs=1e-6)
    for s in det.struts:
        V, _ = volume_props(foil_mesh(s).triangles)
        assert V == pytest.approx(0.685 * s.tc * s.chord**2 * s.length, rel=0.05)
    paths = export(d, Config(), str(tmp_path / "x"))
    assert all((tmp_path / p.split("/")[-1]).stat().st_size > 84 for p in paths)
