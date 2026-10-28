import pytest

cq = pytest.importorskip("cadquery")

from hullopt.params import Config, Design  # noqa: E402
from hullopt.parts import (FrontSpec, PrintSettings, fits_bed, front_checks, front_fairing,  # noqa: E402
                           front_frame, front_geometry, front_strength, keep_out_cones, overhang_report,
                           size_members)

FEASIBLE = dict(loa=1.15, beam=0.13, depth=0.24, spacing_inner=0.395, rear_thr_x_frac=0.17,
                front_thr_cover=0.051, rear_thr_cover=0.051, skid_margin=0.0101, front_thr_gap=-0.03)


@pytest.fixture(scope="module")
def built():
    d, spec, ps = Design(**FEASIBLE), FrontSpec(), PrintSettings()
    cfg, _ = size_members(d, Config(), spec, ps)
    row, det, fp, pod, m = front_geometry(d, cfg)
    inlet, outlet = keep_out_cones(cfg, fp)
    frame, m3, lugs = front_frame(d, cfg, spec, fp, pod, m)
    pod.thr = fp
    fairing = front_fairing(d, cfg, spec, pod, lugs, max(d.pod_wall, ps.wall), inlet)
    checks = front_checks(cfg, spec, fp, frame, fairing, inlet, outlet, m3)
    return d, cfg, spec, ps, fp, pod, m, frame, fairing, checks


def test_parts_valid_and_fit_bed_nose_up(built):
    d, cfg, spec, ps, fp, pod, m, frame, fairing, checks = built
    for p in (frame, fairing):
        assert p.val().isValid()
        ok, dims = fits_bed(p, ps, up_axis=0)
        assert ok, dims


def test_printable_without_supports(built):
    d, cfg, spec, ps, fp, pod, m, frame, fairing, checks = built
    # only small ledges at joints remain (slicer bridges them)
    assert overhang_report(frame, cfg.max_overhang_deg)["overhang_mm2"] < 500
    assert overhang_report(fairing, cfg.max_overhang_deg)["overhang_mm2"] < 100


def test_clearances_cones_removal_access(built):
    *_, checks = built
    for k in ("frame_vs_thruster_mm3", "fairing_vs_thruster_mm3", "frame_vs_fairing_mm3",
              "frame_in_inlet_cone_mm3", "frame_in_outlet_cone_mm3", "fairing_in_outlet_cone_mm3",
              "removal_path_max_mm3"):
        assert checks[k] < 1.0, (k, checks[k])
    assert checks["m3_screws_blocked"] == 0


def test_frame_holds_thrust_and_head_on(built):
    d, cfg, spec, ps, fp, pod, m, *_ = built
    cases, breakaway, _ = front_strength(d, cfg, spec, fp, pod, m, ps)
    for c in cases:
        if not c["primary"]:
            continue
        assert max(x["utilisation"] for x in c["members"]) <= 1.0, c["case"]
        assert max(x["buckling"] for x in c["members"]) <= 0.5, c["case"]
        assert c["screw_shear_max_N"] <= spec.insert_shear, c["case"]
    assert 300 < breakaway < 700
