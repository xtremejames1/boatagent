import numpy as np
import pytest

cq = pytest.importorskip("cadquery")

from hullopt.params import Config, Design  # noqa: E402
from hullopt.parts import (FrontSpec, PrintSettings, fits_bed, front_bracket, front_geometry,  # noqa: E402
                           front_pod, front_pod_bosses, front_strength)

FEASIBLE = dict(loa=1.22, beam=0.13, depth=0.21, spacing_inner=0.395, rear_thr_x_frac=0.17,
                front_thr_cover=0.051, rear_thr_cover=0.051)


@pytest.fixture(scope="module")
def parts():
    d, cfg, spec = Design(**FEASIBLE), Config(), FrontSpec()
    row, det, fp, pod, members = front_geometry(d, cfg)
    br = front_bracket(d, cfg, spec, fp, pod, members)
    pd = front_pod(d, cfg, spec, pod, members, 0.0016).union(front_pod_bosses(cfg, spec, fp, pod, 0.0016))
    c, ax = fp.centre * 1000, fp.axis
    thr = cq.Workplane("XY").add(cq.Solid.makeCylinder(fp.radius * 1000, fp.length * 1000,
                                                       cq.Vector(*(c - ax * fp.length * 500)), cq.Vector(*ax)))
    return d, cfg, spec, fp, pod, members, br, pd, thr


def _vol(a, b):
    try:
        return a.intersect(b).val().Volume()
    except Exception:
        return 0.0


def test_parts_are_valid_and_fit_the_bed(parts):
    *_, br, pd, thr = parts
    for p in (br, pd):
        assert p.val().isValid()
        ok, dims = fits_bed(p, PrintSettings())
        assert ok, dims


def test_no_interference(parts):
    *_, br, pd, thr = parts
    assert _vol(pd, br) < 1.0          # mm^3
    assert _vol(pd, thr) < 1.0
    assert _vol(br, thr) < 1.0


def test_bracket_holds_thrust_and_head_on(parts):
    d, cfg, spec, fp, pod, members, *_ = parts
    cases, breakaway, _ = front_strength(d, cfg, spec, fp, pod, members, PrintSettings())
    for c in cases[:3]:                 # forward, reverse, head-on at breakaway
        assert max(m.utilisation for m in c.members) < 1.0, c.case
        assert max(m.buckling for m in c.members) < 1.0, c.case
        assert c.bolt_tension_max < spec.insert_pullout, c.case
    assert 300 < breakaway < 700
