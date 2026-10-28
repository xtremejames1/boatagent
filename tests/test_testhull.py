import numpy as np
import pytest

from hullopt.params import Design
from hullopt.testhull import INCH, hydro_table, scaled


def test_scaled_keeps_shape_and_sets_gap():
    d = Design()
    s, k = scaled(d, 24 * INCH, 12 * INCH)
    assert s.loa == pytest.approx(0.6096)
    assert s.spacing_inner == pytest.approx(0.3048)
    for f in ("beam", "depth", "rocker_fwd", "rocker_aft"):
        assert getattr(s, f) == pytest.approx(getattr(d, f) * k)
    for f in ("deadrise_deg", "bilge_frac", "entry_half_angle_deg", "transom_frac", "xmax_frac"):
        assert getattr(s, f) == getattr(d, f)


def test_draft_grows_with_load():
    from hullopt.params import Config
    s, _ = scaled(Design(), 24 * INCH, 12 * INCH)
    rows = hydro_table(s, Config(), [1.0, 2.0, 3.0])
    drafts = [r["draft_mm"] for r in rows]
    assert all(a < b for a, b in zip(drafts, drafts[1:]))
