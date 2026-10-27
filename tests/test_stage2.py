import math

import numpy as np
import pytest

from hullopt.michell import _t_grid, michell, spectrum, wigley_offsets, wigley_wetted_area
from hullopt.params import Config, Design
from hullopt.stage2 import Stage2Config, evaluate2, form_factor, main, rank

G = 9.81


def strut(L=1.0, B=0.1, T=0.06, nx=401, nz=8):
    """Wall-sided strut with parabolic waterlines: f = B/2 (1 - (2x/L)^2), depth T."""
    x = np.linspace(-L / 2, L / 2, nx)
    de = np.linspace(0, T, nz + 1)
    F = np.tile((0.5 * B * (1 - (2 * x / L) ** 2))[:, None], (1, nz))
    return x, de, F


def strut_spectrum_exact(lam, U, L, B, T):
    k0 = G / U**2
    a, k = k0 * lam, k0 * lam**2
    # int_{-L/2}^{L/2} f_x e^{iax} dx with f_x = -4 B x / L^2
    ix = 2j * (np.sin(a * L / 2) / a**2 - (L / 2) * np.cos(a * L / 2) / a)
    return (-4 * B / L**2) * ix * (1 - np.exp(-k * T)) / k


@pytest.mark.parametrize("Fr", [0.25, 0.4, 0.6])
def test_spectrum_matches_analytic_strut(Fr):
    L, B, T = 1.0, 0.1, 0.06
    U = Fr * math.sqrt(G * L)
    x, de, F = strut(L, B, T)
    t, lam, S = spectrum(x, de, F, U)
    Se = strut_spectrum_exact(lam, U, L, B, T)
    # piecewise-linear f: error O(dx^2); compare energy-weighted
    err = np.sqrt(np.trapezoid(np.abs(S - Se) ** 2, t) / np.trapezoid(np.abs(Se) ** 2, t))
    assert err < 1e-3


def test_resistance_matches_analytic_strut():
    L, B, T, U = 1.0, 0.1, 0.06, 0.4 * math.sqrt(G)
    x, de, F = strut(L, B, T)
    R = michell(x, de, F, U)
    t = _t_grid(G / U**2, L, 0.0, 10.0, 8)
    lam = 1 + t**2
    Se = strut_spectrum_exact(lam, U, L, B, T)
    Re = 4 * 1000 * G**2 / (math.pi * U**2) * np.trapezoid(np.abs(Se) ** 2 * lam**2 * 2 / np.sqrt(lam + 1), t)
    assert R == pytest.approx(Re, rel=2e-3)


def test_thin_ship_scales_with_beam_squared():
    x, de, F = wigley_offsets()
    U = 0.4 * math.sqrt(G)
    assert michell(x, de, 2 * F, U) == pytest.approx(4 * michell(x, de, F, U), rel=1e-9)


def test_zero_spacing_catamaran_is_double_beam_hull():
    x, de, F = wigley_offsets()
    U = 0.45 * math.sqrt(G)
    assert michell(x, de, F, U, spacing=0.0) == pytest.approx(michell(x, de, 2 * F, U), rel=1e-9)


def test_froude_scaling():
    # geometrically similar hull at the same Fr: R scales with length^3
    x, de, F = wigley_offsets(L=1.0, B=0.1, T=0.0625)
    x2, de2, F2 = wigley_offsets(L=2.0, B=0.2, T=0.125)
    Fr = 0.35
    R1 = michell(x, de, F, Fr * math.sqrt(G))
    R2 = michell(x2, de2, F2, Fr * math.sqrt(2 * G))
    assert R2 == pytest.approx(8 * R1, rel=1e-3)


def test_wide_catamaran_tends_to_two_hulls():
    x, de, F = wigley_offsets()
    U = 0.5 * math.sqrt(G)
    mono = michell(x, de, F, U)
    ratios = [michell(x, de, F, U, spacing=s) / (2 * mono) for s in (2.0, 4.0)]
    assert abs(ratios[1] - 1) < abs(ratios[0] - 1) + 0.02
    assert abs(ratios[1] - 1) < 0.1


def test_lam_max_converged():
    x, de, F = wigley_offsets()
    for Fr in (0.2, 0.5):
        U = Fr * math.sqrt(G)
        assert michell(x, de, F, U, lam_max=10) == pytest.approx(michell(x, de, F, U, lam_max=25), rel=2e-3)


def test_wigley_reference_values():
    """Regression values for the Wigley hull (L/B=10, L/T=16). Wetted area S/L^2 = 0.1488.
    Michell overpredicts the measured hump (Cw ~ 2e-3 at Fr 0.5) by about 2x, as is known."""
    assert wigley_wetted_area() == pytest.approx(0.1488, abs=2e-4)
    x, de, F = wigley_offsets(nx=101, nz=30)
    S = wigley_wetted_area()
    cw = {}
    for Fr in (0.3, 0.35, 0.5):
        U = Fr * math.sqrt(G)
        cw[Fr] = michell(x, de, F, U) / (0.5 * 1000 * U * U * S)
    assert cw[0.3] == pytest.approx(2.14e-3, rel=0.03)
    assert cw[0.35] < cw[0.3]           # hollow between humps
    assert cw[0.5] == pytest.approx(4.52e-3, rel=0.03)


def test_form_factor_couser():
    c2 = Stage2Config()
    assert form_factor(c2, 6.0) == pytest.approx(3.03 * 6.0 ** -0.4)
    assert form_factor(Stage2Config(form_factor=1.2), 6.0) == 1.2


def test_evaluate2_knee_like_design():
    d = Design(loa=1.30, beam=0.155, depth=0.197, spacing_inner=0.395,
               rear_thr_x_frac=0.17)
    row = evaluate2(d, Config(), Stage2Config())
    assert row["stage2_ok"]
    # components positive and add up
    for U in ("0.5", "1", "1.5"):
        parts = sum(row[f"{k}_{U}"] for k in ("R_wave", "R_visc", "R_app", "R_transom"))
        assert row[f"R_total_{U}"] == pytest.approx(parts)
        assert row[f"R_wave_{U}"] >= 0 and row[f"R_visc_{U}"] > 0
    assert row["R_total_1.5"] > row["R_total_0.5"]
    assert row["R_T_weighted"] == pytest.approx(row["R_hull_weighted"] + row["R_app_weighted"])
    assert 1.2 < row["form_factor"] < 1.8
    assert row["v_sway_max_s2"] > 0 and row["yaw_rate_max_dps"] > 0
    assert row["thrust_margin_top"] > 1


def test_rank_and_cli(tmp_path):
    cands = tmp_path / "c.csv"
    designs = [Design(loa=1.30, beam=0.155, depth=0.197, spacing_inner=0.395,
                      rear_thr_x_frac=0.17),
               Design(loa=1.35, beam=0.15, depth=0.20, spacing_inner=0.36,
                      rear_thr_x_frac=0.17)]
    keys = list(designs[0].to_dict())
    cands.write_text(",".join(keys) + "\n" + "\n".join(",".join(str(getattr(d, k)) for k in keys)
                                                      for d in designs) + "\n")
    out = tmp_path / "s2"
    main(["run", "--candidates", str(cands), "--out", str(out), "--finalists", "1"])
    assert (out / "finalists.csv").exists() and (out / "curves.csv").exists()
    rows = [evaluate2(d, Config(), Stage2Config()) for d in designs]
    pareto, finals = rank(rows, 1)
    assert len(finals) == 1
