import numpy as np
import pytest

from hullopt.geometry import Mesh, loft
from hullopt.hydrostatics import Pose, equilibrium, gz_curve, hydro, volume_props

RHO = 1000.0


def box(L, B, D, y0=0.0, n=11):
    """Box hull from x=0..L, |y-y0|<=B/2, z=0..D, built with the same lofter as the hulls."""
    X = np.linspace(0, L, n)
    ys = np.array([0, B / 2, B / 2, 0, -B / 2, -B / 2]) + y0
    zs = np.array([0, 0, D, D, D, 0])
    return loft(X, np.tile(ys, (n, 1)), np.tile(zs, (n, 1)))


def test_closed_volume_and_centroid():
    m = box(2.0, 0.5, 0.3)
    V, c = volume_props(m.triangles)
    assert V == pytest.approx(2.0 * 0.5 * 0.3)
    assert c == pytest.approx([1.0, 0.0, 0.15])


def test_box_hydrostatics_analytic():
    L, B, D, T = 2.0, 0.5, 0.3, 0.1
    m = box(L, B, D)
    hy = hydro(m.triangles, Pose(heave=T))
    assert hy.volume == pytest.approx(L * B * T)
    assert hy.cob == pytest.approx([L / 2, 0, -T / 2])
    assert hy.wp.area == pytest.approx(L * B)
    assert hy.wp.i_t == pytest.approx(L * B**3 / 12)
    assert hy.wp.i_l == pytest.approx(B * L**3 / 12)
    assert hy.wetted_area == pytest.approx(L * B + 2 * L * T + 2 * B * T)
    assert hy.lateral_area == pytest.approx(L * T)


def test_twin_box_parallel_axis():
    L, B, D, T, s = 1.5, 0.2, 0.3, 0.1, 0.6
    m = Mesh.merge(box(L, B, D, y0=s / 2), box(L, B, D, y0=-s / 2))
    hy = hydro(m.triangles, Pose(heave=T))
    expect = 2 * (L * B**3 / 12 + L * B * (s / 2) ** 2)
    assert hy.wp.i_t == pytest.approx(expect)


def test_equilibrium_level_and_trimmed():
    L, B, D = 2.0, 0.5, 0.3
    m = box(L, B, D)
    mass = 80.0
    pose, hy = equilibrium(m.triangles, mass, [L / 2, 0, 0.2], RHO, pivot=[L / 2, 0, 0])
    assert pose.trim == pytest.approx(0, abs=1e-7)
    assert pose.heave == pytest.approx(mass / (RHO * L * B))
    # Shift G forward: boat trims bow down, small-angle tan(theta) = dx / GM_L
    dx = 0.02
    pose2, hy2 = equilibrium(m.triangles, mass, [L / 2 + dx, 0, 0.2], RHO, pivot=[L / 2, 0, 0])
    T = mass / (RHO * L * B)
    gml = T / 2 + L**2 / (12 * T) - (0.2 - T)  # KB + BM_L - KG (KG from keel = 0.2)
    assert pose2.trim > 0
    assert np.tan(pose2.trim) == pytest.approx(dx / gml, rel=0.05)


def test_gz_small_angle_matches_gm():
    L, B, D = 2.0, 0.5, 0.3
    m = box(L, B, D)
    mass, kg = 80.0, 0.12
    pose, _ = equilibrium(m.triangles, mass, [L / 2, 0, kg], RHO, pivot=[L / 2, 0, 0])
    T = pose.heave
    gm = T / 2 + B**2 / (12 * T) - kg
    gz = gz_curve(m.triangles, mass, [L / 2, 0, kg], RHO, pose, [2.0])
    assert gz[0] == pytest.approx(gm * np.sin(np.radians(2.0)), rel=0.02)
