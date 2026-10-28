import numpy as np
import pytest

from hullopt.layout import Strut
from hullopt.strength import Material, Section, analyse


def cantilever(length=0.2):
    return Strut("c", np.array([0.0, 0, 0]), np.array([length, 0, 0]), 0.024, 0.5, aligned=False)


def test_section_properties():
    s = Section(0.012, 0.024, 0.0016)
    assert s.area == pytest.approx(0.012 * 0.024 - 0.0088 * 0.0208)
    assert s.I_h == pytest.approx((0.012 * 0.024**3 - 0.0088 * 0.0208**3) / 12)
    assert s.I_b == pytest.approx((0.024 * 0.012**3 - 0.0208 * 0.0088**3) / 12)


def test_cantilever_axial_and_bending_match_closed_form():
    m = Material()
    L, sec = 0.2, Section(0.012, 0.024, 0.0016)
    st = cantilever(L)
    tip = st.p1.copy()
    # axial pull
    r = analyse([st], [sec], tip, np.array([100.0, 0, 0]), np.zeros(3), m)
    assert r.members[0].axial == pytest.approx(100.0)
    assert r.members[0].sigma == pytest.approx(100.0 / sec.area)
    # vertical tip load: upright member bends about its strong axis, root moment = F L
    r = analyse([st], [sec], tip, np.array([0, 0, 50.0]), np.zeros(3), m)
    assert abs(r.members[0].m_depth) == pytest.approx(50.0 * L, rel=1e-6)
    assert r.members[0].sigma == pytest.approx(50.0 * L * (sec.h / 2) / sec.I_h, rel=1e-6)
    # sideways tip load: weak axis
    r = analyse([st], [sec], tip, np.array([0, 50.0, 0]), np.zeros(3), m)
    assert abs(r.members[0].m_thick) == pytest.approx(50.0 * L, rel=1e-6)


def test_two_parallel_members_share_load():
    m = Material()
    sec = Section(0.012, 0.024, 0.0016)
    a = Strut("a", np.array([0.0, 0.0, 0.05]), np.array([0.2, 0.0, 0.05]), 0.024, 0.5, aligned=False)
    b = Strut("b", np.array([0.0, 0.0, -0.05]), np.array([0.2, 0.0, -0.05]), 0.024, 0.5, aligned=False)
    r = analyse([a, b], [sec, sec], np.array([0.2, 0, 0]), np.array([0, 0, 100.0]), np.zeros(3), m)
    # symmetric: equal and opposite axial forces (a couple) plus member bending
    na, nb = r.members[0].axial, r.members[1].axial
    assert na == pytest.approx(-nb, rel=1e-6)
    # equilibrium at the roots: couple + both bending moments = applied F * L
    couple = abs(na) * 0.10
    bending = abs(r.members[0].m_depth) + abs(r.members[1].m_depth)
    assert couple + bending == pytest.approx(100.0 * 0.2, rel=1e-6)
    assert couple > 0.2 * 100.0 * 0.2      # the couple carries a real share


def test_bolt_group_tension():
    m = Material()
    sec = Section(0.012, 0.024, 0.0016)
    st = Strut("s", np.array([0.0, 0, 0]), np.array([0.1, 0, 0]), 0.024, 0.5, aligned=False)
    bolt = {"origin": np.zeros(3), "normal": np.array([1.0, 0, 0]),
            "axes": (np.array([0, 1.0, 0]), np.array([0, 0, 1.0])),
            "bolts": [(sy * 0.02, sz * 0.03) for sy in (-1, 1) for sz in (-1, 1)]}
    # pure pull along +x: equal share
    r = analyse([st], [sec], st.p1, np.array([400.0, 0, 0]), np.zeros(3), m, bolt)
    assert r.bolt_tension_max == pytest.approx(100.0)
    # downward tip load: moment lifts the upper bolts. M = F L, tension = M z / sum z^2
    F, L = 100.0, 0.1
    r = analyse([st], [sec], st.p1, np.array([0, 0, -F]), np.zeros(3), m, bolt)
    assert r.bolt_tension_max == pytest.approx(F * L * 0.03 / (4 * 0.03**2), rel=1e-6)
