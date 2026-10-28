import numpy as np
import pytest

from hullopt.frame import Frame
from hullopt.strength import Material, Section


def cantilever(L=0.2):
    f = Frame()
    a, b = f.node([0, 0, 0]), f.node([L, 0, 0])
    f.fixed.add(a)
    sec = Section(0.012, 0.024, 0.0016)
    f.add("c", a, b, sec, [0, 0, 1.0])          # depth vertical
    return f, a, b, sec


def test_cantilever_tip_load_matches_beam_theory():
    m = Material()
    L, P = 0.2, 50.0
    f, a, b, sec = cantilever(L)
    u, res = f.solve({b: (np.array([0, 0, -P]), np.zeros(3))}, m)
    assert u[6 * b + 2] == pytest.approx(-P * L**3 / (3 * m.E * sec.I_h), rel=1e-6)
    assert res[0].m_depth == pytest.approx(P * L, rel=1e-6)
    # sideways: weak axis
    u, res = f.solve({b: (np.array([0, P, 0]), np.zeros(3))}, m)
    assert u[6 * b + 1] == pytest.approx(P * L**3 / (3 * m.E * sec.I_b), rel=1e-6)
    assert res[0].m_thick == pytest.approx(P * L, rel=1e-6)
    # axial and torsion
    u, res = f.solve({b: (np.array([100.0, 0, 0]), np.array([2.0, 0, 0]))}, m)
    assert res[0].axial == pytest.approx(100.0)
    assert res[0].torque == pytest.approx(2.0)
    assert u[6 * b] == pytest.approx(100.0 * L / (m.E * sec.area))


def test_reactions_balance_loads():
    m = Material()
    f, a, b, sec = cantilever(0.2)
    F = np.array([10.0, -20.0, 30.0])
    u, _ = f.solve({b: (F, np.zeros(3))}, m)
    R, M = f.reactions(u, m)[a]
    assert R == pytest.approx(-F, abs=1e-6)
    assert M == pytest.approx(-np.cross(f.nodes[b] - f.nodes[a], F), abs=1e-6)


def test_two_member_truss_axial_split():
    # symmetric A-frame: load at the apex shares equally, members mostly axial
    m = Material()
    f = Frame()
    a, b, c = f.node([0, -0.1, 0]), f.node([0, 0.1, 0]), f.node([0.2, 0, 0])
    f.fixed |= {a, b}
    sec = Section(0.012, 0.024, 0.0016)
    f.add("l", a, c, sec, [0, 0, 1]), f.add("r", b, c, sec, [0, 0, 1])
    u, res = f.solve({c: (np.array([-100.0, 0, 0]), np.zeros(3))}, m)
    assert res[0].axial == pytest.approx(res[1].axial, rel=1e-6)
    L = np.hypot(0.2, 0.1)
    assert res[0].axial == pytest.approx(-50.0 * L / 0.2, rel=0.02)   # compression
