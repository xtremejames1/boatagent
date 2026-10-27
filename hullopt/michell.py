"""Michell thin-ship wave resistance, monohull or catamaran.

    R = 4 rho g^2 / (pi U^2) * int_1^inf |S(lam)|^2 lam^2 / sqrt(lam^2 - 1) * I(lam) dlam
    S(lam) = int int f_x(x, z) exp(k0 lam^2 z + i k0 lam x) dx dz,   k0 = g / U^2

f(x, z) is the demihull half-breadth on the centre plane (z <= 0 below the
free surface). For two identical demihulls with centreline spacing s and no
stagger, the free waves add with the phase k_y s, k_y = k0 lam sqrt(lam^2-1):
    I(lam) = 2 [1 + cos(k0 lam sqrt(lam^2 - 1) s)]   (I = 1 for one hull)

Panels: f is linear in x between stations and constant in z within a depth
cell, so S is integrated exactly over each panel (no aliasing at high lam).
The singularity at lam = 1 is removed with lam = 1 + t^2.

The hull is not closed at a transom: the integral runs over the hull only,
which treats the transom flow as a constant-section wake. Add a transom
pressure term (see stage2) for wet transoms.
"""
from __future__ import annotations

import math

import numpy as np


def _t_grid(k0: float, length: float, spacing: float, lam_max: float, pts_per_period: int):
    """Uniform grid in t that resolves the fastest phase in the integrand."""
    t_max = math.sqrt(lam_max - 1.0)
    # d(phase)/dlam <= k0 (L + 2 lam s); dlam/dt = 2t
    rate = k0 * (length + 2.0 * lam_max * spacing) * 2.0 * t_max
    n = int(math.ceil(t_max * rate * pts_per_period / (2 * math.pi))) + 400
    return np.linspace(0.0, t_max, n + 1)


def spectrum(x, depth_edges, F, U, g=9.81, t=None, lam_max=10.0, spacing=0.0, pts_per_period=8):
    """Return (t, lam, S) for offsets F (nx, nz) at stations x and depth cells depth_edges (nz+1)."""
    x = np.asarray(x, float)
    de = np.asarray(depth_edges, float)
    F = np.asarray(F, float)
    k0 = g / U**2
    L = x[-1] - x[0]
    if t is None:
        t = _t_grid(k0, L, spacing, lam_max, pts_per_period)
    lam = 1.0 + t**2
    fx = np.diff(F, axis=0) / np.diff(x)[:, None]          # (nx-1, nz)
    xm = x - 0.5 * (x[0] + x[-1])
    a = k0 * lam
    k = k0 * lam**2
    z_top, z_bot = -de[:-1], -de[1:]
    S = np.zeros(len(t), complex)
    # chunk over lam to bound memory
    step = max(1, 2_000_000 // max(len(x), 1))
    for i0 in range(0, len(t), step):
        sl = slice(i0, i0 + step)
        E = np.exp(1j * a[sl, None] * xm[None, :])
        X = (E[:, 1:] - E[:, :-1]) / (1j * a[sl, None])     # exact int of exp(i a x) per x cell
        Z = (np.exp(k[sl, None] * z_top) - np.exp(k[sl, None] * z_bot)) / k[sl, None]
        S[sl] = ((X @ fx) * Z).sum(1)
    return t, lam, S


def michell(x, depth_edges, F, U, spacing=None, rho=1000.0, g=9.81, lam_max=10.0, pts_per_period=8) -> float:
    """Wave resistance (N). spacing=None: one hull. Otherwise two hulls, centreline spacing in m."""
    s = 0.0 if spacing is None else float(spacing)
    t, lam, S = spectrum(x, depth_edges, F, U, g=g, lam_max=lam_max, spacing=s,
                         pts_per_period=pts_per_period)
    w = np.abs(S) ** 2 * lam**2 * 2.0 / np.sqrt(lam + 1.0)   # lam^2/sqrt(lam^2-1) dlam = this dt
    if spacing is not None:
        w = w * 2.0 * (1.0 + np.cos(k0_of(U, g) * lam * np.sqrt(lam**2 - 1.0) * s))
    return float(4.0 * rho * g**2 / (math.pi * U**2) * np.trapezoid(w, t))


def k0_of(U: float, g: float = 9.81) -> float:
    return g / U**2


def wigley_offsets(L=1.0, B=0.1, T=0.0625, nx=81, nz=24):
    """Wigley parabolic hull: f = B/2 (1 - (2x/L)^2)(1 - (z/T)^2)."""
    x = np.linspace(-L / 2, L / 2, nx)
    de = np.linspace(0, T, nz + 1)
    dm = 0.5 * (de[1:] + de[:-1])
    F = 0.5 * B * (1 - (2 * x[:, None] / L) ** 2) * (1 - (dm[None, :] / T) ** 2)
    return x, de, F


def wigley_wetted_area(L=1.0, B=0.1, T=0.0625, n=400) -> float:
    """Wetted area of the Wigley hull by direct surface integration."""
    x = np.linspace(-L / 2, L / 2, n)
    z = np.linspace(0, T, n)
    X, Z = np.meshgrid(x, z, indexing="ij")
    fx = 0.5 * B * (-8 * X / L**2) * (1 - (Z / T) ** 2)
    fz = 0.5 * B * (1 - (2 * X / L) ** 2) * (-2 * Z / T**2)
    dA = np.sqrt(1 + fx**2 + fz**2)
    return float(2 * np.trapezoid(np.trapezoid(dA, z, axis=1), x))
