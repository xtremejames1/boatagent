"""Thrusters, shrouds, struts, skids, and the geometric checks between them and the hulls.

Thrusters are cylinders with horizontal axes, each on a 3D-printed mount.

Front pair: one ahead of each hull's bow, on that hull's centreline, inside an
impact shroud. Axis at `front_thr_angle_deg` from +x, pointing forward-outboard,
so the jet runs aft-inboard toward the other hull or forward-outboard, depending
on thrust sign. The mount is an arm forward from the stem plus a faired strut
down to the shroud.

Rear pair: unshrouded, axis fore-aft, wholly in the gap between the hulls
(never under a hull). The mount is a faired strut down from the cross structure.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .geometry import halfbeam, halfbreadth_at, hull_centre_y, keel_z
from .params import Config, Design


@dataclass
class Thruster:
    name: str
    centre: np.ndarray   # body frame
    heading: float       # axis angle from +x in plan (rad), positive toward port
    radius: float        # envelope radius (thruster + shroud if any)
    length: float        # envelope length

    @property
    def axis(self) -> np.ndarray:
        return np.array([np.cos(self.heading), np.sin(self.heading), 0.0])

    def envelope_points(self, n_ax: int = 5, n_circ: int = 16) -> np.ndarray:
        u = self.axis
        n = np.array([-u[1], u[0], 0.0])
        e = np.array([0.0, 0.0, 1.0])
        s = np.linspace(-self.length / 2, self.length / 2, n_ax)[:, None, None]
        a = np.linspace(0, 2 * np.pi, n_circ, endpoint=False)[None, :, None]
        ring = self.radius * (np.cos(a) * n + np.sin(a) * e)
        return (self.centre + s * u + ring).reshape(-1, 3)

    @property
    def lowest_z(self) -> float:
        return float(self.centre[2] - self.radius)

    @property
    def top_z(self) -> float:
        return float(self.centre[2] + self.radius)


def hull_gap(points: np.ndarray, d: Design) -> float:
    """Smallest lateral gap from any point to either demihull side (negative = inside)."""
    points = np.atleast_2d(points)
    yc = hull_centre_y(d)
    b = halfbreadth_at(points[:, 0], points[:, 2], d)
    best = np.inf
    for yh in (yc, -yc):
        g = np.abs(points[:, 1] - yh) - b
        g = np.where(np.isnan(g), np.inf, g)
        best = min(best, float(g.min()))
    return best


def place_front(d: Design, z: float, heading: float, radius: float, length: float) -> Thruster:
    """Port front thruster: ahead of the port bow, on the hull centreline, with
    `front_thr_gap` between the stem and the nearest point of the shroud."""
    half_x = 0.5 * length * abs(np.cos(heading)) + radius * abs(np.sin(heading))
    x = d.loa + d.front_thr_gap + half_x
    return Thruster("front_port", np.array([x, hull_centre_y(d), z]), heading, radius, length)


def place_rear(d: Design, x: float, z: float, radius: float, length: float, clearance: float) -> Thruster:
    """Port rear thruster: axis fore-aft, as far outboard as the gap between the
    hulls allows (`clearance` from the widest hull side alongside it, at any depth)."""
    xs = np.clip(np.linspace(x - length / 2, x + length / 2, 9) / d.loa, 0, 1)
    inner_face = hull_centre_y(d) - float(halfbeam(xs, d).max())
    y = inner_face - clearance - 1e-3 - radius
    return Thruster("rear_port", np.array([x, y, z]), 0.0, radius, length)


def pair_gap(t: Thruster) -> float:
    """Gap between a port thruster envelope and its starboard mirror image."""
    return 2.0 * float(np.abs(t.envelope_points()[:, 1]).min()) if t.centre[1] > 0 else np.nan


def jet_free_length(t: Thruster, direction: np.ndarray, d: Design, jet_radius: float,
                    max_len: float = 2.0, step: float = 0.005, n_circ: int = 12) -> float:
    """Distance from the thruster face along `direction` before any point of the
    jet's circular edge (radius `jet_radius`, no spreading) enters a hull."""
    u = direction / np.linalg.norm(direction)
    n = np.cross(u, [0.0, 0.0, 1.0])
    n = n / np.linalg.norm(n) if np.linalg.norm(n) > 0 else np.array([0.0, 1.0, 0.0])
    e = np.cross(n, u)
    a = np.linspace(0, 2 * np.pi, n_circ, endpoint=False)
    ring = jet_radius * (np.cos(a)[:, None] * n + np.sin(a)[:, None] * e)
    start = t.centre + u * t.length / 2
    s = np.arange(0.0, max_len, step)
    pts = (start + s[:, None, None] * u + ring[None]).reshape(-1, 3)
    yc = hull_centre_y(d)
    b = halfbreadth_at(pts[:, 0], pts[:, 2], d)
    hit = np.zeros(len(pts), bool)
    for yh in (yc, -yc):
        hit |= np.nan_to_num(np.abs(pts[:, 1] - yh) - b, nan=np.inf) < 0
    hit = hit.reshape(len(s), n_circ).any(1)
    idx = np.nonzero(hit)[0]
    return float(s[idx[0]]) if len(idx) else max_len


def point_to_ray(p: np.ndarray, origin: np.ndarray, u: np.ndarray) -> float:
    """Distance from p to the half-line origin + s u (s >= 0)."""
    u = u / np.linalg.norm(u)
    s = max(0.0, float(np.dot(p - origin, u)))
    return float(np.linalg.norm(p - (origin + s * u)))


@dataclass
class Strut:
    """Vertical faired strut from deck level down to a thruster or shroud top."""

    name: str
    x: float
    y: float
    z_top: float
    z_bot: float
    chord: float
    tc: float

    @property
    def length(self) -> float:
        return max(self.z_top - self.z_bot, 0.0)

    def submerged_length(self, z_wl: float) -> float:
        return max(0.0, min(self.z_top, z_wl) - self.z_bot)

    def mass(self, cfg: Config) -> float:
        # hollow printed foil: perimeter ~2.04 c (plus the arm/clamp, see mount_mass)
        return self.length * 2.04 * self.chord * cfg.strut_wall * cfg.fil_density * cfg.mass_calibration

    def volume_submerged(self, z_wl: float) -> float:
        return self.submerged_length(z_wl) * 0.685 * self.tc * self.chord ** 2

    def wetted_area(self, z_wl: float) -> float:
        return self.submerged_length(z_wl) * 2.04 * self.chord

    @property
    def cg(self) -> np.ndarray:
        return np.array([self.x, self.y, 0.5 * (self.z_top + self.z_bot)])


def mount_mass(strut: Strut, arm_length: float, cfg: Config) -> float:
    """Printed mount: faired strut + horizontal arm + clamp/bracket."""
    arm = arm_length * cfg.mount_arm_kg_per_m * cfg.mass_calibration
    return strut.mass(cfg) + arm + cfg.mount_fixed_kg


def shroud_ring(cfg: Config, d: Design):
    """Mass, enclosed material volume and wetted area of one front shroud ring."""
    r = cfg.thr_radius + cfg.shroud_gap + 0.5 * d.shroud_guard
    L = cfg.thr_length + 0.01
    vol = 2 * np.pi * r * L * d.shroud_guard
    mass = vol * cfg.fil_density * cfg.mass_calibration
    wet = 2 * (2 * np.pi * r * L)
    return mass, vol, wet


@dataclass
class Skid:
    x0: float
    x1: float
    z_bottom: float
    width: float
    heights: np.ndarray
    xs: np.ndarray

    @property
    def length(self) -> float:
        return self.x1 - self.x0

    @property
    def max_height(self) -> float:
        return float(self.heights.max())

    @property
    def volume(self) -> float:
        # rounded rail: ~0.9 of the bounding rectangle
        return float(0.9 * self.width * np.trapezoid(self.heights, self.xs))

    @property
    def centroid(self) -> np.ndarray:
        a = self.heights
        A = np.trapezoid(a, self.xs)
        xc = np.trapezoid(a * self.xs, self.xs) / A
        zc = np.trapezoid(a * (self.z_bottom + a / 2), self.xs) / A
        return np.array([xc, 0.0, zc])

    @property
    def wetted_area(self) -> float:
        return float(np.trapezoid(2 * self.heights + self.width, self.xs))


def make_skid(d: Design, cfg: Config, lowest_other_z: float) -> Skid:
    """Flat-bottomed rail: bottom sits `skid_margin` below the lowest thruster/shroud
    point, and at least `skid_min_height` below the keel."""
    zb = min(lowest_other_z - d.skid_margin, -cfg.skid_min_height)
    half = 0.5 * d.skid_len_frac * d.loa
    xc = d.skid_x_frac * d.loa
    x0, x1 = max(0.0, xc - half), min(d.loa, xc + half)
    xs = np.linspace(x0, x1, 41)
    h = keel_z(xs / d.loa, d) - zb
    return Skid(x0, x1, zb, cfg.skid_width, h, xs)
