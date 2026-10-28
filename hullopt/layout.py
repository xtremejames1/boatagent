"""Thrusters, shrouds, struts, skids, and the geometric checks between them and the hulls.

Thrusters are cylinders with horizontal axes, each on a 3D-printed mount.

Front pair: one ahead of each hull's bow, on that hull's centreline, inside a
faired pod (see Pod) with a load-taking nose and side windows. Thruster axis at
`front_thr_angle_deg` from +x. With toe-in (default) the axis points
forward-inboard, so the jet runs aft-outboard past the outside of the hull, or
forward-inboard ahead of the boat, depending on thrust sign.
The mount is two diagonal members from the bow stem (upper: stem top to pod
top; lower: low on the stem to the pod tail), a triangle bracket that shares
the bow interface with the pod.

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


@dataclass
class Pod:
    """Faired, impact-absorbing pod around a front thruster, aligned with +x.

    Section: ellipse with half-width `half_w` and half-height `half_h` (outer).
    From aft to fore: tail fairing, the thruster bay (thruster at 45 deg inside,
    side windows let it draw and discharge water sideways), then a crushable nose.
    The pod is flooded: only its printed material gives buoyancy.
    """

    name: str
    x0: float            # tail tip
    x1: float            # nose tip
    y: float
    z: float             # section centre (= thruster axis height)
    half_w: float
    half_h: float
    nose_len: float      # crush zone ahead of the thruster bay
    wall: float
    window_area: float   # open area per side
    bay_side_area: float  # side area of the thruster bay (per side)
    bay_len: float = 0.0  # fore-aft length of the thruster bay (inside the walls)

    @property
    def length(self) -> float:
        return self.x1 - self.x0

    @property
    def centre(self) -> np.ndarray:
        return np.array([0.5 * (self.x0 + self.x1), self.y, self.z])

    @property
    def lowest_z(self) -> float:
        return self.z - self.half_h

    @property
    def top_z(self) -> float:
        return self.z + self.half_h

    @property
    def frontal_area(self) -> float:
        return float(np.pi * self.half_w * self.half_h)

    @property
    def d_eq(self) -> float:
        return float(2 * np.sqrt(self.half_w * self.half_h))

    @property
    def outer_area(self) -> float:
        a, b = self.half_w, self.half_h
        perim = np.pi * (3 * (a + b) - np.sqrt((3 * a + b) * (a + 3 * b)))  # Ramanujan
        return float(0.8 * perim * self.length)  # 0.8: nose and tail taper

    def material_volume(self, cfg: Config) -> float:
        shell = max(self.outer_area - 2 * self.window_area, 0.0) * self.wall
        crush = 0.5 * self.frontal_area * self.nose_len * cfg.pod_crush_fill
        return float(shell + crush)

    def mass(self, cfg: Config) -> float:
        return self.material_volume(cfg) * cfg.fil_density * cfg.mass_calibration

    def drag_area(self, cf: float, cfg: Config) -> float:
        """Cd x A (m^2): Hoerner streamlined-body drag on frontal area, plus window drag."""
        ld = self.length / self.d_eq
        cd_body = cf * (3 * ld + 4.5 * ld ** -0.5 + 21 * ld ** -2)
        return cd_body * self.frontal_area + 2 * self.window_area * cfg.pod_window_cd

    def envelope_points(self, n_ax: int = 9, n_circ: int = 16) -> np.ndarray:
        xs = np.linspace(self.x0, self.x1, n_ax)[:, None]
        a = np.linspace(0, 2 * np.pi, n_circ, endpoint=False)[None, :]
        y = self.y + self.half_w * np.cos(a) + 0 * xs
        z = self.z + self.half_h * np.sin(a) + 0 * xs
        return np.stack([xs + 0 * a, y, z], -1).reshape(-1, 3)

    def print_ratio(self, bed) -> float:
        """Largest pod dimension / matching bed dimension (<= 1 fits in one piece)."""
        dims = sorted([self.length, 2 * self.half_w, 2 * self.half_h])
        return float(max(dd / b for dd, b in zip(dims, sorted(bed))))


def place_front(d: Design, cfg: Config, z: float, heading: float) -> tuple[Thruster, Pod]:
    """Port front thruster (axis at height z) and its pod: ahead of the port bow on the
    hull centreline, `front_thr_gap` between the bow face and the pod's tail."""
    R, l, g, w = cfg.thr_radius, cfg.thr_length, cfg.pod_gap, d.pod_wall
    sn, cs = abs(np.sin(heading)), abs(np.cos(heading))
    xe = 0.5 * l * cs + R * sn + g          # half x-extent of the bay
    half_w = 0.5 * l * sn + R * cs + g + w   # thruster at 45 deg sets the width
    top = pod_top_extent(d, cfg)       # above the axis
    bot = R + g + w                    # below the axis
    half_h = 0.5 * (top + bot)
    x0 = d.loa + d.front_thr_gap
    xt = x0 + d.pod_tail_len + w + xe
    x1 = xt + xe + w + d.pod_nose_len
    yc = hull_centre_y(d)
    win = cfg.pod_window_ratio * np.pi * R * R
    zc = z + 0.5 * (top - bot)          # section centre (the mount plate sits above the thruster)
    pod = Pod("front_port", x0, x1, yc, zc, half_w, half_h, d.pod_nose_len, w, win,
              2 * xe * 2 * (R + g), 2 * xe)
    return Thruster("front_port", np.array([xt, yc, z]), heading, R, l), pod


def pod_top_extent(d: Design, cfg: Config) -> float:
    """Thruster axis to the pod's outer top: T200 base, mount plate, gap, wall."""
    return cfg.t200_base_offset + cfg.mount_plate_t + cfg.pod_gap + d.pod_wall


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
    """Straight printed member from p0 to p1 (body frame), hollow faired section:
    `chord` x `tc * chord` thick.

    aligned=True: chord along +x (a foil facing the flow; rear struts).
    aligned=False: chord in the member's vertical plane, thickness across y
    (front bow members: 12 mm faces the flow, 24+ mm upright).
    """

    name: str
    p0: np.ndarray
    p1: np.ndarray
    chord: float
    tc: float
    aligned: bool = True

    # convenience views
    @property
    def x(self) -> float:
        return float(self.p0[0])

    @property
    def y(self) -> float:
        return float(self.p0[1])

    @property
    def y1(self) -> float:
        return float(self.p1[1])

    @property
    def z_bot(self) -> float:
        return float(self.p0[2])

    @property
    def z_top(self) -> float:
        return float(self.p1[2])

    @property
    def length(self) -> float:
        return float(np.linalg.norm(self.p1 - self.p0))

    @property
    def thickness(self) -> float:
        return self.tc * self.chord

    def _sub(self, z_wl: float) -> tuple[float, float]:
        """Submerged share of the length, and the parameter t of its midpoint."""
        z0, z1 = self.p0[2], self.p1[2]
        if abs(z1 - z0) < 1e-12:
            return (1.0, 0.5) if z0 < z_wl else (0.0, 0.5)
        t = (z_wl - z0) / (z1 - z0)
        if z1 > z0:
            f = float(np.clip(t, 0, 1))
            return f, f / 2
        f = float(np.clip(1 - t, 0, 1))
        return f, 1 - f / 2

    def submerged_length(self, z_wl: float) -> float:
        return self.length * self._sub(z_wl)[0]

    def submerged_height(self, z_wl: float) -> float:
        """Vertical extent below the waterline."""
        return abs(self.p1[2] - self.p0[2]) * self._sub(z_wl)[0]

    def submerged_centre(self, z_wl: float) -> np.ndarray:
        return self.p0 + self._sub(z_wl)[1] * (self.p1 - self.p0)

    def side_area(self, z_wl: float) -> float:
        """Submerged area seen from the side (xz projection)."""
        d = self.p1 - self.p0
        if self.aligned:
            return self.chord * self.submerged_height(z_wl)
        return self.chord * self._sub(z_wl)[0] * float(np.hypot(d[0], d[2]))

    def frontal_area(self, z_wl: float) -> float:
        """Submerged area facing the flow for a non-aligned member (thickness x height)."""
        return 0.0 if self.aligned else self.thickness * self.submerged_height(z_wl)

    def mass(self, cfg: Config) -> float:
        # hollow printed section: perimeter ~2.04 c
        return self.length * 2.04 * self.chord * cfg.strut_wall * cfg.fil_density * cfg.mass_calibration

    def volume_submerged(self, z_wl: float) -> float:
        return self.submerged_length(z_wl) * 0.685 * self.tc * self.chord ** 2

    def wetted_area(self, z_wl: float) -> float:
        return self.submerged_length(z_wl) * 2.04 * self.chord

    @property
    def cg(self) -> np.ndarray:
        return 0.5 * (self.p0 + self.p1)


def rear_strut(d: Design, cfg: Config, t: Thruster) -> Strut:
    """Diagonal strut from the top of a rear thruster, outboard and up at
    `rear_strut_angle_deg` from vertical, ending where it meets its own hull
    (inner side or bottom). Falls back to deck height if it never does."""
    sign = 1.0 if t.centre[1] >= 0 else -1.0
    a = np.radians(cfg.rear_strut_angle_deg)
    x, y0, z0 = t.centre[0], t.centre[1], t.top_z
    s = np.arange(0.0, 0.8, 0.002)
    ys, zs = y0 + sign * s * np.sin(a), z0 + s * np.cos(a)
    yc = sign * hull_centre_y(d)
    b = halfbreadth_at(np.full_like(s, x), zs, d)
    inside = np.nan_to_num(b - np.abs(ys - yc), nan=-1.0) >= 0
    idx = np.nonzero(inside | (zs >= d.depth))[0]
    k = idx[0] if len(idx) else len(s) - 1
    return Strut(t.name, np.array([x, y0, z0]), np.array([x, float(ys[k]), float(zs[k])]),
                 cfg.rear_strut_chord, cfg.rear_strut_tc, aligned=True)


def bow_face(d: Design) -> tuple[float, float]:
    """Bottom and top z of the flat bow face (keel at the bow to the deck)."""
    return float(keel_z(np.array([1.0]), d)[0]), d.depth


def mount_plate_centre(cfg: Config, t: Thruster) -> np.ndarray:
    """Centre of the plate the T200 base bolts to (on top of the thruster)."""
    return t.centre + np.array([0.0, 0.0, cfg.t200_base_offset + cfg.mount_plate_t / 2])


def front_members(d: Design, cfg: Config, t: Thruster, pod: Pod) -> list[Strut]:
    """Front bracket: a back plate bolted to the flat bow face, and two members
    from it to the T200 mount plate (a triangle). Upper: from near the top of the
    bow face down-forward; lower: from near the bottom, roughly level."""
    zb, zt = bow_face(d)
    m = cfg.bow_plate_margin
    x0 = d.loa + cfg.bow_plate_t
    pc = mount_plate_centre(cfg, t)
    half = cfg.mount_plate_size / 2
    upper = Strut(t.name + "_upper", np.array([x0, pod.y, zt - m]), pc + np.array([half, 0, 0]),
                  cfg.front_strut_chord, cfg.front_strut_tc, aligned=False)
    lower = Strut(t.name + "_lower", np.array([x0, pod.y, zb + m]), pc - np.array([half, 0, 0]),
                  cfg.front_strut_chord, cfg.front_strut_tc, aligned=False)
    return [upper, lower]


def mount_mass(members: list, cfg: Config) -> float:
    """Printed mount: its members plus one clamp/bracket."""
    return sum(m.mass(cfg) for m in members) + cfg.mount_fixed_kg


@dataclass
class SkidPad:
    """One trapezoidal pad (side view): flat bottom of length `flat` at z_bottom,
    ends sloping up at `slope_deg` until they meet the keel; top follows the keel."""

    xc: float
    flat: float
    xs: np.ndarray
    heights: np.ndarray

    @property
    def flat_x0(self) -> float:
        return self.xc - self.flat / 2

    @property
    def flat_x1(self) -> float:
        return self.xc + self.flat / 2

    @property
    def top_length(self) -> float:
        return float(self.xs[-1] - self.xs[0])


@dataclass
class Skid:
    """Two printed pads per hull, fore and aft (same on both hulls).

    The boat rests on the four flat pad bottoms; it sits stable on flat ground
    when the LCG falls between the aft and forward pads (x0..x1).
    """

    pads: list
    z_bottom: float
    width: float

    @property
    def x0(self) -> float:          # aft end of the support span
        return min(p.flat_x0 for p in self.pads)

    @property
    def x1(self) -> float:          # forward end of the support span
        return max(p.flat_x1 for p in self.pads)

    @property
    def xs(self) -> np.ndarray:     # combined profile (zero between pads)
        return np.concatenate([p.xs for p in self.pads])

    @property
    def heights(self) -> np.ndarray:
        return np.concatenate([p.heights for p in self.pads])

    @property
    def length(self) -> float:      # reference length for friction (mean pad length)
        return float(np.mean([p.top_length for p in self.pads]))

    @property
    def max_height(self) -> float:
        return float(max(p.heights.max() for p in self.pads))

    def _int(self, f) -> float:
        return float(sum(np.trapezoid(f(p), p.xs) for p in self.pads))

    @property
    def volume(self) -> float:
        # rounded section: ~0.9 of the bounding rectangle
        return 0.9 * self.width * self._int(lambda p: p.heights)

    @property
    def centroid(self) -> np.ndarray:
        A = self._int(lambda p: p.heights)
        xc = self._int(lambda p: p.heights * p.xs) / A
        zc = self._int(lambda p: p.heights * (self.z_bottom + p.heights / 2)) / A
        return np.array([xc, 0.0, zc])

    @property
    def wetted_area(self) -> float:
        return self._int(lambda p: 2 * p.heights + self.width)

    @property
    def lateral_area(self) -> float:
        return self._int(lambda p: p.heights)


def make_pad(d: Design, xc: float, flat: float, zb: float, slope_deg: float) -> SkidPad:
    t = np.tan(np.radians(slope_deg))
    # widest possible span: flat plus slopes up to the highest keel point
    reach = (max(d.rocker_fwd, d.rocker_aft) - zb) / t
    xs = np.linspace(max(0.0, xc - flat / 2 - reach), min(d.loa, xc + flat / 2 + reach), 161)
    bottom = zb + np.maximum(np.abs(xs - xc) - flat / 2, 0.0) * t
    h = keel_z(xs / d.loa, d) - bottom
    keep = h > 0
    idx = np.nonzero(keep)[0]
    xs, h = xs[idx[0]: idx[-1] + 1], np.maximum(h[idx[0]: idx[-1] + 1], 0.0)
    return SkidPad(xc, flat, xs, h)


def make_skid(d: Design, cfg: Config, lowest_other_z: float) -> Skid:
    """Fore and aft pads. Their flat bottoms sit `skid_margin` below the lowest
    thruster/pod point, and at least `skid_min_height` below the keel."""
    zb = min(lowest_other_z - d.skid_margin, -cfg.skid_min_height)
    L, f = d.loa, d.skid_pad_len
    pads = [make_pad(d, float(np.clip(x * L, f / 2, L - f / 2)), f, zb, cfg.skid_end_slope_deg)
            for x in (d.skid_aft_x_frac, d.skid_fwd_x_frac)]
    return Skid(pads, zb, cfg.skid_width)
