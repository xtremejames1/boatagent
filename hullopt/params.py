"""Design variables, their bounds, and fixed configuration.

All values are SI (m, kg, s, N) unless the name says otherwise.

Body frame used everywhere: x forward from the transom (x=0) to the bow
(x=LOA), y to port, z up, z=0 at the lowest keel point of the bare hull.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, fields, replace

INCH = 0.0254
LB = 0.45359237


@dataclass
class Design:
    """One candidate design. Every field here is an optimizer variable."""

    # --- demihull ---
    loa: float = 1.45                # hull length overall (m)
    beam: float = 0.20               # demihull max beam (m)
    depth: float = 0.24              # lowest keel point to deck (m)
    deadrise_deg: float = 10.0       # bottom deadrise angle
    bilge_frac: float = 0.6          # bilge radius / local half-beam (0 = hard chine)
    rocker_fwd: float = 0.03         # keel rise at the bow (m)
    rocker_aft: float = 0.01         # keel rise at the transom (m)
    entry_half_angle_deg: float = 12.0  # target half-angle of entry (plan view)
    transom_frac: float = 0.55       # transom half-beam / max half-beam
    xmax_frac: float = 0.50          # max beam and lowest keel, fraction of LOA from transom

    # --- layout ---
    spacing_inner: float = 12 * INCH  # gap between demihull inner faces (m)
    front_thr_x_frac: float = 0.92    # front thruster centre, fraction of LOA from transom
    rear_thr_x_frac: float = 0.06     # rear thruster centre, fraction of LOA from transom
    front_thr_cover: float = 0.06     # water above the top of the front shroud at the design waterline (m)
    rear_thr_cover: float = 0.06      # water above the top of the rear thrusters at the design waterline (m)

    # --- skid (replaceable printed rail under each hull) ---
    skid_margin: float = 0.015       # skid bottom below lowest thruster/shroud point (m)
    skid_len_frac: float = 0.55      # skid length / LOA
    skid_x_frac: float = 0.45        # skid centre, fraction of LOA from transom

    # --- front shroud ---
    shroud_guard: float = 0.008      # shroud wall thickness around/below thruster (m)
    shroud_chord: float = 0.12       # fairing strut chord (m)
    shroud_tc: float = 0.15          # fairing strut thickness/chord

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "Design":
        names = {f.name for f in fields(cls)}
        return cls(**{k: float(v) for k, v in d.items() if k in names})


# Optimizer bounds. Keys must match Design fields.
BOUNDS: dict[str, tuple[float, float]] = {
    "loa": (1.20, 1.52),
    "beam": (0.12, 0.30),
    "depth": (0.15, 0.35),
    "deadrise_deg": (0.0, 25.0),
    "bilge_frac": (0.0, 1.0),
    "rocker_fwd": (0.0, 0.08),
    "rocker_aft": (0.0, 0.04),
    "entry_half_angle_deg": (6.0, 25.0),
    "transom_frac": (0.2, 0.9),
    "xmax_frac": (0.35, 0.65),
    "spacing_inner": (10 * INCH, 0.50),
    "front_thr_x_frac": (0.75, 1.02),
    "rear_thr_x_frac": (0.0, 0.25),
    "front_thr_cover": (0.05, 0.15),
    "rear_thr_cover": (0.05, 0.15),
    "skid_margin": (0.010, 0.020),
    "skid_len_frac": (0.30, 0.80),
    "skid_x_frac": (0.30, 0.60),
    "shroud_guard": (0.004, 0.015),
    "shroud_chord": (0.08, 0.20),
    "shroud_tc": (0.10, 0.25),
}


@dataclass
class Config:
    """Fixed inputs. Items marked PLACEHOLDER need a real number from the team."""

    # --- water ---
    rho: float = 1000.0              # fresh water (kg/m^3)
    g: float = 9.81
    nu: float = 1.0e-6               # kinematic viscosity (m^2/s), ~20 C fresh water

    # --- masses and centres (z measured above the deck) ---
    payload_kg: float = 35 * LB      # everything carried except hulls/thrusters/frame
    lidar_kg: float = 1.0            # part of payload, on the mast
    lidar_x_frac: float = 0.60
    lidar_z_above_deck: float = 0.60  # PLACEHOLDER: mast height to lidar CG
    payload_z_above_deck: float = 0.10  # CG of electronics/battery box
    payload_x_range: tuple = (0.30, 0.65)  # movable payload LCG range (fraction of LOA)
    frame_fixed_kg: float = 2.0      # deck plate, brackets
    frame_kg_per_m_width: float = 1.2  # cross beams, scales with overall width
    frame_z_above_deck: float = 0.03

    # --- thrusters (Blue Robotics T200) ---
    thr_mass: float = 0.344          # in air, with cable
    thr_volume: float = 1.88e-4      # (344 g air - 156 g water) / rho
    thr_length: float = 0.113
    thr_radius: float = 0.050        # duct outer radius
    thr_force_kgf: float = 5.25      # forward thrust at 16 V
    thr_rev_force_kgf: float = 4.1   # reverse thrust at 16 V
    front_thr_angle_deg: float = 45.0  # axis angle from +x; axis points forward-outboard
    shroud_gap: float = 0.005        # radial gap thruster -> shroud
    rear_strut_chord: float = 0.08
    rear_strut_tc: float = 0.20

    # --- printing ---
    fil_density: float = 1270.0      # PETG (kg/m^3)
    wall_thickness: float = 0.0016   # 4 perimeters x 0.4 mm
    rib_vol_frac: float = 0.012      # ribs/infill as fraction of enclosed volume
    sealer_kg_per_m2: float = 0.12
    strut_wall: float = 0.0012
    skid_width: float = 0.020
    skid_min_height: float = 0.005
    skid_fill: float = 0.6           # skid printed fairly dense
    mass_calibration: float = 1.0    # measured / predicted printed mass, from a test print
    bed: tuple = (0.256, 0.256, 0.256)  # printer build volume (m)
    joint_kg_per_m: float = 0.04     # flange/bond mass per metre of joint girth

    # --- constraints ---
    loa_max: float = 1.52
    width_max: float = 3 * 12 * INCH  # PLACEHOLDER: overall width limit
    spacing_min: float = 10 * INCH
    gm_t_min: float = 0.30           # PLACEHOLDER (X in the brief)
    freeboard_min: float = 0.08      # PLACEHOLDER (Y in the brief)
    trim_max_deg: float = 1.5
    heel_sway_max_deg: float = 5.0
    thr_cover_min: float = 0.05      # water above top of thruster/shroud; verify in CFD
    clearance_min: float = 0.010     # thruster envelope to hull, and between thrusters
    jet_clear_min_D: float = 2.0     # free jet length before hitting a hull, in duct diameters
    skid_max_height: float = 0.12    # PLACEHOLDER: tallest acceptable rail. This limit drives the design.
    skid_lcg_margin: float = 0.05    # LCG must sit this far inside the skid ends
    reserve_buoyancy_min: float = 1.0  # (enclosed - displaced) / displaced

    # --- resistance proxy ---
    surge_speeds: tuple = (0.5, 0.75, 1.0, 1.25, 1.5, 2.0)
    surge_weights: tuple = (0.15, 0.20, 0.25, 0.20, 0.15, 0.05)
    sway_speed: float = 0.4
    form_factor_k: float = 0.10
    roughness_dcf: float = 4.0e-4    # printed + sealed surface allowance
    sway_cd: float = 0.9

    # --- mesh resolution ---
    n_stations: int = 61

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "Config":
        names = {f.name for f in fields(cls)}
        clean = {}
        for k, v in d.items():
            if k in names:
                clean[k] = tuple(v) if isinstance(v, list) else v
        return replace(cls(), **clean)

    @classmethod
    def load(cls, path: str | None) -> "Config":
        if not path:
            return cls()
        with open(path) as f:
            return cls.from_dict(json.load(f))


def clip_to_bounds(d: Design) -> Design:
    vals = d.to_dict()
    for k, (lo, hi) in BOUNDS.items():
        vals[k] = min(max(vals[k], lo), hi)
    return Design.from_dict(vals)

