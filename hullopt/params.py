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
    loa: float = 1.15                # hull length overall (m), front pods excluded
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
    front_thr_gap: float = -0.03      # bow tip to the aft end of the front thruster bay (m);
                                      # negative: the narrow cap tip tucks into the bay's empty corner
    rear_thr_x_frac: float = 0.06     # rear thruster centre, fraction of LOA from transom
    front_thr_cover: float = 0.06     # water above the top of the front pod at the design waterline (m)
    rear_thr_cover: float = 0.06      # water above the top of the rear thrusters at the design waterline (m)

    # --- skid: two replaceable printed trapezoid pads under each hull ---
    skid_margin: float = 0.015       # skid bottom below lowest thruster/pod point (m)
    skid_aft_x_frac: float = 0.25    # aft pad centre, fraction of LOA from transom
    skid_fwd_x_frac: float = 0.70    # forward pad centre, fraction of LOA from transom
    skid_pad_len: float = 0.08       # flat bottom length of each pad (m)

    # --- front pod (faired, impact-absorbing, around the front thruster) ---
    pod_wall: float = 0.003          # pod shell thickness (m)
    pod_nose_len: float = 0.006      # pod skin ahead of the thruster bay, up to the front hoop (m)

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "Design":
        names = {f.name for f in fields(cls)}
        return cls(**{k: float(v) for k, v in d.items() if k in names})


# Optimizer bounds. Keys must match Design fields.
BOUNDS: dict[str, tuple[float, float]] = {
    "loa": (1.00, 1.50),  # hull only; front pods sit ahead of it within loa_max
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
    "front_thr_gap": (-0.06, 0.06),
    "rear_thr_x_frac": (0.0, 0.25),
    "front_thr_cover": (0.05, 0.15),
    "rear_thr_cover": (0.05, 0.15),
    "skid_margin": (0.010, 0.020),
    "skid_aft_x_frac": (0.08, 0.45),
    "skid_fwd_x_frac": (0.55, 0.90),
    "skid_pad_len": (0.03, 0.15),
    "pod_wall": (0.0016, 0.006),
    "pod_nose_len": (0.005, 0.05),
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
    front_thr_angle_deg: float = 45.0  # axis angle from +x
    front_thr_toe_in: bool = True    # axis points forward-inboard (jets run aft-outboard)
    pod_gap: float = 0.005           # gap between thruster/mount plate and pod inner wall
    t200_base_offset: float = 0.052  # PLACEHOLDER: thruster axis to T200 mounting face (check BR 3D model)
    t200_bolt_square: float = 0.019  # 4x M3 on a 19 mm square (M6 centre hole unused)
    mount_plate_t: float = 0.006     # front bracket plate the T200 bolts to
    mount_plate_size: float = 0.044  # square (T200 pattern; the upper strut lands on its outboard corner)
    hull_bolt_pattern: tuple = (0.040, 0.060)  # 4x M4 into the hull: width (y) x height (z)
    # bow cap: the hull is printed blunt where it is cap_min_half*2 wide; the mount's
    # cap slides over a spigot on it and narrows at cap_entry_deg to the bow tip (loa)
    cap_entry_deg: float = 25.0      # cap's entry half-angle (plan view) from the cut to the tip
    cap_min_half: float = 0.012      # cut no further forward than where the hull is 24 mm wide
    cap_sleeve_len: float = 0.025    # cap overlaps the hull spigot by this much
    cap_wall: float = 0.003          # cap sleeve wall (hull spigot is inset by this)
    cap_margin: float = 0.015        # strut roots this far inside the cap top
    forefoot_chord: float = 0.050    # forefoot fin below the cap: fore-aft length (max; stops short of the thruster)
    frame_height_max: float = 0.245  # frame's top-to-bottom size (it lies on the bed nose-up)
    forefoot_thick: float = 0.012
    rail_clear: float = 0.004        # lower rail below the thruster
    hoop_len: float = 0.008          # screw lugs at the strut ends: the nose fairing bolts on
    bumper_tip_r: float = 0.010      # small flat at the bumper tip
    max_overhang_deg: float = 35.0   # printable overhang from vertical (no supports)
    inlet_cone_deg: float = 45.0     # keep-out cone at the duct inlet: half-angle,
    inlet_cone_len_D: float = 0.5    #   and length in duct diameters
    outlet_cone_deg: float = 15.0    # keep-out cone for the jet: half-angle,
    outlet_cone_len_D: float = 3.0   #   and length in duct diameters
    pod_window_ratio: float = 1.0    # side window area per side / thruster disk area
    pod_window_max_frac: float = 0.7  # windows may use at most this share of the bay side
    pod_window_cd: float = 0.2       # PLACEHOLDER: drag of the open windows, on window area
    pod_crush_fill: float = 0.0      # nose infill; 0 = hollow nose, printed as a plain shell
    pod_crush_len_min: float = 0.03  # PLACEHOLDER: minimum crush length ahead of the thruster
    pod_bed: tuple = (0.25, 0.25, 0.25)  # each pod must print in one piece on this bed
    impact_speed: float = 2.0        # head-on impact speed for the energy estimate
    rear_strut_chord: float = 0.048  # diagonal rear strut, chord along the flow (m)
    rear_strut_tc: float = 0.25      # -> 12 mm thick
    rear_strut_angle_deg: float = 45.0  # from vertical, outboard-up into the hull side
    strut_min_thickness: float = 0.012  # every mount member: at least 12 mm x 24 mm
    strut_min_chord: float = 0.024
    mount_fixed_kg: float = 0.06     # printed clamp/bracket per thruster mount
    front_upper_depth: float = 0.024  # front upper member, upright: depth (m) x thickness across y
    front_upper_thick: float = 0.012
    front_lower_width: float = 0.040  # front lower member, flat plate: width across y x thickness
    front_lower_thick: float = 0.012
    stem_embed: float = 0.02         # members run this far into the stem
    stem_attach_margin: float = 0.015  # member ends this far below deck / above keel
    member_cd: float = 1.0           # pressure drag of front members on thickness x height

    # --- printing ---
    fil_density: float = 1240.0      # PLA (kg/m^3), prototype material
    wall_thickness: float = 0.0016   # 4 perimeters x 0.4 mm
    rib_vol_frac: float = 0.012      # ribs/infill as fraction of enclosed volume
    sealer_kg_per_m2: float = 0.12
    strut_wall: float = 0.0012
    skid_width: float = 0.020
    skid_min_height: float = 0.005
    skid_fill: float = 0.6           # skid printed fairly dense
    skid_end_slope_deg: float = 45.0  # pad ends rise at this angle to meet the keel
    mass_calibration: float = 1.0    # measured / predicted printed mass, from a test print
    bed: tuple = (0.25, 0.25, 0.25)  # printer build volume (m), standard 250 mm bed
    joint_kg_per_m: float = 0.04     # flange/bond mass per metre of joint girth

    # --- constraints ---
    loa_max: float = 1.52            # whole boat, front pods included
    width_max: float = 3 * 12 * INCH  # PLACEHOLDER: overall width limit
    spacing_min: float = 10 * INCH
    gm_t_min: float = 0.30           # PLACEHOLDER (X in the brief)
    freeboard_min: float = 0.08      # PLACEHOLDER (Y in the brief)
    trim_max_deg: float = 1.5
    heel_sway_max_deg: float = 5.0
    thr_cover_min: float = 0.05      # water above top of rear thruster / front pod; verify in CFD
    clearance_min: float = 0.010     # thruster envelope to hull, and between thrusters
    jet_clear_min_D: float = 2.0     # free jet length before hitting a hull, in duct diameters
    skid_max_height: float = 0.12    # PLACEHOLDER: tallest acceptable rail. This limit drives the design.
    skid_lcg_margin: float = 0.05    # LCG must sit this far inside the pad support span
    reserve_buoyancy_min: float = 1.0  # (enclosed - displaced) / displaced

    # --- resistance proxy ---
    surge_speeds: tuple = (0.5, 0.75, 1.0, 1.25, 1.5, 2.0)
    surge_weights: tuple = (0.15, 0.20, 0.25, 0.20, 0.15, 0.05)
    sway_speed: float = 0.4
    form_factor_k: float = 0.10
    roughness_dcf: float = 4.0e-4    # printed + sealed surface allowance
    sway_cd: float = 0.9

    stage1_wave: bool = True         # add a coarse Michell wave term to the stage 1 drag objective

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

