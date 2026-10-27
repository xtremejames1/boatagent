"""Stage 1: algebraic evaluation and screening of catamaran hull designs.

For each design:
  1. build the demihull and catamaran meshes;
  2. compute printed hull mass and the mass budget;
  3. float the boat level to get the design waterline, then place thrusters,
     pods, struts and skids relative to it;
  4. move the movable payload to trim the boat level (within its deck range);
  5. solve free heave + trim and compute hydrostatics (GM_T, GM_L, freeboard...);
  6. run clearance, submergence, jet, skid and size checks;
  7. compute cheap proxies (friction drag, slenderness, sway drag, heel in sway).

Objectives (passed to stage 2 for ranking on real resistance):
  minimise total mass, minimise weighted friction drag, maximise slenderness.

Usage:
  python -m hullopt.stage1 run --trials 3000 --out results/stage1
  python -m hullopt.stage1 eval [--design d.json] [--stl boat.stl] [--gz]
  python -m hullopt.stage1 dump-config > my_config.json
  python -m hullopt.stage1 calibrate --measured-kg 0.41 --design test_piece.json
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
import time
from dataclasses import dataclass, replace

import numpy as np

from .geometry import (
    actual_entry_angle_deg,
    catamaran_mesh,
    demihull_mesh,
    hull_centre_y,
    write_stl,
)
from .hydrostatics import Appendages, Pose, equilibrium, gz_curve, hydro, solve_heave
from .layout import (
    Strut,
    Thruster,
    hull_gap,
    jet_free_length,
    make_skid,
    pair_gap,
    mount_mass,
    Pod,
    place_front,
    place_rear,
    rear_strut,
    point_to_ray,
)
from .mass import MassBudget, demihull_mass
from .params import BOUNDS, Config, Design

CONSTRAINTS = [
    "c_loa", "c_width", "c_spacing", "c_gm_t", "c_freeboard", "c_trim", "c_heel_sway",
    "c_thr_cover", "c_clearance", "c_jet", "c_skid_height", "c_skid_lcg", "c_reserve",
    "c_pod_print", "c_pod_crush", "c_pod_window", "c_strut_size",
]
OBJECTIVES = [("obj_mass", "min"), ("obj_drag", "min"), ("obj_slenderness", "max")]
FAIL = 10.0


def ittc57(re: float) -> float:
    return 0.075 / (math.log10(max(re, 1e3)) - 2.0) ** 2


@dataclass
class Detail:
    """Geometry kept for export and plotting."""

    design: Design
    cat_tri: np.ndarray
    pose: Pose | None
    thrusters: list
    skid: object
    budget: MassBudget
    struts: list
    z_wl: float           # body z of the design waterline (level float)
    pods: list            # front pods (port, starboard)


def _failed(d: Design, reason: str) -> dict:
    row = d.to_dict()
    row.update({k: FAIL for k in CONSTRAINTS})
    row.update({"feasible": False, "fail_reason": reason,
                "obj_mass": 1e3, "obj_drag": 1e3, "obj_slenderness": 0.0})
    return row


def evaluate(d: Design, cfg: Config | None = None, detail: bool = False):
    cfg = cfg or Config()
    L = d.loa
    rho, g = cfg.rho, cfg.g
    pivot = np.array([L / 2, 0.0, 0.0])
    yc = hull_centre_y(d)
    deck = d.depth
    width = d.spacing_inner + 2 * d.beam

    demi = demihull_mesh(d, cfg.n_stations)
    cat = catamaran_mesh(d, cfg.n_stations)
    tri = cat.triangles
    hm = demihull_mass(demi.triangles, d, cfg)

    # ---- fixed part of the mass budget
    bud = MassBudget()
    bud.add("hulls", hm.mass, hm.cg + [0, yc, 0])
    bud.add("hulls", hm.mass, hm.cg + [0, -yc, 0])
    frame = cfg.frame_fixed_kg + cfg.frame_kg_per_m_width * width
    bud.add("frame", frame, [L / 2, 0, deck + cfg.frame_z_above_deck])
    bud.add("lidar", cfg.lidar_kg, [cfg.lidar_x_frac * L, 0, deck + cfg.lidar_z_above_deck])
    m_mov = cfg.payload_kg - cfg.lidar_kg

    fixed = list(bud.items)
    psi = math.radians(cfg.front_thr_angle_deg)
    pod_half_h = cfg.thr_radius + cfg.pod_gap + d.pod_wall

    def mirror(t: Thruster, name: str) -> Thruster:
        c = t.centre.copy()
        c[1] = -c[1]
        return Thruster(name, c, -t.heading, t.radius, t.length)

    def mirror_pod(p: Pod) -> Pod:
        return replace(p, name="front_stbd", y=-p.y)

    # Thrusters hang a set cover below the design waterline, but the waterline
    # depends on the appendage masses. Pass 1 uses an estimate; pass 2 the result.
    m_est = bud.total + m_mov + 4 * cfg.thr_mass + 0.6
    z_wl = solve_heave(tri, m_est, rho, Pose(0.0, pivot=pivot), h_lo=-0.05, h_hi=deck)
    for _pass in range(2):
        if z_wl is None:
            return (_failed(d, "sinks"), None) if detail else _failed(d, "sinks")
        bud = MassBudget(list(fixed))

        # ---- thrusters (front: in pods ahead of the bows; rear: between the hulls)
        fz = z_wl - d.front_thr_cover - pod_half_h
        rz = z_wl - d.rear_thr_cover - cfg.thr_radius
        # port heading: toe-in points the axis forward-inboard (toward -y)
        fp, pp = place_front(d, cfg, fz, -psi if cfg.front_thr_toe_in else psi)
        pods = [pp, mirror_pod(pp)]
        rp = place_rear(d, d.rear_thr_x_frac * L, rz, cfg.thr_radius, cfg.thr_length, cfg.clearance_min)
        fs, rs = mirror(fp, "front_stbd"), mirror(rp, "rear_stbd")
        thrusters = [fp, fs, rp, rs]

        struts = []
        for t in (fp, fs):   # front: vertical pylon from the deck-level arm down to the pod
            struts.append(Strut(t.name, t.centre[0], t.centre[1], deck, pp.top_z, d.pylon_chord, d.pylon_tc))
        for t in (rp, rs):   # rear: diagonal from the thruster up into the hull side
            struts.append(rear_strut(d, cfg, t))

        lowest = min(b.lowest_z for b in pods + [rp, rs])
        skid = make_skid(d, cfg, lowest)
        skid_mass = skid.volume * cfg.fil_density * cfg.skid_fill * cfg.mass_calibration

        app = Appendages()
        for t in thrusters:
            bud.add("thrusters", cfg.thr_mass, t.centre)
            app.add(cfg.thr_volume, t.centre)
        for p_ in pods:
            bud.add("pods", p_.mass(cfg), p_.centre)
            app.add(p_.material_volume(cfg), p_.centre)
        for s_ in struts:
            arm = max(s_.x - L, 0.0) if s_.name.startswith("front") else 0.0
            bud.add("mounts", mount_mass(s_, arm, cfg), s_.cg)
            vs = s_.volume_submerged(z_wl)
            if vs > 0:
                app.add(vs, s_.submerged_centre(z_wl))
        sc = skid.centroid
        for sy in (yc, -yc):
            bud.add("skids", skid_mass, [sc[0], sy, sc[2]])
            app.add(skid.volume, [sc[0], sy, sc[2]])

        m_tot = bud.total + m_mov
        h1 = solve_heave(tri, m_tot, rho, Pose(0.0, pivot=pivot), app, h_lo=-0.05, h_hi=deck)
        if _pass == 0:
            z_wl = h1
    if h1 is None:
        return (_failed(d, "sinks"), None) if detail else _failed(d, "sinks")

    # ---- movable payload: put LCG over the level-trim LCB, within its deck range
    lcb = hydro(tri, Pose(h1, pivot=pivot), app).cob[0] + pivot[0]
    others = sum(i.mass * i.pos[0] for i in bud.items)
    x_mov = (lcb * m_tot - others) / m_mov
    lo, hi = cfg.payload_x_range
    x_mov = float(np.clip(x_mov, lo * L, hi * L))
    bud.add("payload", m_mov, [x_mov, 0, deck + cfg.payload_z_above_deck])

    M, G = bud.total, bud.cg
    pose, hy = equilibrium(tri, M, G, rho, pivot, app, h_hi=deck)
    if pose is None:
        return (_failed(d, "no_equilibrium"), None) if detail else _failed(d, "no_equilibrium")

    # ---- hydrostatics
    Ge = pose.to_earth(G[None])[0]
    V = hy.volume
    gm_t = hy.cob[2] + hy.wp.i_t / V - Ge[2]
    gm_l = hy.cob[2] + hy.wp.i_l / V - Ge[2]
    ve = pose.to_earth(cat.vertices)
    deck_v = cat.vertices[:, 2] >= deck - 1e-9
    freeboard = float(ve[deck_v, 2].min())
    hull_draft = float(-ve[:, 2].min())
    skid_bottom_e = pose.to_earth(np.array([[sc[0], yc, skid.z_bottom]]))[0, 2]
    total_draft = max(hull_draft, -skid_bottom_e)

    port_tri = demi.triangles + np.array([0, yc, 0])
    hy_p = hydro(port_tri, pose)
    lwl = hy_p.wp.x_max - hy_p.wp.x_min
    bwl = hy_p.wp.y_max - hy_p.wp.y_min
    t_hull = hull_draft
    vol_hull = hy_p.mesh_volume
    cb = vol_hull / max(lwl * bwl * t_hull, 1e-9)
    cwp = hy_p.wp.area / max(lwl * bwl, 1e-9)
    disp_per_hull = M / (2 * rho)
    slender = lwl / disp_per_hull ** (1 / 3)
    reserve = (2 * hm.enclosed_volume - hy.mesh_volume) / hy.mesh_volume

    # ---- resistance proxies
    s_hull = hy.wetted_area
    app_wet = []  # (wetted area, reference length)
    for s in struts:
        a = s.wetted_area(z_wl)
        if a > 0:
            app_wet.append((a, s.chord))
    app_wet += [(p_.outer_area, p_.length) for p_ in pods] + [(skid.wetted_area, skid.length)] * 2
    k = cfg.form_factor_k
    r_f = []
    for v in cfg.surge_speeds:
        q = 0.5 * rho * v * v
        r = q * s_hull * (1 + k) * (ittc57(v * lwl / cfg.nu) + cfg.roughness_dcf)
        r += sum(q * a * (1 + k) * (ittc57(v * lr / cfg.nu) + cfg.roughness_dcf) for a, lr in app_wet)
        r_f.append(r)
    w = np.asarray(cfg.surge_weights, float)
    drag_w = float(np.dot(w, r_f) / w.sum())
    fr = [v / math.sqrt(g * lwl) for v in cfg.surge_speeds]

    a_lat = hy.lateral_area + 2 * float(np.trapezoid(skid.heights, skid.xs))
    a_lat += sum(s.submerged_height(z_wl) * s.chord for s in struts)
    a_lat += sum(p_.length * 2 * p_.half_h for p_ in pods)
    sway_drag = 0.5 * rho * cfg.sway_speed ** 2 * cfg.sway_cd * a_lat
    f_sway = 2 * cfg.thr_force_kgf * g * math.sin(psi)
    v_sway_max = math.sqrt(2 * f_sway / (rho * cfg.sway_cd * a_lat))
    f_surge = 2 * cfg.thr_force_kgf * g * (math.cos(psi) + 1.0)

    z_thr_e = pose.to_earth(fp.centre[None])[0, 2]
    lever = abs(z_thr_e - hy.lateral_centroid[1])
    heel_sway = math.degrees(math.atan(f_sway * lever / (M * g * max(gm_t, 1e-6))))
    clr_offset = (hy.lateral_centroid[0] - Ge[0]) / lwl  # + means lateral area forward of G

    # ---- checks
    bodies = pods + [rp, rs]   # outer envelopes in the water
    covers = []
    for b in bodies:
        top = pose.to_earth(np.array([[b.centre[0], b.centre[1], b.top_z]]))[0, 2]
        covers.append(-top)
    cover = float(min(covers))
    gaps = [hull_gap(b.envelope_points(), d) for b in bodies]
    hull_clear = float(min(gaps))
    pair_f = 2 * (pp.y - pp.half_w)
    pair_r = pair_gap(rp)
    clearance = min(hull_clear, pair_f, pair_r, d.front_thr_gap)

    jr = cfg.thr_radius
    jet_aft = jet_free_length(fp, -fp.axis, d, jr)   # aft-outboard with toe-in
    jet_fwd = jet_free_length(fp, fp.axis, d, jr)
    jet_rear_fwd = jet_free_length(rp, rp.axis, d, jr)
    start = fp.centre - fp.axis * fp.length / 2
    fjet_to_rear = min(point_to_ray(t.centre, start, -fp.axis) for t in (rp, rs)) - rp.radius - jr
    rstart = rp.centre + rp.axis * rp.length / 2
    rjet_to_front = min(point_to_ray(p_.centre, rstart, rp.axis) for p_ in pods) - pp.half_w - jr
    jet_min = cfg.jet_clear_min_D * 2 * cfg.thr_radius
    jet_len = min(jet_aft, jet_fwd)

    env = np.vstack([b.envelope_points() for b in bodies])
    x_lo, x_hi = min(0.0, env[:, 0].min()), max(L, env[:, 0].max())
    loa_overall = x_hi - x_lo
    lcg = G[0]
    skid_lcg = min(lcg - skid.x0, skid.x1 - lcg)

    # ---- front pod: one-piece print, crush zone, side windows, impact energy
    pod_print = pp.print_ratio(cfg.pod_bed)
    win_frac = pp.window_area / pp.bay_side_area
    rear_st = [s_ for s_ in struts if s_.name.startswith("rear")]
    strut_size = max(max(cfg.strut_min_thickness / s_.thickness, cfg.strut_min_chord / s_.chord)
                     for s_ in rear_st)
    impact_e = 0.5 * M * cfg.impact_speed ** 2

    c = {
        "c_loa": (loa_overall - cfg.loa_max) / cfg.loa_max,
        "c_width": (width - cfg.width_max) / cfg.width_max,
        "c_spacing": (cfg.spacing_min - d.spacing_inner) / cfg.spacing_min,
        "c_gm_t": (cfg.gm_t_min - gm_t) / cfg.gm_t_min,
        "c_freeboard": (cfg.freeboard_min - freeboard) / cfg.freeboard_min,
        "c_trim": (abs(math.degrees(pose.trim)) - cfg.trim_max_deg) / cfg.trim_max_deg,
        "c_heel_sway": (heel_sway - cfg.heel_sway_max_deg) / cfg.heel_sway_max_deg,
        "c_thr_cover": (cfg.thr_cover_min - cover) / cfg.thr_cover_min,
        "c_clearance": (cfg.clearance_min - clearance) / cfg.clearance_min,
        "c_jet": (jet_min - jet_len) / jet_min,
        "c_skid_height": (skid.max_height - cfg.skid_max_height) / cfg.skid_max_height,
        "c_skid_lcg": (cfg.skid_lcg_margin - skid_lcg) / cfg.skid_lcg_margin,
        "c_reserve": (cfg.reserve_buoyancy_min - reserve) / cfg.reserve_buoyancy_min,
        "c_pod_print": pod_print - 1.0,
        "c_pod_crush": (cfg.pod_crush_len_min - pp.nose_len) / cfg.pod_crush_len_min,
        "c_pod_window": (win_frac - cfg.pod_window_max_frac) / cfg.pod_window_max_frac,
        "c_strut_size": strut_size - 1.0,
    }
    c = {k: float(np.clip(v, -FAIL, FAIL)) for k, v in c.items()}
    feasible = all(v <= 0 for v in c.values())

    row = d.to_dict()
    row.update({
        "feasible": feasible,
        "fail_reason": "" if feasible else ",".join(k[2:] for k, v in c.items() if v > 0),
        "obj_mass": M,
        "obj_drag": drag_w,
        "obj_slenderness": slender,
        "mass_total": M,
        **{f"mass_{k}": v for k, v in bud.by_name().items()},
        "hull_printed_each": hm.shell + hm.ribs + hm.joints,
        "print_pieces_per_hull": hm.pieces,
        "lcg": lcg, "vcg": G[2],
        "payload_x": x_mov,
        "heave": pose.heave, "trim_deg": math.degrees(pose.trim),
        "draft_hull": hull_draft, "draft_total": total_draft,
        "freeboard": freeboard,
        "gm_t": gm_t, "gm_l": gm_l,
        "lwl": lwl, "bwl": bwl, "cb": cb, "cwp": cwp,
        "loa_overall": loa_overall, "width": width,
        "slenderness": slender,
        "reserve_buoyancy": reserve,
        "wetted_hulls": s_hull,
        "entry_angle_actual_deg": actual_entry_angle_deg(d),
        **{f"fr_{v:g}": f for v, f in zip(cfg.surge_speeds, fr)},
        **{f"rf_{v:g}": r for v, r in zip(cfg.surge_speeds, r_f)},
        "drag_fric_weighted": drag_w,
        "lateral_area": a_lat,
        "sway_drag": sway_drag,
        "f_sway": f_sway, "f_surge": f_surge,
        "v_sway_max_est": v_sway_max,
        "heel_sway_deg": heel_sway,
        "clr_offset": clr_offset,
        "thr_cover_min": cover,
        "thr_hull_clearance": hull_clear,
        "front_pair_gap": pair_f, "rear_pair_gap": pair_r,
        "front_thr_y": fp.centre[1], "rear_thr_y": rp.centre[1],
        "front_thr_z": fp.centre[2], "rear_thr_z": rp.centre[2],
        "jet_front_aft": jet_aft,
        "jet_front_fwd": jet_fwd,
        "jet_rear_fwd": jet_rear_fwd,
        "jet_front_to_rear_thr": fjet_to_rear,
        "jet_rear_to_front_thr": rjet_to_front,
        "skid_bottom_z": skid.z_bottom,
        "skid_max_height": skid.max_height,
        "skid_lcg_margin": skid_lcg,
        "pod_length": pp.length, "pod_width": 2 * pp.half_w, "pod_height": 2 * pp.half_h,
        "pod_mass_each": pp.mass(cfg),
        "pod_print_ratio": pod_print,
        "pod_window_frac": win_frac,
        "rear_strut_length": rear_st[0].length,
        "rear_strut_top_z": rear_st[0].z_top,
        "impact_energy": impact_e,
        "impact_crush_force": impact_e / pp.nose_len,
        **c,
    })
    if detail:
        return row, Detail(d, tri, pose, thrusters, skid, bud, struts, z_wl, pods)
    return row


# ---------------------------------------------------------------- selection
def objective_matrix(rows: list[dict]) -> np.ndarray:
    """Objectives as a minimisation matrix."""
    F = np.array([[r[k] for k, _ in OBJECTIVES] for r in rows], float)
    for j, (_, sense) in enumerate(OBJECTIVES):
        if sense == "max":
            F[:, j] = -F[:, j]
    return F


def nondominated_fronts(F: np.ndarray) -> list[np.ndarray]:
    n = len(F)
    if n == 0:
        return []
    le = (F[:, None, :] <= F[None, :, :]).all(-1)
    lt = (F[:, None, :] < F[None, :, :]).any(-1)
    dom = le & lt                       # dom[i, j]: i dominates j
    count = dom.sum(0)
    fronts, remaining = [], np.ones(n, bool)
    while remaining.any():
        f = np.nonzero(remaining & (count == 0))[0]
        if len(f) == 0:
            break
        fronts.append(f)
        remaining[f] = False
        count = count - dom[f].sum(0)
        count[~remaining] = -1
    return fronts


def crowding(F: np.ndarray) -> np.ndarray:
    n, m = F.shape
    cd = np.zeros(n)
    if n <= 2:
        return np.full(n, np.inf)
    for j in range(m):
        o = np.argsort(F[:, j])
        span = F[o[-1], j] - F[o[0], j] or 1.0
        cd[o[0]] = cd[o[-1]] = np.inf
        cd[o[1:-1]] += (F[o[2:], j] - F[o[:-2], j]) / span
    return cd


def select_candidates(rows: list[dict], n: int) -> tuple[list[dict], list[dict]]:
    """Pareto front of the feasible rows, and up to n diverse candidates for stage 2."""
    feas = [r for r in rows if r["feasible"]]
    if not feas:
        return [], []
    F = objective_matrix(feas)
    fronts = nondominated_fronts(F)
    pareto = [feas[i] for i in fronts[0]]
    chosen: list[int] = []
    for f in fronts:
        if len(chosen) + len(f) <= n:
            chosen += list(f)
        else:
            cd = crowding(F[f])
            chosen += list(f[np.argsort(-cd)][: n - len(chosen)])
            break
    for rank, f in enumerate(fronts):
        for i in f:
            feas[i]["pareto_rank"] = rank
    return pareto, [feas[i] for i in chosen]


# ---------------------------------------------------------------- drivers
def lhs(n: int, rng: np.random.Generator) -> list[Design]:
    keys = list(BOUNDS)
    u = (rng.permuted(np.tile(np.arange(n), (len(keys), 1)), axis=1).T + rng.random((n, len(keys)))) / n
    out = []
    for row in u:
        out.append(Design.from_dict({k: lo + x * (hi - lo) for x, (k, (lo, hi)) in zip(row, BOUNDS.items())}))
    return out


def run(cfg: Config, trials: int, seed: int, sampler: str, pop: int, log_every: int = 200) -> list[dict]:
    rows: list[dict] = []
    t0 = time.time()

    def log():
        if len(rows) % log_every == 0:
            nf = sum(r["feasible"] for r in rows)
            print(f"  {len(rows):5d} evals, {nf} feasible, {time.time() - t0:.0f} s", file=sys.stderr)

    if sampler == "nsga2":
        try:
            import optuna
        except ImportError:
            print("optuna not installed; using Latin hypercube sampling", file=sys.stderr)
            sampler = "lhs"
    if sampler == "lhs":
        for d in [Design()] + lhs(trials - 1, np.random.default_rng(seed)):
            rows.append(evaluate(d, cfg))
            log()
        return rows

    optuna.logging.set_verbosity(optuna.logging.WARNING)
    # optuna >= 5 takes constraints per trial; older versions need constraints_func.
    per_trial = hasattr(optuna.trial.Trial, "set_constraint")
    s = optuna.samplers.NSGAIISampler(
        population_size=pop, seed=seed,
        constraints_func=None if per_trial else (lambda t: t.user_attrs["constraints"]),
    )
    study = optuna.create_study(directions=[sense + "imize" for _, sense in OBJECTIVES], sampler=s)
    study.enqueue_trial(Design().to_dict())

    def objective(trial):
        vals = {k: trial.suggest_float(k, lo, hi) for k, (lo, hi) in BOUNDS.items()}
        row = evaluate(Design.from_dict(vals), cfg)
        row["trial"] = trial.number
        if per_trial:
            for k in CONSTRAINTS:
                trial.set_constraint(k, row[k])
        else:
            trial.set_user_attr("constraints", [row[k] for k in CONSTRAINTS])
        rows.append(row)
        log()
        return tuple(row[k] for k, _ in OBJECTIVES)

    study.optimize(objective, n_trials=trials)
    return rows


def write_csv(path: str, rows: list[dict]) -> None:
    if not rows:
        open(path, "w").close()
        return
    keys: list[str] = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})


def export_design(row: dict, cfg: Config, folder: str, name: str) -> None:
    d = Design.from_dict(row)
    with open(os.path.join(folder, f"{name}.json"), "w") as f:
        json.dump(d.to_dict(), f, indent=2)
    write_stl(os.path.join(folder, f"{name}_catamaran.stl"), catamaran_mesh(d, cfg.n_stations).triangles, name)
    write_stl(os.path.join(folder, f"{name}_demihull.stl"), demihull_mesh(d, cfg.n_stations).triangles, name)


def representatives(pareto: list[dict]) -> dict[str, dict]:
    if not pareto:
        return {}
    F = objective_matrix(pareto)
    lo, hi = F.min(0), F.max(0)
    N = (F - lo) / np.where(hi > lo, hi - lo, 1)
    return {
        "min_mass": pareto[int(np.argmin(F[:, 0]))],
        "min_drag": pareto[int(np.argmin(F[:, 1]))],
        "max_slenderness": pareto[int(np.argmin(F[:, 2]))],
        "knee": pareto[int(np.argmin(np.linalg.norm(N, axis=1)))],
    }


SUMMARY_COLS = ["loa", "beam", "depth", "spacing_inner", "lwl", "draft_hull", "mass_total",
                "drag_fric_weighted", "slenderness", "gm_t", "freeboard", "trim_deg",
                "heel_sway_deg", "print_pieces_per_hull"]


def write_summary(path: str, rows, pareto, cands, reps, cfg: Config, seconds: float) -> str:
    nf = sum(r["feasible"] for r in rows)
    fails: dict[str, int] = {}
    for r in rows:
        for k in (r["fail_reason"] or "").split(","):
            if k:
                fails[k] = fails.get(k, 0) + 1
    lines = [
        "# Stage 1 summary", "",
        f"- evaluations: {len(rows)} in {seconds:.0f} s",
        f"- feasible: {nf}",
        f"- Pareto front (mass, friction drag, slenderness): {len(pareto)}",
        f"- stage 2 candidates: {len(cands)}", "",
        "## Constraint violations (count over all evaluations)", "",
    ]
    lines += [f"- {k}: {v}" for k, v in sorted(fails.items(), key=lambda kv: -kv[1])]
    lines += ["", "## Representative Pareto designs", "",
              "| design | " + " | ".join(SUMMARY_COLS) + " |",
              "|---" * (len(SUMMARY_COLS) + 1) + "|"]
    for name, r in reps.items():
        lines.append(f"| {name} | " + " | ".join(f"{r[k]:.3g}" for k in SUMMARY_COLS) + " |")
    text = "\n".join(lines) + "\n"
    with open(path, "w") as f:
        f.write(text)
    return text


def cmd_run(a) -> None:
    cfg = Config.load(a.config)
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    rows = run(cfg, a.trials, a.seed, a.sampler, a.pop)
    dt = time.time() - t0
    pareto, cands = select_candidates(rows, a.select)
    write_csv(os.path.join(a.out, "all.csv"), rows)
    write_csv(os.path.join(a.out, "feasible.csv"), [r for r in rows if r["feasible"]])
    write_csv(os.path.join(a.out, "pareto.csv"), pareto)
    write_csv(os.path.join(a.out, "stage2_candidates.csv"), cands)
    with open(os.path.join(a.out, "config.json"), "w") as f:
        json.dump(cfg.to_dict(), f, indent=2)
    reps = representatives(pareto)
    dd = os.path.join(a.out, "designs")
    os.makedirs(dd, exist_ok=True)
    for name, r in reps.items():
        export_design(r, cfg, dd, name)
    print(write_summary(os.path.join(a.out, "summary.md"), rows, pareto, cands, reps, cfg, dt))
    print(f"results in {a.out}/")


def cmd_eval(a) -> None:
    cfg = Config.load(a.config)
    d = Design()
    if a.design:
        with open(a.design) as f:
            d = Design.from_dict(json.load(f))
    row, det = evaluate(d, cfg, detail=True)
    for k, v in row.items():
        print(f"{k:28s} {v:.5g}" if isinstance(v, float) else f"{k:28s} {v}")
    if det is None:
        return
    if a.stl:
        write_stl(a.stl, det.cat_tri, "catamaran")
        print(f"wrote {a.stl}")
    if a.gz:
        angles = np.arange(0, 61, 5)
        gz = gz_curve(det.cat_tri, det.budget.total, det.budget.cg, cfg.rho, det.pose, angles)
        print("\nheel_deg  GZ_m")
        for ang, v in zip(angles, gz):
            print(f"{ang:8.0f}  {v:.4f}")


def cmd_calibrate(a) -> None:
    cfg = Config.load(a.config)
    with open(a.design) as f:
        d = Design.from_dict(json.load(f))
    hm = demihull_mass(demihull_mesh(d, cfg.n_stations).triangles, d, cfg)
    printed = (hm.shell + hm.ribs + hm.joints) / cfg.mass_calibration
    print(f"predicted printed mass (uncalibrated): {printed:.4f} kg")
    print(f"mass_calibration = {a.measured_kg / printed:.4f}")


def main(argv=None) -> None:
    p = argparse.ArgumentParser(prog="hullopt.stage1", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="screen designs and pick stage 2 candidates")
    r.add_argument("--trials", type=int, default=3000)
    r.add_argument("--seed", type=int, default=0)
    r.add_argument("--sampler", choices=["nsga2", "lhs"], default="nsga2")
    r.add_argument("--pop", type=int, default=80, help="NSGA-II population size")
    r.add_argument("--select", type=int, default=200, help="number of stage 2 candidates")
    r.add_argument("--out", default="results/stage1")
    r.add_argument("--config")
    r.set_defaults(func=cmd_run)

    e = sub.add_parser("eval", help="evaluate one design and print every metric")
    e.add_argument("--design", help="JSON with Design fields (defaults for missing ones)")
    e.add_argument("--config")
    e.add_argument("--stl", help="write the catamaran hull mesh to this STL")
    e.add_argument("--gz", action="store_true", help="print the GZ curve")
    e.set_defaults(func=cmd_eval)

    c = sub.add_parser("dump-config", help="print the default config as JSON")
    c.set_defaults(func=lambda a: print(json.dumps(Config().to_dict(), indent=2)))

    k = sub.add_parser("calibrate", help="mass_calibration from a weighed test print")
    k.add_argument("--design", required=True, help="JSON of the printed test hull")
    k.add_argument("--measured-kg", type=float, required=True)
    k.add_argument("--config")
    k.set_defaults(func=cmd_calibrate)

    a = p.parse_args(argv)
    a.func(a)


if __name__ == "__main__":
    main()
