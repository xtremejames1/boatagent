"""Stage 2: semi-empirical resistance, sway and yaw for the stage 1 candidates.

Per design and speed:
  R_T = R_wave + R_visc + R_app + R_transom
  R_wave    Michell thin-ship integral with demihull interference (hullopt.michell),
            hull at its stage 1 sinkage and trim. Times `wave_factor`.
  R_visc    ITTC-57 friction x catamaran form factor (1 + beta k) + roughness allowance.
            Form factor from Couser et al. (1997): 1 + beta k = 3.03 (L / Vol^(1/3))^-0.4.
  R_app     rear struts (foil friction, Hoerner thickness factor), front bow
            members (friction + pressure on their 12 mm face), skids,
            front pods (Hoerner streamlined-body drag + open side windows), and
            bare rear thrusters (pressure drag).
  R_transom Holtrop (1984) wet-transom pressure term.

Sway: crossflow drag of the hull lateral profile (Cd from bilge radius / draft,
leeward hull partly shielded) plus appendages. Gives top sway speed on the
front thrusters' lateral thrust.

Yaw: top turn-in-place rate. Thruster yaw torque (each thruster at full forward or
reverse thrust) against crossflow yaw damping, int Cd h(x) |x - xG|^3 dx, plus
hull friction from the hulls moving fore-aft at r * yc.

Ranking: Pareto on (weighted R_T over 0.5-1.5 m/s, total mass, sway drag), then
finalists for stage 3 by Pareto rank and crowding.

Usage:
  python -m hullopt.stage2 run [--candidates results/stage1/stage2_candidates.csv] [--out results/stage2]
  python -m hullopt.stage2 eval [--design d.json]
  python -m hullopt.stage2 wigley
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
import time
from dataclasses import asdict, dataclass, fields, replace

import numpy as np

from .export import export as export_stl
from .geometry import catamaran_mesh, hull_centre_y, keel_z, offsets_grid, write_stl
from .michell import michell, wigley_offsets, wigley_wetted_area
from .params import Config, Design
from .stage1 import crowding, evaluate, nondominated_fronts, write_csv


@dataclass
class Stage2Config:
    nx: int = 101                    # stations for the Michell offsets
    nz: int = 24                     # depth cells
    lam_max: float = 10.0
    wave_factor: float = 1.0         # scale on Michell R_wave (calibrate from stage 3)
    form_factor: str = "couser"      # "couser" or a number for (1 + beta k)
    score_max_speed: float = 1.5     # weighted score uses speeds up to this
    top_speed: float = 2.0           # report thrust margin here
    # appendages
    skid_form_k: float = 0.3
    pod_cd_rear: float = 0.3         # aligned, flow-through duct
    # sway / yaw crossflow
    cd_sharp: float = 1.2            # hard chine section
    cd_round: float = 0.7            # full bilge radius
    lee_factor: float = 0.6          # leeward hull loading relative to windward
    cd_appendage: float = 1.2
    # selection
    n_finalists: int = 8
    sweep: tuple = tuple(np.round(np.arange(0.3, 2.51, 0.05), 3))

    def to_dict(self):
        return asdict(self)

    @classmethod
    def from_dict(cls, d):
        names = {f.name for f in fields(cls)}
        return replace(cls(), **{k: (tuple(v) if isinstance(v, list) else v) for k, v in d.items() if k in names})


def ittc57(re):
    return 0.075 / (np.log10(np.maximum(re, 1e3)) - 2.0) ** 2


def form_factor(c2: Stage2Config, slenderness: float) -> float:
    if c2.form_factor == "couser":
        return 3.03 * slenderness ** -0.40
    return float(c2.form_factor)


@dataclass
class Prepared:
    """Everything stage 2 needs from one design, computed once."""

    d: Design
    row1: dict
    x: np.ndarray
    depth_edges: np.ndarray
    F: np.ndarray
    spacing: float
    t_local: np.ndarray      # local draft at each station (centreline profile)
    s_hull: float            # wetted area, both hulls
    lwl: float
    k_form: float
    transom_area: float      # per hull, immersed at rest
    app_fric: list           # (wetted area, reference length, form factor)
    pod_area_cd: float       # sum Cd * frontal area of the bare rear thrusters
    pods: list               # front pods (drag from Pod.drag_area)
    app_lat_cd: float        # sum Cd * lateral area of appendages
    skid_profile: tuple      # (xs, heights)
    xg: float
    yaw_torque: float
    f_sway: float
    f_surge: float


def prepare(d: Design, cfg: Config, c2: Stage2Config) -> Prepared | None:
    r1, det = evaluate(d, cfg, detail=True)
    if det is None:
        return None
    pose = det.pose
    L = d.loa
    xs = np.linspace(0.0, L, c2.nx)
    z_wl = (pose.heave + math.sin(pose.trim) * (xs - L / 2)) / math.cos(pose.trim)
    t_loc = np.clip(z_wl - keel_z(xs / L, d), 0.0, None)
    t_max = max(float(t_loc.max()), 1e-4)
    de = np.linspace(0.0, t_max, c2.nz + 1)
    dm = 0.5 * (de[1:] + de[:-1])
    F = offsets_grid(xs, z_wl[:, None] - dm[None, :], d)
    transom = 2.0 * float(np.sum(F[0] * np.diff(de)))

    app = []
    pod = 0.0
    lat = 0.0
    for s in det.struts:
        sub = s.submerged_length(det.z_wl)
        if sub <= 0:
            continue
        if s.aligned:   # foil facing the flow: friction with Hoerner thickness factor
            app.append((2.04 * s.chord * sub, s.chord, 1 + 2 * s.tc + 60 * s.tc**4))
        else:           # bow member: friction along it plus pressure on its 12 mm face
            app.append((2.04 * s.chord * sub, sub, 1.0))
            pod += cfg.member_cd * s.frontal_area(det.z_wl)
        lat += c2.cd_appendage * s.side_area(det.z_wl)
    for t in det.thrusters:
        if t.name.startswith("front"):
            continue  # inside its pod; the pod carries the drag
        r, ln = t.radius, t.length
        app.append((2 * (2 * math.pi * r * ln), ln, 1.0))
        pod += c2.pod_cd_rear * math.pi * r * r
        lat += c2.cd_appendage * 2 * r * ln
    for pd in det.pods:
        lat += c2.cd_appendage * pd.length * 2 * pd.half_h
    sk = det.skid
    app += [(sk.wetted_area, sk.length, 1 + c2.skid_form_k)] * 2
    lat += c2.cd_appendage * 2 * float(np.trapezoid(sk.heights, sk.xs))

    # yaw torque: every thruster at full thrust, one of each pair forward, one reverse
    G = det.budget.cg
    fp = next(t for t in det.thrusters if t.name == "front_port")
    rp = next(t for t in det.thrusters if t.name == "rear_port")
    # z-moment of a unit force along the port front axis, about G
    tau_front = abs((fp.centre[0] - G[0]) * math.sin(fp.heading) - (fp.centre[1] - G[1]) * math.cos(fp.heading))
    pair = (cfg.thr_force_kgf + cfg.thr_rev_force_kgf) * cfg.g
    yaw_torque = pair * (abs(rp.centre[1]) + tau_front)

    return Prepared(
        d=d, row1=r1, x=xs, depth_edges=de, F=F, spacing=2 * hull_centre_y(d), t_local=t_loc,
        s_hull=r1["wetted_hulls"], lwl=r1["lwl"],
        k_form=form_factor(c2, r1["slenderness"]), transom_area=transom,
        app_fric=app, pod_area_cd=pod, pods=det.pods, app_lat_cd=lat, skid_profile=(sk.xs, sk.heights),
        xg=float(G[0]), yaw_torque=yaw_torque,
        f_sway=r1["f_sway"], f_surge=r1["f_surge"],
    )


def resistance(p: Prepared, U: float, cfg: Config, c2: Stage2Config, mono: bool = False) -> dict:
    rho, g, nu = cfg.rho, cfg.g, cfg.nu
    q = 0.5 * rho * U * U
    rw = c2.wave_factor * michell(p.x, p.depth_edges, p.F, U, spacing=p.spacing, rho=rho, g=g,
                                  lam_max=c2.lam_max)
    rv = q * p.s_hull * (p.k_form * ittc57(U * p.lwl / nu) + cfg.roughness_dcf)
    ra = sum(q * a * (k * ittc57(U * lr / nu) + cfg.roughness_dcf) for a, lr, k in p.app_fric)
    ra += q * p.pod_area_cd
    ra += sum(q * pd.drag_area(ittc57(U * pd.length / nu) + cfg.roughness_dcf, cfg) for pd in p.pods)
    rt_ = 0.0
    if p.transom_area > 1e-6:
        bwl, cwp = p.row1["bwl"], p.row1["cwp"]
        fnt = U / math.sqrt(2 * g * p.transom_area / (bwl * (1 + cwp)))
        c6 = 0.2 * (1 - 0.2 * fnt) if fnt < 5 else 0.0
        rt_ = 2 * q * p.transom_area * c6
    out = {"U": U, "Fr": U / math.sqrt(g * p.lwl), "R_wave": rw, "R_visc": rv, "R_app": ra,
           "R_transom": rt_, "R_total": rw + rv + ra + rt_}
    if mono:
        r1 = michell(p.x, p.depth_edges, p.F, U, spacing=None, rho=rho, g=g, lam_max=c2.lam_max)
        out["interference"] = rw / (2 * c2.wave_factor * r1) if r1 > 0 else float("nan")
    return out


def sway_coeff(p: Prepared, c2: Stage2Config) -> float:
    """C such that sway drag = 0.5 rho V^2 C (m^2)."""
    d = p.d
    a_lat = float(np.trapezoid(p.t_local, p.x))           # one hull
    i = int(np.argmax(p.t_local))
    ratio = min(d.bilge_frac * d.beam / 2 / max(p.t_local[i], 1e-6), 1.0)
    cd = c2.cd_sharp - (c2.cd_sharp - c2.cd_round) * ratio
    return cd * a_lat * (1 + c2.lee_factor) + p.app_lat_cd, cd


def yaw_rate_max(p: Prepared, cfg: Config, c2: Stage2Config, cd_hull: float) -> tuple[float, float]:
    """Top turn-in-place rate (rad/s) and the share of lateral profile aft of G."""
    rho, nu = cfg.rho, cfg.nu
    xs, h = p.x, p.t_local.copy()
    sx, sh = p.skid_profile
    h += np.interp(xs, sx, sh, left=0.0, right=0.0)
    arm = np.abs(xs - p.xg)
    i3 = float(np.trapezoid(h * arm**3, xs))                 # one hull
    yc = p.spacing / 2

    def damping(r):
        m_cf = 0.5 * rho * cd_hull * r * r * i3 * (1 + c2.lee_factor)
        u = r * yc
        q = 0.5 * rho * u * u
        m_f = 2 * yc * q * (p.s_hull / 2) * (p.k_form * ittc57(max(u, 1e-3) * p.lwl / nu) + cfg.roughness_dcf)
        return m_cf + m_f

    lo, hi = 0.0, 20.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if damping(mid) < p.yaw_torque:
            lo = mid
        else:
            hi = mid
    aft = float(np.trapezoid(np.where(xs < p.xg, h, 0.0), xs) / np.trapezoid(h, xs))
    return lo, aft


def evaluate2(d: Design, cfg: Config, c2: Stage2Config) -> dict:
    p = prepare(d, cfg, c2)
    if p is None or not p.row1["feasible"]:
        row = (p.row1 if p else evaluate(d, cfg)).copy()
        row["stage2_ok"] = False
        return row
    row = dict(p.row1)
    speeds = np.asarray(cfg.surge_speeds, float)
    w = np.asarray(cfg.surge_weights, float) * (speeds <= c2.score_max_speed + 1e-9)
    res = [resistance(p, U, cfg, c2) for U in speeds]
    rt = np.array([r["R_total"] for r in res])
    for r in res:
        tag = f"{r['U']:g}"
        for k in ("R_wave", "R_visc", "R_app", "R_transom", "R_total"):
            row[f"{k}_{tag}"] = r[k]
    top = resistance(p, c2.top_speed, cfg, c2)
    c_sway, cd_h = sway_coeff(p, c2)
    q = 0.5 * cfg.rho
    r_yaw, aft = yaw_rate_max(p, cfg, c2, cd_h)
    row.update({
        "stage2_ok": True,
        "R_T_weighted": float(np.dot(w, rt) / w.sum()),
        "R_wave_weighted": float(np.dot(w, [r["R_wave"] for r in res]) / w.sum()),
        "R_hull_weighted": float(np.dot(w, [r["R_total"] - r["R_app"] for r in res]) / w.sum()),
        "R_app_weighted": float(np.dot(w, [r["R_app"] for r in res]) / w.sum()),
        "form_factor": p.k_form,
        "transom_area": p.transom_area,
        "R_T_top": top["R_total"],
        "thrust_margin_top": p.f_surge / top["R_total"],
        "sway_cd_hull": cd_h,
        "sway_drag_s2": q * cfg.sway_speed**2 * c_sway,
        "v_sway_max_s2": math.sqrt(p.f_sway / (q * c_sway)),
        "yaw_torque": p.yaw_torque,
        "yaw_rate_max_dps": math.degrees(r_yaw),
        "lat_profile_aft_frac": aft,
    })
    return row


OBJ2 = ["R_T_weighted", "mass_total", "sway_drag_s2"]


def rank(rows: list[dict], n_final: int):
    ok = [r for r in rows if r.get("stage2_ok")]
    if not ok:
        return [], []
    F = np.array([[r[k] for k in OBJ2] for r in ok], float)
    fronts = nondominated_fronts(F)
    for rnk, f in enumerate(fronts):
        for i in f:
            ok[i]["s2_rank"] = rnk
    pareto = sorted((ok[i] for i in fronts[0]), key=lambda r: r["R_T_weighted"])
    chosen: list[int] = []
    for f in fronts:
        if len(chosen) + len(f) <= n_final:
            chosen += list(f)
        else:
            cd = crowding(F[f])
            chosen += list(f[np.argsort(-cd)][: n_final - len(chosen)])
            break
    finals = sorted((ok[i] for i in chosen), key=lambda r: r["R_T_weighted"])
    return pareto, finals


def curves(d: Design, cfg: Config, c2: Stage2Config, name: str) -> list[dict]:
    p = prepare(d, cfg, c2)
    out = []
    for U in c2.sweep:
        r = resistance(p, float(U), cfg, c2, mono=True)
        r["design"] = name
        out.append(r)
    return out


def plot_curves(rows: list[dict], path: str) -> bool:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return False
    names = list(dict.fromkeys(r["design"] for r in rows))
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for n in names:
        rr = [r for r in rows if r["design"] == n]
        U = [r["U"] for r in rr]
        ax[0].plot(U, [r["R_total"] for r in rr], label=n)
        ax[1].plot(U, [r["R_wave"] for r in rr], label=n)
    for a, t in zip(ax, ["Total resistance (hull + appendages)", "Wave resistance (Michell, with interference)"]):
        a.axvspan(0.5, 1.5, color="0.92", zorder=0)
        a.set_xlabel("speed (m/s)")
        a.set_ylabel("N")
        a.set_title(t)
        a.grid(alpha=0.3)
    ax[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return True


SUMMARY = ["loa", "beam", "depth", "spacing_inner", "mass_total", "slenderness",
           "R_T_weighted", "R_hull_weighted", "R_wave_weighted", "R_app_weighted",
           "R_T_top", "thrust_margin_top", "sway_drag_s2", "v_sway_max_s2",
           "yaw_rate_max_dps", "gm_t", "skid_max_height"]


def read_candidates(path: str) -> list[Design]:
    with open(path) as f:
        return [Design.from_dict({k: v for k, v in r.items() if v != ""}) for r in csv.DictReader(f)]


def cmd_run(a):
    cfg = Config.load(a.config)
    c2 = Stage2Config()
    if a.stage2_config:
        with open(a.stage2_config) as f:
            c2 = Stage2Config.from_dict(json.load(f))
    c2 = replace(c2, n_finalists=a.finalists)
    designs = read_candidates(a.candidates)
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    rows = []
    for i, d in enumerate(designs):
        rows.append(evaluate2(d, cfg, c2))
        rows[-1]["candidate"] = i
        if (i + 1) % 25 == 0:
            print(f"  {i + 1}/{len(designs)} designs, {time.time() - t0:.0f} s", file=sys.stderr)
    pareto, finals = rank(rows, c2.n_finalists)
    write_csv(os.path.join(a.out, "all.csv"), rows)
    write_csv(os.path.join(a.out, "pareto.csv"), pareto)
    write_csv(os.path.join(a.out, "finalists.csv"), finals)
    with open(os.path.join(a.out, "stage2_config.json"), "w") as f:
        json.dump(c2.to_dict(), f, indent=2, default=float)
    dd = os.path.join(a.out, "designs")
    os.makedirs(dd, exist_ok=True)
    curve_rows = []
    for k, r in enumerate(finals):
        name = f"finalist_{k + 1}"
        r["finalist"] = name
        d = Design.from_dict(r)
        with open(os.path.join(dd, f"{name}.json"), "w") as f:
            json.dump(d.to_dict(), f, indent=2)
        write_stl(os.path.join(dd, f"{name}_catamaran.stl"), catamaran_mesh(d, cfg.n_stations).triangles, name)
        export_stl(d, cfg, os.path.join(dd, name))   # _pod.stl and _assembly.stl for viewing
        curve_rows += curves(d, cfg, c2, name)
    write_csv(os.path.join(a.out, "finalists.csv"), finals)
    write_csv(os.path.join(a.out, "curves.csv"), curve_rows)
    plotted = plot_curves(curve_rows, os.path.join(a.out, "resistance.png"))
    dt = time.time() - t0
    lines = [
        "# Stage 2 summary", "",
        f"- designs: {len(rows)} ({sum(bool(r.get('stage2_ok')) for r in rows)} evaluated) in {dt:.0f} s",
        f"- Pareto front on (weighted R_T 0.5-{c2.score_max_speed:g} m/s, mass, sway drag): {len(pareto)}",
        f"- finalists for stage 3: {len(finals)}", "",
        "## Finalists", "",
        "| name | " + " | ".join(SUMMARY) + " |",
        "|---" * (len(SUMMARY) + 1) + "|",
    ]
    for r in finals:
        lines.append(f"| {r['finalist']} | " + " | ".join(f"{r[k]:.3g}" for k in SUMMARY) + " |")
    if plotted:
        lines += ["", "![resistance curves](resistance.png)"]
    text = "\n".join(lines) + "\n"
    with open(os.path.join(a.out, "summary.md"), "w") as f:
        f.write(text)
    print(text)
    print(f"results in {a.out}/")


def cmd_eval(a):
    cfg = Config.load(a.config)
    d = Design()
    if a.design:
        with open(a.design) as f:
            d = Design.from_dict(json.load(f))
    c2 = Stage2Config()
    row = evaluate2(d, cfg, c2)
    keys = [k for k in row if k.startswith(("R_", "stage2", "form_", "transom", "thrust_margin", "sway_",
                                            "v_sway", "yaw", "lat_profile", "mass_total", "feasible",
                                            "fail_reason"))]
    for k in keys:
        v = row[k]
        print(f"{k:24s} {v:.5g}" if isinstance(v, float) else f"{k:24s} {v}")
    if row.get("stage2_ok"):
        print("\n   U     Fr   R_wave  R_visc   R_app  R_tran  R_total  interf")
        for r in curves(d, cfg, c2, "eval")[::4]:
            print(f"{r['U']:5.2f} {r['Fr']:5.2f} {r['R_wave']:7.2f} {r['R_visc']:7.2f} {r['R_app']:7.2f} "
                  f"{r['R_transom']:7.2f} {r['R_total']:8.2f}  {r['interference']:6.2f}")


def cmd_wigley(a):
    """Michell Cw of the standard Wigley hull (L/B=10, L/T=16), mono and catamaran."""
    L = 1.0
    x, de, F = wigley_offsets(L=L, nx=101, nz=30)
    S = wigley_wetted_area(L=L)
    print("  Fr   Cw_mono   Cw_cat s/L=0.2 (per demihull)   s/L=0.3")
    for Fr in np.arange(0.2, 0.81, 0.05):
        U = Fr * math.sqrt(9.81 * L)
        q = 0.5 * 1000 * U * U
        m = michell(x, de, F, U) / (q * S)
        c2_ = michell(x, de, F, U, spacing=0.2 * L) / (2 * q * S)
        c3 = michell(x, de, F, U, spacing=0.3 * L) / (2 * q * S)
        print(f"{Fr:5.2f}  {m:.3e}   {c2_:.3e}                     {c3:.3e}")


def main(argv=None):
    p = argparse.ArgumentParser(prog="hullopt.stage2", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="score stage 1 candidates and pick stage 3 finalists")
    r.add_argument("--candidates", default="results/stage1/stage2_candidates.csv")
    r.add_argument("--out", default="results/stage2")
    r.add_argument("--finalists", type=int, default=8)
    r.add_argument("--config")
    r.add_argument("--stage2-config")
    r.set_defaults(func=cmd_run)
    e = sub.add_parser("eval", help="stage 2 metrics and resistance curve for one design")
    e.add_argument("--design")
    e.add_argument("--config")
    e.set_defaults(func=cmd_eval)
    w = sub.add_parser("wigley", help="Michell Cw for the Wigley hull (validation reference)")
    w.set_defaults(func=cmd_wigley)
    a = p.parse_args(argv)
    a.func(a)


if __name__ == "__main__":
    main()
