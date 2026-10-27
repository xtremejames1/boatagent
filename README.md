# boatagent

Multi-fidelity hull optimization for a small holonomic autonomous catamaran
(RoboBoat scale, 4 × T200, 3D-printed hulls). See `CLAUDE.md` for the design brief.

Status: **stages 1 (algebraic screening) and 2 (semi-empirical resistance) are
done.** Stage 3 (OpenFOAM) is next, on the 8 finalists in `results/stage2/`.

## Quick start

```bash
pip install -r requirements.txt
python -m hullopt.stage1 eval --gz                # one design, every metric + GZ curve
python -m hullopt.stage1 run --trials 3000        # optimizer -> results/stage1/
python -m hullopt.stage2 run                      # score candidates -> results/stage2/
python -m pytest -q                               # tests
```

Other commands:

```bash
python -m hullopt.stage1 eval --design results/stage1/designs/knee.json --stl knee.stl
python -m hullopt.stage1 dump-config > my.json    # edit, then pass --config my.json
python -m hullopt.stage1 calibrate --design test_piece.json --measured-kg 0.41
```

## What stage 1 does

For each design (about 0.1 s each):

1. **Geometry.** A parametric demihull: deadrise bottom, bilge radius, vertical
   side, flat deck. Plan view has a power-law forebody set by the entry angle and
   a transom width. The keel has rocker (separate bow and stern rise). The script
   meshes it, then mirrors it to a catamaran with the given inner-face spacing.
2. **Mass.** Printed hull = shell area × wall × filament density + ribs (fraction
   of enclosed volume) + sealer + joint flanges. The print-bed split sets the
   joint count. `mass_calibration` scales the printed part. Set it from a weighed
   test print with the `calibrate` command. Frame mass grows with overall width.
   The payload, LIDAR mast, thrusters, pods, mounts and skids all have positions,
   so the code gets a real KG and LCG.
3. **Thrusters, pods and skids, placed against the design waterline.**
   - Each thruster hangs a set water cover (an optimizer variable) below the
     waterline.
   - Front thrusters sit one ahead of each hull's bow, on the hull centreline,
     `front_thr_gap` clear of the stem. Each is angled 45° inside a faired pod
     aligned with the flow. The pod has an elliptical section, side windows so the
     thruster can draw water sideways, and a solid nose that takes head-on hits
     instead of the thruster. The thrusters toe in: the forward end of each points
     toward the centreline, so the jets run aft-outboard past the hulls. It
     covers the thruster mount and must print in one piece on a 250 mm bed.
   - Rear thrusters point fore-aft and have no pod. They sit wholly in the gap
     between the hulls, as far outboard as clearance allows.
   - Every thruster sits on a 3D-printed mount: a faired strut down from deck
     level plus a clamp. The front mounts add an arm forward from the stem. Mount
     mass is counted.
   - Each hull has two small printed skid pads, one fore and one aft. Each is a
     trapezoid in side view: a flat bottom `skid_margin` below the lowest thruster
     or pod point, with 45° ends rising to the keel. Positions and flat length are
     optimizer variables.
   - The front mount is an arm along the deck from the bow stem to the pod's
     pylon. The rear struts hang from a crossbeam between the hulls.
   - Placement runs twice, so the waterline includes the appendage masses.
4. **Trim and float.**
   - The movable payload (batteries, electronics) slides within its deck range to
     trim the boat level.
   - The code then solves free heave and trim with mesh hydrostatics: exact surface
     integrals on the clipped mesh, waterplane by Green's theorem.
5. **Checks** (constraint ≤ 0 means pass):
   - overall length (hull + front pods ≤ 1.52 m) and overall width
   - front pod: one-piece fit on the 250 mm bed, nose length, window area
   - inner-face spacing
   - GM_T and freeboard
   - trim
   - heel under full sway thrust
   - water cover over the thrusters
   - thruster clearance to the hulls and to each other
   - jet path: how far the front-thruster jet travels before hitting a hull,
     in both directions
   - skid pad height
   - LCG between the fore and aft pads (so the boat rests stable on the ground)
   - reserve buoyancy
6. **Proxies.**
   - ITTC-57 friction with form factor and roughness allowance, weighted over
     0.5–2 m/s. Appendages use their own Reynolds numbers.
   - slenderness L/∇^(1/3) and Froude numbers
   - lateral area, sway crossflow drag, estimated top sway speed
   - lateral-area offset from G (a yaw proxy)
   - GZ curve (`eval --gz`)

The optimizer is NSGA-II (Optuna) on three objectives: minimize total mass,
minimize weighted friction drag, maximize slenderness. Friction alone favours
short hulls and wave drag favours long ones, so slenderness stands in for wave
drag until stage 2 computes it. `--sampler lhs` runs without Optuna.

### Outputs (`results/stage1/`)

| file | content |
|---|---|
| `summary.md` | counts, which constraints bind, representative designs |
| `all.csv` | every evaluation (not committed; regenerate) |
| `feasible.csv`, `pareto.csv` | feasible designs (not committed), and the non-dominated ones |
| `stage2_candidates.csv` | up to `--select` (200) diverse feasible designs, by Pareto rank then crowding |
| `designs/*.json`, `*.stl` | min-mass, min-drag, max-slenderness and knee designs |
| `config.json` | the config used |

## What stage 2 does

`python -m hullopt.stage2 run` reads `results/stage1/stage2_candidates.csv`. For
each design it reruns stage 1 to get sinkage and trim, then computes, at each
speed:

| term | model |
|---|---|
| wave | Michell thin-ship integral with demihull interference, `2[1 + cos(k0 λ√(λ²−1) s)]`. Exact per-panel integration of the hull offsets at the floating waterline (`hullopt/michell.py`). Scaled by `wave_factor` (1.0 until stage 3 calibrates it). |
| viscous | ITTC-57 × catamaran form factor `1+βk = 3.03 (L/∇^⅓)^−0.4` (Couser et al. 1997) + roughness ΔCf |
| appendages | struts and pylons (foil friction, Hoerner thickness factor), skids, front pods (Hoerner streamlined-body drag + side-window drag `pod_window_cd` 0.2), bare rear thrusters (`pod_cd_rear` 0.3) |
| transom | Holtrop (1984) wet-transom term |

It also computes:

- **Sway:** crossflow drag of the lateral profile. Cd runs from 1.2 for a hard
  chine to 0.7 for a full bilge radius, and the leeward hull takes 0.6 of the
  load. Appendages are added. The output is a top sway speed on the front
  thrusters' lateral thrust.
- **Yaw:** the top turn-in-place rate. Full thrust (forward on one side, reverse on
  the other) is balanced against crossflow yaw damping `∫ Cd h(x) |x−x_G|³ dx`
  and hull friction. Deep lateral area far aft costs yaw rate.
- **Top speed:** the margin between full surge thrust and R_T at 2 m/s.

The score is weighted R_T over 0.5–1.5 m/s (the 2 m/s weight is dropped). The
code finds the Pareto front on (weighted R_T, total mass, sway drag) and picks
8 finalists by rank and crowding. Outputs: `all.csv`, `pareto.csv`,
`finalists.csv`, `curves.csv` (resistance components from 0.3 to 2.5 m/s),
`resistance.png`, and `designs/finalist_*.json` plus STLs. Each finalist gets
`_catamaran.stl` (hulls only, for CFD), `_pod.stl` (front pod with the windows
cut out, T200 and pylon) and `_assembly.stl` (whole boat). GitHub shows STLs in
3D. `python -m hullopt.export --design d.json --out name` writes them for any
design.

`python -m hullopt.stage2 wigley` prints Michell Cw for the Wigley hull, mono and
catamaran. Use it as the reference when validating the stage 3 CFD setup.
`python -m hullopt.stage2 eval --design d.json` prints one design's breakdown.

### Validation in the tests

- The integral matches a closed-form answer for a parabolic strut to better
  than 0.2 %.
- R scales exactly as beam² (thin-ship) and as L³ at a fixed Froude number.
- At zero spacing the catamaran equals one hull of double beam; at wide spacing
  it tends to two separate hulls.
- Results don't depend on the integration cutoff.
- The Wigley hull gives the known wetted area (0.1488 L²) and regression Cw
  values: 2.14e-3 at Fr 0.3 and 4.52e-3 at Fr 0.5. Michell theory overpredicts
  the measured hump at Fr 0.5 by about 2×, which is a known limitation.

### What stage 2 shows (toe-in pods, fore/aft skid pads)

- **Best weighted resistance at 0.5–1.5 m/s is 13.9 N.** Of that, 4.5 N is the
  hull and 9.4 N the appendages. At 2 m/s the thrusters have 3.4× the thrust
  needed. The 1.26 m hulls (slenderness about 5.4) have the least drag. The
  1.0–1.2 m hulls are lighter and turn faster but carry 7–15 % more drag.
- **Toe-in helps turning and clears the jets.** The front jets now run
  aft-outboard and clear every hull in both directions. The top turn-in-place
  rate rose from about 180 to 240 °/s, because both front thrusters' lines of
  action now sit further from G.
- **Skid pads are small.** Finalist 1 has a 3 cm flat with 45° ends, pads at 35 %
  and 60 % of hull length, 7.9 cm tall at most, and 0.46 kg for all four pads.
- **Thruster height is at its limit.** The optimizer already puts the thrusters
  as high as `thr_cover_min` allows, because deeper thrusters need taller, heavier
  pads. Each centimetre less cover makes the pads about 1 cm shorter. On
  finalist 1: 5 cm cover gives 7.2 cm pads, and 3 cm cover gives 5.1 cm pads.
  Check ventilation (stage 3) before lowering it.
- **What's left of the appendage drag**, at 1 m/s:

  | source | drag |
  |---|---|
  | pod side windows | ≈3 N |
  | bare rear thrusters | ≈3 N |
  | pod bodies | ≈1.8 N |

  The window and rear-thruster numbers rest on guessed drag coefficients, so stage
  3 should measure them.
- **The rear jets still graze the front pods** (`jet_rear_to_front_thr` ≈ −2 cm).
  This is reported but doesn't fail a design.

## Coordinates and units

SI throughout. Body frame: x forward from the transom, y to port, z up, z = 0
at the lowest keel point. Front thruster axes point forward-inboard at 45° (toe-in).
`loa` is the hull length only; the front pods sit ahead of it.

## Placeholders to confirm

These values live in `hullopt/params.py` (`Config`). Change them there or in a
JSON passed with `--config`:

| item | default | note |
|---|---|---|
| `gm_t_min` (X in the brief) | 0.30 m | rarely binding: cats here have GM_T ≈ 0.5–0.9 m |
| `freeboard_min` (Y in the brief) | 0.08 m | |
| `width_max` | 0.914 m (3 ft) | check against current RoboBoat rules |
| `skid_max_height` | 0.12 m | **the main design driver**, see below |
| `thr_cover_min` | 0.05 m | water above the rear thruster / front pod top; verify ventilation in CFD |
| `pod_window_ratio`, `pod_window_cd` | 1.0, 0.2 | side-window area per side ÷ thruster disk area; window drag |
| `pod_crush_len_min`, `impact_speed` | 3 cm, 2 m/s | load-taking nose length; impact energy estimate |
| `lidar_z_above_deck`, `lidar_kg` | 0.60 m, 1.0 kg | sets KG |
| `frame_*`, `payload_z_above_deck` | | mass and KG of the frame and electronics |
| `wall_thickness`, `rib_vol_frac`, `mass_calibration` | 1.6 mm, 1.2 %, 1.0 | calibrate with a test print |
| `bed` | 256 mm cube | printer build volume |
| `roughness_dcf`, `form_factor_k` | 4e-4, 0.10 | friction proxy |

## What the first run shows

A T200 in a pod with 5 cm of water above it reaches about 17 cm below the
waterline. A 25–30 kg catamaran with slender hulls draws only about 8–11 cm.
So the skid must be about 8–12 cm tall. With an 8 cm cap almost no slender design
was feasible. With 12 cm, many are. The skid's mass, wetted area and lateral area
count against every design, so the optimizer still keeps the skid small.
Shallower thruster cover, a slimmer pod, or a deeper hull all relax this.

## Layout

```
hullopt/params.py        design variables, bounds, fixed config
hullopt/geometry.py      demihull/catamaran mesh, half-breadth lookup, STL writer
hullopt/hydrostatics.py  clipped-mesh hydrostatics, equilibrium, GZ
hullopt/mass.py          printed mass, print-bed segmentation, mass budget
hullopt/layout.py        thrusters, front pods, mounts, skids, clearance and jet checks
hullopt/stage1.py        evaluation, NSGA-II driver, candidate selection, CLI
hullopt/michell.py       Michell wave resistance (mono / catamaran), Wigley offsets
hullopt/stage2.py        resistance build-up, sway, yaw, finalist selection, CLI
hullopt/export.py        viewing STLs: front pod and full assembly
tests/                   analytic box checks, mesh checks, end-to-end checks
```
