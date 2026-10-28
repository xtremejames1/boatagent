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
python -m hullopt.parts --design results/stage2/designs/finalist_1.json   # printable parts
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
     `front_thr_gap` from the bow tip. That gap may be negative: the narrow bow cap
     tucks into the empty corner of the angled thruster's bay. Each thruster is
     angled 45° with toe-in: the forward end, which is the inlet, points toward
     the centreline, so the jets run aft-outboard past the hulls. A separate
     faired nose fairing sits ahead of each one (see "Printable parts").
   - Each front thruster hangs from a printed frame merged with the bow cap: an
     upper strut in two segments and a level rail under the thruster, fixed to a
     forefoot fin. The frame's lowest point counts toward the skid pad depth.
   - Rear thrusters point fore-aft and have no pod. They sit wholly in the gap
     between the hulls, as far outboard as clearance allows. Each hangs from one
     diagonal strut that runs outboard and up at 45° from the thruster top into
     its own hull's side (no crossbeam). The strut section is 48 mm chord × 12 mm
     thick, against a minimum of 24 × 12 mm.
   - Mount mass is counted.
   - Each hull has two small printed skid pads, one fore and one aft. Each is a
     trapezoid in side view: a flat bottom `skid_margin` below the lowest thruster,
     pod or frame point, with 45° ends rising to the keel. Positions and flat
     length are optimizer variables.
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
| appendages | rear struts (foil friction, Hoerner thickness factor), front bow members (friction + pressure on their 12 mm face, `member_cd` 1.0), skids, front pods (Hoerner streamlined-body drag + side-window drag `pod_window_cd` 0.2), bare rear thrusters (`pod_cd_rear` 0.3) |
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
cut out, T200 and its two bow members) and `_assembly.stl` (whole boat). GitHub shows STLs in
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

### What stage 2 shows (current run: frame v2, nose fairing)

- **Best weighted resistance at 0.5–1.5 m/s is 16.7 N**, on a 1.22 m hull: 4.8 N
  hull and 11.8 N appendages. At 2 m/s the thrusters have 3.0× the thrust needed.
- **Skid pads limit the design.** The frame's rail runs under the front thruster,
  so it is the boat's lowest point, and the pads must reach below it (about
  11 cm). Most designs that fail in stage 1 fail on skid height. Narrow, deep
  hulls do best (finalists have a 122 mm beam).
- **The bow cap tucks into the thruster bay.** The angled thruster leaves an
  empty corner on the centreline, so the narrow cap tip sits about 3 cm into the
  bay. That keeps the frame short enough to print nose-up.
- **Guessed drag coefficients.** The window, rear-thruster and frame drag all
  rest on guessed coefficients; stage 3 should measure them.

## Printable parts (front frame and nose fairing)

`python -m hullopt.parts --design results/stage2/designs/finalist_1.json --out results/parts`
builds the port-side parts with CadQuery; starboard is the mirror image. Each
part is exported as STEP (edit in CAD), as STL in the boat frame, and as STL
oriented as printed (`*_print.stl`). It also writes a report.

- **front_frame** (one part per side):
  - **Bow cap:** the hull is printed blunt where it is about 24 mm wide. The cap
    slides over a spigot on it and is fixed with 4 × M4 sideways into heat-set
    inserts in the hull. Ahead of the cut, the cap narrows at 25° to the bow tip,
    so it supplies the knife edge.
  - **Forefoot:** a 12 mm fin below the cap, down to the rail.
  - **Upper strut:** two segments, cap → T200 plate → top fairing lug. It runs
    outboard of the 4 M3 thruster screws, so a tool can reach them from above.
  - **Rail:** level, under the thruster, to the bottom fairing lug.

  The thruster comes off by undoing the 4 M3 screws, lowering it slightly and
  sliding it out sideways under the upper strut.
- **front_fairing:** the nose and bumper.
  - It starts at the thruster's inlet.
  - Tall slots open the skin where the inlet keep-out cone passes through it.
  - It closes at ≤ 35° to a small flat tip.
  - Two nylon M3 screws hold it to the lugs. They shear at about 420 N, so a hard
    hit breaks the fairing off, and the frame and thruster survive.
- **hull_bow_spigot:** a reference model of the hull end the cap fits over, with
  the insert holes.
- **Printing:** both parts print nose-up (boat +x up) with no supports.
  - Every member leans at most 35° from the print direction.
  - The cap's cavity ends in a wedge.
  - The T200 plate is a kite shape.
  - Sideways holes are teardrops.
  - The frame's height is limited so it lies on a 250 mm bed.
- **Keep-out cones** (starting values; CFD should refine them):
  - inlet: 45° half-angle, 0.5 D long. Inflow arrives from a wide region.
  - jet: 15° half-angle, 3 D long. A turbulent jet spreads at about 10–12°.
- **Checks on every build:**
  - overhang area over the limit
  - frame or fairing material inside either cone
  - interference with the T200
  - the thruster's removal path
  - tool access to each M3
  - bed fit
- **Strength** (`hullopt/frame.py`, a 3D frame solver): printed PLA, 25 MPa
  allowable, walls only (infill ignored). Load cases:
  - 3× forward and reverse thrust
  - a head-on hit at 1.25× the fairing's breakaway load
  - a side hit at 0.5× breakaway (secondary, reported as a capacity)

  Member sizes (minimum 12 × 24 mm) are picked automatically as the smallest
  that pass. The report also gives the shear on each cap screw.
- **Slicer mass:** PrusaSlicer, headless, with PLA, 4 walls, 10 % gyroid and
  0.2 mm layers.

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
hullopt/parts.py         printable CAD parts (CadQuery), slicer mass, report
hullopt/strength.py      sections, material, simple strength checks
hullopt/frame.py         3D frame solver used for the front mount
tests/                   analytic box checks, mesh checks, end-to-end checks
```
