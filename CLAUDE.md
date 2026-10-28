# Catamaran Hull Optimization — Project Brief

Context for Claude Code. This project builds a multi-fidelity hull optimization
pipeline for a small holonomic autonomous catamaran (RoboBoat-scale).

## Design decisions (settled)

- **Configuration:** catamaran, holonomic.
- **Thrusters:** 4 × Blue Robotics T200 (assume for now), all on 3D-printed mounts.
  - Front pair: one in front of each hull's bow, angled ~45° toe-in (forward end
    toward the centreline; the forward end is the inlet, jets aft-outboard).
  - Front mount = one printed frame per side, merged with the bow: a cap that
    slides over a spigot on the blunt-cut hull (M4 into heat-set inserts in the
    hull; the cap carries the knife edge), a forefoot fin, an upper strut
    (cap -> T200 plate -> fairing lug) and a level rail under the thruster.
    Members >= 12 x 24 mm. Impacts must go into the struts, not the thruster.
    Thruster attaches with its 4 M3 only (M6 unused), cable exits the top and is
    routed by the mount, and it must come off without removing the fairing.
  - Separate nose fairing/bumper ahead of the thruster, its aft edge at the inlet,
    open where the inlet needs free space, takes head-on hits and breaks off
    (nylon M3); hollow, printed normally. Side hits matter less.
  - Printing: PLA prototype, min 4 walls and 10 % infill, no supports, 35° max
    overhang, every part fits a 250 x 250 x 250 mm bed.
  - Keep a free cone at the inlet and the outlet (defaults 45° x 0.5 D inlet,
    15° x 3 D jet; refine with CFD).
  - Rear pair: between the hulls, no pod. Only the front pair has pods. Each hangs
    from one diagonal strut running outboard-up (45°) straight into its hull's side
    (no crossbeam); strut section at least 12 mm thick and 24 mm wide.
  - Overall length ≤ 1.52 m includes the front pods (confirmed; hull `loa` is
    shorter).
- **Grounding:** the boat must rest on the hulls, not the thrusters. Achieved with a
  **separate, replaceable 3D-printed skid** under each hull: two small trapezoid
  pads (fore and aft), flat bottoms below the lowest thruster point. The hull draft is therefore a free variable,
  not set by thruster clearance.
- **Hull spacing:** 10–14 in measured between hull **inner faces**; wider is allowed.
- **Length:** overall length ≤ 1.52 m (5 ft). Optimizer range ~1.2–1.52 m.
  Compute waterline length from geometry; do not assume LWL = LOA.
- **Payload:** ~35 lb (~16 kg).
- **Fabrication:** 3D printed, sprayed with sealer.
- **Speed objective:** weighted multi-speed. Most weight on 0.5–1.5 m/s
  (surge). Check 2 m/s as top speed. Include sway at a lower speed.

## Objectives and constraints

- Maximize thrust-to-mass → in practice: minimize hull mass and weighted resistance
  with thrusters fixed.
- Forward (surge) drag is the priority; sway drag also matters (sway thrust is only
  the front pair's lateral components, ~1/3 of surge thrust).
- Must be "pretty stable": check transverse GM with realistic KG (LIDAR mast is
  high) and a heeling moment from sway thrust.
- Must be reasonably maneuverable (yaw). Rocker helps; avoid deep aft lateral area.

## Key physics notes (estimates from planning, to be verified)

- At 1.35–1.52 m, 0.5–1.5 m/s is Fr ≈ 0.13–0.41, below the main wave-resistance
  hump. Shorter hulls move the priority band onto the hump.
- Slenderness L/∇^(1/3) per hull ≈ 5.9–6.6 at 1.52 m for 25–35 kg total.
  Wave drag drops sharply above ~6. Mass is as important as length.
- Demihull wave interference is significant at this Fr; spacing is a real
  optimization variable (no longer stability-limited — rough BM_T ≈ 0.5 m).
- Front thrusters near the bow risk ventilation and jet impingement on the
  opposite hull / rear thrusters. Check in CFD.

## Parameters

| Group | Parameters |
|---|---|
| Demihull | length, beam, draft, section shape (bilge radius/deadrise), rocker, bow entry angle, transom width |
| Layout | hull spacing (inner-face, ≥10 in), fore-aft thruster positions |
| Skid | two trapezoid pads per hull (fore/aft): depth below lowest thruster (small margin, ~1–2 cm), flat length, fore-aft positions. |
| Front pod | depth (cover), wall, fairing straight length (`pod_nose_len`), bow gap (may be negative: cap tip tucks into the bay). Frame member sizes picked by the strength check (min 12 x 24 mm). |

Hard constraints: skid bottom below lowest thruster, displacement ≥ mass, GM_T ≥ X,
freeboard ≥ Y, LOA ≤ 1.52 m, overall width limit, thrusters submerged enough to
avoid ventilation.

## Pipeline (multi-fidelity)

1. **Stage 1 — algebraic (thousands of evals):** parametric demihull generator,
   hydrostatics from mesh (displacement, waterline, LWL, GM_T, GM_L), clearance
   checks, printed-mass model (shell area × wall thickness × filament density +
   ribs/infill; calibrate against a test print), print-bed segmentation penalty.
2. **Stage 2 — semi-empirical (hundreds):** Michell thin-ship wave resistance with
   demihull interference + ITTC-57 friction + roughness allowance (printed +
   sealed surface, Re ~3e6); crossflow-drag estimate for sway; empirical yaw proxy.
   Score on the weighted 0.5–1.5 m/s objective.
3. **Stage 3 — OpenFOAM interFoam (5–10 finalists):** free sinkage/trim, surge and
   sway cases, actuator-disk thrusters for interaction and ventilation, front
   pod window cases. Validate setup first on a Wigley catamaran benchmark; do a mesh
   independence study.

Optimizer: Optuna or pymoo (NSGA-II for multi-objective).

## Code

- `hullopt/` package; stage 1 is `python -m hullopt.stage1` (see README.md).
- Body frame: x forward from the transom, y to port, z up, z=0 at lowest keel. SI units.
- Fixed inputs and placeholders (X, Y, width limit, skid height cap, KG items) are in
  `hullopt/params.py` `Config`; optimizer variables and bounds are `Design` / `BOUNDS`.
- Run `python -m pytest -q` after changes.

## Status

- Stage 1 done: geometry, mesh hydrostatics + free trim, mass model, print
  segmentation, clearance/jet/skid checks, friction proxy, NSGA-II (Optuna).
  Output `results/stage1/stage2_candidates.csv`.
- Stage 2 done (`hullopt/michell.py`, `hullopt/stage2.py`): Michell + interference,
  ITTC-57 with Couser form factor, appendage and transom drag, sway crossflow,
  yaw-rate proxy. Output `results/stage2/finalists.csv` (8 designs, STLs in
  `results/stage2/designs/`).
- Findings (current run): best weighted R_T 16.7 N (hull 4.8, appendages 11.8; all
  appendage Cd are placeholders), 1.22 m hull. Skid height binds (the frame rail under
  the front thruster is the lowest point). Front parts for finalist 1 in results/parts:
  frame 168 g, fairing 89 g (slicer), all checks clean except small ledges.
- Mount layout (front ahead of bows, rear between hulls) is set; mount details will
  still change. Keep thruster/pod/mount geometry parametric (`hullopt/layout.py`).

- Front pod model (`Pod` in `hullopt/layout.py`): elliptical section aligned with +x,
  sized around the 45° T200 (+5 mm gap + wall); tail fairing, thruster bay with side
  windows (area per side = `pod_window_ratio` x thruster disk), load-taking nose
  (`pod_nose_len`, infill `pod_crush_fill`). Flooded. Drag: Hoerner streamlined-body
  formula + window drag (`pod_window_cd`, placeholder). Constraints: one-piece fit on
  `pod_bed`, nose ≥ `pod_crush_len_min`, windows ≤ `pod_window_max_frac` of bay side.
  Reports impact energy at `impact_speed` and average force over the nose length.
- Thruster height: optimizer already drives cover to `thr_cover_min` (5 cm,
  placeholder) because deeper thrusters need taller, heavier skid pads. Lowering
  that limit is the only way higher; check ventilation in CFD first.
- Viewing STLs: `python -m hullopt.export`; stage 2 writes `_pod.stl` and
  `_assembly.stl` for each finalist.
- Printable parts: `python -m hullopt.parts` (CadQuery; STEP + STL + print-oriented STL;
  PrusaSlicer mass; checks; frame strength via `hullopt/frame.py`). Hull interface: bow
  cap with 4x M4 sideways (staggered, 2 per side) into heat-set inserts. T200: 4x M3 on
  a 19 mm square (base offset from the axis is a placeholder, `t200_base_offset`).
  Fairing breakaway ~420 N (2 nylon M3); frame holds 3x thrust and 1.25x breakaway
  head-on; side hits reported as a capacity.
- Workflow: user reviews parts (STL links / renders) and gives feedback; iterate
  the parametric parts, then feed slicer masses back into the optimizer.

## Next step

Iterate the front pod and bracket with the user, then the rear mount and skid pads
as parts. Then:

Stage 3: OpenFOAM interFoam. Validate on a Wigley catamaran first (reference Michell
values from `python -m hullopt.stage2 wigley`), mesh independence study, then the
finalists: free sinkage/trim surge cases, sway cases, actuator-disk thrusters,
pod window cases. Use CFD to calibrate `wave_factor`, `pod_window_cd` and the
rear-thruster Cd.
