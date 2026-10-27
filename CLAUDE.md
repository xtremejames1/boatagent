# Catamaran Hull Optimization — Project Brief

Context for Claude Code. This project builds a multi-fidelity hull optimization
pipeline for a small holonomic autonomous catamaran (RoboBoat-scale).

## Design decisions (settled)

- **Configuration:** catamaran, holonomic.
- **Thrusters:** 4 × Blue Robotics T200 (assume for now), all on 3D-printed mounts.
  - Front pair: one in front of each hull's bow, angled ~45°, inside a shroud that
    protects it from impacts and must be reasonably hydrodynamic.
  - Rear pair: between the hulls, unshrouded. Only the front pair has shrouds.
  - Overall length ≤ 1.52 m includes the front shrouds (assumption; hull `loa` is
    shorter).
- **Grounding:** the boat must rest on the hulls, not the thrusters. Achieved with a
  **separate, replaceable 3D-printed skid/rail** on each hull bottom that sits
  below the lowest thruster point. The hull draft is therefore a free variable,
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
| Skid | depth below lowest thruster (small margin, ~1–2 cm), length, fore-aft position. Low rounded rail, not a fin. |
| Shroud | depth, fairing chord/thickness. Fairing alignment tested as 2–3 discrete cases at stage 3 only (e.g. aligned with thruster axis vs. one intermediate). |

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
   sway cases, actuator-disk thrusters for interaction and ventilation, shroud
   fairing cases. Validate setup first on a Wigley catamaran benchmark; do a mesh
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
- Findings: skid height binds in stage 1 (`skid_max_height` = 0.12 m placeholder).
  Appendages (mostly the 45° front pods) are ~80 % of weighted resistance; long
  hulls (1.45–1.5 m) have ~30 % less hull resistance than short ones.
- Mount layout (front ahead of bows, rear between hulls) is set; mount details will
  still change. Keep thruster/shroud/mount geometry parametric (`hullopt/layout.py`).

## Next step

Stage 3: OpenFOAM interFoam. Validate on a Wigley catamaran first (reference Michell
values from `python -m hullopt.stage2 wigley`), mesh independence study, then the
finalists: free sinkage/trim surge cases, sway cases, actuator-disk thrusters,
shroud fairing cases. Use CFD to calibrate `wave_factor` and pod Cd values.
