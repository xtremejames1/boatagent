# Test hull, 24 in (610 mm)

The full-size hull (loa 1217 mm) scaled by 0.501. Hull gap (inner faces) 304.8 mm = 12.0 in; hull centrelines 366.1 mm apart; overall width 427 mm.

| | test hull | full size |
|---|---|---|
| length | 610 mm | 1217 mm |
| beam (each hull) | 61.3 mm | 122.4 mm |
| depth (keel to deck) | 112.1 mm | 223.7 mm |
| rocker bow / stern | 38.8 / 17.2 mm | 77.4 / 34.4 mm |
| deadrise, bilge, entry | 5.0°, 0.90, 21.0° | same |

Angles and proportions are unchanged; the hull gap is not scaled (fixed at 12 in), so the two hulls sit relatively further apart than on the full-size boat.

Model-scale (Froude) testing: speed scales by √0.501 = 0.708, so a full-size 1.0 m/s is 0.71 m/s here. The matching total mass is about 3.23 kg (full size 25.7 kg × 0.501³).

## Draft against load (both hulls, level)

| total mass (kg) | draft (mm) | freeboard (mm) | waterline length (mm) | KB (mm) | BM_T (mm) | GM_T if KG at deck (mm) |
|---|---|---|---|---|---|---|
| 1.616 | 38.7 | 73.4 | 609 | 24.1 | 1299 | 1211 |
| 2.425 | 51.7 | 60.4 | 610 | 31.1 | 869 | 788 |
| 3.233 | 64.6 | 47.5 | 610 | 37.9 | 652 | 578 |
| 4.041 | 77.6 | 34.5 | 610 | 44.5 | 521 | 454 |
| 4.849 | 90.5 | 21.6 | 610 | 51.1 | 435 | 374 |
| 6.466 | sinks | | | | | |

## Print pieces (PLA, 4 walls, 10 % gyroid, 250 mm bed)

Each hull prints in pieces, standing on the aft cut face (no supports). Print two sets (port and starboard are identical). The cut faces are flat; glue them, or add pins.

| piece | x range (mm) | size (mm) | overhang > 35° (mm²) | mass, 10 % infill (g) | print time | mass, hollow (g) | print time |
|---|---|---|---|---|---|---|---|
| demihull_piece1 | 0.0–203.2 | 203.2 x 60.7 x 111.7 | 0.0 | 249.11 | 14h 48m 57s | 128.08 | 8h 9m 16s |
| demihull_piece2 | 203.2–406.4 | 203.2 x 61.3 x 112.1 | 0.0 | 293.36 | 17h 18m 35s | 141.72 | 9h 2m 10s |
| demihull_piece3 | 406.4–609.6 | 202.8 x 59.8 x 104.2 | 0.0 | 192.64 | 11h 23m 4s | 105.08 | 6h 29m 45s |

One hull: 735 g with 10 % infill, 375 g hollow (4 walls, closed at the cut faces). Two hulls: 1470 g / 750 g, against the 3.23 kg model-scale total.

Files: `demihull.step` (whole hull, solid), `demihull_pieceN.step` / `_print.stl` (print pieces, standing on the aft cut), `demihull.stl`, `catamaran.stl` (both hulls at the 12 in gap, for viewing or CFD), `design.json`.
