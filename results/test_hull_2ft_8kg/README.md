# Test hull, 2 ft, for 8 kg payload

610 mm long, 123 mm beam, 140 mm keel to deck, 12 in (304.8 mm) between the hulls' inner faces; hull centrelines 427.5 mm apart; overall width 550 mm.

Based on the 2 ft scaled hull, with twice the beam, 140 mm depth, a wider transom (75 % of max beam), a fuller bow (25° entry) and less rocker (25 mm bow, 10 mm stern).

Each hull is about 575 g printed hollow (PLA, 4 walls); estimated from the sliced 2 ft hull, scaled by surface area.

## Draft against payload (fresh water, level)

| payload (kg) | total incl. hulls (kg) | draft (mm) | freeboard (mm) | GM_T if KG at deck (mm) |
|---|---|---|---|---|
| 0 | 1.15 | 21.7 | 118.3 | 3587 |
| 2 | 3.15 | 40.4 | 99.6 | 1609 |
| 4 | 5.15 | 57.1 | 82.9 | 990 |
| 6 | 7.15 | 73.7 | 66.3 | 692 |
| 8 | 9.15 | 90.2 | 49.8 | 528 |
| 9 | 10.15 | 98.5 | 41.5 | 472 |
| 10 | 11.15 | 106.8 | 33.2 | 426 |

## Files

- `demihull_smooth.step`: one smooth solid (B-spline sides and bottom, flat deck, flat transom, a 1–2 mm flat at the bow tip), mm. x forward from the transom, z = 0 at the lowest keel point, centred on y = 0. Copy it to the other side with centrelines 427.5 mm apart.
- `demihull.stl`, `catamaran.stl`: meshes for viewing.
- `design.json`: the hull parameters; `hullopt.cadhull.smooth_demihull` rebuilds the STEP from them.
