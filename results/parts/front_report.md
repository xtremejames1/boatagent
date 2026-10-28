# Front pod and bracket

Print settings: PLA, 4 walls (1.6 mm), 10 % gyroid, 0.2 mm layers, 250 x 250 x 250 mm bed.

## Parts

| part | size (mm) | fits bed | slicer mass (g) | print time | CAD volume (cm³) |
|---|---|---|---|---|---|
| front_bracket | 177.6 x 63.9 x 167.4 | yes | 114.3 | 9h 15m 57s | 180.8 |
| front_pod | 205.1 x 186.7 x 140.6 | yes | 295.12 | 22h 43m 6s | 604.4 |

Optimizer's estimate: pod 215 g; all four mounts 324 g.

## Strength

Allowable stress 25 MPa (printed PLA, across layers). Utilisation = combined stress / allowable; buckling = compression / Euler load. Both must stay under 1.

Pod breakaway: 2 nylon M3 screws in shear ≈ 420 N. Hull inserts: assumed 800 N pull-out each (test this).

| load case | member | axial (N) | stress (MPa) | utilisation | buckling | max bolt tension (N) | max bolt shear (N) |
|---|---|---|---|---|---|---|---|
| forward thrust x3 (155 N) | front_port_upper | -95.0 | 5.97 | 0.25 | 0.025 | 341.0 | 124.2 |
|  | front_port_lower | 197.5 | 15.95 | 0.672 | 0.0 |  |  |
| reverse thrust x3 (121 N) | front_port_upper | 74.2 | 4.66 | 0.195 | 0.0 | 223.6 | 97.0 |
|  | front_port_lower | -154.2 | 12.46 | 0.525 | 0.01 |  |  |
| head-on impact at breakaway x1.25 (525 N) | front_port_upper | 429.3 | 12.21 | 0.488 | 0.0 | 410.5 | 0.0 |
|  | front_port_lower | -923.2 | 16.99 | 0.68 | 0.059 |  |  |
| side hit at breakaway x1.25 (525 N) | front_port_upper | 0.0 | 4.59 | 0.257 | -0.0 | 1604.4 | 582.0 |
|  | front_port_lower | 0.0 | 93.69 | 3.855 | -0.0 |  |  |
