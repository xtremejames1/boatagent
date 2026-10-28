# Front pod and bracket

Print settings: PLA, 4 walls (1.6 mm), 10 % gyroid, 0.2 mm layers, 250 x 250 x 250 mm bed.

## Parts

| part | size (mm) | fits bed | slicer mass (g) | print time | CAD volume (cm³) |
|---|---|---|---|---|---|
| front_bracket | 177.6 x 63.9 x 167.4 | yes | 112.32 | 8h 48m 34s | 182.1 |
| front_pod | 205.3 x 186.7 x 140.6 | yes | 174.29 | 16h 7m 25s | 137.5 |

Optimizer's estimate: pod 146 g; all four mounts 324 g.

## Strength

Allowable stress 25 MPa (printed PLA, across layers). Utilisation = combined stress / allowable; buckling = compression / Euler load. Both must stay under 1.

Pod breakaway: 2 nylon M3 screws in shear ≈ 420 N. Hull inserts: assumed 800 N pull-out each (test this).

Side hits are secondary: the bracket holds about 136 N sideways at the pod nose before a member or hull insert reaches its limit.

| load case | member | axial (N) | stress (MPa) | utilisation | buckling | max bolt tension (N) | max bolt shear (N) |
|---|---|---|---|---|---|---|---|
| forward thrust x3 (155 N) | front_port_upper | -95.1 | 5.97 | 0.25 | 0.025 | 341.0 | 124.2 |
|  | front_port_lower | 197.6 | 15.95 | 0.672 | 0.0 |  |  |
| reverse thrust x3 (121 N) | front_port_upper | 74.3 | 4.66 | 0.195 | 0.0 | 223.7 | 97.0 |
|  | front_port_lower | -154.3 | 12.46 | 0.525 | 0.01 |  |  |
| head-on impact at breakaway x1.25 (525 N) | front_port_upper | 429.9 | 12.21 | 0.488 | 0.0 | 411.0 | 0.0 |
|  | front_port_lower | -923.7 | 16.98 | 0.679 | 0.059 |  |  |
| side hit, secondary, x0.5 breakaway (210 N) | front_port_upper | 0.0 | 1.83 | 0.102 | -0.0 | 641.8 | 233.0 |
|  | front_port_lower | 0.0 | 37.48 | 1.542 | -0.0 |  |  |
