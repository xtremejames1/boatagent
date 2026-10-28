# Front frame and nose fairing

Print: PLA, 4 walls (1.6 mm), 10 % gyroid, 0.2 mm layers, 250 x 250 x 250 mm bed, nose-up (boat +x up), no supports, overhang limit 35° from vertical.

Sized members: upper strut 12 x 32 mm (thick x deep); rail 24 x 12 mm (wide x thick); forefoot 20 mm thick.

## Parts

| part | size x,y,z (mm) | fits bed nose-up | overhang > limit (mm²) | slicer mass (g) | print time |
|---|---|---|---|---|---|
| front_frame | 208.8 x 68.4 x 247.3 | yes | 368.2 | 168.0 | 12h 45m 13s |
| front_fairing | 124.6 x 187.4 x 141.3 | yes | 0.0 | 89.04 | 6h 13m 41s |

Optimizer's estimate: pod 172 g; all four mounts 365 g.

## Checks

- frame overlaps thruster: 0.0 mm³
- fairing overlaps thruster: 0.0 mm³
- frame overlaps fairing: 0.0 mm³
- frame inside inlet keep-out: 0.0 mm³
- frame inside jet keep-out: 17.9 mm³
- fairing inside jet keep-out: 0.0 mm³
- thruster removal path blocked by: 0.0 mm³
- M3 screws without tool access: 0

## Strength

Allowable 25 MPa (printed PLA across layers, walls only). Fairing breaks off at about 420 N (2 nylon M3). Cap screws: M4 into heat-set inserts, 1000 N each assumed in shear (test this). Side hits are secondary: the frame takes about 78 N sideways at the nose.

| load case | member | axial (N) | stress (MPa) | utilisation | buckling | max cap screw shear (N) |
|---|---|---|---|---|---|---|
| forward thrust x3 (155 N) | upper_a | 18.9 | 21.67 | 0.937 | 0.0 | 205.7 |
|  | upper_b | -66.1 | 12.99 | 0.573 | 0.006 |  |
|  | forefoot | -8.4 | 9.51 | 0.412 | 0.001 |  |
|  | rail | 69.3 | 11.41 | 0.458 | 0.0 |  |
| reverse thrust x3 (121 N) | upper_a | -14.7 | 16.92 | 0.731 | 0.002 | 160.7 |
|  | upper_b | 51.7 | 10.15 | 0.448 | 0.0 |  |
|  | forefoot | 6.6 | 7.43 | 0.322 | 0.0 |  |
|  | rail | -54.1 | 8.91 | 0.357 | 0.028 |  |
| head-on hit x1.25 breakaway (525 N) | upper_a | -201.6 | 21.87 | 0.879 | 0.03 | 964.3 |
|  | upper_b | -233.6 | 14.06 | 0.562 | 0.02 |  |
|  | forefoot | 7.4 | 24.59 | 0.986 | 0.0 |  |
|  | rail | -287.9 | 11.2 | 0.449 | 0.147 |  |
| side hit, secondary, x0.5 breakaway (210 N) | upper_a | 22.1 | 66.29 | 2.695 | 0.0 | 0.0 |
|  | upper_b | 52.3 | 41.64 | 1.675 | 0.0 |  |
|  | forefoot | -0.1 | 16.03 | 1.285 | 0.0 |  |
|  | rail | -53.4 | 62.42 | 2.498 | 0.027 |  |
