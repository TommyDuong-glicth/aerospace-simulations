Re = 1e6, AR = 8.0, e = 0.85; window -2° to 8°; v1 C_D0 = 0.01 (assumed).

**NACA 2412**

| Metric | v1 inviscid theory | XFOIL viscous | Error / divergence |
|---|---|---|---|
| 2D lift slope [1/rad] | 6.283 | 6.224 | +0.9 % |
| α_L0 [deg] | -2.08 | -2.30 | +0.22° |
| RMSE of c_l, linear window | n/a | n/a | 0.0256 (4.4 % of mean abs c_l) |
| c_l at α_stall = 16.5° | 2.037 | 1.582 (max) | 0.455 |
| (L/D)max, 3D | 23.11 | 31.24 | -8.13 |
| α at (L/D)max [deg] | 3.38 | 1.93 | +1.45° |

**NACA 0012**

| Metric | v1 inviscid theory | XFOIL viscous | Error / divergence |
|---|---|---|---|
| 2D lift slope [1/rad] | 6.283 | 6.566 | -4.3 % |
| α_L0 [deg] | 0.00 | 0.04 | -0.04° |
| RMSE of c_l, linear window | n/a | n/a | 0.0237 (6.2 % of mean abs c_l) |
| c_l at α_stall = 15.5° | 1.700 | 1.389 (max) | 0.311 |
| (L/D)max, 3D | 23.11 | 28.61 | -5.50 |
| α at (L/D)max [deg] | 5.45 | 3.54 | +1.91° |
