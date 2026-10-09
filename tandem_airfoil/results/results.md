# Tandem airfoil results: S1223 + NACA 0009

| item | value |
|---|---|
| configuration | S1223 (c=0.30 m) + NACA 0009 (c=0.10 m), gap 0.20 m |
| flow | air ISA sea level, Re = 300000, U = 14.607 m/s, AoA = 0.0 deg, Tu = 0.1 % |
| model | kOmegaSSTLM (k-omega SST + gamma-ReTheta transition) |
| mesh | 66413 cells (66407 quads) |
| first layer | S1223 8.270e-06 m, NACA 0009 6.385e-06 m |
| iterations | 761 (stopped on user request at iteration 761 (8 MPI ranks); no divergence) |

## Aerodynamic coefficients

| body | ref. chord | Cl | Cd | Cm (c/4, +nose-up) | Cl mean +- std (last 100 it.) |
|---|---|---|---|---|---|
| S1223 | 0.30 m | 0.9930 | 0.03005 | -0.2487 | 0.9968 +- 0.0036 |
| NACA 0009 | 0.10 m | -0.4057 | -0.01515 | 0.0024 | -0.3838 +- 0.0124 |
| both (total) | 0.30 m | 0.8577 | 0.02500 | -0.0456 | 0.8688 +- 0.0075 |

## Transition (x/c from each leading edge)

| airfoil | surface | transition onset | method | laminar separation (Cf<0) | reattachment |
|---|---|---|---|---|---|
| S1223 | upper | 0.490 | gamma>0.5 & nut/nu>1 | 0.323 | 0.538 |
| S1223 | lower | 0.060 | gamma>0.5 & nut/nu>1 | 0.023 | 0.085 |
| NACA 0009 | upper | 0.996 | bubble: Cf minimum | 0.975 | 0.998 |
| NACA 0009 | lower | 0.549 | bubble: Cf minimum | 0.382 | 0.682 |

Transition onset: first station where, within 0.5 mm of the wall, the intermittency > 0.5 and nut/nu > 1; when a laminar separation bubble closes before that, the Cf minimum inside the bubble (separation-induced transition).

## y+

| airfoil | max | mean |
|---|---|---|
| S1223 | 0.699 | 0.185 |
| NACA 0009 | 0.714 | 0.133 |
