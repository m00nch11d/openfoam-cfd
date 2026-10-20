# Plunging NACA 0009 - results

| item | value |
|---|---|
| airfoil | NACA 0009, c = 0.30 m, AoA = 0 deg |
| flow | air ISA sea level, Re = 300000, U = 14.607 m/s, Tu = 0.1 % |
| motion | y = 0.300 m sin(2 pi 3.00 Hz t); max plunge speed 5.65 m/s, max effective AoA 21.2 deg, k = 0.194, St = 0.123 |
| model | kOmegaSSTLM (k-omega SST + gamma-ReTheta transition), pimpleFoam (rigid mesh motion) |
| mesh | 52403 cells (52395 quads), first cell 4.50e-06 m |
| time step | 1.667e-04 s (2000 per cycle), 5.00 cycles computed |
| run | max_cycles - reached the limit of 5 cycles without meeting the 10 % cycle criterion |

## Coefficients over the last cycle

| quantity | min | max | mean |
|---|---|---|---|
| Cl | -2.682 | 2.157 | -0.083 (rms 1.347) |
| Cd | -0.1672 | 0.1842 | -0.0322 |
| Cm (c/4, +nose-up) | -0.1654 | 0.5857 | 0.0461 |

Thrust coefficient Ct = -mean(Cd) = 0.0322; input power coefficient mean(-Cl dy/dt / U) = 0.3544; propulsive efficiency Ct/Cp = 0.091.

## Cycle-to-cycle convergence

| cycle | Cl min | Cl max | Cl rms | mean Cd | largest Cl change vs previous cycle | vs two cycles before |
|---|---|---|---|---|---|---|
| 1 | -2.147 | 2.537 | 1.202 | -0.0477 | - | - |
| 2 | -2.309 | 2.041 | 1.213 | -0.0429 | 24.3 % | - |
| 3 | -2.390 | 2.325 | 1.424 | -0.0162 | 14.8 % | 15.5 % |
| 4 | -2.226 | 2.015 | 1.220 | -0.0362 | 16.7 % | 3.7 % |
| 5 | -2.682 | 2.157 | 1.347 | -0.0322 | 17.0 % | 10.9 % |

Convergence criterion: Cl max, Cl min and Cl rms of a cycle each within 10 % of their values one cycle earlier (period 1) or two cycles earlier (period 2, a pattern dynamic stall can lock into), checked from cycle 3 on.

Averages over the last two cycles: Cl rms 1.286, mean Cd -0.0342 (Ct = 0.0342), power coefficient 0.3400, propulsive efficiency 0.100.

## Four phases of the last cycle

| phase | y [m] | dy/dt [m/s] | alpha_eff | Cl | Cd | Cm | transition upper | transition lower |
|---|---|---|---|---|---|---|---|---|
| 0 deg | +0.000 | +5.65 | -21.2 deg | -1.759 | -0.0883 | +0.0471 | - (gamma>0.5 & nut/nu>1) | 0.022 (gamma>0.5 & nut/nu>1) |
| 90 deg | +0.300 | -0.00 | +0.0 deg | +0.363 | +0.0156 | -0.0745 | 0.997 (bubble: Cf minimum) | 0.190 (bubble: Cf minimum) |
| 180 deg | -0.000 | -5.65 | +21.2 deg | +1.905 | -0.1070 | -0.1237 | 0.022 (gamma>0.5 & nut/nu>1) | - (gamma>0.5 & nut/nu>1) |
| 270 deg | -0.300 | +0.00 | -0.0 deg | -0.479 | +0.0362 | +0.0685 | 0.125 (gamma>0.5 & nut/nu>1) | - (gamma>0.5 & nut/nu>1) |

Transition x/c from the leading edge: inner-BL intermittency > 0.5 and nut/nu > 1 within 0.5 mm of the wall, or the Cf minimum inside a closed laminar separation bubble. The full phase dependence is in 11_transition_vs_phase.png.

## y+ (maximum over the last cycle)

| surface | max y+ |
|---|---|
| upper | 0.528 |
| lower | 0.725 |
