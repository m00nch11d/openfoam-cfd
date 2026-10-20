# Plunging NACA 0009 - results

| item | value |
|---|---|
| airfoil | NACA 0009, c = 0.30 m, AoA = 0 deg |
| flow | air ISA sea level, Re = 300000, U = 14.607 m/s, Tu = 0.1 % |
| motion | y = 0.300 m sin(2 pi 3.00 Hz t); max plunge speed 5.65 m/s, max effective AoA 21.2 deg, k = 0.194, St = 0.123 |
| model | kOmegaSSTLM (k-omega SST + gamma-ReTheta transition), pimpleFoam (rigid mesh motion) |
| mesh | 52403 cells (52395 quads), first cell 4.50e-06 m |
| time step | 1.667e-04 s (2000 per cycle), 3.00 cycles computed |
| run | converged - cycle statistics changed by 7.9 % (< 10 %) from cycle 2 to 3 |

## Coefficients over the last cycle

| quantity | min | max | mean |
|---|---|---|---|
| Cl | -2.390 | 2.325 | -0.010 (rms 1.422) |
| Cd | -0.1750 | 0.1631 | -0.0163 |
| Cm (c/4, +nose-up) | -0.6184 | 0.6215 | 0.0014 |

Thrust coefficient Ct = -mean(Cd) = 0.0163; input power coefficient mean(-Cl dy/dt / U) = 0.3808; propulsive efficiency Ct/Cp = 0.043.

## Cycle-to-cycle convergence

| cycle | Cl min | Cl max | Cl rms | mean Cd |
|---|---|---|---|---|
| 1 | -2.147 | 2.537 | 1.202 | -0.0477 |
| 2 | -2.309 | 2.041 | 1.213 | -0.0429 |
| 3 | -2.390 | 2.325 | 1.424 | -0.0162 |

## Four phases of the last cycle

| phase | y [m] | dy/dt [m/s] | alpha_eff | Cl | Cd | Cm | transition upper | transition lower |
|---|---|---|---|---|---|---|---|---|
| 0 deg | +0.000 | +5.65 | -21.2 deg | -1.675 | -0.1175 | +0.0860 | 0.022 (bubble: Cf minimum) | 0.022 (gamma>0.5 & nut/nu>1) |
| 90 deg | +0.300 | -0.00 | +0.0 deg | +0.289 | +0.0238 | -0.0771 | - (gamma>0.5 & nut/nu>1) | 0.167 (gamma>0.5 & nut/nu>1) |
| 180 deg | -0.000 | -5.65 | +21.2 deg | +2.199 | -0.0576 | -0.1559 | 0.022 (gamma>0.5 & nut/nu>1) | 0.022 (bubble: Cf minimum) |
| 270 deg | -0.300 | +0.00 | -0.0 deg | -0.308 | +0.0201 | +0.0748 | 0.220 (bubble: Cf minimum) | - (gamma>0.5 & nut/nu>1) |

Transition x/c from the leading edge: inner-BL intermittency > 0.5 and nut/nu > 1 within 0.5 mm of the wall, or the Cf minimum inside a closed laminar separation bubble. The full phase dependence is in 11_transition_vs_phase.png.

## y+ (maximum over the last cycle)

| surface | max y+ |
|---|---|
| upper | 0.563 |
| lower | 0.616 |
