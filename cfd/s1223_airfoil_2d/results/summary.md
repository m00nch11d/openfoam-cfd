# S1223, Re = 300000, alpha = 0 deg, kOmegaSSTLM

Run status: MAX_ITERATIONS - reached nIterations without meeting the Cl criterion  
Iterations: 4000

Coefficients averaged over iterations 3000-4000; fields: averaged from iteration 3000.

| Coefficient | Mean | Std | Min .. max in window | Final iteration |
|---|---|---|---|---|
| Cl | 1.0656 | 0.0108 | 1.0366 .. 1.0850 | 1.0757 |
| Cd | 0.0200 | 0.0020 | 0.0163 .. 0.0247 | 0.0199 |
| Cm (c/4, nose-up +) | -0.2517 | 0.0015 | -0.2545 .. -0.2476 | -0.2543 |

Elements: 103238 (75900 quadrilaterals, 27338 triangles)

| Transition | Upper surface | Lower surface |
|---|---|---|
| transition (BL nut/nu > 2), x/c | 0.473 | 0.169 |
| laminar separation (Cf < 0), x/c | 0.326 | 0.017 |
| reattachment, x/c | 0.519 | 0.285 |

y+ max / mean: upper 0.66 / 0.29, lower 0.88 / 0.21

Mesh quality:

```
cells                                      103238  (quads 75900, triangles 27338)
airfoil faces                              242  (120 upper, 120 lower, 2 base)
first cell height                          1.158e-05 m
equiangle quality min / 1st pct / mean     0.574 / 0.743 / 0.960
  quads min, triangles min                 0.603, 0.574
orthogonal quality min / 1st pct / mean    0.782 / 0.923 / 0.994
  quads min, triangles min                 0.812, 0.782
far-field (C-boundary) edge length min / max 0.581 / 0.589 m
```

The literal rule 'stop when |dCl/Cl| between two consecutive iterations < 5 %' would have stopped at iteration 8 with Cl = 0.0864.
