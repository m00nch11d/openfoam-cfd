# S1223 half wing with sine ("chainsaw") leading edge - OpenFOAM RANS + transition

Rectangular wing, S1223 section, mean chord c = 0.30 m, half span b = 1.00 m,
leading edge x_LE(y) = -0.006 cos(2 pi y / 0.125) m (+-2 % c, wavelength
125 mm = 0.417 c, 8 waves; a peak sits on the root so the LE meets the
symmetry plane with zero slope). Sections are the S1223 scaled about the
straight trailing edge. Re = 3e5 on the mean chord, ISA sea-level air
(rho = 1.225 kg/m3, mu = 1.7894e-5 Pa s, U = 14.61 m/s), **AoA = 0 deg**
(not given in the request).

## Results (iteration 1494, run stopped by the user with `stopAt writeNow`)

| Quantity | Final iteration | Mean of last 100 iterations |
|---|---|---|
| CL | 0.8558 | 0.8542 |
| CD | 0.06218 | 0.06188 |
| CM (about the mean-chord c/4, nose-up +) | -0.2633 | -0.2632 |
| L/D | 13.76 | |
| Total cells | 2,944,680 | 1,996,800 hexahedra + 947,880 prisms |

Reference area S = c b = 0.30 m2 (half wing), reference length c = 0.30 m.
Independent integration of wall pressure + shear on the patch gives
CL = 0.8554 (check against forceCoeffs).

Convergence: residuals flat at Ux 3.5e-5, Uz 2e-4, p 1.9e-3, k 9e-5,
ReThetat 5e-5, gammaInt 4e-5. Peak-to-peak variation over the last 300
iterations: CL 0.48 %, CD 2.1 % (unsteadiness of the laminar separation
bubbles in a steady solver). No divergence. The run was stopped by the user
at iteration 1493 of the 4000 allowed; the |dCL|/CL < 1e-5 convergence stop
had not triggered yet.

### Spanwise stations (peaks and troughs of the LE wave)

Transition onset = first x/c where the boundary-layer maximum of nut/nu rises
above 2 (laminar 0.1-0.4, turbulent > 10). LSB = laminar separation bubble
from the sign of the streamwise wall shear. All x/c are on the local chord.

| Station | y [mm] | LE | local c [mm] | cl | cd | cm (c/4) | x_tr/c upper | x_tr/c lower | LSB upper (sep-reatt) | LSB lower (sep-reatt) |
|---|---|---|---|---|---|---|---|---|---|---|
| S00 | 0.0 | peak | 306.0 | 0.9352 | 0.06120 | -0.2732 | 0.904 | 0.260 | 0.817-0.880 | 0.011-0.031 |
| S01 | 62.5 | trough | 294.0 | 0.9687 | 0.04145 | -0.2687 | 0.397 | 0.035 | 0.305-0.392 | 0.009-0.016 |
| S02 | 125.0 | peak | 306.0 | 0.9363 | 0.06182 | -0.2741 | 0.887 | 0.259 | 0.799-0.912 | 0.011-0.069 |
| S03 | 187.5 | trough | 294.0 | 0.9546 | 0.03899 | -0.2602 | 0.430 | 0.039 | 0.312-0.481 | 0.009-0.014 |
| S04 | 250.0 | peak | 306.0 | 0.9291 | 0.06503 | -0.2751 | 0.891 | 0.227 | 0.794-0.921 | 0.011-0.031 |
| S05 | 312.5 | trough | 294.0 | 0.9618 | 0.04264 | -0.2679 | 0.414 | 0.033 | 0.305-0.363 | 0.007-0.016 |
| S06 | 375.0 | peak | 306.0 | 0.9091 | 0.06260 | -0.2692 | 0.897 | 0.226 | 0.811-0.912 | 0.011-0.073 |
| S07 | 437.5 | trough | 294.0 | 0.9429 | 0.04443 | -0.2661 | 0.399 | 0.035 | 0.305-0.370 | 0.007-0.014 |
| S08 | 500.0 | peak | 306.0 | 0.8902 | 0.06774 | -0.2726 | 0.910 | 0.241 | 0.811-0.871 | 0.011-0.019 |
| S09 | 562.5 | trough | 294.0 | 0.9253 | 0.05309 | -0.2784 | 0.446 | 0.028 | 0.312-0.370 | 0.007-0.025 |
| S10 | 625.0 | peak | 306.0 | 0.8476 | 0.06817 | -0.2695 | 0.905 | 0.262 | 0.823-0.921 | 0.011-0.019 |
| S11 | 687.5 | trough | 294.0 | 0.8641 | 0.05961 | -0.2628 | 0.441 | 0.015 | 0.312-0.384 | 0.007-0.016 |
| S12 | 750.0 | peak | 306.0 | 0.7920 | 0.07290 | -0.2635 | 0.917 | 0.154 | 0.817-0.895 | 0.009-0.028 |
| S13 | 812.5 | trough | 294.0 | 0.7609 | 0.06747 | -0.2506 | 0.434 | 0.018 | 0.312-0.384 | 0.007-0.016 |
| S14 | 875.0 | peak | 306.0 | 0.6702 | 0.08023 | -0.2448 | 0.892 | 0.053 | 0.860-- | 0.009-0.022 |
| S15 | 937.5 | trough | 294.0 | 0.6790 | 0.08905 | -0.2521 | 0.975 | 0.019 | 0.341-0.548 | 0.007-0.022 |
| S16 | 1000.0 | peak | 306.0 | 0.7198 | 0.15579 | -0.3155 | 0.497 | 0.035 | 0.925-- | 0.009-0.019 |

Pattern: behind the LE **troughs** a mid-chord bubble (x/c 0.31-0.39) trips
the upper surface at x/c ~0.40-0.45, while behind the **peaks** the upper
boundary layer stays laminar to a trailing-edge bubble (transition x/c ~0.89-0.92).
Troughs therefore carry ~0.03 more cl and ~35 % less cd than peaks inboard.
Stations 937.5 and 1000 mm are dominated by the tip vortex.

### y+

Airfoil surfaces: max 1.41, mean 0.295,
99th percentile 0.94, 99.2 % of faces below 1.
First layer height 1.456e-5 m from the flat-plate estimate for y+ = 0.7
(Cf = 0.0576 Re^-0.2 = 4.62e-3, u_tau = 0.702 m/s, y = 0.7 nu / u_tau).
The flat tip cap is resolved only by the 15.6 mm spanwise cell (mean y+ 81).

## Images (`results/`)

* `wing_3d_pressure.png`, `wing_3d_velocity.png`, `wing_3d_turbulence_intensity.png`
  (log scale), `wing_3d_yplus.png` - isometric from the top / leading-edge side
* `residuals.png`, `convergence_CL_CD_CM.png`
* `stations/Sxx_*_Cp.png` - upper and lower Cp with transition marks
* `stations/Sxx_*_pressure.png`, `_velocity.png`, `_turb_intensity.png` - section contours
* `stations.csv`, `summary.json` - all numbers
* `mesh_2d_section.png` - section mesh

## Mesh

* 2D template: 200 cells on each of the upper and lower surfaces
  (0.7 sine + 0.3 uniform spacing: fine at LE and TE, coarse mid-chord);
  structured quads in a 5 cm wall layer (40 layers, first 1.456e-5 m) and
  along the flat 30 c wake cut (120 half-cosine cells, fine at the TE);
  triangles elsewhere growing to 0.5 m (> c) at the C-domain boundary
  (radius 20 c around the TE, outlet 30 c behind the TE).
* Extruded spanwise: 64 layers over the wing (15.6 mm, 8 per LE wave), 14
  layers growing x1.28 to y = 2.71 m beyond the tip. Each plane is morphed so
  the section matches the local chord (morph fades out 0.05 -> 0.5 m from the wall).
* Patches: wing (wall, incl. flat tip cap), root (symmetryPlane), side,
  farfield (freestream), outlet.

## Solver set-up

OpenFOAM v1912 (Ubuntu 24.04 package), simpleFoam, `kOmegaSSTLM`
(Langtry-Menter k-omega SST gamma-ReTheta transition model), 4 MPI ranks.
Free stream Tu = 0.2 %, nut/nu = 5, ReThetat = 1061. SIMPLE with U 0.6,
p 0.3, turbulence 0.5 relaxation; linearUpwind for U, limitedLinear for
turbulence. (SIMPLEC with 0.8/0.7 diverged around the square tip at
iteration ~60.)

The monitor (`monitor.py`) stops the run on NaN/fatal errors, unphysical
coefficients or residuals > 1 (divergence), and on convergence when
|dCL|/CL < 1e-5 for 100 consecutive iterations. The requested "< 5 % between
two consecutive iterations" test is already met at iteration 16, long before
the flow has developed (CL was still moving by 10 % over the next 500
iterations), so it was tightened to 1e-5.

The Ubuntu v1912 package aborts any run with function objects
(`OSHA1stream` error in `dictionary::digest()`); `shim/digest_shim.c` is a
small LD_PRELOAD library that replaces that digest with a constant.

## Re-running

```
./run.sh            # mesh, case, decompose, solve (4 ranks), reconstruct, post-process
RESUME=1 ./run.sh   # continue from the latest decomposed time
```
Parameters live in `params.py` (AoA, LE wave, mesh sizes).
