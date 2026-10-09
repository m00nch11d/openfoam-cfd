# S1223 half wing with a sine ("chainsaw") leading edge - OpenFOAM RANS + transition

Rectangular wing, S1223 section, mean chord c = 0.30 m, half span b = 1.00 m,
root = symmetry plane. Leading edge x_LE(y) = -0.006 cos(2 pi y / 0.125) m
(+-2 % c, wavelength 125 mm = 0.417 c, 8 waves; a peak sits on the root so
the LE meets the symmetry plane with zero slope); sections are the S1223
scaled about the straight trailing edge. Re = 3e5 on the mean chord, ISA
sea-level air (rho 1.225 kg/m3, mu 1.7894e-5 Pa s, U = 14.61 m/s),
AoA = 0 deg. k-omega SST gamma-ReTheta transition model (`kOmegaSSTLM`).

## Run it on your computer (Ubuntu 24.04)

```
cd wing_cfd
./run_all.sh            # installs everything (asks for sudo once), meshes, runs, post-processes
```

* **Edit `params.py`** - the single setup file: geometry, flow (AoA, Re, Tu),
  mesh sizes, number of processes (`N_PROCS = 24`: 12 cores x 2 threads),
  maximum iterations (`N_ITER_MAX = 1000`), stop rules, pathline length.
* Options: `--clean` (delete `case/` + `results/` and start over),
  `--resume` (continue a run), `--post-only` (redo the images/tables),
  `--no-install` (skip apt/pip).
* What it installs: Ubuntu's `openfoam` package (v1912), OpenMPI, build tools,
  OpenGL/OSMesa libraries (apt) and a Python venv `.venv` with gmsh, pyvista,
  numpy, scipy, matplotlib (`requirements.txt`).
* Expected run time: ~4.2 M cells. On the 4-core test machine 7.7 s per
  iteration (~2 h for 1000 iterations); scales roughly with core count.
* Memory: ~7 GB for the solver, ~6 GB for post-processing.
* Output: `results/` (see below); logs in `case/log.*`; run status in
  `case/run_status.txt`.

Files: `params.py` (settings) -> `mesh2d.py` (2D hybrid C-mesh) ->
`extrude.py` (3D polyMesh) -> `setup_case.py` (0/, constant/, system/) ->
`monitor.py` (stop rules during the run) -> `postprocess.py` (results);
`shim/digest_shim.c` (work-around, see below); `s1223.dat` (airfoil).

## Results of the test run (iteration 932, stopped by the user before the 1000 limit; no divergence)

| Quantity | Final iteration | Mean of last 100 iterations |
|---|---|---|
| CL | 0.8361 | 0.8354 |
| CD | 0.06162 | 0.06087 |
| CM (about the mean-chord c/4, nose-up +) | -0.2601 | -0.2586 |
| L/D | 13.57 | |
| Total cells | 4,236,896 | 2,867,200 hexahedra + 1,369,696 prisms |

Reference area S = c b = 0.30 m2 (half wing), reference length c = 0.30 m.
Independent integration of wall pressure + shear gives CL = 0.8356.
Wing surface: 200 x 80 = 16,000 faces on the upper and on the lower surface, 726 on the tip cap.

Convergence: final initial residuals Ux 5e-5, Uy 1.7e-4, Uz 2.2e-4, p 1.4e-3,
k 9e-5, ReThetat 6e-5, gammaInt 5e-5. Peak-to-peak over the last 100
iterations: CL 0.30 %, CD 2.9 %; over the last 300: CL 1.6 %, CD 4.1 % -
the solution was still settling slowly (laminar separation bubbles in a
steady solver) when stopped; let it run to 1000 locally.

### Spanwise stations (peaks and troughs of the LE wave)

Transition onset = first x/c where the boundary-layer maximum of nut/nu rises
above 2 (laminar 0.1-0.4, turbulent > 10). LSB = laminar separation bubble
from the sign of the streamwise wall shear. x/c on the local chord. The tip
peak (y = 1000 mm) is sampled 7.8 mm inboard: exactly at the square tip edge
the section integrates the edge suction and gives a meaningless negative cd.

| Station | y [mm] | LE | local c [mm] | cl | cd | cm (c/4) | x_tr/c upper | x_tr/c lower | LSB upper (sep-reatt) | LSB lower (sep-reatt) |
|---|---|---|---|---|---|---|---|---|---|---|
| S00 | 0.0 | peak | 306.0 | 0.9293 | 0.06174 | -0.2728 | 0.891 | 0.252 | 0.794-- | 0.011-0.028 |
| S01 | 62.5 | trough | 294.0 | 0.9685 | 0.03929 | -0.2629 | 0.427 | 0.039 | 0.305-0.488 | 0.009-0.014 |
| S02 | 125.0 | peak | 306.0 | 0.9194 | 0.06132 | -0.2694 | 0.902 | 0.254 | 0.805-0.876 | 0.012-0.037 |
| S03 | 187.5 | trough | 294.0 | 0.9454 | 0.03884 | -0.2521 | 0.406 | 0.025 | 0.305-0.362 | 0.006-0.018 |
| S04 | 250.0 | peak | 306.0 | 0.9104 | 0.06279 | -0.2694 | 0.901 | 0.227 | 0.794-0.895 | 0.012-0.023 |
| S05 | 312.5 | trough | 294.0 | 0.9264 | 0.03994 | -0.2491 | 0.396 | 0.030 | 0.297-0.384 | 0.006-0.015 |
| S06 | 375.0 | peak | 306.0 | 0.8945 | 0.06517 | -0.2681 | 0.906 | 0.221 | 0.823-0.895 | 0.012-0.020 |
| S07 | 437.5 | trough | 294.0 | 0.9080 | 0.04241 | -0.2452 | 0.394 | 0.020 | 0.297-0.340 | 0.005-0.015 |
| S08 | 500.0 | peak | 306.0 | 0.8673 | 0.06749 | -0.2663 | 0.905 | 0.225 | 0.811-0.890 | 0.013-0.023 |
| S09 | 562.5 | trough | 294.0 | 0.8862 | 0.04863 | -0.2515 | 0.416 | 0.024 | 0.297-0.362 | 0.005-0.012 |
| S10 | 625.0 | peak | 306.0 | 0.8269 | 0.06979 | -0.2645 | 0.907 | 0.251 | 0.800-0.908 | 0.013-0.029 |
| S11 | 687.5 | trough | 294.0 | 0.8367 | 0.05519 | -0.2439 | 0.450 | 0.026 | 0.311-0.503 | 0.005-0.015 |
| S12 | 750.0 | peak | 306.0 | 0.7562 | 0.07110 | -0.2542 | 0.901 | 0.240 | 0.834-- | 0.010-0.018 |
| S13 | 812.5 | trough | 294.0 | 0.7523 | 0.05935 | -0.2380 | 0.473 | 0.018 | 0.318-0.547 | 0.005-0.010 |
| S14 | 875.0 | peak | 306.0 | 0.6316 | 0.07347 | -0.2334 | 0.560 | 0.054 | 0.526-0.541 | 0.010-0.020 |
| S15 | 937.5 | trough | 294.0 | 0.6040 | 0.07888 | -0.2198 | 0.687 | 0.012 | 0.333-0.391 | 0.005-0.010 |
| S16 | 1000.0 | peak | 306.0 | 0.6619 | 0.17781 | -0.3322 | 0.321 | 0.539 | 0.963-0.994 | 0.011-0.037 |

Behind the LE **troughs** a laminar bubble starting at x/c ~0.30-0.32 trips
the upper surface at x/c ~0.39-0.47; behind the **peaks** (root to 750 mm) the
upper boundary layer stays laminar to a bubble near the trailing edge
(transition x/c ~0.89-0.91). Inboard of 600 mm the troughs carry 0.01-0.04
more cl and ~30-40 % less cd (0.039-0.049 vs 0.061-0.067) than the peaks.
The outer ~150 mm are dominated by the tip vortex.

### y+

Airfoil surfaces: mean 0.35, 97.8 % of faces below 1,
max 2.10. Inboard of y = 0.94 m: 99.2 % below 1, max 1.38 (at the
LE suction peaks). Most faces above 1 sit within 0.5 mm of the sharp square
tip edge, where the flow turns round the corner. Tip cap: mean y+ 4.1
(first spanwise cell 0.12 mm).
First layer height 1.456e-5 m = flat-plate estimate for y+ = 0.7
(Cf = 0.0576 Re^-0.2 = 4.62e-3, u_tau = 0.702 m/s, y = 0.7 nu / u_tau).

## Images (`results/`)

* `wing_3d_pressure.png`, `wing_3d_velocity.png`,
  `wing_3d_turbulence_intensity.png` (log scale), `wing_3d_yplus.png` -
  isometric view from the top / leading-edge side
* `wing_3d_tip_vortex_pathlines.png` (isometric) and `..._rear.png`
  (looking upstream) - particle paths through the tip region to TE + 3c
* `residuals.png`, `convergence_CL_CD_CM.png`
* `stations/Sxx_*_Cp.png` - upper and lower Cp with the transition marked
* `stations/Sxx_*_pressure.png`, `_velocity.png`, `_turb_intensity.png` -
  section contours at each station
* `stations.csv`, `summary.json` - all numbers; `mesh_2d_section.png`

## Mesh

* 2D template: 200 cells on each of the upper and lower surfaces
  (0.7 sine + 0.3 uniform spacing: fine at LE and TE, coarse mid-chord);
  structured quads in a 5 cm wall layer (40 layers, first 1.456e-5 m) and
  along the flat 30 c wake cut (120 half-cosine cells, fine at the TE);
  triangles elsewhere, growing to 0.5 m (> c) at the C-domain boundary
  (radius 20 c around the TE, outlet 30 c behind the TE).
* Spanwise: 61 uniform cells (15.4 mm, 8 per LE wave) from the root, then
  clustered (x1.3) to 0.12 mm on both sides of the tip plane, growing to
  0.41 m outboard (side boundary at y = 2.77 m). Every plane is morphed so
  the section matches the local chord (morph fades 0.05 -> 0.5 m from the wall).
* Patches: wing (wall, incl. flat tip cap), root (symmetryPlane), side,
  farfield (freestream), outlet. checkMesh: max non-orthogonality 77 deg
  (128 faces > 70 deg at the sharp-TE slivers), high aspect ratio in the wall
  layers only (expected).

## Solver set-up and stop rules

simpleFoam, `kOmegaSSTLM`, Tu = 0.2 %, nut/nu = 5, ReThetat = 1061.
SIMPLE with U 0.6 / p 0.3 / turbulence 0.5 relaxation (SIMPLEC 0.8/0.7
diverged at the square tip), linearUpwind for U, limitedLinear for turbulence,
GAMG (DICGaussSeidel, relTol 0.1) for p. Decomposition: hierarchical
(2 x N/2 for N >= 8 processes).

`monitor.py` stops the run with a final write on divergence (NaN/fatal error,
|CL| > 10, CD < 0 or > 5 after iteration 100, or an initial residual > 1) and
on convergence when |CL_n - CL_n-1| / |CL_n| < `CL_TOL` (1e-5) for 100
consecutive iterations after iteration 300. The requested "< 5 % between two
consecutive iterations" is already met at iteration ~16, while CL still moves
by > 10 %; set `CL_TOL = 0.05, CL_WINDOW = 1, MIN_ITER = 1` in `params.py` to
use it literally.

**Ubuntu OpenFOAM v1912 work-around:** the package aborts any run with
function objects (`OSHA1stream` error in `dictionary::digest()`).
`shim/digest_shim.c` is a tiny LD_PRELOAD library (built automatically) that
replaces that digest with a constant; it only affects detecting edits to
function-object settings during a run.

**Pathlines:** VTK cannot locate points reliably inside the micron-thin
wall/wake cells, so the tip region is resampled onto a uniform 5 mm grid
before the particle paths are traced.
