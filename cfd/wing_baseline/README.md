# wing_baseline — external aerodynamics of the S1223 wing

Steady RANS (`simpleFoam`, `kOmegaSST`) of the CadQuery wing in `wing.py`,
meshed with `blockMesh` + `snappyHexMesh` from
`output/wing_outer_mold_line.stl`.

Written for **OpenFOAM v2412 (ESI / openfoam.com)**. It uses ESI-only dictionary
features (`#eval`, `#ifeq`), so it will not run as-is on the Foundation releases
(openfoam.org v11+ have no `simpleFoam`).

## Geometry (from `wing.py`)

| | |
|---|---|
| Airfoil | Selig S1223, blunt TE (0.3 % c) |
| Span (tip to tip) | 2.000 m |
| Root / tip chord | 0.300 / 0.060 m (taper 0.2) |
| Quarter-chord sweep | 25° |
| Flaperons | outer 10 % chord, y = 0.1 m to the tip, both sides, undeflected |
| Axes | x aft, y right wing, z up; origin at the root leading edge |

The STL is in **millimetres** (bounding box 586 × 2000 × 45). `Allrun`
scales it by 0.001 to metres with `surfaceTransformPoints`, then checks that
the scaled span is ~2 m and stops if it is not.

The flow is symmetric, so only the right half (y ≥ 0) is meshed, with a
symmetry plane at the root.

## Setting the conditions: `caseParameters`

Everything you would normally change is in this one file:

```
Umag            20;     // freestream speed [m/s]
alpha           5;      // angle of attack [deg]
nIterations     2000;
meshResolution  fine;   // fine | coarse
```

The angle of attack rotates the inflow vector, not the wing: the inlet velocity
and the lift/drag directions are derived from `Umag` and `alpha`. Nothing else
needs editing for an alpha sweep.

## Force coefficients

The `forceCoeffs` function object (`system/forceCoeffs`) reports Cl, Cd and Cm
every iteration to `postProcessing/forceCoeffs/0/coefficient.dat`. The
reference values are computed in `caseParameters` from the wing parameters:

| Quantity | Value | Meaning |
|---|---|---|
| `Aref` | 0.18 m² | half of the 0.36 m² trapezoidal planform (half model), so the coefficients are full-wing values |
| `lRef` | 0.2067 m | mean aerodynamic chord, 2/3·c_r·(1+λ+λ²)/(1+λ) |
| `CofR` | (0.2563, 0, 0) m | quarter-chord point of the MAC (at y = 0.389 m); moment reference |
| `liftDir` / `dragDir` | ⊥ / ∥ to the freestream | |
| `pitchAxis` | (0 1 0) | Cm is positive nose-up |

Dimensional forces use ρ = 1.225 kg/m³ (`rhoInf`). For the full wing, double
them.

## Running

```bash
source /path/to/openfoam2412/etc/bashrc   # or: openfoam2412
cd cfd/wing_baseline
./Allrun              # all cores
./Allrun -np 8        # or choose the number of processes
python3 plot_convergence.py               # replot at any time
./Allclean            # remove mesh and results
```

`Allrun` runs: STL scaling and check → `surfaceFeatureExtract` → `blockMesh`
→ `decomposePar` → `snappyHexMesh` → `checkMesh` → `potentialFoam` (initial
field) → `simpleFoam` → reconstruct → `plot_convergence.py`. Each step writes a
`log.<application>` file and the script stops at the first one that fails.
To rerun a step, delete its log (or run `./Allclean`).

The solver stops when `nIterations` is reached or when all residuals pass
`residualControl` in `system/fvSolution`. `plot_convergence.py` writes
`convergence.png` (residuals, Cl, Cd and Cm against iteration) and prints Cl, Cd,
Cm and L/D averaged over the last 100 iterations (`--avg N` changes that). Use
those averages, not the last iteration, because the values usually oscillate
slightly in a steady run. To view the flow fields, open `wing.foam` in ParaView
(`touch wing.foam`).

### Mesh resolution

| `meshResolution` | background cell | wing surface cells | prism layers | approx. size |
|---|---|---|---|---|
| `coarse` | 0.5 m | 7.8 / 3.9 mm | 3 | 0.19 M cells, smoke test only |
| `fine` | 0.25 m | 3.9 / 2.0 mm | 5 | a few million cells |

`coarse` was used to check that the whole chain runs without errors (see
below). Do not use its coefficients for anything. For a mesh study, change
the `#else` (fine) branch in `caseParameters`: `surfaceLevel`, `featureLevel`, `nLayers` and `baseCell`.

### Coarse-mesh check (done before committing)

OpenFOAM v2412, `meshResolution coarse`, `nIterations 400`, `./Allrun -np 4`:

* Total run time was about 2 minutes. Every step finished without errors.
* snappyHexMesh: 189 k cells, layers added, no faces failing the snappy
  quality checks. checkMesh reports one face with skewness 4.3 (limit 4).
* Residuals reached roughly 1e-5 (p 5e-5) after 400 iterations. Cl, Cd and Cm
  levelled off after about 200 iterations.
* On that coarse mesh: Cl ≈ 1.16, Cd ≈ 0.086, Cm ≈ −0.25, y+ average 31.
  These numbers are not mesh-converged. They only show that the values are
  in a sensible range.
* There are occasional small `bounding omega` / `bounding k` messages in cells
  next to the coarse wing surface. That is harmless at this resolution, but
  check whether they persist on the fine mesh.

## Modelling notes and limitations

* **Reynolds number.** About 2.8 × 10⁵ on the MAC at 20 m/s. The S1223 is a
  low-Re airfoil and at this Re its performance depends on laminar separation
  bubbles and transition. `kOmegaSST` treats the boundary layer as fully
  turbulent, so expect Cl to come out somewhat low and Cd somewhat high
  compared with wind-tunnel or XFOIL data. A transition model (e.g.
  `kOmegaSSTLM`) is the next step if that matters.
* **Wall treatment.** `nutUSpaldingWallFunction` is valid at any y+. The
  `yPlus` function object writes the y+ field at each write time; check it in
  ParaView. Snapped layers on the fine mesh give y+ of roughly 20–50.
* **Flaperon gaps.** The STL keeps the 1 mm hinge and end gaps between the
  flaperons and the main wing. At 2 mm surface cells the mesh does not resolve
  them, and snappy closes them, which is the intended "sealed gap" behaviour.
  If you refine below ~1 mm, snappy will start to mesh into the gaps.
  Either seal them in the geometry or check `checkMesh` there carefully.
* **Domain.** x ∈ [−5, 10] m, y ∈ [0, 5] m, z ∈ [−5, 5] m (about 25 MAC
  upstream, 50 downstream and 5 semi-spans to the side). The inlet uses a fixed
  velocity, the outlet a fixed pressure, and top, bottom and side use
  `freestreamVelocity`/`freestreamPressure`.
