# Tandem airfoil CFD: S1223 + NACA 0009 (OpenFOAM, k-ω SST γ-Reθ transition)

2-D steady RANS analysis of an S1223 airfoil (chord 30 cm) with a NACA 0009
(chord 10 cm) placed on the same horizontal line, 20 cm behind the S1223 trailing
edge. Air at ISA sea level, Re = 3·10⁵ (based on the S1223 chord), AoA = 0°.

## Run it (one command)

```bash
cd tandem_airfoil
./Allrun.sh                      # 8 MPI ranks, max 1000 iterations
./Allrun.sh -np 16 -maxiter 1500 # options
```

The script:

1. finds OpenFOAM (uses the environment if `simpleFoam` is on the PATH, otherwise
   sources the newest `/usr/lib/openfoam/openfoam*/etc/bashrc`,
   `/opt/openfoam*/etc/bashrc` or `~/OpenFOAM/OpenFOAM-v*/etc/bashrc`),
2. creates `./.venv` with `gmsh numpy scipy matplotlib` if the system Python
   lacks them,
3. builds the mesh, converts it (`gmshToFoam`), runs `checkMesh`,
4. writes the case, decomposes it, runs `simpleFoam` in parallel under a monitor
   that stops and reports on divergence and stops early on Cl convergence,
5. reconstructs the case and writes all figures and tables to `./results`.

The OpenFOAM case is left in `./run`. Open `run/case.foam` in ParaView.

**Requirements:** an openfoam.com (ESI) release, v2006 or newer (tested with
**v2412**). The Foundation releases (openfoam.org v11+) have no `simpleFoam`. On
Ubuntu, Python needs `python3-venv`. If `import gmsh` fails, install
`libglu1-mesa libxcursor1 libxinerama1 libxft2`.

> Tested here with OpenFOAM v2412 on a 4-core machine. The Ubuntu 24.04
> `openfoam` apt package (v1912) is **not** usable: every function object
> (forces, y+) crashes with an `IOstream "sha1"` error in that build.

## Set-up

| item | value |
|---|---|
| fluid | air, ISA sea level: ρ = 1.225 kg/m³, μ = 1.7894·10⁻⁵ Pa·s, ν = 1.4607·10⁻⁵ m²/s |
| free stream | U∞ = Re·ν/c = **14.607 m/s**, AoA = 0°, Tu = 0.1 %, ν_t/ν = 1 |
| model | `kOmegaSSTLM`: k-ω SST + Langtry-Menter γ-Reθt transition (4 equations), with ambient decay control |
| solver | `simpleFoam`, SIMPLEC, second-order upwind for U, limitedLinear for turbulence |
| BCs | farfield: freestream velocity/pressure; outlet: p = 0; walls: no-slip, k = 0, γ/Reθ zero gradient |
| forces | Cl, Cd, Cm (about c/4, positive nose-up) for each airfoil on its own chord, total on the S1223 chord |

### First layer thickness (y⁺)

Flat-plate estimate: C_f = 0.058·Re_c^-0.2, u_τ = U∞·√(C_f/2), Δy = y⁺·ν/u_τ, with y⁺ = 0.7.

| airfoil | Re_c | C_f | u_τ [m/s] | estimate Δy | used Δy |
|---|---|---|---|---|---|
| S1223 | 3.0·10⁵ | 0.00466 | 0.705 | 1.45·10⁻⁵ m | **8.27·10⁻⁶ m** |
| NACA 0009 | 1.0·10⁵ | 0.00580 | 0.787 | 1.30·10⁻⁵ m | **6.39·10⁻⁶ m** |

With the flat-plate value, y⁺ peaked at 1.23 and 1.42 at the leading-edge
suction peaks. The spacing was therefore scaled by 0.7/peak. The final peak
y⁺ is about 0.7 on both airfoils (mean 0.14 to 0.19).

### Mesh (`scripts/mesh.py`)

* **C-domain:** a semicircle of radius 20c (6 m) centred on the S1223 trailing
  edge, with a flat wake section 30c (9 m) long behind it. The NACA 0009 sits
  inside that wake.
* **Surfaces:** 200 + 200 elements on the S1223 and 100 + 100 on the NACA 0009,
  with a sine distribution. Elements are small at the LE and TE and largest
  mid-chord (90 % cosine law, 10 % uniform).
* **Quad layer:** a 4 cm structured quad layer around each airfoil, built as a
  local C-grid with a short wake cut along the TE bisector. Wall-normal growth is
  1.10, giving 65 (S1223) and 68 (NACA) layers.
* **Wake:** a flat wake line from the NACA block to the outlet. Node spacing
  follows a cosine law, x' = L(1 − cos θ) with θ uniform: small near the airfoils
  and growing downstream.
* **Outer region:** quad-dominant (gmsh frontal-quad + blossom). Cell size grows
  smoothly from the quad-layer edge (growth 0.08 m per m) to 0.40 m (> c) on the
  far boundary.
* **Total:** 66,413 cells (66,407 hexahedra, 6 prisms). checkMesh: max
  non-orthogonality 60°, max skewness 1.0, "Mesh OK".

### Run control (`scripts/run_solver.py`)

* **Divergence:** a NaN, a floating-point exception, a fatal error, a residual
  above 10, or |Cl| above 50 stops the run. The status is reported, the exit
  code is 2, and no post-processing is done.
* **Convergence:** the run stops when the relative Cl change between two
  consecutive iterations stays below `CL_TOL` for 50 consecutive iterations,
  counted after iteration 200. The default `CL_TOL` is 1·10⁻⁴; change it with
  `-cltol`. When triggered, `stopAt writeNow` is written into `controlDict` so
  the solver writes the current state and exits cleanly.
* **Iteration limit:** `-maxiter`, default 1000.

### Transition detection (`scripts/post.py`)

The transition onset is the first station from the leading edge where both
conditions hold in the cells within 0.5 mm of the wall:

* the model intermittency γ is above 0.5,
* ν_t/ν is above 1.

In the laminar boundary layer γ ≈ 0.02 and ν_t/ν < 0.1. The thin band keeps the
free stream (γ = 1) and the S1223 wake (which impinges on the NACA) out of the
measure. Laminar separation and reattachment come from the sign of C_f.

## Outputs (`results/`)

| file | content |
|---|---|
| `01_pressure_distribution_Cp.png` | Cp on upper and lower surfaces of both airfoils, with transition marked |
| `02_pressure_contours.png` | Cp contours: whole domain, TE/gap region, near field |
| `03_velocity_contours.png` | velocity magnitude contours |
| `04_turbulence_intensity_log.png` | turbulence intensity, log colour scale |
| `05_yplus_distribution.png` | y⁺ along both airfoils |
| `06_residuals.png` | initial residuals of all equations, plus Cl history |
| `07_pathlines.png` | 26 pathlines from x = −0.2 m to 60 cm behind the S1223 TE |
| `08_mesh.png` | mesh overview and LE/TE close-ups |
| `results.md`, `results.json` | Cl, Cd, Cm, transition, separation, y⁺ |

Every figure is stamped with the airfoils, element counts, flow conditions,
AoA and turbulence model.

## Files

```
Allrun.sh             one-command driver
s1223.dat             S1223 coordinates (Selig format)
scripts/params.py     all parameters (flow, mesh, solver, run control)
scripts/mesh.py       mesh generator (gmsh + structured quad layers -> msh2)
scripts/setup_case.py OpenFOAM dictionaries
scripts/run_solver.py parallel run + divergence / Cl-convergence monitor
scripts/post.py       coefficients, transition, figures
results/              results of the reference run
```
