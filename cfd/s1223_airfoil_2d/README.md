# s1223_airfoil_2d — S1223 airfoil, Re = 3×10⁵, α = 0°, γ–Reθ transition

2D steady RANS of the Selig S1223 (`../../s1223.dat`) with OpenFOAM v2412
(openfoam.com): `simpleFoam` with **`kOmegaSSTLM`**, the k-ω SST model with
the Langtry–Menter γ–Reθt transition model. The mesh is a hybrid C-mesh
built by `tools/make_mesh.py` and written directly as an OpenFOAM polyMesh.

## Results (`results/`)

| | |
|---|---|
| **Cl** | **1.066** (± 0.011) |
| **Cd** | **0.0200** (± 0.0020) |
| **Cm** (about c/4, nose-up positive) | **−0.252** (± 0.0015) |
| Elements | 103 238 (75 900 quadrilaterals, 27 338 triangles) |
| Transition, upper surface | x/c = **0.47**, inside a laminar separation bubble from x/c = 0.33 to 0.52 |
| Transition, lower surface | x/c = **0.17**, after a laminar separation bubble from x/c = 0.02 to 0.29 |
| y+ (first cell centre) | max 0.66 upper / 0.88 lower; mean 0.29 / 0.21 |

Values are means (± one standard deviation) over iterations 3000–4000. The
laminar separation bubbles make the steady solution oscillate slightly, so a
single-iteration value would be misleading (see *Convergence* below).
`results/summary.md` and `summary.json` hold the full numbers.

| Image | Content |
|---|---|
| `cp_distribution.png` | Cp over the upper and lower surfaces, with the transition points |
| `pressure_contours.png` | Cp field, whole domain and near the airfoil |
| `velocity_contours.png` | \|U\|/U∞, whole domain and near the airfoil |
| `turbulence_intensity.png` | Tu = √(2k/3)/U∞ on a logarithmic colour scale |
| `yplus_distribution.png` | y+ over the upper and lower surfaces |
| `residuals.png` | initial residuals of all 7 equations |
| `mesh_domain.png` | the whole C-mesh |
| `mesh_airfoil.png` | 2 × 2 chords around the airfoil |
| `mesh_details.png` | leading edge, trailing edge, blunt base close-ups |
| `cf_distribution.png`, `transition_indicator.png` | skin friction and boundary-layer νt/ν: how the transition points were found |
| `force_coefficients.png` | Cl, Cd, Cm against iteration, with the averaging window |

## Flow conditions (`caseParameters`)

Everything you are likely to change is in `caseParameters`: flow, mesh and run control.

* Air, ISA sea level: ρ = 1.225 kg/m³, μ = 1.7894×10⁻⁵ Pa s, so
  ν = 1.4607×10⁻⁵ m²/s.
* Chord c = 0.3 m and Re = 3×10⁵, so **U∞ = Re·ν/c = 14.607 m/s**. α = 0°
  (the inflow is rotated, not the airfoil).
* Freestream turbulence Tu = 0.1 %, νt/ν = 1, held constant up to the airfoil by
  `decayControl`. Reθt∞ = 1137 from the Langtry–Menter correlation.
* Boundary conditions:
  * C-boundary: `freestreamVelocity` / `freestreamPressure`.
  * Outlet: fixed p, `inletOutlet` U.
  * Airfoil: wall-resolved (k = 0, `omegaWallFunction` in its viscous-sublayer
    branch, `nutLowReWallFunction`).
* Reference values: chord 0.3 m, area chord × depth, moments about the
  quarter chord.

## Mesh (`tools/make_mesh.py`)

* **Airfoil:** the spline through `s1223.dat` is cut straight at x/c = 0.999,
  giving a **blunt base 0.0315 mm** thick (2 elements, treated as wall).
  **120 elements on each surface** with sine (cosine) spacing: 0.055 mm at
  the leading and trailing edges, 4.2 mm at mid-chord.
* **First layer:**
  * Flat-plate estimate: Cf = (2 log₁₀Re − 0.65)^−2.3 = 0.00468, so
    u_τ = U∞√(Cf/2) = 0.7065 m/s. For y+ = 0.7 the first cell **centre** sits
    at y = 1.45×10⁻⁵ m, so the cell height is 2.9×10⁻⁵ m.
  * The first run with that height gave a mean y+ of 0.6 (on target), but
    peaks of 1.3–2.3 at the nose and at the TE corners. So the height is
    scaled by `firstCellFactor 0.4` to **1.16×10⁻⁵ m**, which keeps y+ < 1
    everywhere (max 0.88).
  * The height has to be the same all around the wall. Any variation tapers
    these very thin cells and lowers their orthogonal quality.
* **Quadrilaterals:** a structured C-band covers the wall and 10 cm away from
  it, plus the wake, with 80 layers at a growth ratio of 1.086.
  * The **wake extends 30 chords** behind the TE, with one-sided cosine
    spacing (0.09 mm at the TE, 40 mm at the outlet).
  * The wake cut leaves the TE along the TE bisector (−35°) and curves to
    horizontal within a quarter chord. A flat cut would shear the TE cells.
  * The thin base sits at the head of a 2-cell slit block that runs along
    the wake cut.
* **Triangles:** gmsh Frontal-Delaunay fills the rest. Triangle size grows
  from the band edge size at 10 % per unit distance up to 0.6 m.
  * The **far-field edges are 0.58–0.59 m**, about 2 chords.
  * The C-boundary is a semicircle of radius 30 chords around the leading
    edge, plus straight lines to the outlet.
* **Quality (target ≥ 0.7):**
  * **Orthogonal quality** (Fluent definition: cosine between each face
    normal and the cell-to-face and cell-to-neighbour vectors) is
    **min 0.78**, mean 0.99. It is the metric used for the 0.7 target
    because it is the one that drives the solver's non-orthogonality error
    and it does not penalise the long, thin cells that y+ < 1 requires.
  * Equiangle quality (1 − equiangle skewness) is min 0.57; 99 % of cells
    are ≥ 0.74. The lowest values are:
    * the first slit cells, which join the vertical base to the −35° wake
      cut;
    * a few triangles.
  * `checkMesh`: max non-orthogonality 36°. It warns about aspect ratio
    (up to about 3500) in the thin wake-cut cells near the outlet; those
    cells are aligned with the flow.

## Running

```bash
source /path/to/openfoam2412/etc/bashrc
pip install numpy scipy matplotlib gmsh       # mesh + plots
cd cfd/s1223_airfoil_2d
./Allrun              # nProcs from caseParameters (24)
./Allrun -np 8        # or choose the number of MPI ranks
./Allclean
```

`Allrun` builds the mesh, runs `checkMesh`, `renumberMesh`,
`decomposePar` and `potentialFoam`, then runs `simpleFoam` with
`tools/monitor.py` watching it, then `reconstructPar` and
`tools/postprocess.py`. If your Python is not `python3`, set it with
`PYTHON=/path/to/python ./Allrun`.

### Threads

The case is set up for **24 MPI ranks** (`nProcs 24`; decomposition into 24
subdomains was tested). The run in `results/` used **4 ranks**, because the
machine it ran on has only 4 cores. On 4 cores, 24 oversubscribed ranks ran
at about 38 s per iteration (MPI ranks spin while waiting), against
0.1 s per iteration with 4 ranks. The partitioning does not change the
converged result.

### Stopping rules (`tools/monitor.py`)

* **Divergence:** a NaN or inf in a residual or coefficient, |Cl| > 20 or
  |Cd| > 10, or a solver crash. The run stops (`stopAt writeNow`), the
  reason is written to `results/run_status.txt`, and `Allrun` exits with
  code 2. None occurred.
* **Maximum iterations:** 4000.
* **Cl convergence:** after 500 iterations, stop once |ΔCl/Cl| between
  consecutive iterations stays below `clTol` (1×10⁻⁵) for `clWindow` (200)
  iterations in a row.
  * The requested threshold, 5 % between two consecutive iterations, is met
    almost at once: at **iteration 8, with Cl = 0.086**, while the flow is
    still developing from the potential-flow start.
  * The check is therefore kept but with a tight tolerance over a window.
    Both values are parameters in `caseParameters`.
  * In this run Cl kept a small oscillation (±1 %, from the separation
    bubbles), so the run went to 4000 iterations.

### Convergence

* Residuals level off at about 10⁻⁴–10⁻⁵ for U<sub>x</sub>, k, ω, γ and Reθ,
  and at about 2×10⁻³ for p and U<sub>y</sub>.
* The forces settle into a limit cycle: Cl 1.04–1.09, Cd 0.016–0.025.
* This is the usual behaviour of steady RANS with unsteady laminar
  separation bubbles. That is why fields and coefficients are averaged
  (`fieldAverage` from `averageStart` = 3000).
* An unsteady (pimpleFoam) run would be the next step if the bubble
  dynamics matter.

### How the transition points are found

Near-wall γ stays small (about 0.02) everywhere, because the wall cell sits
in the viscous sublayer, so it cannot mark transition. Instead,
`postprocess.py` takes, at each wall face, the peak eddy-viscosity ratio
νt/ν across the boundary layer (0.01–5 mm from the wall). That ratio is
below 0.6 while laminar and 10–50 once turbulent. Transition is the first
x/c after which it stays above 2. The laminar separation and reattachment
points are where Cf changes sign. Both are plotted in
`transition_indicator.png` and `cf_distribution.png`.

## Notes

* The fully-turbulent `kOmegaSST` model on the same mesh gives Cl ≈ 1.03,
  Cd ≈ 0.022, Cm ≈ −0.245. The transition model adds lift and removes drag,
  as expected at this Reynolds number.
* As a check of the mesher, boundary conditions and schemes, NACA 0012 at
  α = 5° (same tools, fully turbulent) gave Cl = 0.50 and Cd = 0.0146.
