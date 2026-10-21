# Plunging NACA 0009 (OpenFOAM, URANS, k-ω SST γ-Reθ transition)

A 2-D unsteady RANS analysis of a NACA 0009 airfoil oscillating up and down in a
uniform stream of air.

| | |
|---|---|
| airfoil | NACA 0009, chord c = 0.30 m, AoA = 0° |
| flow | air at ISA sea level, Re = 3·10⁵ → U∞ = 14.607 m/s, Tu = 0.1 % |
| motion | plunge y(t) = A sin(2πft) with A = 1 c = 0.30 m and f = 3 Hz (T = 0.333 s) |
| derived | peak plunge speed Aω = 5.65 m/s; peak effective AoA ±21.2°; reduced frequency k = πfc/U∞ = 0.19; Strouhal number 2fA/U∞ = 0.12 |
| model | `kOmegaSSTLM`: k-ω SST + Langtry-Menter γ-Reθt transition, ambient decay control |
| solver | `pimpleFoam` on a rigidly moving mesh (`solidBody` motion), initialised from a steady `simpleFoam` solution |

## Run it (one command)

```bash
cd plunging_naca0009
./Allrun.sh                       # 8 MPI ranks, up to 5 plunge cycles
./Allrun.sh -np 16 -maxcycles 8   # options
./Allrun.sh -cycletol 0.05        # periodic-convergence tolerance (default 0.10)
./Allrun.sh -noanim               # skip the animation
./Allrun.sh -postonly             # redo figures/animation from an existing ./run
./Allrun.sh -smoketest            # ~5 min check that your OpenFOAM version accepts
                                  # everything (writes to smoketest_run/ and
                                  # smoketest_results/, never touches run/ or results/)
```

> On an OpenFOAM version other than v2412 (e.g. v2606), run `./Allrun.sh -smoketest`
> first. It meshes, runs the steady initialisation, 60 moving-mesh time steps and the
> post-processing, and ends with `SMOKE TEST PASSED`.

The script runs these steps:

1. Finds OpenFOAM. If `pimpleFoam` is not on the PATH, it sources the newest
   `/usr/lib/openfoam/openfoam*/etc/bashrc`, `/opt/openfoam*/etc/bashrc` or
   `~/OpenFOAM/OpenFOAM-v*/etc/bashrc`.
2. Creates `./.venv` with gmsh, numpy, scipy and matplotlib if the system
   Python lacks them.
3. Meshes the domain, converts the mesh with `gmshToFoam` and checks it.
4. Runs the steady initialisation (`run/steady`, simpleFoam, 400 iterations,
   static airfoil).
5. Maps that solution onto the plunging case (`run/unsteady`) and runs
   `pimpleFoam` in parallel under a monitor. The monitor stops and reports on
   divergence, and stops early once the flow is periodic.
6. Reconstructs the last cycle and writes all figures, the animation and the
   tables to `./results`.

**Requirements:** an openfoam.com (ESI) release, v2006 or newer; tested with
**v2412**. The Ubuntu 24.04 apt package (v1912) is not usable, because its
function objects crash. On Ubuntu, `python3-venv` is needed; if `import gmsh`
fails, install `libglu1-mesa libxcursor1 libxinerama1 libxft2`. Install
`ffmpeg` to also get an MP4 of the animation.

**Disk:** the case keeps 40 field snapshots (the last cycle) in ascii, about
2–3 GB including the processor directories.

## How the requirements were applied

### Mesh (`scripts/mesh.py`)

* **Domain:** a C-domain with an arc of radius 20c centred on the trailing edge
  (TE) and a flat wake section extending 30c (9 m) behind the TE.
* **Surface elements:** 200 on the upper and 200 on the lower surface. They
  follow a sine distribution (80 % cosine law, 20 % uniform): small at the
  leading edge (LE) and TE, largest mid-chord. The share of the cosine law was
  lowered from 90 % to 80 % so the TE cells are about 0.3 mm. The smaller TE
  cells made the time-accurate solver unstable.
* **Quad layer:** the wall and the first 4 cm are a structured quad layer, a
  local C-grid with a short wake cut along the TE bisector. The wall-normal
  growth ratio is 1.10.
* **Wake:** the wake line uses a cosine node distribution, small at the TE and
  growing downstream.
* **Wake band:** the whole mesh moves with the airfoil, so in the mesh frame the
  wake leaves the TE at y = 0 and swings to ±A downstream. A refined band
  therefore opens as a wedge from the TE: 5 mm cells at the TE, growing
  slowly over 8 chords, with a half-height of up to A + 0.1 m.
* **Outer region:** quad-dominant, with cell size growing gradually from the
  quad-layer edge to 0.40 m (> c) on the far boundary.

### First layer thickness (y⁺ < 1, target 0.7)

The flat-plate estimate uses the **peak relative velocity**,
√(U∞² + (Aω)²) = 15.66 m/s:
C_f = 0.058 Re^-0.2 = 0.00459, u_τ = 0.751 m/s, so Δy(y⁺ = 0.7) = 1.36·10⁻⁵ m.
At effective angles of attack around 20°, the leading-edge suction peak raises
the local wall shear well above the flat-plate value. The spacing used is
therefore 0.33 × the estimate, **4.50·10⁻⁶ m**. The peak y⁺ over the whole
cycle is checked in post-processing and reported in `results.md` and
`05_yplus_distribution.png`.

### Time stepping and stability

What it took to get a stable run, in the order the problems appeared:

1. **Impulsive start.** `oscillatingLinearMotion` starts at full plunge speed
   (5.65 m/s) against 4.5 µm wall cells, and the run diverged immediately. The
   motion is now a table (`tabulated6DoFMotion`, one entry per time step) with
   the amplitude ramped smoothly over the first half cycle.
2. **Trailing-edge cells.** The sine distribution was eased from 90 % to 80 %
   cosine weight, giving TE cells of about 0.3 mm, and the wake cut starts at
   twice the TE spacing. This lowers the Courant number in the thin
   TE/wake-cut cells.
3. **Dynamic stall.** With 3 PIMPLE outer iterations per step, the run went
   through the first downstroke and then diverged at t ≈ 0.205 s. That is when
   the leading-edge (dynamic-stall) vortex reaches the trailing edge
   (Cl ≈ 2, effective AoA ≈ +20°). Restarts from t = 0.200 s compared four
   remedies:

| variant | result |
|---|---|
| velocity limiter (`limitVelocity`, 3 U∞) + upwind turbulence convection | survives, but the limiter clipped up to ~300 cells continuously and Cl was noisy: rejected |
| 12 outer iterations (2 correctors, relaxation 0.3) | stable, smooth Cl, 4 s/step |
| 8 outer × 1 corrector | crashed |
| **6 outer × 2 correctors, relaxation p 0.5 / equations 0.8** | **same Cl as the 12-outer run, 2.3 s/step: used** |
| 4 outer, same relaxation | noisy Cl: rejected |

With converged outer iterations, Δt = T/2000 = 1.67·10⁻⁴ s (about 120 steps
per chord-transit time) reproduces the T/3000 result within about 3 %, so
2000 steps per cycle are used. The maximum Courant number is 20–45, all of it
in the thin near-wall cells.

> **Reference run note:** the run in `results/` covered the first 0.2 s
> (60 % of the start-up cycle) with the earlier settings (3000 steps per cycle,
> 3 outer iterations). It was then continued from t = 0.2 s with the final
> settings above. All reported cycles after the first used the final settings.
> `./Allrun.sh` uses the final settings from t = 0.

### Run control (`scripts/run_solver.py`)

* **Divergence:** a NaN, a floating-point exception, a fatal error, a residual
  above 10, a Courant number above 10⁴, or |Cl| above 50 stops the run with a
  report, and no post-processing is done.
* **Convergence:** "Cl change < 10 %" is applied the way a periodic flow
  requires: cycle to cycle. After each completed cycle, the monitor compares
  that cycle's Cl max, Cl min and Cl rms with the previous cycle's. The run
  stops when each has changed by less than 10 % of its own value, checked from
  cycle 3 on.
* **Two-cycle periodicity:** dynamic stall can lock into a pattern that repeats
  every *second* cycle, with the leading-edge vortex alternating between two
  shedding patterns. So the cycle two back is also compared, and agreement
  within 10 % there counts as convergence too.
* **Mean Cd:** reported but not used for stopping. It is a small difference
  of large thrust and drag phases, and it fluctuates much more from cycle to
  cycle than Cl.
* **Limit:** 5 cycles (10,000 time steps). Each time step is itself iterated
  with 6 PIMPLE outer correctors.

### Post-processing (`scripts/post.py`, `scripts/animate.py`)

* **Snapshots:** the last cycle is stored as 40 snapshots. Figures show four
  plunge phases (y = A sin φ):
  * φ = 0°: mid-stroke moving up
  * φ = 90°: top of the stroke
  * φ = 180°: mid-stroke moving down
  * φ = 270°: bottom of the stroke
* **Moment:** Cm is taken about the **moving** quarter chord. The forces
  function object uses a fixed reference point, so post-processing corrects
  it: Cm_c/4 = Cm_fixed − (y/c)·Cd.
* **Thrust and efficiency:**
  * thrust coefficient: Ct = −mean(Cd)
  * input power coefficient: mean(−Cl·ẏ/U∞)
  * propulsive efficiency: their ratio
* **Transition:** in cells within 0.5 mm of the wall, intermittency above 0.5
  and turbulent viscosity above the molecular viscosity. When a laminar
  bubble closes first, the Cf minimum inside the bubble is used instead.
  Transition is evaluated at all 40 snapshots.
* **Pathlines:** particles are tracked in physical time through the moving,
  time-dependent flow. Snapshots are interpolated in phase, and the mesh
  position is known exactly from y(t).

## Reference results (5 cycles, 4 MPI ranks)

Full tables are in `results/results.md`.

| quantity (last cycle) | value |
|---|---|
| Cl | −2.68 … +2.16 (rms 1.35, mean −0.08) |
| Cd | −0.167 … +0.184, **mean −0.032 (net thrust, Ct = 0.032)** |
| Cm (c/4, +nose-up) | −0.17 … +0.59 |
| propulsive efficiency Ct / Cp_in | 0.09 (0.10 averaged over the last two cycles) |
| max y⁺ over the cycle | 0.73 (lower surface), 0.53 (upper surface) |

* **Phase behaviour:**
  * At mid-stroke (effective AoA ±21°) the suction side transitions at
    x/c ≈ 0.02, just behind the suction peak, and a leading-edge (dynamic
    stall) vortex forms and travels over the airfoil. It shows as the moving
    Cp bump near x/c ≈ 0.2.
  * The pressure side stays laminar.
  * At the top and bottom of the stroke, transition retreats to
    x/c ≈ 0.1–0.2, inside laminar separation bubbles.
* **Convergence:** the run reached the 5-cycle limit without meeting the 10 %
  criterion. From cycle 2 on, Cl rms stays between 1.21 and 1.42, but the peak
  values vary by 7–17 % from cycle to cycle (cycle 4 vs 2: 3.7 %; cycle 5
  vs 3: 10.9 %). The dynamic-stall vortex shedding is not strictly periodic at
  this effective AoA, so expect cycle-to-cycle scatter of this size. The mean
  Cd varies from −0.016 to −0.043 between cycles (about −0.032 averaged over
  cycles 2–5). The multi-cycle averages are the robust numbers.

## Outputs (`results/`)

| file | content |
|---|---|
| `01_pressure_distribution_Cp.png` | Cp on the upper and lower surface at four phases, with transition marked |
| `02_pressure_contours.png` | Cp contours: whole domain, wake, near field at four phases |
| `03_velocity_contours.png` | velocity magnitude, same layout |
| `04_turbulence_intensity_log.png` | turbulence intensity, logarithmic colour scale |
| `05_yplus_distribution.png` | y⁺ along the airfoil at four phases, plus the cycle maximum |
| `06_residuals.png` | residuals per time step and the maximum Courant number |
| `07_pathlines.png` | 25 particles released at x = −0.2 m at four phases, tracked to 60 cm behind the TE |
| `08_coefficients_history.png` | Cl, Cd, Cm history, last cycle vs phase, Cl–α_eff and Cd–α_eff loops |
| `09_mesh.png`, `10_mesh_3d.png` | mesh views (2-D and 3-D) |
| `11_transition_vs_phase.png` | transition and laminar separation over the cycle |
| `12_plunging_animation.gif` / `.mp4` | vorticity, moving airfoil and particles over the last cycle |
| `results.md`, `results.json` | all numbers |

Every figure is stamped with the airfoil, element counts, flow conditions,
AoA and the plunge motion.
