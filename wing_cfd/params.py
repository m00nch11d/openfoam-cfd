"""
=============================================================================
 SETUP FILE - edit this file only, then run   ./run_all.sh
=============================================================================
S1223 rectangular half wing with a sine ("chainsaw") leading edge,
OpenFOAM simpleFoam + k-omega SST gamma-ReTheta transition (kOmegaSSTLM).

All values are SI. Axes: x streamwise (aft +), y spanwise (root y=0 ->
tip y=+B), z up. Mean-chord leading edge at x=0, straight trailing edge at
x=C. The root (y=0) is a symmetry plane.
"""
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))

# ============================================================== geometry
AIRFOIL_FILE = os.path.join(HERE, "s1223.dat")   # Selig format, unit chord
C = 0.30                 # mean chord [m]
B = 1.00                 # half span [m]
LE_AMP = 0.02 * C        # leading-edge sine amplitude (+-2 % chord) [m]
LE_WAVELENGTH = 0.125    # leading-edge sine wavelength [m]; B / (lambda/2)
                         # must be an integer (stations = peaks + troughs)
# x_LE(y) = -LE_AMP*cos(2 pi y / LE_WAVELENGTH): a peak (LE most forward,
# longest chord) sits on the symmetry plane so the LE meets it with zero
# slope. Sections are the S1223 scaled about the fixed trailing edge.

# ============================================================== flow
RE = 3.0e5               # Reynolds number on the mean chord
RHO = 1.225              # ISA sea level density [kg/m3]
MU = 1.7894e-5           # ISA sea level viscosity [Pa s]
NU = MU / RHO            # 1.4607e-5 m2/s
U_INF = RE * NU / C      # 14.61 m/s
AOA_DEG = 0.0            # angle of attack [deg]
TU = 0.002               # free-stream turbulence intensity (0.2 %)
VISC_RATIO = 5.0         # free-stream nut/nu

# ============================================================== 2D mesh
YPLUS_TARGET = 0.7       # sets the first-layer height (see H1 below)
N_SURF = 200             # cells on EACH of the upper and lower surfaces
SINE_BLEND = 0.7         # surface spacing = SINE_BLEND*sine + rest uniform
                         # (1.0 = pure sine: 10 um cells at LE/TE)
N_WAKE = 120             # cells along the flat 30c wake cut (half-cosine)
WAKE_LEN = 30.0 * C      # wake / domain length behind the TE
QUAD_THICK = 0.05        # structured quad layer around wall and wake [m]
N_LAYERS = 40            # quad layers across QUAD_THICK
R_FAR = 20.0 * C         # C-domain radius (centred at the TE)
FAR_SIZE = 0.50          # far-field triangle size [m] (>= chord)
GROWTH = 0.15            # triangle size gradient (size grows 0.15 m per m)

# ============================================================== span mesh
DY_WING = 0.015625       # spanwise cell size on the wing [m] (8 per LE wave)
TIP_DY_FIRST = 1.2e-4    # spanwise cell next to the flat tip cap, both
                         # sides [m] (gives y+ < 1 on the tip cap)
TIP_RATIO = 1.3          # spanwise growth away from the tip plane
TIP_DY_MAX = 0.45        # outboard growth stops at this cell size [m]

# ============================================================== solver / run
N_PROCS = 4              # MPI processes (threads); use <= physical cores
N_ITER_MAX = 1000        # maximum iterations
WRITE_INTERVAL = 250     # write fields every N iterations
RELAX_U = 0.6            # SIMPLE relaxation (SIMPLEC 0.8/0.7 diverges at the tip)
RELAX_P = 0.3
RELAX_TURB = 0.5

# stop rules (monitor.py)
CL_TOL = 1.0e-5          # stop when |CL_n - CL_n-1| / |CL_n| < CL_TOL ...
CL_WINDOW = 100          # ... for this many consecutive iterations ...
MIN_ITER = 300           # ... but not before this iteration.
# NOTE: CL_TOL = 0.05, CL_WINDOW = 1, MIN_ITER = 1 is the literal "5 % between
# two consecutive iterations" rule; it is met at iteration ~16, long before
# the flow has developed, so a tighter value is the default.
DIVERGE_CL = 10.0        # |CL| above this (after iteration 100) = divergence

# ============================================================== post-processing
PATHLINE_X_BEHIND = 3.0 * C   # tip-vortex pathlines run to TE + this [m]


# ---------------------------------------------------------- derived values
def first_layer_height(yplus=YPLUS_TARGET):
    """Flat-plate estimate (Schlichting: Cf = 0.0576 Re^-0.2)."""
    cf = 0.0576 * RE ** -0.2
    tau_w = 0.5 * RHO * U_INF ** 2 * cf
    u_tau = math.sqrt(tau_w / RHO)
    return yplus * NU / u_tau


H1 = first_layer_height()   # 1.456e-5 m for y+ = 0.7


def x_le(y):
    return -LE_AMP * math.cos(2 * math.pi * y / LE_WAVELENGTH)


def local_chord(y):
    return C - x_le(min(y, B))


def stations():
    """Peaks (LE most forward) and troughs of the LE wave, root to tip."""
    n = int(round(2 * B / LE_WAVELENGTH))
    return [(i * LE_WAVELENGTH / 2, "peak" if i % 2 == 0 else "trough")
            for i in range(n + 1)]
