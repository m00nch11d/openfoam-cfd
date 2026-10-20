"""
Single source of truth for the plunging NACA 0009 case.

The airfoil (c = 0.30 m, AoA = 0) plunges as  y(t) = A sin(2 pi f t)
with A = 1 chord and f = 3 Hz in a uniform stream (Re = 3e5).
"""
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

# ------------------------------------------------------------------ air
RHO = 1.225            # ISA sea level
MU = 1.7894e-5
NU = MU / RHO          # 1.4607e-5 m2/s

# ------------------------------------------------------------------ flow / motion
RE = 3.0e5
AOA_DEG = 0.0
C = 0.30               # chord [m]
U_INF = RE * NU / C    # 14.607 m/s
AMP = 1.0 * C          # plunge amplitude [m]
FREQ = 3.0             # plunge frequency [Hz]
OMEGA = 2.0 * math.pi * FREQ
PERIOD = 1.0 / FREQ
V_PLUNGE = AMP * OMEGA                                   # peak plunge velocity
ALPHA_EFF_MAX = math.degrees(math.atan(V_PLUNGE / U_INF))  # peak effective AoA
K_RED = math.pi * FREQ * C / U_INF                       # reduced frequency
STROUHAL = 2.0 * FREQ * AMP / U_INF
TU = 0.001             # free-stream turbulence intensity (0.1 %)
VISC_RATIO = 1.0

# ------------------------------------------------------------------ geometry
X_LE = 0.0
X_TE = C
R_FAR = 20.0 * C       # C-domain arc radius (centred on the TE)
L_WAKE = 30.0 * C      # wake length behind the TE
X_OUT = X_TE + L_WAKE

# ------------------------------------------------------------------ mesh
N_SURF = 200           # elements on EACH of upper / lower surface
SINE_WEIGHT = 0.80     # 1 = pure cosine (sine) clustering, 0 = uniform
YPLUS_TARGET = 0.7
BL_THICK = 0.04        # structured quad layer thickness [m]
BL_RATIO = 1.10        # wall-normal growth ratio
BL_SPREAD = 12.0       # smoothing of the outer-row spacing of the quad layer
WAKE_CUT_START = 2.0   # first wake-cut cell behind the TE / TE surface element
WAKE_CUT_RATIO = 1.15  # growth along the wake cut
FAR_SIZE = 0.40        # far-boundary element size (>= chord) [m]
GROWTH = 0.08          # size increase per metre of distance (outer region)
N_WAKE = 45            # cosine-distributed wake nodes over 30 c
# the airfoil sweeps +-A, so the near wake is refined in a band around it
WAKE_BAND = AMP + 0.10         # max half-height of the refined wake band [m]
WAKE_BAND_SLOPE = 0.5          # band half-height growth per metre behind the TE
WAKE_BAND_LEN = 8.0 * C        # length of the refined band behind the TE [m]
WAKE_BAND_H0 = 0.005           # cell size in the band at the TE [m]
WAKE_BAND_GROWTH = 0.012       # size growth per metre inside the band
SPAN = 0.01            # extrusion depth (1 cell) [m]


def first_cell_height(u_ref, yplus=YPLUS_TARGET):
    """Flat-plate estimate of the wall spacing giving the requested y+."""
    re_c = u_ref * C / NU
    cf = 0.058 * re_c ** -0.2
    u_tau = u_ref * math.sqrt(cf / 2.0)
    return yplus * NU / u_tau, re_c, cf, u_tau


# The relevant velocity is the peak relative velocity (free stream + plunge).
U_REL_MAX = math.hypot(U_INF, V_PLUNGE)
H1_EST, RE_REL, CF_EST, UTAU_EST = first_cell_height(U_REL_MAX)
# Leading-edge suction peaks raise the local wall shear far above the flat
# plate value. The factor below was set from the peak y+ over a full cycle
# in the first run (see README): it brings the cycle maximum to ~0.7.
H1_FACTOR = 0.33
H1 = H1_EST * H1_FACTOR

K_INF = 1.5 * (TU * U_INF) ** 2
OMEGA_INF = K_INF / (NU * VISC_RATIO)
NUT_INF = NU * VISC_RATIO
_tu = TU * 100.0
RETHETA_INF = (1173.51 - 589.428 * _tu + 0.2196 / _tu ** 2
               if _tu <= 1.3 else 331.5 * (_tu - 0.5658) ** -0.671)

# ------------------------------------------------------------------ solver / run control
N_PROCS = 8
STEADY_ITER = 400                 # simpleFoam initialisation (static airfoil)
STEPS_PER_CYCLE = 2000            # time steps per plunge period
DT = PERIOD / STEPS_PER_CYCLE
MAX_CYCLES = 5                    # hard limit
MIN_CYCLES = 3                    # the first cycle is a start-up transient
# PIMPLE: 6 outer correctors per time step make each step (nearly) converged;
# with only 3 the run diverged when the dynamic-stall vortex reached the TE
N_OUTER = 6
RELAX_P = 0.5                     # under-relaxation of the non-final outer iterations
RELAX_EQ = 0.8
SNAPSHOTS_PER_CYCLE = 40          # fields written per cycle (last cycle kept)
# periodic convergence: the cycle statistics (Cl rms, Cl max, Cl min, mean Cd)
# change by less than CYCLE_TOL between two consecutive cycles
CYCLE_TOL = float(os.environ.get("CYCLE_TOL", "0.10"))
NUT_RATIO_TR = 1.0                # transition detection (post-processing)
TR_BAND = 5.0e-4


RAMP_TIME = 0.5 * PERIOD          # smooth start of the motion (first half cycle)


def ramp(t):
    """0 -> 1 over RAMP_TIME with zero slope at both ends."""
    if t >= RAMP_TIME:
        return 1.0, 0.0
    s = math.pi * t / RAMP_TIME
    return 0.5 * (1.0 - math.cos(s)), 0.5 * math.pi / RAMP_TIME * math.sin(s)


def plunge(t):
    """Displacement and velocity of the airfoil at time t:
    y = A r(t) sin(omega t), r = smooth ramp (r = 1 after the first half cycle)."""
    r, dr = ramp(t)
    return (AMP * r * math.sin(OMEGA * t),
            AMP * (dr * math.sin(OMEGA * t) + r * OMEGA * math.cos(OMEGA * t)))


def summary():
    return "\n".join([
        "Air (ISA sea level): rho=%.4f kg/m3, nu=%.4e m2/s" % (RHO, NU),
        "Re=%.0f (c=%.2f m) -> U_inf=%.4f m/s, AoA=%.1f deg" % (RE, C, U_INF, AOA_DEG),
        "Plunge: A=%.3f m (%.2f c), f=%.2f Hz, T=%.4f s, V_max=%.3f m/s, alpha_eff,max=%.1f deg,"
        " k=%.3f, St=%.3f" % (AMP, AMP / C, FREQ, PERIOD, V_PLUNGE, ALPHA_EFF_MAX, K_RED, STROUHAL),
        "First cell: flat plate at U_rel,max=%.3f m/s (Re=%.0f, Cf=%.5f, u_tau=%.4f) -> %.3e m;"
        " used %.3e m (x%.2f for LE suction peaks)" % (U_REL_MAX, RE_REL, CF_EST, UTAU_EST,
                                                     H1_EST, H1, H1_FACTOR),
        "Time step %.3e s (%d per cycle), %d-%d cycles" % (DT, STEPS_PER_CYCLE, MIN_CYCLES,
                                                         MAX_CYCLES),
    ])


def label(n_cells=None):
    cells = "%d cells" % n_cells if n_cells else ""
    return ("NACA 0009, c=%.0f cm, %d+%d surface elements | %s | plunging y = %.2f m sin(2pi %.0f Hz t)"
            " (alpha_eff = +-%.1f deg)\nAir ISA sea level, Re=%.0f, U=%.2f m/s, AoA=%.0f deg, Tu=%.1f %%"
            " | k-omega SST + gamma-ReTheta (kOmegaSSTLM), pimpleFoam"
            % (C * 100, N_SURF, N_SURF, cells, AMP, FREQ, ALPHA_EFF_MAX, RE, U_INF, AOA_DEG,
               TU * 100))


if __name__ == "__main__":
    print(summary())
