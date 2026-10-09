"""
Single source of truth for the tandem-airfoil case (S1223 + NACA 0009).

Everything else (mesh, OpenFOAM dictionaries, post-processing) imports
these values, so a change here propagates through the whole pipeline.
"""
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

# ------------------------------------------------------------------ air
# ISA standard atmosphere, sea level
RHO = 1.225            # kg/m3
MU = 1.7894e-5         # Pa s
NU = MU / RHO          # m2/s  (= 1.4607e-5)
T_INF = 288.15         # K
P_INF = 101325.0       # Pa

# ------------------------------------------------------------------ flow
RE = 3.0e5             # based on the main (S1223) chord
AOA_DEG = 0.0
C1 = 0.30              # S1223 chord [m]
C2 = 0.10              # NACA 0009 chord [m]
GAP = 0.20             # S1223 TE -> NACA LE [m]
U_INF = RE * NU / C1   # 14.607 m/s
TU = 0.001             # free-stream turbulence intensity (0.1 %, low-turbulence tunnel)
VISC_RATIO = 1.0       # nut/nu in the free stream

# ------------------------------------------------------------------ geometry
X_LE1 = 0.0
X_TE1 = X_LE1 + C1
X_LE2 = X_TE1 + GAP
X_TE2 = X_LE2 + C2
R_FAR = 20.0 * C1      # radius of the C-domain arc (centred on the S1223 TE)
L_WAKE = 30.0 * C1     # wake length behind the S1223 trailing edge
X_OUT = X_TE1 + L_WAKE

# ------------------------------------------------------------------ mesh
N_SURF1 = 200          # elements on EACH of upper / lower surface, S1223
N_SURF2 = 100          # elements on EACH of upper / lower surface, NACA 0009
SINE_WEIGHT = 0.90     # 1 = pure cosine (sine) clustering, 0 = uniform
YPLUS_TARGET = 0.7
BL_THICK = 0.04        # quad (structured) layer thickness around each airfoil [m]
BL_RATIO = 1.10        # wall-normal growth ratio in the quad layer
BL_SPREAD = 12.0       # smoothing (in nodes) of the outer-row spacing of the quad layer
FAR_SIZE = 0.40        # element size on the far boundary (>= chord) [m]
GROWTH = 0.08          # size increase per metre of distance from the walls
N_WAKE = 45            # cosine-distributed nodes along the 30c wake line
SPAN = 0.01            # extrusion depth of the 2-D mesh (1 cell) [m]

# ------------------------------------------------------------------ solver / run control
N_PROCS = 8            # MPI ranks
MAX_ITER = 1000
RELAX_P = 1.0          # SIMPLEC (consistent) -> no pressure under-relaxation
RELAX_U = 0.9
RELAX_TURB = 0.8
RELAX_TRANS = 0.8
# convergence: relative Cl change between two consecutive iterations must stay
# below CL_TOL for CL_WINDOW consecutive iterations (after MIN_ITER)
CL_TOL = float(os.environ.get("CL_TOL", "1.0e-4"))   # override: ./Allrun.sh -cltol X
CL_WINDOW = 50
MIN_ITER = 200
NUT_RATIO_TR = 1.0     # transition onset criterion in post-processing ...
TR_BAND = 5.0e-4       # ... evaluated in cells closer than this to the wall [m]


def first_cell_height(chord, yplus=YPLUS_TARGET):
    """Flat-plate estimate of the wall spacing giving the requested y+."""
    re_c = U_INF * chord / NU
    cf = 0.058 * re_c ** -0.2          # Schlichting turbulent skin friction
    tau_w = 0.5 * RHO * U_INF ** 2 * cf
    u_tau = math.sqrt(tau_w / RHO)
    return yplus * NU / u_tau, re_c, cf, u_tau


H1_1_EST, RE_1, CF_1, UTAU_1 = first_cell_height(C1)
H1_2_EST, RE_2, CF_2, UTAU_2 = first_cell_height(C2)
# The flat-plate estimate (y+ = 0.7) gave peak y+ of 1.23 (S1223) and 1.42
# (NACA 0009) at the leading-edge suction peaks in the first run. The wall
# spacing is therefore scaled by 0.7 / peak so that y+ < 1 everywhere with
# the maximum close to the 0.7 target.
YPLUS_PEAK_RUN1 = (1.228, 1.425)
H1_1 = H1_1_EST * YPLUS_TARGET / YPLUS_PEAK_RUN1[0]
H1_2 = H1_2_EST * YPLUS_TARGET / YPLUS_PEAK_RUN1[1]

# free-stream turbulence quantities
K_INF = 1.5 * (TU * U_INF) ** 2
OMEGA_INF = K_INF / (NU * VISC_RATIO)
NUT_INF = NU * VISC_RATIO
_tu_pct = TU * 100.0
RETHETA_INF = (1173.51 - 589.428 * _tu_pct + 0.2196 / _tu_pct ** 2
               if _tu_pct <= 1.3 else 331.5 * (_tu_pct - 0.5658) ** -0.671)


def summary():
    lines = [
        "Air (ISA sea level): rho=%.4f kg/m3, mu=%.4e Pa s, nu=%.4e m2/s"
        % (RHO, MU, NU),
        "Re=%.0f (c=%.2f m)  ->  U_inf=%.4f m/s, AoA=%.1f deg"
        % (RE, C1, U_INF, AOA_DEG),
        "Tu=%.2f %%, k=%.3e m2/s2, omega=%.3e 1/s, ReTheta_t=%.1f"
        % (TU * 100, K_INF, OMEGA_INF, RETHETA_INF),
        "S1223   : Re_c=%.0f Cf=%.5f u_tau=%.4f m/s -> estimate %.3e m (y+=%.2f), used %.3e m"
        % (RE_1, CF_1, UTAU_1, H1_1_EST, YPLUS_TARGET, H1_1),
        "NACA0009: Re_c=%.0f Cf=%.5f u_tau=%.4f m/s -> estimate %.3e m (y+=%.2f), used %.3e m"
        % (RE_2, CF_2, UTAU_2, H1_2_EST, YPLUS_TARGET, H1_2),
    ]
    return "\n".join(lines)


def label(n_cells=None):
    """Annotation text stamped on every figure."""
    cells = "%d cells" % n_cells if n_cells else ""
    return (
        "S1223 (c=%.0f cm, %d+%d surf. el.) + NACA 0009 (c=%.0f cm, %d+%d surf. el.), "
        "gap %.0f cm | %s\n"
        "Air ISA sea level, Re=%.0f, U=%.2f m/s, AoA=%.0f deg, Tu=%.1f %% | "
        "k-omega SST + gamma-ReTheta (kOmegaSSTLM)"
        % (C1 * 100, N_SURF1, N_SURF1, C2 * 100, N_SURF2, N_SURF2, GAP * 100,
           cells, RE, U_INF, AOA_DEG, TU * 100)
    )


if __name__ == "__main__":
    print(summary())
