"""Shared case parameters (SI units).

Axes: x streamwise (aft +), y spanwise (root y=0 -> tip y=+B), z up.
Mean-chord leading edge at x=0, trailing edge (straight) at x=C.
"""
import math

# ---------------------------------------------------------------- geometry
AIRFOIL_FILE = "../s1223.dat"
C = 0.30                 # mean chord [m]
B = 1.00                 # half span [m]
LE_AMP = 0.02 * C        # leading-edge sine amplitude (+-2 % chord) [m]
LE_WAVELENGTH = 0.125    # leading-edge sine wavelength [m] (0.417 c, 8 waves)
# x_LE(y) = -LE_AMP*cos(2 pi y / lambda): a peak (LE most forward, longest
# chord) sits on the symmetry plane so the LE meets it with zero slope.
# Sections are the S1223 scaled about the fixed trailing edge.

# ---------------------------------------------------------------- flow
RE = 3.0e5
RHO = 1.225              # ISA sea level [kg/m3]
MU = 1.7894e-5           # ISA sea level [Pa s]
NU = MU / RHO            # 1.4607e-5 m2/s
U_INF = RE * NU / C      # 14.61 m/s
AOA_DEG = 0.0            # angle of attack (not specified -> 0 deg)
TU = 0.002               # free-stream turbulence intensity at inlet (0.2 %)
VISC_RATIO = 5.0         # free-stream nut/nu

# ---------------------------------------------------------------- mesh
YPLUS_TARGET = 0.7
N_SURF = 200             # cells on each of upper / lower surface
SINE_BLEND = 0.7         # surface spacing = 0.7*sine + 0.3*uniform
N_WAKE = 120             # cells along the 30c wake cut (half-cosine)
WAKE_LEN = 30.0 * C      # wake / domain length behind the TE
QUAD_THICK = 0.05        # structured quad layer thickness [m]
N_LAYERS = 40           # quad layers across QUAD_THICK
R_FAR = 20.0 * C         # C-domain radius (centred at TE)
FAR_SIZE = 0.50          # far-field element size [m] (>= chord)
GROWTH = 0.15            # triangle size gradient

DY_WING = 0.015625      # spanwise cell size on the wing [m] (8 per wave)
N_TIP = 14              # spanwise layers outboard of the tip
TIP_RATIO = 1.28        # growth of those layers


def first_layer_height(yplus=YPLUS_TARGET):
    """Flat-plate estimate (Schlichting: Cf = 0.0576 Re^-0.2)."""
    cf = 0.0576 * RE ** -0.2
    tau_w = 0.5 * RHO * U_INF ** 2 * cf
    u_tau = math.sqrt(tau_w / RHO)
    return yplus * NU / u_tau


H1 = first_layer_height()


def x_le(y):
    return -LE_AMP * math.cos(2 * math.pi * y / LE_WAVELENGTH)


def local_chord(y):
    return C - x_le(min(y, B))


def stations():
    """Peaks (LE most forward) and troughs of the LE wave, root to tip."""
    n = int(round(2 * B / LE_WAVELENGTH))
    out = []
    for i in range(n + 1):
        y = i * LE_WAVELENGTH / 2
        out.append((y, "peak" if i % 2 == 0 else "trough"))
    return out
