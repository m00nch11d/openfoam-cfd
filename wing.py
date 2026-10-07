"""
CadQuery model of a tapered, swept S1223 wing with carbon-fibre skin,
foam core and a flaperon on each side.

Units: mm.  Axes: x = chordwise (aft +), y = spanwise (right wing +), z = up.
Origin: leading edge of the root chord.
"""
import math

import cadquery as cq
import numpy as np
from scipy.interpolate import splev, splprep
from shapely.geometry import LineString, Polygon, box
from shapely.geometry.polygon import orient

# ---------------------------------------------------------------- parameters
AIRFOIL_FILE = "s1223.dat"
ROOT_CHORD = 300.0           # mm
SPAN = 2000.0                # mm (tip to tip)
TAPER = 0.2                  # tip chord / root chord
SWEEP_C4_DEG = 25.0          # quarter-chord sweep
SKIN_T = 0.30                # carbon skin thickness (2 plies)
TE_GAP = 0.003               # blunt trailing edge, fraction of local chord
MIN_CORE_T = 1.0             # foam core truncated where thinner than this (TE)

FLAP_Y_START = 100.0         # flaperon starts 10 cm from centreline
FLAP_Y_END = SPAN / 2        # runs to the tip
FLAP_CHORD_FRAC = 0.10       # flaperon chord / local chord
FLAP_GAP = 1.0               # hinge-line and inboard-end clearance

HALF_SPAN = SPAN / 2
TIP_CHORD = ROOT_CHORD * TAPER
TAN_C4 = math.tan(math.radians(SWEEP_C4_DEG))


def chord(y):
    return ROOT_CHORD + (TIP_CHORD - ROOT_CHORD) * abs(y) / HALF_SPAN


def x_le(y):
    """Leading-edge x so that the quarter-chord line is swept SWEEP_C4_DEG."""
    return 0.25 * ROOT_CHORD + abs(y) * TAN_C4 - 0.25 * chord(y)


def x_hinge(y):
    return x_le(y) + (1.0 - FLAP_CHORD_FRAC) * chord(y)


# ------------------------------------------------------------------- airfoil
def load_airfoil(path=AIRFOIL_FILE):
    """Selig-format points (TE -> upper -> LE -> lower -> TE), TE opened."""
    pts = np.loadtxt(path, skiprows=1)
    i_le = int(np.argmin(pts[:, 0]))
    pts = pts.copy()
    pts[:i_le, 1] += pts[:i_le, 0] * TE_GAP / 2      # upper surface
    pts[i_le + 1:, 1] -= pts[i_le + 1:, 0] * TE_GAP / 2  # lower surface
    return pts


AIRFOIL = load_airfoil()


def dense_airfoil(n=1200):
    tck, _ = splprep([AIRFOIL[:, 0], AIRFOIL[:, 1]], s=0, k=3)
    u = (1 - np.cos(np.linspace(0, np.pi, n))) / 2
    x, z = splev(u, tck)
    return np.column_stack([x, z])


AIRFOIL_DENSE = dense_airfoil()


def section_wire(pts2d, y):
    """Closed wire: spline through (x, z) points + straight blunt TE, at span y."""
    vecs = [cq.Vector(float(x), float(y), float(z)) for x, z in pts2d]
    spline = cq.Edge.makeSpline(vecs)
    te = cq.Edge.makeLine(vecs[-1], vecs[0])
    return cq.Wire.assembleEdges([spline, te])


def outer_section(y):
    c = chord(y)
    pts = AIRFOIL * c
    pts[:, 0] += x_le(y)
    return section_wire(pts, y)


def _resample(line, n):
    """n points along a polyline, clustered at both ends (cosine spacing)."""
    ls = LineString(line)
    s = (1 - np.cos(np.linspace(0, np.pi, n))) / 2 * ls.length
    return np.array([ls.interpolate(d).coords[0] for d in s])


def inner_points(y, t=SKIN_T, n=70):
    """Foam-core profile: outer airfoil offset inward by t, TE truncated."""
    c = chord(y)
    poly = Polygon(AIRFOIL_DENSE * c)
    inner = poly.buffer(-t, join_style="round", quad_segs=32)
    if inner.geom_type == "MultiPolygon":
        inner = max(inner.geoms, key=lambda g: g.area)

    # truncate the core where it gets thinner than MIN_CORE_T
    xs = np.linspace(inner.bounds[0] + 0.5 * c, inner.bounds[2], 400)
    x_cut = xs[0]
    for x in xs:
        seg = inner.intersection(LineString([(x, -c), (x, c)]))
        if seg.is_empty or seg.length < MIN_CORE_T:
            break
        x_cut = x
    core = inner.intersection(box(-c, -c, x_cut, c))
    if core.geom_type == "MultiPolygon":
        core = max(core.geoms, key=lambda g: g.area)

    # CCW ring starting at the upper TE corner: upper TE -> LE -> lower TE
    ring = np.array(orient(core, 1.0).exterior.coords)[:-1]
    te = np.where(np.abs(ring[:, 0] - x_cut) < 1e-6)[0]
    ring = np.roll(ring, -te[np.argmax(ring[te, 1])], axis=0)
    i_le = int(np.argmin(ring[:, 0]))
    upper = ring[: i_le + 1]
    lower = ring[i_le:]
    up = _resample(upper, n)
    lo = _resample(lower, n)
    pts = np.vstack([up, lo[1:]])
    pts[:, 0] += x_le(y)
    return pts


def inner_section(y):
    return section_wire(inner_points(y), y)


# ------------------------------------------------------------------- solids
def loft(wires):
    return cq.Solid.makeLoft(wires, ruled=True)


# Spanwise loft stations. The outer surface is exact with root+tip alone, but
# shorter faces mesh far more reliably (BRepMesh leaves slits in the STL of a
# single 1 m long trimmed spline face).
STATIONS = [0.0, 50.0, FLAP_Y_START, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0,
            800.0, 900.0, HALF_SPAN]


def half_outer():
    # linear taper + straight sweep => ruled lofts between stations are exact
    return loft([outer_section(y) for y in STATIONS])


def half_core():
    ys = STATIONS[:-1] + [HALF_SPAN - SKIN_T]   # leave a carbon tip closure
    return loft([inner_section(y) for y in ys])


def flap_region(y0, hinge_shift):
    """Prism aft of the (shifted) hinge line, from span y0 past the tip."""
    y1 = HALF_SPAN + 50.0
    pts = [
        (x_hinge(y0) + hinge_shift, y0),
        (x_hinge(y1) + hinge_shift, y1),
        (x_le(y1) + 2 * ROOT_CHORD, y1),
        (x_le(y0) + 2 * ROOT_CHORD, y0),
    ]
    return (
        cq.Workplane("XY").workplane(offset=-100)
        .polyline(pts).close().extrude(200).val()
    )


def mirror(s):
    return s.mirror("XZ")


def both(s):
    return s.fuse(mirror(s)).clean()


def build():
    t, g = SKIN_T, FLAP_GAP
    outer_h = half_outer()
    core_h = half_core()

    # main wing: remove the flaperon bay (with clearance gap)
    cut_main_outer = flap_region(FLAP_Y_START - g, -g)
    cut_main_core = flap_region(FLAP_Y_START - g - t, -g - t)
    # flaperon: region aft of the hinge line, core inset by skin thickness
    flap_outer_reg = flap_region(FLAP_Y_START, 0.0)
    flap_core_reg = flap_region(FLAP_Y_START + t, t)

    main_outer_h = outer_h.cut(cut_main_outer)
    main_core_h = core_h.cut(cut_main_core)
    flap_outer_h = outer_h.intersect(flap_outer_reg)
    flap_core_h = core_h.intersect(flap_core_reg)

    main_outer = both(main_outer_h)
    main_core = both(main_core_h)
    parts = {
        "main_skin": main_outer.cut(main_core),
        "main_core": main_core,
        "flaperon_R_skin": flap_outer_h.cut(flap_core_h),
        "flaperon_R_core": flap_core_h,
        "flaperon_L_skin": mirror(flap_outer_h.cut(flap_core_h)),
        "flaperon_L_core": mirror(flap_core_h),
    }
    outer_shapes = {
        "main": main_outer,
        "flaperon_R": flap_outer_h,
        "flaperon_L": mirror(flap_outer_h),
    }
    return parts, outer_shapes


def volume(shape):
    """Accurate volume (mm^3). Shape.Volume()'s default integration is too
    coarse for the many-knot airfoil splines and under-reports by ~35 %."""
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps

    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape.wrapped, props, 1e-6, False, False)
    return props.Mass()
