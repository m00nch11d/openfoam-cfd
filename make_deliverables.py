"""
Builds the wing (wing.py) and writes the deliverables to ./output:
  * STL (full assembly + one file per body) and a STEP assembly
  * wing_3view_drawing.png  - third-angle style 3-view with dimensions
  * wing_materials.png      - material specification / mass breakdown
"""
import math
import os

import cadquery as cq
import matplotlib
import numpy as np
import trimesh

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Arc, Rectangle  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402
from OCP.BRepLib import BRepLib  # noqa: E402
from OCP.BRepTools import BRepTools  # noqa: E402
from OCP.GCPnts import GCPnts_QuasiUniformDeflection  # noqa: E402
from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt  # noqa: E402
from OCP.HLRAlgo import HLRAlgo_Projector  # noqa: E402
from OCP.HLRBRep import HLRBRep_Algo, HLRBRep_HLRToShape  # noqa: E402

import wing as W  # noqa: E402

OUT = "output"
os.makedirs(OUT, exist_ok=True)

# --------------------------------------------------------------- materials
# Typical values - confirm against the supplier datasheet before sizing.
CARBON = dict(
    name="Carbon fibre / epoxy (woven)",
    spec="2 plies 3K 2x2 twill, 160 g/m2, [0/90 , +-45]",
    density=1.50,  # g/cm3, Vf ~ 50 %
    rows=[
        ("Density", "1500 kg/m3"),
        ("Fibre volume fraction", "~50 %"),
        ("Tensile modulus E1 = E2", "60 GPa"),
        ("Shear modulus G12", "5 GPa"),
        ("Tensile strength", "600 MPa"),
        ("Compressive strength", "500 MPa"),
        ("Poisson's ratio v12", "0.05"),
        ("Cured ply thickness", "0.15 mm"),
        ("Laminate thickness", f"{W.SKIN_T:.2f} mm"),
        ("Max service temp.", "~80 C (RT-cure epoxy)"),
    ],
)
FOAM = dict(
    name="PMI structural foam (Rohacell 51 IG-F)",
    spec="CNC / hot-wire cut core, closed cell",
    density=0.052,  # g/cm3
    rows=[
        ("Density", "52 kg/m3"),
        ("Compressive strength", "0.9 MPa"),
        ("Tensile strength", "1.9 MPa"),
        ("Tensile modulus", "70 MPa"),
        ("Shear strength", "0.8 MPa"),
        ("Shear modulus", "19 MPa"),
        ("Max service temp.", "~130 C"),
        ("Alternative", "XPS foam, 33 kg/m3"),
    ],
)
XPS_DENSITY = 0.033

C_CARBON = "#2b2b2b"
C_FOAM = "#f2c14e"
C_DIM = "#1f4e9c"


# ------------------------------------------------------------------ export
def open_edges(path):
    """Number of mesh edges used by only one triangle (0 = watertight)."""
    m = trimesh.load(path, process=False)
    _, inv = np.unique(m.vertices.round(9), axis=0, return_inverse=True)
    f = inv.reshape(-1)[m.faces]
    e = np.sort(np.vstack([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]]), axis=1)
    _, n = np.unique(e, axis=0, return_counts=True)
    return int(np.sum(n == 1))


def export_stl(shape, path):
    """BRepMesh occasionally leaves slits on the trimmed spline skin for some
    tolerance combinations - retry until the STL is watertight and its volume
    matches the CAD solid."""
    v_cad = W.volume(shape)
    for tol, ang in ((0.01, 0.1), (0.02, 0.1), (0.005, 0.1), (0.03, 0.15)):
        BRepTools.Clean_s(shape.wrapped)   # drop cached triangulation
        cq.exporters.export(shape, path, tolerance=tol, angularTolerance=ang)
        v_stl = trimesh.load(path, process=False).volume
        if open_edges(path) == 0 and abs(v_stl / v_cad - 1) < 0.02:
            return tol
    raise RuntimeError(f"could not mesh {path} watertight")


def export(parts, outer):
    for name, shp in parts.items():
        tol = export_stl(shp, f"{OUT}/wing_{name}.stl")
        print(f"  wing_{name}.stl  watertight (tol {tol} mm)")
    for name, shp in outer.items():
        export_stl(shp, f"{OUT}/_oml_{name}.stl")
    combine = lambda files, dst: trimesh.util.concatenate(  # noqa: E731
        [trimesh.load(f, process=False) for f in files]).export(dst)
    combine([f"{OUT}/wing_{n}.stl" for n in parts], f"{OUT}/wing_assembly_all_bodies.stl")
    oml = [f"{OUT}/_oml_{n}.stl" for n in outer]
    combine(oml, f"{OUT}/wing_outer_mold_line.stl")
    for f in oml:
        os.remove(f)

    assy = cq.Assembly(name="S1223_wing")
    for name, shp in parts.items():
        col = cq.Color(0.17, 0.17, 0.17) if "skin" in name else cq.Color(0.95, 0.76, 0.31)
        assy.add(shp, name=name, color=col)
    assy.export(f"{OUT}/wing_assembly.step")


# --------------------------------------------------------------------- HLR
def _polyline(edge, defl=0.05):
    curve = edge._geomAdaptor()
    pts = GCPnts_QuasiUniformDeflection(curve, defl, curve.FirstParameter(),
                                        curve.LastParameter())
    return np.array([(pts.Value(i).X(), pts.Value(i).Y())
                     for i in range(1, pts.NbPoints() + 1)])


def hlr(shape, view_dir, x_dir):
    """Visible / hidden 2D polylines of `shape` seen from `view_dir`."""
    algo = HLRBRep_Algo()
    algo.Add(shape.wrapped)
    algo.Projector(HLRAlgo_Projector(gp_Ax2(gp_Pnt(), gp_Dir(*view_dir), gp_Dir(*x_dir))))
    algo.Update()
    algo.Hide()
    h = HLRBRep_HLRToShape(algo)

    def collect(comps):
        lines = []
        for c in comps:
            if c.IsNull():
                continue
            BRepLib.BuildCurves3d_s(c, 1e-6)
            lines += [_polyline(e) for e in cq.Shape.cast(c).Edges()]
        return lines

    vis = collect([h.VCompound(), h.Rg1LineVCompound(), h.OutLineVCompound()])
    hid = collect([h.HCompound(), h.OutLineHCompound()])
    return vis, hid


def draw_lines(ax, lines, origin, hidden=False):
    ox, oy = origin
    for p in lines:
        if hidden:
            ax.plot(p[:, 0] + ox, p[:, 1] + oy, color="0.55", lw=0.45, ls=(0, (4, 3)))
        else:
            ax.plot(p[:, 0] + ox, p[:, 1] + oy, color="black", lw=0.8)


# -------------------------------------------------------------- dimensions
def dim(ax, a, b, orient, at, text, fs=8.5, gap=4, over=8, color=C_DIM):
    """Linear dimension between points a and b.
    orient 'h': horizontal dim line at y=at ; 'v': vertical dim line at x=at."""
    (ax_, ay), (bx, by) = a, b
    kw = dict(color=color, lw=0.5)
    if orient == "h":
        for x, y in ((ax_, ay), (bx, by)):
            s = 1 if at > y else -1
            ax.plot([x, x], [y + s * gap, at + s * over], **kw)
        q1, q2 = (ax_, at), (bx, at)
        rot = 0
    else:
        for x, y in ((ax_, ay), (bx, by)):
            s = 1 if at > x else -1
            ax.plot([x + s * gap, at + s * over], [y, y], **kw)
        q1, q2 = (at, ay), (at, by)
        rot = 90
    ax.annotate("", q1, q2, arrowprops=dict(arrowstyle="<|-|>", color=color, lw=0.6,
                                            shrinkA=0, shrinkB=0, mutation_scale=7))
    mid = ((q1[0] + q2[0]) / 2, (q1[1] + q2[1]) / 2)
    ax.text(*mid, text, rotation=rot, ha="center", va="center", fontsize=fs, color=color,
            bbox=dict(fc="white", ec="none", pad=0.6))


def note(ax, xy, xytext, text, fs=8, color=C_DIM, ha="left"):
    ax.annotate(text, xy, xytext, fontsize=fs, color=color, ha=ha, va="center",
                arrowprops=dict(arrowstyle="-|>", color=color, lw=0.5, mutation_scale=7))


def angle(ax, center, r, th0, th1, text, text_r=None, fs=8.5):
    ax.add_patch(Arc(center, 2 * r, 2 * r, theta1=th0, theta2=th1, color=C_DIM, lw=0.6))
    tr = text_r or r + 25
    tm = math.radians((th0 + th1) / 2)
    ax.text(center[0] + tr * math.cos(tm), center[1] + tr * math.sin(tm), text,
            fontsize=fs, color=C_DIM, ha="left", va="center")


def view_label(ax, x, y, text, sub=None):
    ax.text(x, y, text, fontsize=11, fontweight="bold", ha="center")
    if sub:
        ax.text(x, y - 32, sub, fontsize=8, ha="center", color="0.3")


# ------------------------------------------------------------- geometry data
def wing_data():
    cr, ct, b = W.ROOT_CHORD, W.TIP_CHORD, W.SPAN
    lam = W.TAPER
    S = b * (cr + ct) / 2
    mac = 2 / 3 * cr * (1 + lam + lam ** 2) / (1 + lam)
    y_mac = b / 6 * (1 + 2 * lam) / (1 + lam)
    le_sweep = math.degrees(math.atan((W.x_le(W.HALF_SPAN) - W.x_le(0)) / W.HALF_SPAN))
    hinge_sweep = math.degrees(math.atan((W.x_hinge(W.FLAP_Y_END) - W.x_hinge(W.FLAP_Y_START))
                                         / (W.FLAP_Y_END - W.FLAP_Y_START)))
    tc = 0.1214
    return dict(S=S, AR=b ** 2 / S, mac=mac, y_mac=y_mac, le_sweep=le_sweep,
                hinge_sweep=hinge_sweep, tc=tc)


def section(parts, y):
    """Cut every body with the plane y=const -> list of (name, outer, holes)."""
    plane = cq.Face.makePlane(4000, 4000, basePnt=(0, y, 0), dir=(0, 1, 0))
    out = []
    for name, shp in parts.items():
        cut = shp.intersect(plane)
        for f in cut.Faces():
            def pts(w):
                p = w.positions(np.linspace(0, 1, 2500).tolist())
                return np.array([(v.x, v.z) for v in p])
            out.append((name, pts(f.outerWire()), [pts(w) for w in f.innerWires()]))
    return out


# ----------------------------------------------------------------- drawing
def make_drawing(parts, outer, masses):
    D = wing_data()
    shape = cq.Compound.makeCompound(list(outer.values()))
    hb = W.HALF_SPAN
    xle_t, ct, cr = W.x_le(hb), W.TIP_CHORD, W.ROOT_CHORD
    zmax = float(W.AIRFOIL[:, 1].max() * cr)
    zmin = float(W.AIRFOIL[:, 1].min() * cr)

    fig = plt.figure(figsize=(23.4, 16.5), dpi=150)       # A2 landscape
    ax = fig.add_axes([0.01, 0.01, 0.98, 0.98])
    ax.set_xlim(-170, 3080)
    ax.set_ylim(-640, 1650)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.add_patch(Rectangle((-150, -620), 3210, 2250, fill=False, lw=1.6))

    # ---------------- TOP VIEW  (looking down, leading edge up): 2D = (y, -x)
    TO = (1000, 1380)
    vis, hid = hlr(shape, (0, 0, 1), (0, 1, 0))
    draw_lines(ax, hid, TO, hidden=True)
    draw_lines(ax, vis, TO)
    T = lambda y, x: (TO[0] + y, TO[1] - x)  # noqa: E731
    # centre line, quarter-chord line, MAC
    ax.plot([TO[0]] * 2, [TO[1] + 60, TO[1] - cr - 60], color="0.3", lw=0.5, ls="-.")
    for s in (-1, 1):
        ax.plot(*zip(T(0, 0.25 * cr), T(s * hb, W.x_le(hb) + 0.25 * ct)),
                color="tab:red", lw=0.6, ls="-.")
        ym = s * D["y_mac"]
        ax.plot(*zip(T(ym, W.x_le(ym)), T(ym, W.x_le(ym) + W.chord(ym))),
                color="tab:green", lw=1.0)
    ax.text(*T(-hb - 40, 0.25 * cr + W.HALF_SPAN * W.TAN_C4 - 40), "c/4 line", fontsize=8,
            color="tab:red", ha="right")
    note(ax, T(-D["y_mac"], W.x_le(D["y_mac"]) + 0.5 * D["mac"]),
         T(-D["y_mac"] - 330, W.x_le(D["y_mac"]) + 0.5 * D["mac"] + 210),
         f"MAC = {D['mac']:.1f}\nat y = {D['y_mac']:.1f}", color="tab:green")
    # overall span / semi-span
    dim(ax, T(-hb, xle_t), T(hb, xle_t), "h", TO[1] + 110, f"{W.SPAN:.0f}  (span)")
    dim(ax, T(0, 0), T(hb, xle_t), "h", TO[1] + 50, f"{hb:.0f}")
    # chords and sweep offsets
    dim(ax, T(0, 0), T(0, cr), "v", TO[0] - 70, f"{cr:.0f} (root chord)")
    dim(ax, T(hb, xle_t), T(hb, xle_t + ct), "v", TO[0] + hb + 60, f"{ct:.0f}")
    ax.text(TO[0] + hb + 105, TO[1] - xle_t - ct / 2, "tip chord", fontsize=8, color=C_DIM,
            va="center")
    dim(ax, T(0, 0), T(hb, xle_t), "v", TO[0] + hb + 230, f"{xle_t:.1f}")
    dim(ax, T(0, 0), T(hb, xle_t + ct), "v", TO[0] + hb + 330, f"{xle_t + ct:.1f}")
    # sweep angles
    c4 = T(0, 0.25 * cr)
    ax.plot([c4[0], c4[0] + 420], [c4[1]] * 2, color=C_DIM, lw=0.5)
    angle(ax, c4, 380, -W.SWEEP_C4_DEG, 0, f"$\\Lambda_{{c/4}}$ = {W.SWEEP_C4_DEG:.0f}°")
    le = T(0, 0)
    ax.plot([le[0], le[0] + 620], [le[1]] * 2, color=C_DIM, lw=0.5)
    angle(ax, le, 580, -D["le_sweep"], 0, f"$\\Lambda_{{LE}}$ = {D['le_sweep']:.2f}°")
    # flaperon
    yf0, yf1 = W.FLAP_Y_START, W.FLAP_Y_END
    te0 = W.x_le(yf0) + W.chord(yf0)
    yd = TO[1] - (xle_t + ct) - 60
    dim(ax, T(0, cr), T(yf0, te0), "h", yd, f"{yf0:.0f}")
    dim(ax, T(yf0, te0), T(yf1, xle_t + ct), "h", yd,
        f"{yf1 - yf0:.0f}  (flaperon span)")
    note(ax, T(-(yf0 + yf1) / 2, W.x_hinge((yf0 + yf1) / 2) + 0.05 * W.chord((yf0 + yf1) / 2)),
         T(-(yf0 + yf1) / 2 + 80, xle_t + ct + 150),
         "LEFT FLAPERON", ha="right")
    note(ax, T((yf0 + yf1) / 2, W.x_hinge((yf0 + yf1) / 2) + 0.05 * W.chord((yf0 + yf1) / 2)),
         T(150, xle_t + ct + 170),
         f"RIGHT FLAPERON\nchord = {W.FLAP_CHORD_FRAC:.0%} local chord\n"
         f"= {W.FLAP_CHORD_FRAC * W.chord(yf0):.1f} (inbd) to "
         f"{W.FLAP_CHORD_FRAC * ct:.1f} (tip)\n"
         f"hinge line at {1 - W.FLAP_CHORD_FRAC:.0%} c, sweep {D['hinge_sweep']:.2f}°\n"
         f"gap {W.FLAP_GAP:.0f} mm (hinge + inbd end)")
    # section marker A-A
    ys = 200.0
    a0, a1 = T(ys, W.x_le(ys) - 40), T(ys, W.x_le(ys) + W.chord(ys) + 40)
    ax.plot(*zip(a0, a1), color="k", lw=1.2, ls=(0, (10, 3, 2, 3)))
    for p, dy in ((a0, 22), (a1, -22)):
        ax.text(p[0] + 18, p[1] + dy, "A", fontsize=10, fontweight="bold")
    view_label(ax, TO[0] - 300, TO[1] - 760, "TOP VIEW", "planform")

    # ---------------- FRONT VIEW  (looking aft): 2D = (-y, z)
    FO = (1000, 300)
    vis, hid = hlr(shape, (-1, 0, 0), (0, -1, 0))
    draw_lines(ax, hid, FO, hidden=True)
    draw_lines(ax, vis, FO)
    ax.plot([FO[0]] * 2, [FO[1] - 60, FO[1] + 100], color="0.3", lw=0.5, ls="-.")
    dim(ax, (FO[0] - hb, FO[1]), (FO[0] + hb, FO[1]), "h", FO[1] - 90, f"{W.SPAN:.0f}")
    dim(ax, (FO[0], FO[1] + zmin), (FO[0], FO[1] + zmax), "v", FO[0] - 180,
        f"{zmax - zmin:.1f}")
    ax.text(FO[0] - 300, FO[1] + 75, "root height (camber + thickness)", fontsize=8,
            color=C_DIM, ha="center")
    ax.text(FO[0] + 650, FO[1] + 75, "dihedral 0°,  geometric twist 0°", fontsize=8.5,
            color=C_DIM)
    ax.text(FO[0] - hb - 10, FO[1] - 30, "R", fontsize=9, fontweight="bold", ha="right")
    ax.text(FO[0] + hb + 10, FO[1] - 30, "L", fontsize=9, fontweight="bold")
    view_label(ax, FO[0], FO[1] - 175, "FRONT VIEW", "looking aft")

    # ---------------- SIDE VIEW  (from left): 2D = (x, z)
    SO = (2330, 300)
    vis, hid = hlr(shape, (0, -1, 0), (1, 0, 0))
    draw_lines(ax, hid, SO, hidden=True)
    draw_lines(ax, vis, SO)
    dim(ax, (SO[0], SO[1]), (SO[0] + xle_t + ct, SO[1]), "h", SO[1] - 90,
        f"{xle_t + ct:.1f}")
    dim(ax, (SO[0], SO[1]), (SO[0] + cr, SO[1]), "h", SO[1] - 45, f"{cr:.0f}")
    dim(ax, (SO[0], SO[1] + zmin), (SO[0], SO[1] + zmax), "v", SO[0] - 50,
        f"{zmax - zmin:.1f}")
    view_label(ax, SO[0] + 290, SO[1] - 175, "SIDE VIEW", "from left")

    # ---------------- SECTION A-A  (scale 2:1 relative to views)
    k = 2.0
    AO = (2330, 1150)
    xl = W.x_le(ys)
    for name, ext, holes in section(parts, ys):
        col = C_CARBON if "skin" in name else C_FOAM
        ax.fill((ext[:, 0] - xl) * k + AO[0], ext[:, 1] * k + AO[1], color=col, lw=0)
        for h in holes:   # skin holes are filled by the core drawn on top
            ax.fill((h[:, 0] - xl) * k + AO[0], h[:, 1] * k + AO[1], color="white", lw=0)
    for name, ext, holes in section(parts, ys):
        if "core" in name:
            ax.fill((ext[:, 0] - xl) * k + AO[0], ext[:, 1] * k + AO[1], color=C_FOAM, lw=0)
    c = W.chord(ys)
    S = lambda x, z: (AO[0] + k * x, AO[1] + k * z)  # noqa: E731
    up = W.AIRFOIL[W.AIRFOIL[:, 0].argmin()::-1]
    lo = W.AIRFOIL[W.AIRFOIL[:, 0].argmin():]
    xt = 0.1985
    zt_u = np.interp(xt, up[:, 0], up[:, 1])
    zt_l = np.interp(xt, lo[:, 0], lo[:, 1])
    dim(ax, S(0, 0), S(c, 0), "h", AO[1] - 70, f"{c:.1f}  (local chord)")
    dim(ax, S(0.9 * c, 0), S(c, 0), "h", AO[1] - 35, f"{0.1 * c:.1f}")
    dim(ax, S(xt * c, zt_l * c), S(xt * c, zt_u * c), "v", AO[0] + k * xt * c,
        f"{(zt_u - zt_l) * c:.1f}", fs=7.5)
    ax.text(*S(xt * c + 8, (zt_u + zt_l) / 2 * c), f"  t/c = {D['tc']:.1%}\n  at x/c = {xt:.2f}",
            fontsize=7.5, color=C_DIM, va="center")
    note(ax, S(0.45 * c, float(np.interp(0.45, up[:, 0], up[:, 1])) * c),
         S(0.45 * c + 40, 0.20 * c), f"CFRP skin {W.SKIN_T:.1f} mm (2 plies)")
    note(ax, S(0.6 * c, 0.05 * c), S(0.6 * c + 35, -0.10 * c), "foam core")
    note(ax, S(0.9 * c, 0.03 * c), S(0.9 * c + 30, 0.13 * c),
         f"hinge gap {W.FLAP_GAP:.0f} mm")
    view_label(ax, AO[0] + 280, AO[1] - 140, "SECTION A-A  (2:1)",
               f"y = {ys:.0f} mm,  airfoil S1223")

    # ---------------- isometric pictorial
    iso = fig.add_axes([0.03, 0.035, 0.36, 0.25], projection="3d")
    iso.patch.set_alpha(0)
    for name, shp in outer.items():
        v, t = shp.tessellate(0.5, 0.3)
        v = np.array([p.toTuple() for p in v])
        tris = v[np.array(t)]
        col = "#3a3a3a" if name == "main" else "#c0392b"
        iso.add_collection3d(Poly3DCollection(tris[:, :, [1, 0, 2]] * [1, -1, 1], facecolors=col,
                                              lw=0, alpha=1, shade=True))
    iso.set_xlim(-1000, 1000)
    iso.set_ylim(-600, 0)
    iso.set_zlim(-100, 100)
    iso.set_box_aspect((2000, 600, 200), zoom=1.3)
    iso.view_init(elev=30, azim=-115)
    iso.axis("off")
    fig.text(0.17, 0.275, "ISOMETRIC (flaperons in red)", fontsize=10, fontweight="bold",
             ha="center")

    # ---------------- data table + title block
    rows = [
        ("Airfoil", "Selig S1223 (UIUC), TE thickness 0.3 % c"),
        ("Span  b", f"{W.SPAN:.0f} mm"),
        ("Root / tip chord", f"{cr:.0f} / {ct:.0f} mm"),
        ("Taper ratio", f"{W.TAPER}"),
        ("Sweep  c/4  /  LE", f"{W.SWEEP_C4_DEG:.0f}°  /  {D['le_sweep']:.2f}°"),
        ("Wing area  S", f"{D['S'] / 1e6:.3f} m²"),
        ("Aspect ratio", f"{D['AR']:.2f}"),
        ("MAC  /  y_MAC", f"{D['mac']:.1f} mm  /  {D['y_mac']:.1f} mm"),
        ("Flaperons", f"y = {W.FLAP_Y_START:.0f} to {W.FLAP_Y_END:.0f} mm, "
                      f"{W.FLAP_CHORD_FRAC:.0%} local chord"),
        ("Skin / core", f"CFRP {W.SKIN_T:.1f} mm  /  PMI foam 52 kg/m³"),
        ("Est. mass", f"{masses['total']:.0f} g  (see material sheet)"),
    ]
    x0, y0 = 1150, 0
    ax.text(x0, y0, "WING DATA", fontsize=10, fontweight="bold")
    for i, (k_, v_) in enumerate(rows):
        ax.text(x0, y0 - 40 - i * 40, k_, fontsize=8.5)
        ax.text(x0 + 300, y0 - 40 - i * 40, v_, fontsize=8.5)

    bx, by = 2240, -610
    ax.add_patch(Rectangle((bx, by), 810, 330, fill=False, lw=1.2))
    for yy in (by + 110, by + 190, by + 260):
        ax.plot([bx, bx + 810], [yy, yy], color="k", lw=0.6)
    ax.plot([bx + 405, bx + 405], [by, by + 190], color="k", lw=0.6)
    ax.text(bx + 405, by + 295, "S1223 SWEPT WING WITH FLAPERONS", fontsize=13,
            fontweight="bold", ha="center", va="center")
    ax.text(bx + 405, by + 225, "3-VIEW GENERAL ARRANGEMENT", fontsize=10, ha="center",
            va="center")
    cells = [("UNITS", "mm, degrees"), ("SCALE", "NTS (section A-A 2:1)"),
             ("ORIGIN", "root chord leading edge"), ("DATE", "2026-10-03"),
             ("CAD", f"CadQuery {cq.__version__}"), ("SHEET", "1 / 1")]
    for i, (k_, v_) in enumerate(cells):
        cx = bx + 12 + (i % 2) * 405
        cy = by + 165 - (i // 2) * 55
        ax.text(cx, cy, k_, fontsize=7, color="0.35")
        ax.text(cx, cy - 25, v_, fontsize=9)
    ax.text(-130, -600, "Hidden lines dashed. Views are true hidden-line projections of the "
            "CAD solid.", fontsize=8, color="0.3")

    fig.savefig(f"{OUT}/wing_3view_drawing.png", dpi=150)
    plt.close(fig)


# --------------------------------------------------------------- materials
def make_materials(parts, masses):
    fig = plt.figure(figsize=(16.5, 11.7), dpi=150)   # A3 landscape
    fig.text(0.04, 0.95, "MATERIAL SPECIFICATION  -  S1223 SWEPT WING", fontsize=20,
             fontweight="bold")
    fig.text(0.04, 0.925, "Carbon-fibre skin over foam core, flaperon each side.  "
             "Masses computed from CAD volumes.", fontsize=11, color="0.3")

    # section A-A with LE zoom
    ys = 200.0
    sec = section(parts, ys)
    xl = W.x_le(ys)
    ax = fig.add_axes([0.04, 0.6, 0.6, 0.27])
    ax2 = fig.add_axes([0.66, 0.56, 0.3, 0.32])
    for a in (ax, ax2):
        for name, ext, holes in sec:
            if "skin" in name:
                a.fill(ext[:, 0] - xl, ext[:, 1], color=C_CARBON, lw=0)
                for h in holes:
                    a.fill(h[:, 0] - xl, h[:, 1], color="white", lw=0)
        for name, ext, holes in sec:
            if "core" in name:
                a.fill(ext[:, 0] - xl, ext[:, 1], color=C_FOAM, lw=0)
        a.set_aspect("equal")
    ax.set_title(f"Section A-A at y = {ys:.0f} mm  (local chord {W.chord(ys):.1f} mm)",
                 fontsize=12, loc="left")
    ax.set_xlabel("x from local LE [mm]")
    ax.set_ylabel("z [mm]")
    ax.grid(alpha=0.25)
    ax.add_patch(Rectangle((-2, -6), 22, 22, fill=False, ec="tab:red", lw=1))
    ax.text(22, 17, "detail", color="tab:red", fontsize=9)
    ax2.set_xlim(-2, 20)
    ax2.set_ylim(-6, 16)
    ax2.set_title("Detail: leading edge (skin 0.3 mm)", fontsize=12, loc="left")
    ax2.grid(alpha=0.25)
    ax2.set_xlabel("x [mm]")
    ax.text(0.9 * W.chord(ys), 30, "flaperon", ha="center", fontsize=10)
    ax.annotate("", (0.9 * W.chord(ys), 4), (0.9 * W.chord(ys), 27),
                arrowprops=dict(arrowstyle="-|>", lw=0.6))
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=C_CARBON, label="CFRP skin"),
                       Patch(color=C_FOAM, label="Foam core")], loc="upper left", ncol=2)
    ax.set_ylim(-8, 52)

    # property tables
    def table(x, y, title, mat, color):
        fig.patches.append(Rectangle((x, y - 0.012), 0.012, 0.024, color=color,
                                     transform=fig.transFigure, figure=fig))
        fig.text(x + 0.02, y, title, fontsize=13, fontweight="bold", va="center")
        fig.text(x + 0.02, y - 0.03, mat["name"], fontsize=11)
        fig.text(x + 0.02, y - 0.053, mat["spec"], fontsize=9.5, color="0.35")
        for i, (k, v) in enumerate(mat["rows"]):
            yy = y - 0.09 - i * 0.026
            if i % 2 == 0:
                fig.patches.append(Rectangle((x + 0.015, yy - 0.012), 0.27, 0.024,
                                             color="0.94", transform=fig.transFigure,
                                             figure=fig, zorder=-1))
            fig.text(x + 0.02, yy, k, fontsize=10, va="center")
            fig.text(x + 0.17, yy, v, fontsize=10, va="center")

    table(0.04, 0.49, "SKIN", CARBON, C_CARBON)
    table(0.36, 0.49, "CORE", FOAM, C_FOAM)

    # mass breakdown
    x = 0.68
    fig.text(x, 0.49, "MASS BREAKDOWN", fontsize=13, fontweight="bold", va="center")
    hdr = ("Body", "Material", "Vol [cm³]", "Mass [g]")
    cols = [x, x + 0.105, x + 0.19, x + 0.25]
    for cx, h in zip(cols, hdr):
        fig.text(cx, 0.455, h, fontsize=10, fontweight="bold")
    for i, (name, v, m, mat) in enumerate(masses["rows"]):
        yy = 0.425 - i * 0.026
        for cx, val in zip(cols, (name, mat, f"{v:8.1f}", f"{m:7.1f}")):
            fig.text(cx, yy, val, fontsize=10)
    yy = 0.425 - len(masses["rows"]) * 0.026 - 0.01
    fig.lines.append(plt.Line2D([x, x + 0.29], [yy + 0.017] * 2, color="k", lw=0.8,
                                transform=fig.transFigure))
    fig.text(x, yy, "TOTAL (structure)", fontsize=10.5, fontweight="bold")
    fig.text(cols[3], yy, f"{masses['total']:7.1f}", fontsize=10.5, fontweight="bold")
    fig.text(x, yy - 0.035, f"with XPS 33 kg/m³ core: {masses['total_xps']:.0f} g",
             fontsize=10, color="0.3")
    fig.text(x, yy - 0.06, f"wing loading of structure: "
             f"{masses['total'] / (wing_data()['S'] / 1e6):.0f} g/m²", fontsize=10,
             color="0.3")

    notes = [
        "NOTES",
        "1. Skin: 2 plies, outer ply +-45 deg (torsion), inner ply 0/90 deg; "
        "vacuum-bagged on core, RT-cure epoxy.",
        f"2. Core is inset {W.SKIN_T} mm from the outer mould line and stops where it is "
        f"thinner than {W.MIN_CORE_T:.0f} mm;",
        "    aft of that the trailing edge is solid laminate. Tip, flaperon bay and "
        "hinge faces are closed with 0.3 mm CFRP.",
        "3. Excluded: spar caps, hinges, servo bays, wiring, paint/filler and adhesive "
        "(typically +15-25 % mass).",
        "4. Property values are typical handbook/datasheet values: confirm with the "
        "supplier before structural sizing.",
    ]
    for i, n in enumerate(notes):
        fig.text(0.04, 0.115 - i * 0.022, n, fontsize=10 if i else 11,
                 fontweight="bold" if i == 0 else "normal")
    fig.savefig(f"{OUT}/wing_materials.png", dpi=150)
    plt.close(fig)


def mass_table(parts):
    rows, tot, tot_xps = [], 0.0, 0.0
    labels = {"main_skin": "Main wing skin", "main_core": "Main wing core",
              "flaperon_R_skin": "R flaperon skin", "flaperon_R_core": "R flaperon core",
              "flaperon_L_skin": "L flaperon skin", "flaperon_L_core": "L flaperon core"}
    for name, shp in parts.items():
        v = W.volume(shp) / 1e3
        is_skin = "skin" in name
        rho = CARBON["density"] if is_skin else FOAM["density"]
        m = v * rho
        tot += m
        tot_xps += m if is_skin else v * XPS_DENSITY
        rows.append((labels[name], v, m, "CFRP" if is_skin else "PMI foam"))
    return dict(rows=rows, total=tot, total_xps=tot_xps)


if __name__ == "__main__":
    parts, outer = W.build()
    for n, s in {**parts, **outer}.items():
        assert s.isValid(), n
    masses = mass_table(parts)
    for r in masses["rows"]:
        print(f"{r[0]:18s} {r[1]:9.1f} cm3 {r[2]:8.1f} g")
    print(f"TOTAL {masses['total']:.1f} g  (XPS core: {masses['total_xps']:.1f} g)")
    export(parts, outer)
    make_drawing(parts, outer, masses)
    make_materials(parts, masses)
    print("written:", sorted(os.listdir(OUT)))
