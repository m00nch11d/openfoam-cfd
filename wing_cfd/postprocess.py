"""Post-processing: wing coefficients, station data and PNG images.

Outputs (results/):
  summary.md, summary.json, stations.csv
  residuals.png, convergence_CL_CD_CM.png
  wing_3d_pressure.png, wing_3d_velocity.png, wing_3d_turbulence_intensity.png,
  wing_3d_yplus.png
  stations/<name>_Cp.png and <name>_{pressure,velocity,turb_intensity}.png
"""
import json
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pyvista as pv  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
from matplotlib.tri import Triangulation  # noqa: E402

import params as P  # noqa: E402

pv.OFF_SCREEN = True
CASE = sys.argv[1] if len(sys.argv) > 1 else "case"
OUT = "results"
SOUT = os.path.join(OUT, "stations")
os.makedirs(SOUT, exist_ok=True)
QINF = 0.5 * P.U_INF ** 2          # kinematic dynamic pressure
A = math.radians(P.AOA_DEG)
LIFT = np.array([-math.sin(A), math.cos(A)])     # (x, z)
DRAG = np.array([math.cos(A), math.sin(A)])


# ------------------------------------------------------------ histories
def read_table(path):
    hdr, rows = None, []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                if "Time" in line:
                    hdr = line[1:].strip().split("\t")
                    hdr = [h.strip() for h in hdr if h.strip()]
                continue
            v = line.strip().split("\t")
            v = [x.strip() for x in v if x.strip()]
            if hdr and len(v) == len(hdr):
                rows.append(v)
    return hdr, rows


def read_all(func):
    """Concatenate a function-object table over restart directories."""
    base = os.path.join(CASE, "postProcessing", func)
    hdr, rows = None, {}
    for d in sorted(os.listdir(base), key=float):
        f = [x for x in os.listdir(os.path.join(base, d)) if x.endswith(".dat")][0]
        h, r = read_table(os.path.join(base, d, f))
        hdr = h
        for x in r:
            rows[float(x[0])] = x
    return hdr, [rows[k] for k in sorted(rows)]


def histories():
    h, r = read_all("forceCoeffs")
    co = {k: np.array([float(x[i]) for x in r]) for i, k in enumerate(h)}
    h, r = read_all("residuals")
    res = {}
    it = np.array([float(x[0]) for x in r])
    for i, k in enumerate(h):
        if k.endswith("_initial"):
            vals = []
            for x in r:
                try:
                    vals.append(float(x[i]))
                except ValueError:
                    vals.append(np.nan)
            res[k[:-8]] = np.array(vals)
    return co, it, res


def plot_histories(co, it, res):
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for k, v in res.items():
        ax.semilogy(it, v, lw=1, label=k)
    ax.set_xlabel("Iteration"); ax.set_ylabel("Initial residual")
    ax.set_title("Residuals (simpleFoam, k-omega SST + gamma-ReTheta transition)")
    ax.grid(True, which="both", alpha=0.3); ax.legend(ncol=3)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "residuals.png"), dpi=150)
    plt.close(fig)

    fig, axs = plt.subplots(3, 1, figsize=(9, 8), sharex=True)
    for ax, k, lab in zip(axs, ("Cl", "Cd", "CmPitch"), ("CL", "CD", "CM (c/4)")):
        ax.plot(co["Time"], co[k], lw=1)
        ax.set_ylabel(lab); ax.grid(True, alpha=0.3)
        v = co[k][len(co[k]) // 2:]
        ax.set_ylim(v.min() - 0.5 * np.ptp(v) - 1e-3, v.max() + 0.5 * np.ptp(v) + 1e-3)
    axs[-1].set_xlabel("Iteration")
    axs[0].set_title("Wing force-coefficient convergence")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "convergence_CL_CD_CM.png"), dpi=150)
    plt.close(fig)


# ------------------------------------------------------------ fields
def latest_time():
    ts = []
    for d in os.listdir(CASE):
        try:
            if float(d) > 0:
                ts.append(d)
        except ValueError:
            pass
    return max(ts, key=float)


def load(t):
    foam = os.path.join(CASE, "case.foam")
    open(foam, "w").close()
    rd = pv.OpenFOAMReader(foam)
    rd.set_active_time_value(float(t))
    rd.enable_all_patch_arrays()
    rd.enable_all_cell_arrays()
    mb = rd.read()
    internal = mb["internalMesh"]
    wing = mb["boundary"]["wing"]
    return internal, wing


def add_derived(m):
    if "k" in m.point_data:
        m.point_data["TI"] = np.sqrt(2.0 / 3.0 * np.maximum(m.point_data["k"], 1e-14)) / P.U_INF * 100
    if "U" in m.point_data:
        m.point_data["Umag"] = np.linalg.norm(m.point_data["U"], axis=1)
    if "p" in m.point_data:
        m.point_data["Cp"] = m.point_data["p"] / QINF
    return m


# ------------------------------------------------------------ stations
def section(wing, y, tau_sign):
    """Slice of the wing patch (cell data) at span y -> ordered segments."""
    w = wing.copy()
    w.cell_data["n"] = w.compute_normals(cell_normals=True, point_normals=False,
                                         auto_orient_normals=False,
                                         consistent_normals=False,
                                         split_vertices=False).cell_data["Normals"]
    s = w.slice(normal=(0, 1, 0), origin=(0, y, 0))
    pts = s.points
    lines = s.lines.reshape(-1, 3)[:, 1:]
    a, b = pts[lines[:, 0]], pts[lines[:, 1]]
    mid = 0.5 * (a + b)
    dl = np.linalg.norm((b - a)[:, [0, 2]], axis=1)
    p = s.cell_data["p"]
    tau = s.cell_data["wallShearStress"] * tau_sign
    n = s.cell_data["n"]
    g = s.cell_data["gammaInt"]
    yp = s.cell_data["yPlus"] if "yPlus" in s.cell_data else np.zeros(len(p))
    # drop tip-cap faces (normal along the span)
    keep = np.abs(n[:, 1]) < 0.9
    return dict(mid=mid[keep], dl=dl[keep], p=p[keep], tau=tau[keep], n=n[keep],
                gamma=g[keep], yplus=yp[keep], seg=[tuple(l) for l in lines[keep]])


def station_coeffs(sec, y):
    xle = P.x_le(y)
    c = P.local_chord(y)
    f = (sec["p"][:, None] * sec["n"] + sec["tau"]) * sec["dl"][:, None]   # per span
    fxz = f[:, [0, 2]].sum(0)
    cl = fxz @ LIFT / (QINF * c)
    cd = fxz @ DRAG / (QINF * c)
    r = sec["mid"][:, [0, 2]] - np.array([xle + 0.25 * c, 0.0])
    my = (r[:, 1] * f[:, 0] - r[:, 0] * f[:, 2]).sum()
    cm = my / (QINF * c * c)
    return cl, cd, cm


def split_surfaces(sec, y):
    """Upper / lower by walking the closed section contour from LE to TE."""
    xle = P.x_le(y)
    c = P.local_chord(y)
    xc = (sec["mid"][:, 0] - xle) / c
    seg = sec["seg"]
    # adjacency of segments through shared end points
    from collections import defaultdict
    at = defaultdict(list)
    for i, (a, b) in enumerate(seg):
        at[a].append(i); at[b].append(i)
    order, used = [0], {0}
    cur_pt = seg[0][1]
    while True:
        nxt = [j for j in at[cur_pt] if j not in used]
        if not nxt:
            break
        j = nxt[0]
        order.append(j); used.add(j)
        cur_pt = seg[j][0] if seg[j][1] == cur_pt else seg[j][1]
    order = np.array(order)
    xo = xc[order]
    i_le, i_te = int(np.argmin(xo)), int(np.argmax(xo))
    n = len(order)
    arc1 = [order[(i_le + k) % n] for k in range((i_te - i_le) % n + 1)]
    arc2 = [order[(i_te + k) % n] for k in range((i_le - i_te) % n + 1)]
    z = sec["mid"][:, 2]
    upper = np.zeros(len(xc), bool)
    up = arc1 if z[arc1].mean() > z[arc2].mean() else arc2
    upper[up] = True
    return xc, upper


NUT_TR = 2.0       # transition threshold on the boundary-layer max of nut/nu


def bl_nut_ratio(internal, sec, y):
    """Max nut/nu along a wall-normal probe (0.02-6 mm) at each face."""
    dist = np.geomspace(2e-5, 6e-3, 25)
    pts = sec["mid"][:, None, :] - sec["n"][:, None, :] * dist[None, :, None]
    pts[:, :, 1] = y
    s = pv.PolyData(pts.reshape(-1, 3)).sample(internal)
    return (s.point_data["nut"].reshape(len(sec["mid"]), len(dist)) / P.NU).max(1)


def transition(xc, cfx, nutr):
    """Transition onset: first x/c where the boundary-layer maximum of
    nut/nu rises above NUT_TR (laminar ~0.1-0.4, turbulent > 10) and stays
    above it for the next 3 faces; None if the surface stays laminar.
    Also laminar separation (Cf < 0) and reattachment (Cf > 0 again)."""
    o = np.argsort(xc)
    x, cf, r = xc[o], cfx[o], nutr[o]
    xt = None
    for i in range(1, len(x) - 3):
        if x[i] > 0.005 and r[i] >= NUT_TR and r[i - 1] < NUT_TR and (r[i:i + 4] >= NUT_TR).all():
            f = (NUT_TR - r[i - 1]) / (r[i] - r[i - 1])
            xt = float(x[i - 1] + f * (x[i] - x[i - 1]))
            break
    m = (x > 0.005) & (x < 0.995)
    x, cf = x[m], cf[m]
    sep = np.where(cf < 0)[0]
    xs = float(x[sep[0]]) if len(sep) else None
    xr = None
    if len(sep):
        rest = np.where((cf > 0) & (np.arange(len(cf)) > sep[0]))[0]
        xr = float(x[rest[0]]) if len(rest) else None
    return xt, xs, xr


def plot_cp(name, y, kind, xc, upper, cp, cl, cd, cm, tr):
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for msk, lab, col in ((upper, "Upper surface", "tab:blue"), (~upper, "Lower surface", "tab:red")):
        o = np.argsort(xc[msk])
        ax.plot(xc[msk][o], cp[msk][o], "-", color=col, lw=1.4, label=lab)
    for key, col, ls in (("upper", "tab:blue", "--"), ("lower", "tab:red", "--")):
        if tr[key]["transition"] is not None:
            ax.axvline(tr[key]["transition"], color=col, ls=ls, lw=1,
                       label=f"Transition {key} x/c={tr[key]['transition']:.3f}")
    ax.invert_yaxis()
    ax.set_xlabel("x / c (local)"); ax.set_ylabel("Cp")
    ax.set_title(f"{name}: y = {y*1000:.1f} mm ({kind}), c = {P.local_chord(y)*1000:.1f} mm\n"
                 f"cl = {cl:.4f}   cd = {cd:.5f}   cm(c/4) = {cm:.4f}")
    ax.grid(True, alpha=0.3); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(SOUT, f"{name}_Cp.png"), dpi=140)
    plt.close(fig)


def plot_slice_fields(name, y, kind, internal):
    s = internal.slice(normal=(0, 1, 0), origin=(0, y + 1e-6, 0))
    s = s.triangulate()
    s = add_derived(s)
    x, z = s.points[:, 0], s.points[:, 2]
    tri = s.faces.reshape(-1, 4)[:, 1:] if s.faces.size else s.regular_faces
    T = Triangulation(x, z, tri)
    win = (-0.15, 0.6, -0.2, 0.2)
    specs = [("pressure", "Cp", "Pressure coefficient Cp", "RdBu_r", None),
             ("velocity", "Umag", "Velocity magnitude |U| [m/s]", "viridis", None),
             ("turb_intensity", "TI", "Turbulence intensity [%] (log scale)", "inferno", "log")]
    for tag, f, lab, cmap, scale in specs:
        v = s.point_data[f]
        fig, ax = plt.subplots(figsize=(10, 5.6))
        if scale == "log":
            vmin = max(np.percentile(v, 1), 1e-3)
            levels = np.logspace(np.log10(vmin), np.log10(max(v.max(), vmin * 10)), 60)
            cs = ax.tricontourf(T, np.clip(v, vmin, None), levels=levels, cmap=cmap,
                                norm=LogNorm(levels[0], levels[-1]))
        else:
            lo, hi = (-3.0, 1.0) if f == "Cp" else (0.0, np.percentile(v, 99.9))
            cs = ax.tricontourf(T, np.clip(v, lo, hi), levels=np.linspace(lo, hi, 61), cmap=cmap)
        fig.colorbar(cs, ax=ax, label=lab)
        ax.set_xlim(win[0], win[1]); ax.set_ylim(win[2], win[3]); ax.set_aspect("equal")
        ax.set_xlabel("x [m]"); ax.set_ylabel("z [m]")
        ax.set_title(f"{name}: section y = {y*1000:.1f} mm ({kind}) - {lab}")
        fig.tight_layout(); fig.savefig(os.path.join(SOUT, f"{name}_{tag}.png"), dpi=130)
        plt.close(fig)


# ------------------------------------------------------------ 3D views
def iso_camera(pl, focus=(0.2, 0.55, 0.0), scale=0.55, d=(-1.0, -1.0, 1.0)):
    """Isometric view from upstream (leading-edge side), above, root side."""
    focus = np.array(focus)
    d = np.array(d, float)
    d /= np.linalg.norm(d)
    pl.camera.focal_point = focus
    pl.camera.position = focus + 3.0 * d
    pl.camera.up = (0, 0, 1)
    pl.camera.parallel_projection = True
    pl.camera.parallel_scale = scale


def render_pathlines(internal, wing):
    """Particle paths through the tip region, from upstream of the LE to
    PATHLINE_X_BEHIND behind the trailing edge (tip vortex)."""
    x_end = P.C + P.PATHLINE_X_BEHIND
    # VTK cannot locate points reliably inside the micron-thin wall/wake
    # cells (paths stop "out of domain"), so the tip region is resampled
    # onto a uniform 5 mm grid first (points inside the wing get U = 0)
    x0, x1, y0, y1, z0, z1 = -0.1, x_end, P.B - 0.25, P.B + 0.25, -0.25, 0.25
    h = 0.005
    grid = pv.ImageData(dimensions=(int(round((x1 - x0) / h)) + 1,
                                    int(round((y1 - y0) / h)) + 1,
                                    int(round((z1 - z0) / h)) + 1),
                        spacing=(h, h, h), origin=(x0, y0, z0))
    box = internal.clip_box((x0 - 0.02, x1 + 0.02, y0 - 0.02, y1 + 0.02,
                             z0 - 0.02, z1 + 0.02), invert=False)
    grid = grid.sample(box)
    # seeds in the vortex-formation region at the tip trailing edge,
    # integrated both ways: backwards they show the flow wrapping round the
    # tip from the lower surface, forwards the roll-up into the tip vortex
    ys, zs = np.meshgrid(np.linspace(P.B - 0.06, P.B + 0.02, 9),
                         np.linspace(-0.01, 0.05, 7))
    seeds = pv.PolyData(np.column_stack([np.full(ys.size, P.C + 0.005),
                                         ys.ravel(), zs.ravel()]))
    sl = grid.streamlines_from_source(seeds, vectors="U", integration_direction="both",
                                      max_steps=20000, initial_step_length=0.5,
                                      max_step_length=1.0, terminal_speed=1e-3)
    sl.point_data["Umag"] = np.linalg.norm(sl.point_data["U"], axis=1)
    tubes = sl.tube(radius=0.0012)
    wing_pd = add_derived(wing.cell_data_to_point_data())
    clim = (0.5 * P.U_INF, 1.4 * P.U_INF)
    for tag, d, scale, focus in (
            ("iso", (-1.0, -1.0, 1.0), 0.5, (0.55, 0.85, 0.0)),
            ("rear", (1.0, -0.12, 0.12), 0.16, (P.C, P.B - 0.02, 0.0))):
        pl = pv.Plotter(off_screen=True, window_size=(1800, 1200))
        pl.set_background("white")
        pl.add_mesh(wing_pd, scalars="Cp", cmap="RdBu_r", clim=(-2.5, 1),
                    smooth_shading=True,
                    scalar_bar_args=dict(title="Cp", color="black", vertical=True,
                                         position_x=0.04, position_y=0.2, height=0.5))
        pl.add_mesh(tubes, scalars="Umag", cmap="viridis", clim=clim,
                    scalar_bar_args=dict(title="|U| [m/s]", color="black", vertical=True,
                                         position_x=0.9, position_y=0.2, height=0.6))
        view = "isometric, top / leading-edge side" if tag == "iso" else \
            "looking upstream from the end of the paths"
        pl.add_text(f"Tip-vortex particle paths to x = TE + {P.PATHLINE_X_BEHIND / P.C:g}c "
                    f"({x_end:.2f} m); {view}", font_size=11, color="black")
        iso_camera(pl, focus=focus, scale=scale, d=d)
        pl.screenshot(os.path.join(OUT, f"wing_3d_tip_vortex_pathlines{'' if tag == 'iso' else '_rear'}.png"))
        pl.close()
    return sl.n_lines


def render_3d(internal, wing):
    wing_pd = add_derived(wing.cell_data_to_point_data())
    box = internal.clip_box((-0.12, 0.75, -0.01, P.B + 0.3, -0.16, 0.16), invert=False)
    ys = [0.0625, 0.5, 0.9375, P.B + 0.05]
    slices = pv.MultiBlock([add_derived(box.slice(normal=(0, 1, 0), origin=(0, y + 1e-6, 0)))
                            for y in ys]).combine()
    specs = [("pressure", "Cp", "Cp", "RdBu_r", (-2.5, 1), False),
             ("velocity", "Umag", "|U| [m/s]", "viridis", (0, 1.6 * P.U_INF), False),
             ("turbulence_intensity", "TI", "Turbulence intensity [%] (log)", "inferno", None, True)]
    for tag, f, lab, cmap, clim, log in specs:
        pl = pv.Plotter(off_screen=True, window_size=(1800, 1200))
        pl.set_background("white")
        if clim is None:
            v = slices.point_data[f]
            clim = (max(np.percentile(v, 0.5), 1e-2), max(v.max(), 1.0))
        sargs = dict(title=lab, color="black", vertical=True, position_x=0.9,
                     position_y=0.2, height=0.6, fmt="%.3g")
        pl.add_mesh(slices, scalars=f, cmap=cmap, clim=clim, log_scale=log,
                    lighting=False, scalar_bar_args=sargs)
        if f == "TI":       # k = 0 on the wall: draw the wing in grey
            pl.add_mesh(wing_pd, color="lightgrey", smooth_shading=True)
        else:
            pl.add_mesh(wing_pd, scalars=f, cmap=cmap, clim=clim, log_scale=log,
                        show_scalar_bar=False, smooth_shading=True)
        pl.add_text(f"Half wing, S1223 sine LE, Re 3e5, AoA {P.AOA_DEG:g} deg: {lab}\n"
                    f"section planes y = {', '.join(f'{y:.4g}' for y in ys)} m "
                    f"(last plane outboard of the tip)", font_size=11, color="black")
        iso_camera(pl)
        pl.screenshot(os.path.join(OUT, f"wing_3d_{tag}.png"))
        pl.close()

    nrm = wing.compute_normals(cell_normals=True, point_normals=False,
                               auto_orient_normals=False, consistent_normals=False,
                               split_vertices=False).cell_data["Normals"]
    tip = np.abs(nrm[:, 1]) > 0.9
    yv = wing.cell_data["yPlus"]
    surf = wing.extract_cells(np.where(~tip)[0])
    pl = pv.Plotter(off_screen=True, window_size=(1800, 1200))
    pl.set_background("white")
    pl.add_mesh(surf, scalars="yPlus", cmap="turbo", clim=(0, 1.5), preference="cell",
                above_color="magenta",
                scalar_bar_args=dict(title="y+ (magenta > 1.5)", color="black", vertical=True,
                                     position_x=0.9, position_y=0.2, height=0.6))
    ys_ = yv[~tip]
    pl.add_text(f"y+ on the wing surface: max {ys_.max():.2f}, mean {ys_.mean():.3f}, "
                f"{(ys_ < 1).mean()*100:.1f}% of faces < 1\n"
                f"(flat tip cap, resolved only by the 15.6 mm spanwise cell: "
                f"mean y+ {yv[tip].mean():.0f}, not shown)", font_size=11, color="black")
    iso_camera(pl)
    pl.screenshot(os.path.join(OUT, "wing_3d_yplus.png"))
    pl.close()
    return yv[~tip], yv[tip]


# ------------------------------------------------------------ main
def main():
    co, it, res = histories()
    plot_histories(co, it, res)
    t = latest_time()
    internal, wing = load(t)
    ncells = internal.n_cells
    status = open(os.path.join(CASE, "run_status.txt")).read().strip() \
        if os.path.exists(os.path.join(CASE, "run_status.txt")) else "unknown"

    # sign of the wallShearStress so that patch integration matches forceCoeffs
    nrm = wing.compute_normals(cell_normals=True, point_normals=False,
                               auto_orient_normals=False, consistent_normals=False,
                               split_vertices=False).cell_data["Normals"]
    area = wing.compute_cell_sizes(length=False, volume=False).cell_data["Area"]
    fp = (wing.cell_data["p"][:, None] * nrm * area[:, None]).sum(0)
    fv = (wing.cell_data["wallShearStress"] * area[:, None]).sum(0)
    cd_fo = co["Cd"][-1]
    sref = P.C * P.B
    cand = {s: ((fp + s * fv)[[0, 2]] @ DRAG) / (QINF * sref) for s in (1, -1)}
    tau_sign = min(cand, key=lambda s: abs(cand[s] - cd_fo))
    cl_int = ((fp + tau_sign * fv)[[0, 2]] @ LIFT) / (QINF * sref)

    yplus, yplus_tip = render_3d(internal, wing)
    n_paths = render_pathlines(internal, wing)

    rows = []
    for k, (y, kind) in enumerate(P.stations()):
        ys = min(max(y, 1e-4), P.B - 1e-4)
        name = f"S{k:02d}_y{y*1000:06.1f}mm_{kind}"
        sec = section(wing, ys, tau_sign)
        cl, cd, cm = station_coeffs(sec, ys)
        xc, upper = split_surfaces(sec, ys)
        cp = sec["p"] / QINF
        # streamwise wall shear (flow over both surfaces is roughly +x)
        cfx = sec["tau"][:, 0] / QINF
        nutr = bl_nut_ratio(internal, sec, ys)
        tr = {}
        for key, msk in (("upper", upper), ("lower", ~upper)):
            xt, xs, xr = transition(xc[msk], cfx[msk], nutr[msk])
            tr[key] = dict(transition=xt, lam_sep=xs, reattach=xr)
        plot_cp(name, ys, kind, xc, upper, cp, cl, cd, cm, tr)
        plot_slice_fields(name, ys, kind, internal)
        rows.append(dict(station=name, y_mm=y * 1000, type=kind,
                         chord_mm=P.local_chord(y) * 1000, cl=cl, cd=cd, cm=cm,
                         xtr_upper=tr["upper"]["transition"], xtr_lower=tr["lower"]["transition"],
                         sep_upper=tr["upper"]["lam_sep"], reat_upper=tr["upper"]["reattach"],
                         sep_lower=tr["lower"]["lam_sep"], reat_lower=tr["lower"]["reattach"],
                         yplus_mean=float(sec["yplus"].mean()),
                         yplus_max=float(sec["yplus"].max())))
        print(name, f"cl={cl:.4f} cd={cd:.5f} cm={cm:.4f}", tr, flush=True)

    n_avg = min(100, len(co["Cl"]))
    summary = dict(
        status=status, iterations=int(co["Time"][-1]), total_cells=int(ncells),
        CL=float(co["Cl"][-1]), CD=float(co["Cd"][-1]), CM=float(co["CmPitch"][-1]),
        CL_avg_last100=float(co["Cl"][-n_avg:].mean()),
        CD_avg_last100=float(co["Cd"][-n_avg:].mean()),
        CM_avg_last100=float(co["CmPitch"][-n_avg:].mean()),
        CL_patch_integration_check=float(cl_int),
        L_over_D=float(co["Cl"][-1] / co["Cd"][-1]),
        yplus_max=float(yplus.max()), yplus_mean=float(yplus.mean()),
        yplus_frac_below_1=float((yplus < 1).mean()),
        yplus_p99=float(np.percentile(yplus, 99)),
        yplus_tip_cap_mean=float(yplus_tip.mean()),
        first_layer_height_m=P.H1, U_inf=P.U_INF, Re=P.RE, AoA_deg=P.AOA_DEG,
        Sref_m2=sref, cref_m=P.C, moment_ref="(0.075, 0, 0) m = mean-chord quarter chord",
        stations=rows)
    with open(os.path.join(OUT, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    with open(os.path.join(OUT, "stations.csv"), "w") as f:
        keys = list(rows[0].keys())
        f.write(",".join(keys) + "\n")
        for r in rows:
            f.write(",".join("" if r[k] is None else (f"{r[k]:.6g}" if isinstance(r[k], float) else str(r[k]))
                             for k in keys) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "stations"}, indent=2))


if __name__ == "__main__":
    main()
