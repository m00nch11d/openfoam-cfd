#!/usr/bin/env python3
"""
Post-process the S1223 case: coefficients, transition points, images.

Writes to results/:
  summary.md, summary.json
  cp_distribution.png        Cp over upper and lower surfaces
  cf_distribution.png        skin friction (used for the transition point)
  pressure_contours.png      Cp field, whole domain and near the airfoil
  velocity_contours.png      |U|/U_inf, whole domain and near the airfoil
  turbulence_intensity.png   Tu = sqrt(2k/3)/U_inf, logarithmic colour scale
  yplus_distribution.png     y+ over upper and lower surfaces
  residuals.png              initial residuals vs iteration
  force_coefficients.png     Cl, Cd, Cm vs iteration
  mesh_domain.png            whole C-mesh
  mesh_airfoil.png           2 x 2 chords around the airfoil
  mesh_details.png           leading and trailing edge close-ups

Usage: python3 tools/postprocess.py [case_dir]
"""
import glob
import json
import math
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402
from matplotlib.colors import LogNorm, TwoSlopeNorm  # noqa: E402
from matplotlib.path import Path  # noqa: E402
from matplotlib.tri import Triangulation  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

CASE = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else
                       os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from foamio import Mesh2D, read_field  # noqa: E402
from make_mesh import read_parameters  # noqa: E402

OUT = os.path.join(CASE, "results")
os.makedirs(OUT, exist_ok=True)
P = read_parameters(os.path.join(CASE, "caseParameters"))
C = P["chord"]
NU = P["mu"] / P["rho"]
UINF = P["Re"] * NU / C
Q = 0.5 * UINF ** 2                      # kinematic dynamic pressure

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#7d5fd3"]
INK, MUTED, GRID = "#1f1f1e", "#6b6a63", "#e4e3dc"
plt.rcParams.update({"axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.titlecolor": INK, "font.size": 10,
                     "axes.spines.top": False, "axes.spines.right": False})


def style(ax, grid=True):
    if grid:
        ax.grid(True, color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)


def save(fig, name):
    path = os.path.join(OUT, name)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote results/{name}")


def table(func, name):
    files = sorted(glob.glob(os.path.join(CASE, "postProcessing", func, "*", name)),
                   key=lambda f: float(os.path.basename(os.path.dirname(f))))
    cols, data = None, {}
    for path in files:
        with open(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    w = line[1:].split()
                    if w and w[0] == "Time":
                        cols = w
                    continue
                v = line.split()
                if cols and len(v) == len(cols):
                    data[float(v[0])] = v
    return cols, [data[t] for t in sorted(data)]


def fname(T, name):
    """Averaged field (fieldAverage) when it exists, else the instantaneous one."""
    return name + "Mean" if os.path.exists(os.path.join(CASE, T, name + "Mean")) else name


def latest_time():
    times = [d for d in os.listdir(CASE) if d.replace(".", "", 1).isdigit() and float(d) > 0]
    return max(times, key=float) if times else None


# ------------------------------------------------------------------ surfaces
def surface_data(M, T):
    info = np.load(os.path.join(CASE, "constant", "meshInfo.npz"))
    xy, nrm, own = M.patch("airfoil")
    base = np.abs(nrm[:, 0]) > 0.99                      # vertical TE base faces
    d_up, _ = cKDTree(info["upper"]).query(xy)
    d_lo, _ = cKDTree(info["lower"]).query(xy)
    up = (d_up < d_lo) & ~base
    lo = (d_lo <= d_up) & ~base
    p = read_field(CASE, T, fname(T, "p"))
    tau = read_field(CASE, T, fname(T, "wallShearStress"), patch="airfoil")
    yp = read_field(CASE, T, fname(T, "yPlus"), patch="airfoil")
    nut = read_field(CASE, T, fname(T, "nut"))
    # peak eddy-viscosity ratio across the boundary layer (wall distance
    # 10 um .. 5 mm along the face normal); ~0 laminar, O(10) turbulent
    tree = cKDTree(M.centres)
    dist = np.geomspace(1e-5, 5e-3, 60)
    nut_max = np.empty(len(xy))
    for f in range(len(xy)):
        _, c = tree.query(xy[f] - nrm[f] * dist[:, None])   # nrm points into the wall
        nut_max[f] = (nut[c] / NU).max()
    out = {}
    for name, m in (("upper", up), ("lower", lo)):
        o = np.argsort(xy[m, 0])
        f_xy, f_n, f_own = xy[m][o], nrm[m][o], own[m][o]
        # tangent pointing from LE to TE; OpenFOAM's wallShearStress is the
        # stress on the fluid, i.e. opposite to the near-wall flow
        t = np.column_stack([f_n[:, 1], -f_n[:, 0]])
        t *= np.sign(t[:, 0])[:, None]
        cf = -np.einsum("ij,ij->i", tau[m][o][:, :2], t) / Q
        out[name] = dict(x=f_xy[:, 0] / C, y=f_xy[:, 1] / C, cp=p[f_own] / Q, cf=cf,
                         yplus=yp[m][o], nut_ratio=nut_max[m][o])
    out["base_cp"] = float(np.mean(p[own[base]]) / Q) if base.any() else None
    return out


NUT_TRANSITION = 2.0      # BL peak nut/nu marking turbulent flow


def transition(s):
    """Transition point on one surface: first x/c after which the peak
    eddy-viscosity ratio across the boundary layer stays above
    NUT_TRANSITION (laminar ~0, turbulent O(10)). Laminar separation and
    reattachment of a separation bubble come from the sign of Cf."""
    x, cf, r = s["x"], s["cf"], s["nut_ratio"]
    m = (x > 0.005) & (x < 0.995)
    x, cf, r = x[m], cf[m], r[m]
    lam = np.where(r < NUT_TRANSITION)[0]
    res = {}
    if len(lam) == 0:
        res["transition_x_c"] = float(x[0])
    elif lam[-1] == len(x) - 1:
        res["transition_x_c"] = None                      # laminar to the TE
    else:
        k = lam[-1]                                       # last laminar point
        x0, x1, r0, r1 = x[k], x[k + 1], r[k], r[k + 1]
        res["transition_x_c"] = float(x0 + (NUT_TRANSITION - r0) / (r1 - r0) * (x1 - x0))
    # laminar separation bubble: the Cf < 0 run that ends closest upstream
    # of (or contains) the transition point
    neg = np.where(cf < 0)[0]
    xt = res["transition_x_c"]
    if len(neg) and xt is not None:
        runs = np.split(neg, np.where(np.diff(neg) > 1)[0] + 1)
        runs = [r_ for r_ in runs if x[r_[0]] <= xt + 0.05]
        if runs:
            # merge runs separated by short attached patches (unsteady bubble)
            first = max(r_[0] for r_ in runs if x[r_[0]] <= xt)
            sel = [r_ for r_ in runs if r_[0] >= first]
            res["separation_x_c"] = float(x[sel[0][0]])
            res["reattachment_x_c"] = float(x[min(sel[-1][-1] + 1, len(x) - 1)])
    return res


# --------------------------------------------------------------- plotting
def airfoil_polygon():
    info = np.load(os.path.join(CASE, "constant", "meshInfo.npz"))
    return np.vstack([info["lower"], info["upper"][1:]])


def field_plot(M, f, name, label, cmap, norm=None, levels=None, ticks=None):
    poly = airfoil_polygon()
    cen = M.centres
    tri = Triangulation(cen[:, 0], cen[:, 1])
    t = tri.triangles
    # mask triangles that reach into the airfoil, and slivers along the hull
    # of the cell centres (edges much longer than the cells they connect)
    path = Path(poly)
    mids = [0.5 * (cen[t[:, a]] + cen[t[:, b]]) for a, b in ((0, 1), (1, 2), (2, 0))]
    inside = path.contains_points(cen[t].mean(axis=1))
    for m in mids:
        inside |= path.contains_points(m)
    size = np.array([np.linalg.norm(p - np.roll(p, 1, axis=0), axis=1).max() for p in M.polys])
    edge = np.max(np.linalg.norm(cen[t] - cen[np.roll(t, 1, axis=1)], axis=2), axis=1)
    sliver = edge > 3.0 * size[t].max(axis=1)
    tri.set_mask(inside | sliver)
    fig, axs = plt.subplots(1, 2, figsize=(16, 6.4), gridspec_kw=dict(width_ratios=[1, 1.35]))
    views = [((-P["farfieldChords"] * C, P["wakeChords"] * C + C),
              (-P["farfieldChords"] * C, P["farfieldChords"] * C), "Whole domain"),
             ((-0.5 * C, 1.5 * C), (-0.55 * C, 0.55 * C), "Near the airfoil")]
    for ax, (xl, yl, title) in zip(axs, views):
        cs = ax.tricontourf(tri, f, levels=levels, cmap=cmap, norm=norm, extend="both")
        ax.fill(poly[:, 0], poly[:, 1], color="#d9d8d2", ec=INK, lw=0.6, zorder=3)
        ax.set_xlim(*xl); ax.set_ylim(*yl); ax.set_aspect("equal")
        ax.set_title(title); ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")
    cb = fig.colorbar(cs, ax=axs, shrink=0.85, pad=0.02, ticks=ticks)
    cb.set_label(label)
    fig.suptitle(f"S1223, Re = {P['Re']:.0f}, alpha = {P['alpha']:g} deg - {label}", x=0.06, ha="left")
    save(fig, name)


def mesh_plot(M, name, views, titles, nq):
    global BB
    BB = np.array([[p[:, 0].min(), p[:, 1].min(), p[:, 0].max(), p[:, 1].max()] for p in M.polys])
    fig, axs = plt.subplots(1, len(views), figsize=(8 * len(views), 7.2))
    axs = np.atleast_1d(axs)
    for ax, (x0, x1, y0, y1), title in zip(axs, views, titles):
        sel = np.where((BB[:, 0] < x1) & (BB[:, 2] > x0) & (BB[:, 1] < y1) & (BB[:, 3] > y0))[0]
        quads = [M.polys[i] for i in sel if len(M.polys[i]) == 4]
        tris = [M.polys[i] for i in sel if len(M.polys[i]) == 3]
        ax.add_collection(PolyCollection(tris, facecolors="#eef3fb", edgecolors=BLUE, linewidths=0.25))
        ax.add_collection(PolyCollection(quads, facecolors="#fdf0ea", edgecolors=ORANGE, linewidths=0.25))
        poly = airfoil_polygon()
        ax.fill(poly[:, 0], poly[:, 1], color="#b9b8b0", ec=INK, lw=0.6, zorder=3)
        ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_aspect("equal")
        ax.set_title(title); ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")
    axs[0].legend(handles=[plt.Rectangle((0, 0), 1, 1, fc="#fdf0ea", ec=ORANGE, label="quadrilaterals"),
                           plt.Rectangle((0, 0), 1, 1, fc="#eef3fb", ec=BLUE, label="triangles")],
                  loc="upper right", frameon=True)
    save(fig, name)


def surface_plot(S, key, ylabel, name, title, invert=False, hlines=()):
    fig, ax = plt.subplots(figsize=(9, 5.2))
    for side, color in (("upper", BLUE), ("lower", ORANGE)):
        ax.plot(S[side]["x"], S[side][key], color=color, lw=2, label=f"{side} surface")
    for y, lab in hlines:
        ax.axhline(y, color=MUTED, lw=1, ls="--")
        ax.annotate(lab, xy=(1.0, y), xycoords=("axes fraction", "data"), xytext=(-4, 3),
                    textcoords="offset points", ha="right", color=INK, fontsize=9)
    if invert:
        ax.invert_yaxis()
    ax.set_xlabel("x/c"); ax.set_ylabel(ylabel); ax.set_title(title, loc="left")
    ax.legend(frameon=False)
    style(ax)
    return fig, ax


# ------------------------------------------------------------------- main
def main():
    T = latest_time()
    status = "UNKNOWN"
    st = os.path.join(OUT, "run_status.txt")
    if os.path.exists(st):
        lines = open(st).read().splitlines()
        status = " - ".join(lines[:2])
    M = Mesh2D(CASE)
    info = np.load(os.path.join(CASE, "constant", "meshInfo.npz"))
    nq, nt = int(info["n_quad"]), int(info["n_tri"])
    print(f"Post-processing time {T}, status {status}")

    # ---- mesh images
    mesh_plot(M, "mesh_domain.png", [(-9.3, 9.6, -9.3, 9.3)], ["C-mesh, whole domain"], nq)
    mesh_plot(M, "mesh_airfoil.png", [(-0.5 * C, 1.5 * C, -C, C)],
              ["Around the airfoil (2 x 2 chords)"], nq)
    bx, by = 0.5 * (info["lower"][0] + info["upper"][-1])          # TE base centre
    mesh_plot(M, "mesh_details.png",
              [(-0.004, 0.008, -0.006, 0.006), (bx - 0.004, bx + 0.002, by - 0.003, by + 0.003),
               (bx - 0.00012, bx + 0.00008, by - 0.0001, by + 0.0001)],
              ["Leading edge", "Trailing edge", "Blunt trailing-edge base (x/c = 0.999)"], nq)

    # ---- histories
    rcols, rrows = table("solverInfo", "solverInfo.dat")
    fcols, frows = table("forceCoeffs", "coefficient.dat")
    summary = dict(status=status, iterations=int(float(frows[-1][0])) if frows else 0,
                   cells=nq + nt, quads=nq, triangles=nt)
    if rrows:
        fig, ax = plt.subplots(figsize=(10, 5.5))
        it = np.array([float(r[0]) for r in rrows])
        names = [c for c in rcols if c.endswith("_initial")]
        for k, c in enumerate(names):
            v = np.array([float(r[rcols.index(c)]) for r in rrows])
            ax.semilogy(it, v, lw=1.3, color=SERIES[k % len(SERIES)], label=c.replace("_initial", ""))
        ax.set_xlabel("Iteration"); ax.set_ylabel("Initial residual")
        ax.set_title("Residuals", loc="left")
        ax.legend(ncol=len(names), frameon=False, loc="upper right", fontsize=9)
        style(ax)
        save(fig, "residuals.png")
    if frows:
        it = np.array([float(r[0]) for r in frows])
        co = {k: np.array([float(r[fcols.index(k)]) for r in frows]) for k in ("Cl", "Cd", "CmPitch")}
        w = it >= P["averageStart"] if it[-1] - P["averageStart"] >= 200 else it >= it[-1] - 199
        summary["average_window"] = [int(it[w][0]), int(it[w][-1])]
        for k, lab in (("Cl", "cl"), ("Cd", "cd"), ("CmPitch", "cm")):
            summary[lab + "_final"] = float(co[k][-1])
            summary[lab] = float(co[k][w].mean())
            summary[lab + "_std"] = float(co[k][w].std())
            summary[lab + "_min"] = float(co[k][w].min())
            summary[lab + "_max"] = float(co[k][w].max())
        rel = np.abs(np.diff(co["Cl"])) / np.maximum(np.abs(co["Cl"][1:]), 1e-12)
        k5 = np.where(rel < 0.05)[0]
        if len(k5):
            summary["literal_5pct_rule"] = dict(iteration=int(it[k5[0] + 1]), cl=float(co["Cl"][k5[0] + 1]))
        fig, axs = plt.subplots(3, 1, figsize=(10, 8.5), sharex=True)
        for ax, (k, lab) in zip(axs, (("Cl", "Cl"), ("Cd", "Cd"), ("CmPitch", "Cm (c/4, nose-up +)"))):
            ax.plot(it, co[k], color=BLUE, lw=1.4)
            tail = co[k][len(co[k]) // 4:]
            lo_, hi_ = tail.min(), tail.max()
            pad = 0.15 * (hi_ - lo_) + 1e-4
            ax.set_ylim(lo_ - pad, hi_ + pad)
            ax.set_ylabel(lab)
            mean = co[k][w].mean()
            ax.axvspan(it[w][0], it[w][-1], color=GRID, alpha=0.5, lw=0)
            ax.axhline(mean, color=MUTED, lw=1, ls="--")
            ax.annotate(f"mean {mean:.4f} (shaded window)", xy=(1, mean), xycoords=("axes fraction", "data"),
                        xytext=(-4, 5), textcoords="offset points", ha="right", color=INK, fontsize=9)
            style(ax)
        axs[0].set_title("Force coefficients", loc="left")
        axs[-1].set_xlabel("Iteration")
        save(fig, "force_coefficients.png")

    # ---- fields
    if T is None:
        print("No solution time found")
        return
    S = surface_data(M, T)
    tr = {s: transition(S[s]) for s in ("upper", "lower")}
    summary["transition"] = tr
    summary["yplus"] = {s: dict(max=float(S[s]["yplus"].max()), mean=float(S[s]["yplus"].mean()))
                        for s in ("upper", "lower")}
    summary["base_cp"] = S["base_cp"]

    fig, ax = surface_plot(S, "cp", "Cp", "cp_distribution.png",
                           "Pressure coefficient over the airfoil", invert=True)
    for s, color in (("upper", BLUE), ("lower", ORANGE)):
        xo = tr[s]["transition_x_c"]
        if xo is None:
            continue
        ax.axvline(xo, color=color, lw=1, ls=":")
        ax.annotate(f"{s} transition x/c = {xo:.3f}", xy=(xo, 0.04 if s == "upper" else 0.12),
                    xycoords=("data", "axes fraction"), xytext=(4, 0), textcoords="offset points",
                    color=INK, fontsize=9)
    save(fig, "cp_distribution.png")

    fig, ax = surface_plot(S, "cf", "Cf", "cf_distribution.png",
                           "Skin friction coefficient (dotted: transition points)",
                           hlines=((0.0, ""),))
    for s, color in (("upper", BLUE), ("lower", ORANGE)):
        if tr[s]["transition_x_c"] is not None:
            ax.axvline(tr[s]["transition_x_c"], color=color, lw=1, ls=":")
    lim = np.percentile(np.concatenate([S["upper"]["cf"], S["lower"]["cf"]]), [1, 99.5])
    ax.set_ylim(lim[0] - 0.002, lim[1] + 0.002)
    save(fig, "cf_distribution.png")

    fig, ax = surface_plot(S, "nut_ratio", "peak nut/nu across the boundary layer",
                           "transition_indicator.png",
                           "Boundary-layer eddy viscosity (laminar ~0, turbulent O(10))",
                           hlines=((NUT_TRANSITION, f"transition threshold {NUT_TRANSITION:g}"),))
    ax.set_yscale("log"); ax.set_ylim(1e-3, 300)
    save(fig, "transition_indicator.png")

    fig, ax = surface_plot(S, "yplus", "y+", "yplus_distribution.png",
                           "y+ of the first cell centre", hlines=((1.0, "y+ = 1"), (0.7, "target 0.7")))
    save(fig, "yplus_distribution.png")

    p = read_field(CASE, T, fname(T, "p"))
    U = read_field(CASE, T, fname(T, "U"))
    k = read_field(CASE, T, fname(T, "k"))
    summary["fields"] = ("averaged from iteration %d" % int(P["averageStart"])
                         if fname(T, "p") == "pMean" else "instantaneous (final iteration)")
    cp = p / Q
    field_plot(M, cp, "pressure_contours.png", "Pressure coefficient Cp", "RdBu_r",
               norm=TwoSlopeNorm(vcenter=0.0, vmin=-1.6, vmax=1.0),
               levels=np.linspace(-1.6, 1.0, 53))
    vmag = np.linalg.norm(U[:, :2], axis=1) / UINF
    field_plot(M, vmag, "velocity_contours.png", "Velocity magnitude |U|/U_inf", "viridis",
               levels=np.linspace(0, 1.6, 49))
    tu = np.sqrt(2.0 * np.maximum(k, 1e-14) / 3.0) / UINF
    lv = np.logspace(-4, 0, 41)
    field_plot(M, np.clip(tu, 1.01e-4, 0.99), "turbulence_intensity.png",
               "Turbulence intensity sqrt(2k/3)/U_inf (log scale)", "magma",
               norm=LogNorm(vmin=1e-4, vmax=1.0), levels=lv, ticks=[1e-4, 1e-3, 1e-2, 1e-1, 1])

    # ---- summary
    q = open(os.path.join(CASE, "constant", "meshQuality.txt")).read()
    summary["mesh_quality"] = q
    with open(os.path.join(OUT, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    lines = [
        f"# S1223, Re = {P['Re']:.0f}, alpha = {P['alpha']:g} deg, kOmegaSSTLM",
        "",
        f"Run status: {status}  ",
        f"Iterations: {summary['iterations']}",
        "",
        f"Coefficients averaged over iterations {summary['average_window'][0]}-"
        f"{summary['average_window'][1]}; fields: {summary['fields']}.",
        "",
        "| Coefficient | Mean | Std | Min .. max in window | Final iteration |",
        "|---|---|---|---|---|",
    ]
    for lab, name in (("cl", "Cl"), ("cd", "Cd"), ("cm", "Cm (c/4, nose-up +)")):
        if lab in summary:
            lines.append(f"| {name} | {summary[lab]:.4f} | {summary[lab+'_std']:.4f} | "
                         f"{summary[lab+'_min']:.4f} .. {summary[lab+'_max']:.4f} | {summary[lab+'_final']:.4f} |")
    lines += ["", f"Elements: {nq + nt} ({nq} quadrilaterals, {nt} triangles)", "",
              "| Transition | Upper surface | Lower surface |", "|---|---|---|"]
    f = lambda v: "laminar to TE" if v is None else f"{v:.3f}"
    lines.append(f"| transition (BL nut/nu > {NUT_TRANSITION:g}), x/c | "
                 f"{f(tr['upper']['transition_x_c'])} | {f(tr['lower']['transition_x_c'])} |")
    for key, lab in (("separation_x_c", "laminar separation (Cf < 0), x/c"),
                     ("reattachment_x_c", "reattachment, x/c")):
        u, l = tr["upper"].get(key), tr["lower"].get(key)
        if u is not None or l is not None:
            g = lambda v: "-" if v is None else f"{v:.3f}"
            lines.append(f"| {lab} | {g(u)} | {g(l)} |")
    lines += ["", f"y+ max / mean: upper {summary['yplus']['upper']['max']:.2f} / "
                  f"{summary['yplus']['upper']['mean']:.2f}, lower {summary['yplus']['lower']['max']:.2f} / "
                  f"{summary['yplus']['lower']['mean']:.2f}", "", "Mesh quality:", "", "```", q.rstrip(), "```"]
    if "literal_5pct_rule" in summary:
        r = summary["literal_5pct_rule"]
        lines += ["", f"The literal rule 'stop when |dCl/Cl| between two consecutive iterations < 5 %' "
                      f"would have stopped at iteration {r['iteration']} with Cl = {r['cl']:.4f}."]
    with open(os.path.join(OUT, "summary.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
