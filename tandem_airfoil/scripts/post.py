"""
Post-processing of the tandem-airfoil case: coefficients, transition point and
PNG figures (Cp, contours, y+, residuals, pathlines, mesh).

usage: python3 post.py <caseDir> <resultsDir> <geom.npz>
"""
import glob
import json
import math
import os
import re
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

import params as P  # noqa: E402
from run_solver import read_all_coeffs  # noqa: E402

Q_INF = 0.5 * P.U_INF ** 2          # kinematic dynamic pressure (p is p/rho)


# =================================================================== readers
def _body(path):
    with open(path) as f:
        txt = f.read()
    i = txt.find("FoamFile")
    if i >= 0:
        txt = txt[txt.index("}", i) + 1:]
    return re.sub(r"//.*", "", txt)


def _list_block(txt, start=0):
    """Parse 'N ( ... )' starting at the first count after `start`."""
    m = re.compile(r"(\d+)\s*\(").search(txt, start)
    n, i = int(m.group(1)), m.end()
    depth, j = 1, i
    while depth:
        c = txt[j]
        if c == "(":
            depth += 1
        elif c == ")":
            depth -= 1
        j += 1
    return n, txt[i:j - 1], j


def read_points(path):
    n, s, _ = _list_block(_body(path))
    return np.array(s.replace("(", " ").replace(")", " ").split(), float).reshape(n, 3)


def read_labels(path):
    n, s, _ = _list_block(_body(path))
    return np.array(s.split(), int)


def read_faces(path):
    txt = _body(path)
    with open(path) as f:
        head = f.read(2000)
    if "faceCompactList" in head:
        n1, s1, j = _list_block(txt)
        n2, s2, _ = _list_block(txt, j)
        off, lab = np.array(s1.split(), int), np.array(s2.split(), int)
        return [lab[off[k]:off[k + 1]] for k in range(len(off) - 1)]
    _, s, _ = _list_block(txt)
    return [np.array(g.split(), int) for g in re.findall(r"\d+\s*\(([^)]*)\)", s)]


def read_boundary(path):
    txt = _body(path)
    out = {}
    for m in re.finditer(r"(\w+)\s*\{([^}]*)\}", txt):
        d = dict(re.findall(r"(\w+)\s+([^;]+);", m.group(2)))
        if "nFaces" in d:
            out[m.group(1)] = (int(d["startFace"]), int(d["nFaces"]))
    return out


def _parse_values(txt, ncomp):
    s = txt.replace("(", " ").replace(")", " ")
    a = np.array(s.split(), float)
    return a.reshape(-1, ncomp) if ncomp > 1 else a


def read_field(path, ncells):
    """Return (internal values, {patch: values}) for an ascii vol field."""
    txt = _body(path)
    vec = "volVectorField" in open(path).read(1500)
    nc = 3 if vec else 1
    i = txt.index("internalField")
    seg = txt[i: i + 200]
    if "nonuniform" in seg:
        _, s, _ = _list_block(txt, i)
        internal = _parse_values(s, nc)
    else:
        v = re.search(r"uniform\s+(\([^)]*\)|[-+.\deE]+)", seg).group(1)
        internal = np.tile(_parse_values(v, nc), (ncells, 1)) if vec else np.full(ncells, float(v))
    patches = {}
    bf = txt.index("boundaryField")
    for m in re.finditer(r"\n\s{4}(\w+)\s*\n\s{4}\{", txt[bf:]):
        name = m.group(1)
        k0 = bf + m.end()
        k1 = txt.index("\n    }", k0)
        blk = txt[k0:k1]
        vm = re.search(r"value\s+nonuniform", blk)
        if vm:
            _, s, _ = _list_block(blk, vm.end())
            patches[name] = _parse_values(s, nc)
        else:
            vm = re.search(r"value\s+uniform\s+(\([^)]*\)|[-+.\deE]+)", blk)
            if vm:
                patches[name] = _parse_values(vm.group(1), nc)
    return internal, patches


# =================================================================== mesh
class Mesh2D:
    def __init__(self, case):
        pm = os.path.join(case, "constant", "polyMesh")
        self.points = read_points(os.path.join(pm, "points"))
        self.faces = read_faces(os.path.join(pm, "faces"))
        self.owner = read_labels(os.path.join(pm, "owner"))
        self.boundary = read_boundary(os.path.join(pm, "boundary"))
        self.ncells = self.owner.max() + 1
        nb = os.path.join(pm, "neighbour")
        if os.path.exists(nb):
            self.ncells = max(self.ncells, read_labels(nb).max() + 1)
        s, n = self.boundary["frontAndBack"]
        z0 = self.points[:, 2].min()
        polys, cells = [], []
        for f in range(s, s + n):
            fv = self.faces[f]
            if np.all(np.abs(self.points[fv, 2] - z0) < 1e-9):
                polys.append(fv)
                cells.append(self.owner[f])
        self.poly_cell = np.array(cells)
        used = np.unique(np.concatenate(polys))
        remap = -np.ones(len(self.points), int)
        remap[used] = np.arange(len(used))
        self.xy = self.points[used, :2]
        self.polys = [remap[p] for p in polys]
        tris, tcell = [], []
        for p, c in zip(self.polys, self.poly_cell):
            for k in range(1, len(p) - 1):
                tris.append((p[0], p[k], p[k + 1]))
                tcell.append(c)
        self.tri = mtri.Triangulation(self.xy[:, 0], self.xy[:, 1], np.array(tris))
        # cell centres & cell->point averaging operator
        cc = np.zeros((self.ncells, 2))
        cnt = np.zeros(self.ncells)
        for p, c in zip(self.polys, self.poly_cell):
            cc[c] += self.xy[p].mean(axis=0)
            cnt[c] += 1
        self.cc = cc / np.maximum(cnt, 1)[:, None]
        rows = np.concatenate(self.polys)
        cols = np.concatenate([[c] * len(p) for p, c in zip(self.polys, self.poly_cell)])
        self._pr, self._pc = rows, cols
        self._pn = np.bincount(rows, minlength=len(self.xy))

    def to_points(self, cellvals):
        return np.bincount(self._pr, weights=cellvals[self._pc], minlength=len(self.xy)) / self._pn

    def patch(self, name):
        """Wall-face centres (x, y), outward unit normals (into the wall), owner cells."""
        s, n = self.boundary[name]
        fc, nrm, own = [], [], []
        for f in range(s, s + n):
            fv = self.points[self.faces[f]]
            c = fv.mean(axis=0)
            nv = np.cross(fv[1] - fv[0], fv[2] - fv[0])
            fc.append(c[:2])
            nrm.append(nv[:2] / np.linalg.norm(nv[:2]))
            own.append(self.owner[f])
        return np.array(fc), np.array(nrm), np.array(own)


def latest_time(case):
    ts = []
    for d in os.listdir(case):
        try:
            ts.append((float(d), d))
        except ValueError:
            pass
    ts = [t for t in ts if t[0] > 0]
    return max(ts)[1] if ts else None


# =================================================================== helpers
def annotate(fig, ncells, extra=""):
    fig.text(0.5, 0.005, P.label(ncells) + extra, ha="center", va="bottom", fontsize=8.5,
             bbox=dict(boxstyle="round,pad=0.3", fc="#f4f4f4", ec="#bbbbbb"))


def draw_airfoils(ax, geom, **kw):
    kw = dict(dict(color="k", lw=0.8, zorder=5), **kw)
    for a, b in (("up1", "lo1"), ("up2", "lo2")):
        poly = np.vstack([geom[a], geom[b][::-1]])
        ax.fill(poly[:, 0], poly[:, 1], color="#555555", zorder=4)
        ax.plot(poly[:, 0], poly[:, 1], **kw)


def surfaces(mesh, geom, name, up, lo, x0, chord):
    """Split a wall patch into upper / lower ordered LE->TE."""
    fc, nrm, own = mesh.patch(name)
    tu, tl = cKDTree(geom[up]), cKDTree(geom[lo])
    du, _ = tu.query(fc)
    dl, _ = tl.query(fc)
    is_up = du <= dl
    res = {}
    for key, mask, ref in (("upper", is_up, geom[up]), ("lower", ~is_up, geom[lo])):
        idx = np.where(mask)[0]
        s_ref = np.r_[0, np.cumsum(np.hypot(*np.diff(ref, axis=0).T))]
        _, k = cKDTree(ref).query(fc[idx])
        order = idx[np.argsort(s_ref[k] + 1e-9 * np.arange(len(idx)))]
        # tangent LE -> TE
        pts = fc[order]
        t = np.gradient(pts, axis=0)
        t /= np.linalg.norm(t, axis=1)[:, None]
        res[key] = dict(idx=order, xy=pts, xc=(pts[:, 0] - x0) / chord, t=t,
                        s=np.r_[0, np.cumsum(np.hypot(*np.diff(pts, axis=0).T))])
    return res, fc, nrm, own


# =================================================================== main
def main():
    case, out, geomfile = sys.argv[1], sys.argv[2], sys.argv[3]
    os.makedirs(out, exist_ok=True)
    geom = dict(np.load(geomfile))
    mesh = Mesh2D(case)
    N = mesh.ncells
    t = latest_time(case)
    print("post-processing time", t, "cells", N)
    tdir = os.path.join(case, t)
    Uc, Ub = read_field(os.path.join(tdir, "U"), N)
    pc, _ = read_field(os.path.join(tdir, "p"), N)
    kc, _ = read_field(os.path.join(tdir, "k"), N)
    gc, _ = read_field(os.path.join(tdir, "gammaInt"), N)
    nutc, _ = read_field(os.path.join(tdir, "nut"), N)
    yp = read_field(os.path.join(tdir, "yPlus"), N)[1] if os.path.exists(os.path.join(tdir, "yPlus")) else {}
    wss = read_field(os.path.join(tdir, "wallShearStress"), N)[1] \
        if os.path.exists(os.path.join(tdir, "wallShearStress")) else {}
    status = {}
    sp = os.path.join(case, "run_status.json")
    if os.path.exists(sp):
        status = json.load(open(sp))

    results = dict(case=dict(
        airfoils="S1223 (c=%.2f m) + NACA 0009 (c=%.2f m), gap %.2f m" % (P.C1, P.C2, P.GAP),
        Re=P.RE, U_inf=P.U_INF, AoA_deg=P.AOA_DEG, rho=P.RHO, nu=P.NU, Tu=P.TU,
        turbulence_model="kOmegaSSTLM (k-omega SST + gamma-ReTheta transition)",
        cells=int(N), quads_2d=int(geom["nquad"]), first_layer_S1223=P.H1_1,
        first_layer_NACA0009=P.H1_2, iterations=int(float(t)), run=status))

    # ---------------------------------------------------------- coefficients
    coeffs = {}
    for name in ("coeffs_S1223", "coeffs_NACA0009", "coeffs_total"):
        d = read_all_coeffs(case, name)
        if not d:
            continue
        # drop iterations beyond the saved time (e.g. from an interrupted run)
        keep = [k for k, tt in enumerate(d["Time"]) if tt <= float(t) + 1e-9]
        d = {c: [v[k] for k in keep] for c, v in d.items()}
        cm_key = "CmPitch" if "CmPitch" in d else "Cm"
        n = len(d["Cl"])
        w = min(100, n)
        coeffs[name] = dict(
            Cl=d["Cl"][-1], Cd=d["Cd"][-1], Cm=d[cm_key][-1],
            Cl_mean_last100=float(np.mean(d["Cl"][-w:])), Cl_std_last100=float(np.std(d["Cl"][-w:])),
            Cd_mean_last100=float(np.mean(d["Cd"][-w:])), Cd_std_last100=float(np.std(d["Cd"][-w:])),
            Cm_mean_last100=float(np.mean(d[cm_key][-w:])), history=d)
    results["coefficients"] = {k: {kk: vv for kk, vv in v.items() if kk != "history"}
                               for k, v in coeffs.items()}

    # ---------------------------------------------------------- surface data
    surf = {}
    for name, up, lo, x0, c, lab in (("airfoil1", "up1", "lo1", P.X_LE1, P.C1, "S1223"),
                                     ("airfoil2", "up2", "lo2", P.X_LE2, P.C2, "NACA 0009")):
        sd, fc, nrm, own = surfaces(mesh, geom, name, up, lo, x0, c)
        # maximum of nut/nu and of the intermittency in the inner part of the
        # boundary layer (cells within TR_BAND of the wall, attributed to their
        # nearest wall face). A thin band keeps the free stream / the S1223
        # wake (gamma = 1, nut/nu > 1) out of the measure.
        dw, kw = cKDTree(fc).query(mesh.cc)
        sel = dw < P.TR_BAND
        ratio_max = np.zeros(len(fc))
        np.maximum.at(ratio_max, kw[sel], nutc[sel] / P.NU)
        gam_max = np.zeros(len(fc))
        np.maximum.at(gam_max, kw[sel], gc[sel])
        cp_all = pc[own] / Q_INF
        g_all = gc[own]
        tau = wss.get(name)
        for side in ("upper", "lower"):
            d = sd[side]
            d["cp"] = cp_all[d["idx"]]
            d["gamma"] = g_all[d["idx"]]
            d["nutr"] = ratio_max[d["idx"]]
            d["gmax"] = gam_max[d["idx"]]
            d["yplus"] = yp[name][d["idx"]] if name in yp else None
            if tau is not None:
                # OpenFOAM wallShearStress = force per area exerted by the fluid
                # (kinematic); Cf positive for attached flow LE -> TE
                tv = tau[d["idx"], :2]
                d["cf"] = -np.sum(tv * d["t"], axis=1) / Q_INF
                if np.median(d["cf"]) < 0:
                    d["cf"] = -d["cf"]
            else:
                d["cf"] = None
        surf[lab] = sd

    # ---------------------------------------------------------- transition
    def transition(d):
        x, r, gm = d["xc"], d["nutr"], d["gmax"]
        out = {}
        m = x > 0.02
        # transition onset: inner-BL intermittency > 0.5 AND turbulent viscosity
        # larger than the molecular one (laminar BL: gamma ~ 0.02, nut/nu < 0.1)
        k = np.where(m & (r > P.NUT_RATIO_TR) & (gm > 0.5))[0]
        out["x_transition"] = float(x[k[0]]) if len(k) else None
        out["max_nut_over_nu_in_BL"] = float(r.max())
        if d["cf"] is not None:
            cf = d["cf"]
            neg = np.where(m & (cf < 0))[0]
            if len(neg):
                out["x_laminar_separation"] = float(x[neg[0]])
                # reattachment: Cf back to positive and staying positive for 5 % chord
                out["x_reattachment"] = None
                for j in range(neg[0] + 1, len(cf)):
                    if cf[j] > 0 and np.all(cf[(x >= x[j]) & (x <= x[j] + 0.05)] > 0):
                        out["x_reattachment"] = float(x[j])
                        break
            else:
                out["x_laminar_separation"] = None
                out["x_reattachment"] = None
            # separation-induced transition: the gamma-ReTheta model triggers it
            # through its internal separation intermittency, not the transported
            # gamma field, so if a bubble closes before the gamma onset the
            # transition point is taken at the Cf minimum inside the bubble
            out["method"] = "gamma>0.5 & nut/nu>%g" % P.NUT_RATIO_TR
            xs, xr = out["x_laminar_separation"], out["x_reattachment"]
            if xs is not None and xr is not None and \
                    (out["x_transition"] is None or out["x_transition"] > xr):
                inb = (x >= xs) & (x <= xr)
                out["x_transition"] = float(x[inb][np.argmin(cf[inb])])
                out["method"] = "bubble: Cf minimum"
            # onset of the turbulent Cf rise: minimum Cf before the largest Cf rise
            dcf = np.gradient(cf, x)
            mm = m & (x < 0.98)
            if mm.any():
                j = np.argmax(np.where(mm, dcf, -np.inf))
                jmin = np.argmin(np.where(mm & (np.arange(len(cf)) <= j), cf, np.inf))
                out["x_cf_min_before_rise"] = float(x[jmin])
                out["x_max_cf_rise"] = float(x[j])
        return out

    tr = {lab: {side: transition(surf[lab][side]) for side in ("upper", "lower")} for lab in surf}
    results["transition"] = tr

    # ---------------------------------------------------------- y+
    results["yplus"] = {lab: dict(max=float(max(np.max(surf[lab][s]["yplus"]) for s in surf[lab])),
                                  mean=float(np.mean(np.r_[surf[lab]["upper"]["yplus"],
                                                           surf[lab]["lower"]["yplus"]])))
                        for lab in surf if surf[lab]["upper"]["yplus"] is not None}

    with open(os.path.join(out, "results.json"), "w") as f:
        json.dump(results, f, indent=2, default=lambda o: None)

    # =============================================================== figures
    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.3})
    c_up, c_lo = "#1f5fbf", "#d1495b"

    # ---- Cp
    fig, axs = plt.subplots(1, 2, figsize=(15, 6.2))
    for ax, lab in zip(axs, surf):
        for side, col in (("upper", c_up), ("lower", c_lo)):
            d = surf[lab][side]
            ax.plot(d["xc"], d["cp"], "-", color=col, lw=1.6, label="%s surface" % side)
        ax.invert_yaxis()
        ax.set_xlabel("x/c")
        ax.set_ylabel("Cp")
        key = "coeffs_S1223" if lab == "S1223" else "coeffs_NACA0009"
        cl = coeffs[key]["Cl"] if key in coeffs else float("nan")
        cd = coeffs[key]["Cd"] if key in coeffs else float("nan")
        cm = coeffs[key]["Cm"] if key in coeffs else float("nan")
        ttl = "%s   Cl=%.4f  Cd=%.5f  Cm(c/4)=%.4f" % (lab, cl, cd, cm)
        for side, col in (("upper", c_up), ("lower", c_lo)):
            xt = tr[lab][side].get("x_transition")
            if xt is not None:
                ax.axvline(xt, color=col, ls=":", lw=1)
                yl = ax.get_ylim()
                ax.text(xt, yl[0] + (0.06 if side == "upper" else 0.16) * (yl[1] - yl[0]),
                        " transition %s x/c=%.3f" % (side, xt), color=col, fontsize=8.5)
        ax.set_title(ttl, fontsize=10.5)
        ax.legend(loc="lower right")
    fig.suptitle("Pressure coefficient distribution (dotted: transition onset, inner-BL intermittency > 0.5 and nut/nu > %g)" % P.NUT_RATIO_TR)
    fig.tight_layout(rect=(0, 0.06, 1, 0.97))
    annotate(fig, N)
    fig.savefig(os.path.join(out, "01_pressure_distribution_Cp.png"), dpi=150)
    plt.close(fig)

    # ---- contour helper
    def contour(fname, cellvals, title, cbar, cmap, levels=None, log=False):
        pv = mesh.to_points(cellvals)
        fig = plt.figure(figsize=(15, 12))
        gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.15])
        views = [(gs[0, 0], (-5.8, 9.4), (-6.1, 6.1), "whole C-domain"),
                 (gs[0, 1], (0.2, 0.75), (-0.14, 0.1), "trailing edge, gap and NACA 0009"),
                 (gs[1, :], (-0.15, 1.0), (-0.28, 0.25), "near field")]
        for k, (g, xl, yl, nm) in enumerate(views):
            ax = fig.add_subplot(g)
            if log:
                pv_ = np.clip(pv, levels[0], levels[-1])
                cs = ax.tricontourf(mesh.tri, pv_, levels=levels, cmap=cmap, norm=LogNorm(),
                                    extend="neither")
            else:
                cs = ax.tricontourf(mesh.tri, np.clip(pv, levels[0], levels[-1]), levels=levels,
                                    cmap=cmap, extend="both")
            if k > 0:
                ax.tricontour(mesh.tri, np.clip(pv, levels[0], levels[-1]), levels=levels[::4],
                              colors="k", linewidths=0.25, norm=LogNorm() if log else None)
            draw_airfoils(ax, geom)
            ax.set_xlim(*xl)
            ax.set_ylim(*yl)
            ax.set_aspect("equal")
            ax.grid(False)
            ax.set_title(nm)
            ax.set_xlabel("x [m]")
            ax.set_ylabel("y [m]")
            cb = fig.colorbar(cs, ax=ax, fraction=0.025, pad=0.01)
            cb.set_label(cbar)
        fig.suptitle(title, fontsize=13)
        fig.tight_layout(rect=(0, 0.05, 1, 0.97))
        annotate(fig, N)
        fig.savefig(os.path.join(out, fname), dpi=150)
        plt.close(fig)

    cp_cells = pc / Q_INF
    contour("02_pressure_contours.png", cp_cells, "Pressure coefficient Cp = (p - p_inf)/q_inf",
            "Cp", "RdBu_r", levels=np.linspace(-2.0, 1.0, 31))
    umag = np.linalg.norm(Uc, axis=1)
    contour("03_velocity_contours.png", umag, "Velocity magnitude |U| [m/s]", "|U| [m/s]",
            "viridis", levels=np.linspace(0, 1.8 * P.U_INF, 37))
    ti = np.sqrt(2.0 / 3.0 * np.maximum(kc, 1e-14)) / P.U_INF * 100.0
    contour("04_turbulence_intensity_log.png", ti,
            "Turbulence intensity Tu = sqrt(2k/3)/U_inf (logarithmic scale)",
            "Tu [%] (log scale)", "magma", levels=np.logspace(-2, 2, 33), log=True)

    # ---- y+
    if results.get("yplus"):
        fig, axs = plt.subplots(1, 2, figsize=(15, 5.8))
        for ax, lab in zip(axs, surf):
            for side, col in (("upper", c_up), ("lower", c_lo)):
                d = surf[lab][side]
                ax.plot(d["xc"], d["yplus"], "-", color=col, lw=1.4, label="%s surface" % side)
            ax.axhline(1.0, color="k", ls="--", lw=0.8, label="y+ = 1")
            ax.axhline(P.YPLUS_TARGET, color="gray", ls=":", lw=0.8,
                       label="design target %.1f" % P.YPLUS_TARGET)
            ax.set_xlabel("x/c")
            ax.set_ylabel("y+ (wall-adjacent cell centre)")
            ax.set_title("%s   max y+ = %.3f, mean y+ = %.3f   (first layer %.2e m)"
                         % (lab, results["yplus"][lab]["max"], results["yplus"][lab]["mean"],
                            P.H1_1 if lab == "S1223" else P.H1_2), fontsize=10.5)
            ax.legend()
        fig.suptitle("y+ distribution over the airfoils")
        fig.tight_layout(rect=(0, 0.06, 1, 0.97))
        annotate(fig, N)
        fig.savefig(os.path.join(out, "05_yplus_distribution.png"), dpi=150)
        plt.close(fig)

    # ---- residuals
    log = os.path.join(case, "log.simpleFoam")
    if os.path.exists(log):
        res, it = {}, 0
        seen = set()
        for line in open(log, errors="replace"):
            if line.startswith("Time = "):
                it = int(float(line.split()[2]))
                seen = set()
                continue
            m = re.search(r"Solving for (\w+), Initial residual = ([-+.\deE]+)", line)
            if m and m.group(1) not in seen:
                seen.add(m.group(1))
                res.setdefault(m.group(1), {})[it] = float(m.group(2))
        # a restart repeats iterations: keep the latest value of each
        res = {k: ([i for i in sorted(v) if i <= float(t)],
                   [v[i] for i in sorted(v) if i <= float(t)]) for k, v in res.items()}
        fig, ax = plt.subplots(figsize=(12, 6.2))
        for k, (x, y) in res.items():
            ax.semilogy(x, y, lw=1.1, label=k)
        ax.set_xlabel("Iteration")
        ax.set_ylabel("Initial residual")
        restarts = sorted({int(float(os.path.basename(os.path.dirname(f))))
                           for f in glob.glob(os.path.join(case, "postProcessing", "coeffs_S1223", "*", "*.dat"))
                           if os.path.basename(os.path.dirname(f)).replace(".", "").isdigit()} - {0})
        for rs in restarts:
            ax.axvline(rs, color="gray", ls=":", lw=1)
            ax.text(rs, 2e-7, " restart from\n saved it. %d" % rs, fontsize=8, color="gray")
        st = status.get("status", "")
        ax.set_title("Residuals (simpleFoam, SIMPLEC)  -  run status: %s %s"
                     % (st, ("(" + status.get("reason", "") + ")") if st else ""), fontsize=10)
        ax.legend(ncol=4)
        if "coeffs_S1223" in coeffs:
            ax2 = ax.twinx()
            h = coeffs["coeffs_S1223"]["history"]
            ax2.plot(h["Time"], h["Cl"], "k--", lw=1.0, label="Cl S1223")
            ax2.set_ylabel("Cl S1223 (dashed)")
            ax2.set_ylim(min(0, min(h["Cl"][len(h["Cl"]) // 5:]) * 1.2),
                         max(h["Cl"][len(h["Cl"]) // 5:]) * 1.3)
            ax2.grid(False)
        fig.tight_layout(rect=(0, 0.07, 1, 1))
        annotate(fig, N)
        fig.savefig(os.path.join(out, "06_residuals.png"), dpi=150)
        plt.close(fig)

    # ---- pathlines
    ux = mtri.LinearTriInterpolator(mesh.tri, mesh.to_points(Uc[:, 0]))
    uy = mtri.LinearTriInterpolator(mesh.tri, mesh.to_points(Uc[:, 1]))

    def vel(p):
        a, b = ux(p[0], p[1]), uy(p[0], p[1])
        if np.ma.is_masked(a) or np.ma.is_masked(b):
            return None
        return np.array([float(a), float(b)])

    x_seed, x_end = -0.20, P.X_TE1 + 0.60
    seeds = np.r_[np.linspace(-0.15, -0.02, 7), np.linspace(-0.012, 0.03, 12),
                  np.linspace(0.04, 0.16, 7)]
    lines = []
    for y0 in seeds:
        p = np.array([x_seed, y0])
        pts = [p.copy()]
        ds = 0.0015
        for _ in range(6000):
            v1 = vel(p)
            if v1 is None:
                break
            sp1 = np.linalg.norm(v1)
            if sp1 < 1e-6:
                break
            pm = p + 0.5 * ds * v1 / sp1
            v2 = vel(pm)
            if v2 is None:
                break
            p = p + ds * v2 / np.linalg.norm(v2)
            pts.append(p.copy())
            if p[0] > x_end + 0.02:
                break
        lines.append(np.array(pts))
    fig, ax = plt.subplots(figsize=(15, 6.8))
    pv = mesh.to_points(umag)
    cs = ax.tricontourf(mesh.tri, np.clip(pv, 0, 1.8 * P.U_INF), levels=np.linspace(0, 1.8 * P.U_INF, 37),
                        cmap="Blues", alpha=0.85)
    fig.colorbar(cs, ax=ax, fraction=0.02, pad=0.01).set_label("|U| [m/s]")
    for L in lines:
        ax.plot(L[:, 0], L[:, 1], color="#c0392b", lw=0.9)
        ax.plot(L[0, 0], L[0, 1], "o", color="#c0392b", ms=2.5)
        for xa in (-0.12, 0.75):
            k = np.searchsorted(L[:, 0], xa) if np.all(np.diff(L[:, 0]) > 0) else \
                int(np.argmin(np.abs(L[:, 0] - xa)))
            if 0 < k < len(L) - 1:
                ax.annotate("", xy=L[k + 1], xytext=L[k - 1],
                            arrowprops=dict(arrowstyle="-|>", color="#c0392b", lw=0.8))
    draw_airfoils(ax, geom)
    ax.axvline(x_end, color="k", ls=":", lw=0.8)
    ax.text(x_end - 0.005, 0.18, "60 cm behind S1223 TE ", fontsize=8, ha="right")
    ax.set_xlim(x_seed - 0.03, x_end + 0.03)
    ax.set_ylim(-0.19, 0.2)
    ax.set_aspect("equal")
    ax.set_title("Pathlines of %d massless particles released at x = %.2f m (steady flow)"
                 % (len(lines), x_seed))
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    annotate(fig, N)
    fig.savefig(os.path.join(out, "07_pathlines.png"), dpi=150)
    plt.close(fig)

    # ---- mesh
    segs = []
    for p in mesh.polys:
        q = mesh.xy[np.r_[p, p[0]]]
        segs.append(q)
    fig, axs = plt.subplots(2, 2, figsize=(15, 11))
    views = [((-5.8, 9.4), (-6.1, 6.1), "whole C-domain (wake 30c)"),
             ((-0.1, 0.75), (-0.2, 0.2), "both airfoils"),
             ((-0.012, 0.02), (-0.012, 0.012), "S1223 leading edge"),
             ((0.27, 0.33), (-0.03, 0.025), "S1223 trailing edge")]
    for ax, (xl, yl, nm) in zip(axs.ravel(), views):
        ax.add_collection(LineCollection(segs, colors="k", linewidths=0.15 if "whole" not in nm else 0.1))
        draw_airfoils(ax, geom, lw=0.5)
        ax.set_xlim(*xl)
        ax.set_ylim(*yl)
        ax.set_aspect("equal")
        ax.grid(False)
        ax.set_title(nm)
    fig.suptitle("Hybrid C-mesh: %d cells (%d quads), 4 cm structured quad layer, "
                 "first cell %.2e m / %.2e m" % (N, int(geom["nquad"]), P.H1_1, P.H1_2))
    fig.tight_layout(rect=(0, 0.05, 1, 0.97))
    annotate(fig, N)
    fig.savefig(os.path.join(out, "08_mesh.png"), dpi=150)
    plt.close(fig)

    mesh3d(mesh, geom, N, os.path.join(out, "09_mesh_3d.png"))

    # ---------------------------------------------------------- report
    write_report(out, results)
    print(json.dumps(results["coefficients"], indent=2))
    print(json.dumps(results["transition"], indent=2))
    print(json.dumps(results.get("yplus"), indent=2))


def mesh3d(mesh, geom, ncells, path, span=0.30, nspan=12):
    """3-D view: both wings extruded over an illustrative span with their
    surface mesh, and the 2-D computational mesh on the symmetry plane."""
    from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection

    xl, yl = (-0.06, 0.70), (-0.13, 0.13)
    segs = []
    for p in mesh.polys:
        q = mesh.xy[np.r_[p, p[0]]]
        if q[:, 0].max() < xl[0] or q[:, 0].min() > xl[1] or \
                q[:, 1].max() < yl[0] or q[:, 1].min() > yl[1]:
            continue
        segs.append(np.c_[q, np.zeros(len(q))])
    zs = np.linspace(0.0, span, nspan + 1)
    light = np.array([0.3, 0.8, -0.5])
    light /= np.linalg.norm(light)
    fig = plt.figure(figsize=(16, 10))
    views = [(fig.add_subplot(1, 1, 1, projection="3d"), 24, -58)]
    for ax, elev, azim in views:
        ax.add_collection3d(Line3DCollection(segs, colors="#2b2b2b", linewidths=0.12, alpha=0.8))
        for a, b, base in (("up1", "lo1", "#8fb3e0"), ("up2", "lo2", "#e59a6b")):
            for surf in (geom[a], geom[b]):
                polys, cols = [], []
                for i in range(len(surf) - 1):
                    p0, p1 = surf[i], surf[i + 1]
                    tvec = np.r_[p1 - p0, 0.0]
                    nvec = np.cross(tvec, [0, 0, 1.0])
                    nvec /= max(np.linalg.norm(nvec), 1e-12)
                    shade = 0.45 + 0.55 * abs(nvec @ light)
                    for k in range(nspan):
                        polys.append([(p0[0], p0[1], zs[k]), (p1[0], p1[1], zs[k]),
                                      (p1[0], p1[1], zs[k + 1]), (p0[0], p0[1], zs[k + 1])])
                        c = np.array(matplotlib.colors.to_rgb(base)) * shade
                        cols.append(np.clip(c, 0, 1))
                pc = Poly3DCollection(polys, facecolors=cols, edgecolors="#1a1a1a",
                                      linewidths=0.08)
                ax.add_collection3d(pc)
            # end caps outline at the far tip
            ring = np.vstack([geom[a], geom[b][::-1]])
            ax.plot(ring[:, 0], ring[:, 1], span, color="k", lw=0.8)
        ax.set_xlim(*xl)
        ax.set_ylim(*yl)
        ax.set_zlim(0, span)
        ax.set_box_aspect((xl[1] - xl[0], yl[1] - yl[0], span), zoom=1.45)
        ax.view_init(elev=elev, azim=azim)
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_zlabel("span z [m]")
        ax.grid(False)
    fig.suptitle("3-D view of the wing meshes: S1223 (blue, %d+%d surface elements) and NACA 0009 "
                 "(orange, %d+%d), with the 2-D computational mesh on the z = 0 plane\n"
                 "(the CFD mesh is one cell deep; the wings are drawn over a %.0f cm span for "
                 "illustration)" % (P.N_SURF1, P.N_SURF1, P.N_SURF2, P.N_SURF2, span * 100),
                 fontsize=11)
    fig.subplots_adjust(left=0, right=1, bottom=0.04, top=0.93)
    annotate(fig, ncells)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def write_report(out, r):
    c = r["coefficients"]
    L = ["# Tandem airfoil results: S1223 + NACA 0009", ""]
    cs = r["case"]
    L += ["| item | value |", "|---|---|",
          "| configuration | %s |" % cs["airfoils"],
          "| flow | air ISA sea level, Re = %.0f, U = %.3f m/s, AoA = %.1f deg, Tu = %.1f %% |"
          % (cs["Re"], cs["U_inf"], cs["AoA_deg"], cs["Tu"] * 100),
          "| model | %s |" % cs["turbulence_model"],
          "| mesh | %d cells (%d quads) |" % (cs["cells"], cs["quads_2d"]),
          "| first layer | S1223 %.3e m, NACA 0009 %.3e m |" % (cs["first_layer_S1223"],
                                                             cs["first_layer_NACA0009"]),
          "| iterations | %d (%s) |" % (cs["iterations"], cs["run"].get("reason", "")), ""]
    L += ["## Aerodynamic coefficients", "",
          "| body | ref. chord | Cl | Cd | Cm (c/4, +nose-up) | Cl mean +- std (last 100 it.) |",
          "|---|---|---|---|---|---|"]
    for k, lab, ch in (("coeffs_S1223", "S1223", P.C1), ("coeffs_NACA0009", "NACA 0009", P.C2),
                       ("coeffs_total", "both (total)", P.C1)):
        if k in c:
            L.append("| %s | %.2f m | %.4f | %.5f | %.4f | %.4f +- %.4f |"
                     % (lab, ch, c[k]["Cl"], c[k]["Cd"], c[k]["Cm"], c[k]["Cl_mean_last100"],
                        c[k]["Cl_std_last100"]))
    L += ["", "## Transition (x/c from each leading edge)", "",
          "| airfoil | surface | transition onset | method | laminar separation (Cf<0) | reattachment |",
          "|---|---|---|---|---|---|"]
    f = lambda v: "-" if v is None else "%.3f" % v
    for lab, d in r["transition"].items():
        for side, v in d.items():
            L.append("| %s | %s | %s | %s | %s | %s |" % (lab, side, f(v.get("x_transition")),
                                                    v.get("method", "-"),
                                                    f(v.get("x_laminar_separation")),
                                                    f(v.get("x_reattachment"))))
    L += ["", "Transition onset: first station where, within %.1f mm of the wall, the "
          "intermittency > 0.5 and nut/nu > %g; when a laminar separation bubble closes "
          "before that, the Cf minimum inside the bubble (separation-induced transition)."
          % (P.TR_BAND * 1e3, P.NUT_RATIO_TR)]
    if r.get("yplus"):
        L += ["", "## y+", "", "| airfoil | max | mean |", "|---|---|---|"]
        for lab, v in r["yplus"].items():
            L.append("| %s | %.3f | %.3f |" % (lab, v["max"], v["mean"]))
    with open(os.path.join(out, "results.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
