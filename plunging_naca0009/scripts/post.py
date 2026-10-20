"""
Post-processing of the plunging NACA 0009 (last computed cycle).

Reads the force-coefficient history and the field snapshots of the last
plunge cycle (the case keeps SNAPSHOTS_PER_CYCLE of them) and writes:
coefficients (incl. thrust and propulsive efficiency), transition points vs
phase, y+, and PNG figures at four plunge phases.

Phase convention: y = A sin(phi), phi = omega t;
  phi =   0 deg : mid stroke, moving UP   at max speed (alpha_eff = -max)
  phi =  90 deg : top of the stroke       (alpha_eff = 0)
  phi = 180 deg : mid stroke, moving DOWN at max speed (alpha_eff = +max)
  phi = 270 deg : bottom of the stroke    (alpha_eff = 0)

usage: python3 post.py <unsteadyCase> <resultsDir> <geom.npz>
"""
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
from foamio import Mesh2D, read_field, time_dirs  # noqa: E402
from run_solver import cycle_change, cycle_stats, read_all_coeffs  # noqa: E402

Q_INF = 0.5 * P.U_INF ** 2
PHASES = (0, 90, 180, 270)
PHASE_NAME = {0: "mid-stroke, moving up", 90: "top of stroke", 180: "mid-stroke, moving down",
              270: "bottom of stroke"}
X_SEED, X_END = -0.20, P.X_TE + 0.60


# =================================================================== helpers
def phase_deg(t):
    return math.degrees(P.OMEGA * t) % 360.0


def alpha_eff(t):
    return -math.degrees(math.atan2(P.plunge(t)[1], P.U_INF))


def annotate(fig, ncells, extra=""):
    fig.text(0.5, 0.005, P.label(ncells) + extra, ha="center", va="bottom", fontsize=8.5,
             bbox=dict(boxstyle="round,pad=0.3", fc="#f4f4f4", ec="#bbbbbb"))


def airfoil_ring(geom, h=0.0):
    ring = np.vstack([geom["up"], geom["lo"][::-1]])
    return ring + [0.0, h]


def draw_airfoil(ax, geom, h=0.0, **kw):
    ring = airfoil_ring(geom, h)
    ax.fill(ring[:, 0], ring[:, 1], color=kw.pop("fill", "#555555"), zorder=4,
            alpha=kw.pop("alpha", 1.0))
    ax.plot(ring[:, 0], ring[:, 1], color=kw.pop("color", "k"), lw=kw.pop("lw", 0.8), zorder=5)


def surfaces(mesh, geom):
    """Wall faces split into upper / lower surface, ordered LE -> TE."""
    fc, nrm, own = mesh.patch("airfoil")
    tu, tl = cKDTree(geom["up"]), cKDTree(geom["lo"])
    is_up = tu.query(fc)[0] <= tl.query(fc)[0]
    res = {}
    for key, mask, ref in (("upper", is_up, geom["up"]), ("lower", ~is_up, geom["lo"])):
        idx = np.where(mask)[0]
        s_ref = np.r_[0, np.cumsum(np.hypot(*np.diff(ref, axis=0).T))]
        _, k = cKDTree(ref).query(fc[idx])
        order = idx[np.argsort(s_ref[k] + 1e-9 * np.arange(len(idx)))]
        pts = fc[order]
        t = np.gradient(pts, axis=0)
        t /= np.linalg.norm(t, axis=1)[:, None]
        res[key] = dict(idx=order, xy=pts, xc=(pts[:, 0] - P.X_LE) / P.C, t=t)
    return res, fc, own


def transition(x, cf, nutr, gmax, cp):
    """Transition onset and laminar bubble on one surface (x = x/c, LE -> TE).
    At large effective AoA the stagnation point moves onto the pressure-side
    surface; the reversed wall flow ahead of it is not a separation, so the
    search starts behind the stagnation point (Cp maximum)."""
    out = {}
    x_stag = float(x[np.argmax(cp)]) if cp.max() > 0.8 else 0.0
    out["x_stagnation"] = x_stag
    m = x > max(0.02, x_stag + 0.005)
    k = np.where(m & (nutr > P.NUT_RATIO_TR) & (gmax > 0.5))[0]
    out["x_transition"] = float(x[k[0]]) if len(k) else None
    out["method"] = "gamma>0.5 & nut/nu>%g" % P.NUT_RATIO_TR
    neg = np.where(m & (cf < 0))[0]
    out["x_laminar_separation"] = float(x[neg[0]]) if len(neg) else None
    out["x_reattachment"] = None
    if len(neg):
        for j in range(neg[0] + 1, len(cf)):
            if cf[j] > 0 and np.all(cf[(x >= x[j]) & (x <= x[j] + 0.05)] > 0):
                out["x_reattachment"] = float(x[j])
                break
    xs, xr = out["x_laminar_separation"], out["x_reattachment"]
    if xs is not None and xr is not None and (out["x_transition"] is None
                                              or out["x_transition"] > xr):
        inb = (x >= xs) & (x <= xr)
        out["x_transition"] = float(x[inb][np.argmin(cf[inb])])
        out["method"] = "bubble: Cf minimum"
    if out["x_transition"] is None and xs is not None and xr is None:
        out["method"] = "separated to the TE (no reattachment)"
    return out


# =================================================================== data
class Snapshot:
    def __init__(self, mesh, case, tname, sd, own, band):
        d = os.path.join(case, tname)
        N = mesh.ncells
        self.t = float(tname)
        self.h, self.v = P.plunge(self.t)
        self.phase = phase_deg(self.t)
        self.alpha = alpha_eff(self.t)
        self.U = read_field(os.path.join(d, "U"), N)[0][:, :2]
        self.p = read_field(os.path.join(d, "p"), N)[0]
        self.k = read_field(os.path.join(d, "k"), N)[0]
        nut = read_field(os.path.join(d, "nut"), N)[0]
        gam = read_field(os.path.join(d, "gammaInt"), N)[0]
        yp = read_field(os.path.join(d, "yPlus"), N)[1].get("airfoil") \
            if os.path.exists(os.path.join(d, "yPlus")) else None
        wss = read_field(os.path.join(d, "wallShearStress"), N)[1].get("airfoil") \
            if os.path.exists(os.path.join(d, "wallShearStress")) else None
        sel, kw = band
        nr = np.zeros(len(own))
        np.maximum.at(nr, kw[sel], nut[sel] / P.NU)
        gm = np.zeros(len(own))
        np.maximum.at(gm, kw[sel], gam[sel])
        self.surf = {}
        for side, s in sd.items():
            e = dict(xc=s["xc"], cp=self.p[own[s["idx"]]] / Q_INF,
                     yplus=yp[s["idx"]] if yp is not None else None)
            if wss is not None:
                cf = -np.sum(wss[s["idx"], :2] * s["t"], axis=1) / Q_INF
                if np.median(cf) < 0:
                    cf = -cf
            else:
                cf = np.zeros(len(s["idx"]))
            e["cf"] = cf
            e["tr"] = transition(s["xc"], cf, nr[s["idx"]], gm[s["idx"]], e["cp"])
            self.surf[side] = e


def load(case, geom):
    mesh = Mesh2D(case)
    sd, fc, own = surfaces(mesh, geom)
    dw, kw = cKDTree(fc).query(mesh.cc)
    band = (dw < P.TR_BAND, kw)
    times = time_dirs(case)
    if not times:
        raise SystemExit("no time directories in %s" % case)
    t_end = times[-1][0]
    names = [n for v, n in times if v > t_end - P.PERIOD + 1e-7 and
             os.path.exists(os.path.join(case, n, "U"))]
    print("reading %d snapshots (t = %s ... %s)" % (len(names), names[0], names[-1]))
    snaps = [Snapshot(mesh, case, n, sd, own, band) for n in names]
    return mesh, snaps


def nearest_phase(snaps, ph):
    return min(snaps, key=lambda s: abs((s.phase - ph + 180.0) % 360.0 - 180.0))


def tri_at(mesh, h):
    return mtri.Triangulation(mesh.xy[:, 0], mesh.xy[:, 1] + h, mesh.tri.triangles)


# =================================================================== figures
def fig_cp(out, snaps, N, coeff_at):
    fig, axs = plt.subplots(2, 2, figsize=(15, 11), sharex=True)
    for ax, ph in zip(axs.ravel(), PHASES):
        s = nearest_phase(snaps, ph)
        for side, col in (("upper", "#1f5fbf"), ("lower", "#d1495b")):
            e = s.surf[side]
            ax.plot(e["xc"], e["cp"], color=col, lw=1.5, label="%s surface" % side)
        ax.invert_yaxis()
        yl = ax.get_ylim()
        for side, col, f in (("upper", "#1f5fbf", 0.06), ("lower", "#d1495b", 0.14)):
            xt = s.surf[side]["tr"]["x_transition"]
            if xt is not None:
                ax.axvline(xt, color=col, ls=":", lw=1)
                ax.text(xt, yl[0] + f * (yl[1] - yl[0]), " transition %s x/c=%.3f"
                        % (side, xt), color=col, fontsize=8)
        cl, cd, cm = coeff_at(s.t)
        ax.set_title("phase %.0f deg (%s): y=%+.3f m, dy/dt=%+.2f m/s, alpha_eff=%+.1f deg\n"
                     "Cl=%+.3f  Cd=%+.4f  Cm(c/4)=%+.4f" % (s.phase, PHASE_NAME[ph], s.h, s.v,
                                                         s.alpha, cl, cd, cm), fontsize=9.5)
        ax.set_ylabel("Cp")
        ax.legend(loc="lower right", fontsize=8)
    for ax in axs[1]:
        ax.set_xlabel("x/c")
    fig.suptitle("Pressure coefficient on the upper and lower surface at four plunge phases "
                 "(last cycle; dotted: transition)")
    fig.tight_layout(rect=(0, 0.06, 1, 0.97))
    annotate(fig, N)
    fig.savefig(os.path.join(out, "01_pressure_distribution_Cp.png"), dpi=140)
    plt.close(fig)


def fig_contours(out, fname, mesh, snaps, geom, N, getter, title, cblabel, cmap, levels,
                 log=False):
    fig = plt.figure(figsize=(15, 15))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1, 1])
    s0 = nearest_phase(snaps, 0)
    panels = [(gs[0, 0], s0, (-5.8, 9.4), (-6.1, 6.1), "whole C-domain, phase %.0f deg" % s0.phase),
              (gs[0, 1], s0, (-0.5, 6.0), (-1.6, 1.6), "wake, phase %.0f deg" % s0.phase)]
    for g, ph in zip((gs[1, 0], gs[1, 1], gs[2, 0], gs[2, 1]), PHASES):
        s = nearest_phase(snaps, ph)
        panels.append((g, s, (-0.25, 1.15), (-0.55, 0.55),
                       "phase %.0f deg - %s (alpha_eff %+.1f deg)" % (s.phase, PHASE_NAME[ph],
                                                                     s.alpha)))
    for k, (g, s, xl, yl, nm) in enumerate(panels):
        ax = fig.add_subplot(g)
        tri = tri_at(mesh, s.h)
        v = np.clip(mesh.to_points(getter(s)), levels[0], levels[-1])
        if log:
            cs = ax.tricontourf(tri, v, levels=levels, cmap=cmap, norm=LogNorm())
        else:
            cs = ax.tricontourf(tri, v, levels=levels, cmap=cmap, extend="both")
        if k > 1:
            ax.tricontour(tri, v, levels=levels[::4], colors="k", linewidths=0.2,
                          norm=LogNorm() if log else None)
            for hh in (-P.AMP, P.AMP):
                ax.axhline(hh, color="w", ls=":", lw=0.6)
        draw_airfoil(ax, geom, s.h)
        ax.set_xlim(*xl)
        ax.set_ylim(*yl)
        ax.set_aspect("equal")
        ax.set_title(nm, fontsize=10)
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        fig.colorbar(cs, ax=ax, fraction=0.03, pad=0.01).set_label(cblabel)
    fig.suptitle(title + " (last cycle; dotted lines: extremes of the stroke)", fontsize=12)
    fig.tight_layout(rect=(0, 0.04, 1, 0.98))
    annotate(fig, N)
    fig.savefig(os.path.join(out, fname), dpi=130)
    plt.close(fig)


def fig_yplus(out, snaps, N):
    fig, axs = plt.subplots(1, 2, figsize=(15, 6))
    res = {}
    for ax, side in zip(axs, ("upper", "lower")):
        ymax = np.max([s.surf[side]["yplus"] for s in snaps], axis=0)
        x = snaps[0].surf[side]["xc"]
        for ph, col in zip(PHASES, ("#1f77b4", "#2ca02c", "#d62728", "#9467bd")):
            s = nearest_phase(snaps, ph)
            ax.plot(x, s.surf[side]["yplus"], color=col, lw=1.0,
                    label="phase %.0f deg" % s.phase)
        ax.plot(x, ymax, "k-", lw=1.8, label="maximum over the cycle")
        ax.axhline(1.0, color="k", ls="--", lw=0.8, label="y+ = 1")
        ax.axhline(P.YPLUS_TARGET, color="gray", ls=":", lw=0.8, label="target 0.7")
        ax.set_xlabel("x/c")
        ax.set_ylabel("y+ (wall-adjacent cell centre)")
        ax.set_title("%s surface: max y+ over the cycle = %.3f" % (side, ymax.max()))
        ax.legend(fontsize=8)
        res[side] = float(ymax.max())
    fig.suptitle("y+ distribution (first cell %.2e m)" % P.H1)
    fig.tight_layout(rect=(0, 0.07, 1, 0.96))
    annotate(fig, N)
    fig.savefig(os.path.join(out, "05_yplus_distribution.png"), dpi=140)
    plt.close(fig)
    return res


def fig_residuals(out, case, N, status):
    log = os.path.join(case, "log.pimpleFoam")
    if not os.path.exists(log):
        return
    res, co, step, seen = {}, [], 0.0, set()
    for line in open(log, errors="replace"):
        if line.startswith("Time = "):
            step = float(line.split()[2]) / P.PERIOD       # time in cycles
            seen = set()
            continue
        if line.startswith("Courant Number"):
            try:
                co.append((step, float(line.split()[-1])))
            except ValueError:
                pass
            continue
        m = re.search(r"Solving for (\w+), Initial residual = ([-+.\deE]+)", line)
        if m and m.group(1) not in seen:
            seen.add(m.group(1))
            res.setdefault(m.group(1), ([], []))
            res[m.group(1)][0].append(step)
            res[m.group(1)][1].append(float(m.group(2)))
    fig, axs = plt.subplots(2, 1, figsize=(14, 9), sharex=True,
                            gridspec_kw=dict(height_ratios=[2.2, 1]))
    for k, (x, y) in res.items():
        axs[0].semilogy(x, y, lw=0.6, label=k)
    axs[0].set_ylabel("initial residual (1st outer iteration)")
    axs[0].legend(ncol=7, fontsize=8)
    axs[0].set_title("Residuals per time step (pimpleFoam, %d outer correctors) - run status: %s"
                     % (P.N_OUTER, status.get("reason", "")), fontsize=10)
    if co:
        c = np.array(co)
        axs[1].semilogy(c[:, 0], c[:, 1], "k-", lw=0.4)
    axs[1].set_ylabel("max Courant number")
    axs[1].set_xlabel("time / T  (plunge cycles)")
    for ax in axs:
        for n in range(1, int(step) + 1):
            ax.axvline(n, color="gray", ls=":", lw=0.8)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    annotate(fig, N)
    fig.savefig(os.path.join(out, "06_residuals.png"), dpi=140)
    plt.close(fig)


def fig_coeffs(out, d, N):
    t = d["Time"]
    tt = t / P.PERIOD
    last = t > t[-1] - P.PERIOD + 1e-9
    ph = np.degrees(P.OMEGA * t) % 360.0
    al = np.array([alpha_eff(x) for x in t])
    fig, axs = plt.subplots(2, 2, figsize=(15, 10))
    ax = axs[0, 0]
    ax.plot(tt, d["Cl"], lw=0.8, label="Cl")
    ax.plot(tt, d["Cd"] * 10, lw=0.8, label="Cd x 10")
    ax.plot(tt, d["Cm"] * 10, lw=0.8, label="Cm(c/4) x 10")
    ax.plot(tt, [P.plunge(x)[0] / P.C for x in t], "k:", lw=0.8, label="y/c")
    ax.set_xlabel("time / T")
    ax.set_title("history of all cycles")
    ax.legend(fontsize=8, ncol=2)
    ax = axs[0, 1]
    o = np.argsort(ph[last])
    for key, f, lab in (("Cl", 1, "Cl"), ("Cd", 10, "Cd x 10"), ("Cm", 10, "Cm(c/4) x 10")):
        ax.plot(ph[last][o], d[key][last][o] * f, lw=1.2, label=lab)
    ax.plot(ph[last][o], al[last][o] / 10.0, "k:", lw=1, label="alpha_eff / 10 [deg]")
    ax.set_xlabel("phase [deg]  (y = A sin(phase))")
    ax.set_xticks(range(0, 361, 45))
    ax.set_title("last cycle vs phase")
    ax.legend(fontsize=8)
    ax = axs[1, 0]
    ncyc = int(np.floor(t[-1] / P.PERIOD + 1e-6))      # complete cycles only
    for n in range(ncyc):
        m = (t > n * P.PERIOD) & (t <= (n + 1) * P.PERIOD)
        ax.plot(al[m], d["Cl"][m], lw=0.9 if n == ncyc - 1 else 0.6,
                alpha=1 if n == ncyc - 1 else 0.45, label="cycle %d" % (n + 1))
    a = np.linspace(-P.ALPHA_EFF_MAX, P.ALPHA_EFF_MAX, 3)
    ax.plot(a, 2 * np.pi * np.radians(a), "k--", lw=0.8, label="2 pi alpha (quasi-steady)")
    ax.set_xlabel("effective angle of attack [deg]")
    ax.set_ylabel("Cl")
    ax.set_title("lift loops (hysteresis)")
    ax.legend(fontsize=8)
    ax = axs[1, 1]
    for n in range(ncyc):
        m = (t > n * P.PERIOD) & (t <= (n + 1) * P.PERIOD)
        ax.plot(al[m], d["Cd"][m], lw=0.9 if n == ncyc - 1 else 0.6,
                alpha=1 if n == ncyc - 1 else 0.45, label="cycle %d" % (n + 1))
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel("effective angle of attack [deg]")
    ax.set_ylabel("Cd  (negative = thrust)")
    ax.set_title("drag / thrust loops")
    ax.legend(fontsize=8)
    fig.suptitle("Force and moment coefficients (Cm about the moving quarter chord)")
    fig.tight_layout(rect=(0, 0.06, 1, 0.97))
    annotate(fig, N)
    fig.savefig(os.path.join(out, "08_coefficients_history.png"), dpi=140)
    plt.close(fig)


# ---------------------------------------------------------- pathlines
class FlowInTime:
    """Velocity at absolute (x, y, t), periodic over the stored cycle."""

    def __init__(self, mesh, snaps):
        self.mesh = mesh
        o = np.argsort([s.phase for s in snaps])
        self.snaps = [snaps[i] for i in o]
        self.ph = np.array([s.phase for s in self.snaps])
        self.ip = []
        for s in self.snaps:
            self.ip.append((mtri.LinearTriInterpolator(mesh.tri, mesh.to_points(s.U[:, 0])),
                            mtri.LinearTriInterpolator(mesh.tri, mesh.to_points(s.U[:, 1]))))

    def __call__(self, p, t):
        ph = phase_deg(t)
        h = P.AMP * math.sin(P.OMEGA * t)     # periodic state (after the start-up ramp)
        j = np.searchsorted(self.ph, ph) % len(self.ph)
        i = j - 1
        p0, p1 = self.ph[i], self.ph[j]
        span = (p1 - p0) % 360.0 or 360.0
        w = ((ph - p0) % 360.0) / span
        x, y = p[:, 0], p[:, 1] - h
        out = np.zeros_like(p)
        for k, wk in ((i, 1 - w), (j, w)):
            ux, uy = self.ip[k]
            out[:, 0] += wk * np.ma.filled(ux(x, y), np.nan)
            out[:, 1] += wk * np.ma.filled(uy(x, y), np.nan)
        return out


def trace(flow, seeds, t0, dt=2.0e-4, t_max=0.25):
    p = seeds.copy()
    path = [p.copy()]
    t = t0
    for _ in range(int(t_max / dt)):
        v1 = flow(p, t)
        v2 = flow(p + 0.5 * dt * v1, t + 0.5 * dt)
        p = p + dt * v2
        p[~np.isfinite(p).all(axis=1) | (p[:, 0] > X_END + 0.03)] = np.nan
        path.append(p.copy())
        t += dt
        if np.isnan(p).all():
            break
    return np.array(path), dt


def fig_pathlines(out, mesh, snaps, geom, N):
    flow = FlowInTime(mesh, snaps)
    seeds_y = np.linspace(-0.50, 0.50, 25)
    seeds = np.c_[np.full(len(seeds_y), X_SEED), seeds_y]
    fig, axs = plt.subplots(2, 2, figsize=(16, 10.5))
    t_last = max(s.t for s in snaps)
    t_base = t_last - P.PERIOD
    for ax, ph in zip(axs.ravel(), PHASES):
        t0 = t_base + math.radians(ph) / P.OMEGA
        path, dt = trace(flow, seeds, t0)
        n = len(path)
        # background: |U| when the release reaches the airfoil (~0.25 m transit)
        t_mid = t0 + (0.25 - X_SEED) / P.U_INF
        s_mid = min(snaps, key=lambda s: abs((s.phase - phase_deg(t_mid) + 180) % 360 - 180))
        h_mid = P.AMP * math.sin(P.OMEGA * t_mid)
        cs = ax.tricontourf(tri_at(mesh, h_mid), np.clip(mesh.to_points(
            np.linalg.norm(s_mid.U, axis=1)), 0, 1.8 * P.U_INF),
            levels=np.linspace(0, 1.8 * P.U_INF, 31), cmap="Blues", alpha=0.8)
        for k in range(len(seeds_y)):
            L = path[:, k]
            L = L[np.isfinite(L).all(axis=1)]
            ax.plot(L[:, 0], L[:, 1], color="#c0392b", lw=0.9)
            ax.plot(L[0, 0], L[0, 1], "o", color="#c0392b", ms=2.5)
        for f in np.linspace(0, 1, 5):         # airfoil positions during the transit
            tt = t0 + f * (n - 1) * dt
            draw_airfoil(ax, geom, P.AMP * math.sin(P.OMEGA * tt), alpha=0.25 + 0.75 * (f == 0.5),
                         lw=0.5)
        ax.axvline(X_END, color="k", ls=":", lw=0.8)
        ax.set_xlim(X_SEED - 0.03, X_END + 0.03)
        ax.set_ylim(-0.6, 0.6)
        ax.set_aspect("equal")
        ax.set_title("released at phase %.0f deg (%s)" % (ph, PHASE_NAME[ph]), fontsize=10)
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
    fig.colorbar(cs, ax=axs, fraction=0.015, pad=0.01).set_label(
        "|U| [m/s] when the particles reach the airfoil")
    fig.suptitle("Pathlines of %d particles released at x = %.2f m and tracked through the "
                 "time-dependent flow to %.0f cm behind the trailing edge\n(airfoil drawn at five "
                 "instants of the transit; solid = middle instant)"
                 % (len(seeds_y), X_SEED, (X_END - P.X_TE) * 100), fontsize=11)
    annotate(fig, N)
    fig.savefig(os.path.join(out, "07_pathlines.png"), dpi=130, bbox_inches="tight")
    plt.close(fig)


def fig_mesh(out, mesh, geom, N):
    segs = [mesh.xy[np.r_[p, p[0]]] for p in mesh.polys]
    fig, axs = plt.subplots(2, 2, figsize=(15, 11))
    views = [((-5.8, 9.4), (-6.1, 6.1), "whole C-domain (wake 30c)"),
             ((-0.2, 3.2), (-0.8, 0.8), "near field and refined wake band"),
             ((-0.012, 0.02), (-0.012, 0.012), "leading edge"),
             ((0.27, 0.33), (-0.03, 0.03), "trailing edge")]
    for ax, (xl, yl, nm) in zip(axs.ravel(), views):
        ax.add_collection(LineCollection(segs, colors="k", linewidths=0.1 if "whole" in nm else 0.15))
        draw_airfoil(ax, geom, lw=0.5)
        ax.set_xlim(*xl)
        ax.set_ylim(*yl)
        ax.set_aspect("equal")
        ax.set_title(nm)
    fig.suptitle("Hybrid C-mesh: %d cells (%d quads), 4 cm structured quad layer, first cell %.2e m "
                 "(the whole mesh moves rigidly with the airfoil)" % (N, int(geom["nquad"]), P.H1))
    fig.tight_layout(rect=(0, 0.05, 1, 0.97))
    annotate(fig, N)
    fig.savefig(os.path.join(out, "09_mesh.png"), dpi=140)
    plt.close(fig)


def fig_mesh3d(out, mesh, geom, N, span=0.45, nspan=10):
    from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection
    xl, yl = (-0.08, 0.40), (-0.36, 0.36)
    segs = []
    for p in mesh.polys:
        q = mesh.xy[np.r_[p, p[0]]]
        if q[:, 0].max() < xl[0] or q[:, 0].min() > xl[1] or abs(q[:, 1]).min() > 0.08:
            continue
        segs.append(np.c_[q, np.zeros(len(q))])
    zs = np.linspace(0.0, span, nspan + 1)
    light = np.array([0.3, 0.8, -0.5])
    light /= np.linalg.norm(light)
    fig = plt.figure(figsize=(15, 10))
    ax = fig.add_subplot(projection="3d")
    ax.add_collection3d(Line3DCollection(segs, colors="#2b2b2b", linewidths=0.12, alpha=0.8))
    for h, alpha in ((-P.AMP, 0.18), (P.AMP, 0.18), (0.0, 1.0)):
        for surf in (geom["up"], geom["lo"]):
            polys, cols = [], []
            for i in range(len(surf) - 1):
                p0, p1 = surf[i], surf[i + 1]
                nvec = np.cross(np.r_[p1 - p0, 0.0], [0, 0, 1.0])
                nvec /= max(np.linalg.norm(nvec), 1e-12)
                shade = 0.45 + 0.55 * abs(nvec @ light)
                for k in range(nspan):
                    polys.append([(p0[0], p0[1] + h, zs[k]), (p1[0], p1[1] + h, zs[k]),
                                  (p1[0], p1[1] + h, zs[k + 1]), (p0[0], p0[1] + h, zs[k + 1])])
                    cols.append(np.clip(np.array(matplotlib.colors.to_rgb("#8fb3e0")) * shade, 0, 1))
            ax.add_collection3d(Poly3DCollection(polys, facecolors=cols, alpha=alpha,
                                                 edgecolors="#1a1a1a" if alpha == 1 else "none",
                                                 linewidths=0.08))
    for h in (-P.AMP, P.AMP):
        ax.plot([P.X_LE + 0.25 * P.C] * 2, [0, h], [span, span], "k--", lw=0.8)
    ax.set_xlim(*xl)
    ax.set_ylim(*yl)
    ax.set_zlim(0, span)
    ax.set_box_aspect((xl[1] - xl[0], yl[1] - yl[0], span), zoom=1.05)
    ax.view_init(elev=24, azim=-62)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_zlabel("span z [m]")
    ax.grid(False)
    fig.suptitle("3-D view of the wing mesh (%d+%d surface elements) at mid stroke, ghosts at the top "
                 "and bottom of the +-%.2f m plunge;\n2-D computational mesh on the z = 0 plane "
                 "(the CFD mesh is one cell deep; span drawn %.0f cm for illustration)"
                 % (P.N_SURF, P.N_SURF, P.AMP, span * 100), fontsize=11)
    fig.subplots_adjust(left=0, right=1, bottom=0.10, top=0.93)
    annotate(fig, N)
    fig.savefig(os.path.join(out, "10_mesh_3d.png"), dpi=140, bbox_inches="tight")
    plt.close(fig)


# =================================================================== main
def main():
    case, out, geomfile = sys.argv[1], sys.argv[2], sys.argv[3]
    os.makedirs(out, exist_ok=True)
    geom = dict(np.load(geomfile))
    status = {}
    sp = os.path.join(case, "run_status.json")
    if os.path.exists(sp):
        status = json.load(open(sp))
    d = read_all_coeffs(case)
    cm_key = "CmPitch" if "CmPitch" in d else "Cm"
    # Cm about the MOVING quarter chord (forceCoeffs uses a fixed CofR at y = 0)
    hh = np.array([P.plunge(t)[0] for t in d["Time"]])
    d["Cm"] = d[cm_key] - hh / P.C * d["Cd"]
    vv = np.array([P.plunge(t)[1] for t in d["Time"]])

    mesh, snaps = load(case, geom)
    N = mesh.ncells

    def coeff_at(t):
        k = int(np.argmin(np.abs(d["Time"] - t)))
        return d["Cl"][k], d["Cd"][k], d["Cm"][k]

    # ---------------------------------------------------------- numbers
    last = d["Time"] > d["Time"][-1] - P.PERIOD + 1e-9
    ct = -float(np.mean(d["Cd"][last]))
    cp_in = float(np.mean(-d["Cl"][last] * vv[last] / P.U_INF))
    cyc = cycle_stats(d)
    res = dict(case=dict(
        airfoil="NACA 0009, c = %.2f m" % P.C, Re=P.RE, U_inf=P.U_INF, AoA_deg=P.AOA_DEG,
        plunge="y = %.3f m sin(2 pi %.2f Hz t)" % (P.AMP, P.FREQ), amplitude_m=P.AMP,
        frequency_Hz=P.FREQ, max_plunge_velocity=P.V_PLUNGE, max_alpha_eff_deg=P.ALPHA_EFF_MAX,
        reduced_frequency=P.K_RED, strouhal=P.STROUHAL, Tu=P.TU,
        turbulence_model="kOmegaSSTLM (k-omega SST + gamma-ReTheta transition)",
        cells=int(N), quads_2d=int(geom["nquad"]), first_layer=P.H1, dt=P.DT,
        steps_per_cycle=P.STEPS_PER_CYCLE, cycles_computed=round(d["Time"][-1] / P.PERIOD, 3),
        run=status))
    res["last_cycle"] = dict(
        Cl_max=float(d["Cl"][last].max()), Cl_min=float(d["Cl"][last].min()),
        Cl_mean=float(d["Cl"][last].mean()),
        Cl_rms=float(np.sqrt(np.mean((d["Cl"][last] - d["Cl"][last].mean()) ** 2))),
        Cd_mean=float(d["Cd"][last].mean()), Cd_max=float(d["Cd"][last].max()),
        Cd_min=float(d["Cd"][last].min()), Cm_mean=float(d["Cm"][last].mean()),
        Cm_max=float(d["Cm"][last].max()), Cm_min=float(d["Cm"][last].min()),
        thrust_coefficient_Ct=ct, power_coefficient=cp_in,
        propulsive_efficiency=(ct / cp_in) if cp_in > 1e-9 else None)
    res["cycles"] = cyc
    two = d["Time"] > d["Time"][-1] - 2 * P.PERIOD + 1e-9
    if d["Time"][-1] > 2.5 * P.PERIOD:
        cp2 = float(np.mean(-d["Cl"][two] * vv[two] / P.U_INF))
        cd2 = float(np.mean(d["Cd"][two]))
        res["last_two_cycles"] = dict(
            Cl_rms=float(np.sqrt(np.mean((d["Cl"][two] - d["Cl"][two].mean()) ** 2))),
            Cd_mean=cd2, power_coefficient=cp2,
            propulsive_efficiency=(-cd2 / cp2) if cp2 > 1e-9 else None)
    res["phases"] = {}
    for ph in PHASES:
        s = nearest_phase(snaps, ph)
        cl, cd, cm = coeff_at(s.t)
        res["phases"]["%d" % ph] = dict(
            t=s.t, phase=s.phase, y=s.h, dydt=s.v, alpha_eff=s.alpha, Cl=float(cl), Cd=float(cd),
            Cm=float(cm), transition={k: v["tr"] for k, v in s.surf.items()})
    trs = {side: [s.surf[side]["tr"]["x_transition"] for s in snaps] for side in ("upper", "lower")}
    res["transition_vs_phase"] = {side: [[s.phase, x] for s, x in zip(snaps, trs[side])]
                                  for side in trs}

    # ---------------------------------------------------------- figures
    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.3})
    fig_cp(out, snaps, N, coeff_at)
    fig_contours(out, "02_pressure_contours.png", mesh, snaps, geom, N, lambda s: s.p / Q_INF,
                 "Pressure coefficient Cp", "Cp", "RdBu_r", np.linspace(-3.0, 1.0, 41))
    fig_contours(out, "03_velocity_contours.png", mesh, snaps, geom, N,
                 lambda s: np.linalg.norm(s.U, axis=1), "Velocity magnitude |U| [m/s]",
                 "|U| [m/s]", "viridis", np.linspace(0, 2.0 * P.U_INF, 41))
    fig_contours(out, "04_turbulence_intensity_log.png", mesh, snaps, geom, N,
                 lambda s: np.sqrt(2.0 / 3.0 * np.maximum(s.k, 1e-14)) / P.U_INF * 100.0,
                 "Turbulence intensity Tu = sqrt(2k/3)/U_inf (logarithmic scale)",
                 "Tu [%] (log scale)", "magma", np.logspace(-2, 2, 33), log=True)
    res["yplus_max"] = fig_yplus(out, snaps, N)
    fig_residuals(out, case, N, status)
    fig_pathlines(out, mesh, snaps, geom, N)
    fig_coeffs(out, d, N)
    fig_mesh(out, mesh, geom, N)
    fig_mesh3d(out, mesh, geom, N)
    fig_transition(out, snaps, N)

    with open(os.path.join(out, "results.json"), "w") as f:
        json.dump(res, f, indent=2, default=lambda o: None)
    write_report(out, res)
    print(open(os.path.join(out, "results.md")).read())


def fig_transition(out, snaps, N):
    fig, ax = plt.subplots(figsize=(12, 5.5))
    for side, col, mk in (("upper", "#1f5fbf", "o"), ("lower", "#d1495b", "s")):
        ph = [s.phase for s in snaps]
        xt = [s.surf[side]["tr"]["x_transition"] for s in snaps]
        xs = [s.surf[side]["tr"]["x_laminar_separation"] for s in snaps]
        ax.plot(ph, [np.nan if v is None else v for v in xt], mk, color=col, ms=5,
                label="%s: transition onset" % side)
        ax.plot(ph, [np.nan if v is None else v for v in xs], mk, color=col, ms=4, mfc="none",
                label="%s: laminar separation" % side)
    ax2 = ax.twinx()
    phs = np.linspace(0, 360, 361)
    ax2.plot(phs, -np.degrees(np.arctan(P.V_PLUNGE * np.cos(np.radians(phs)) / P.U_INF)), "k:",
             lw=1)
    ax2.set_ylabel("alpha_eff [deg] (dotted)")
    ax2.grid(False)
    ax.set_xlim(0, 360)
    ax.set_xticks(range(0, 361, 45))
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("phase [deg]  (y = A sin(phase))")
    ax.set_ylabel("x/c")
    ax.set_title("Transition point and laminar separation over the last cycle "
                 "(missing marker: none found on that surface)")
    ax.legend(fontsize=8, ncol=2, loc="lower left")
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    annotate(fig, N)
    fig.savefig(os.path.join(out, "11_transition_vs_phase.png"), dpi=140)
    plt.close(fig)


def write_report(out, r):
    c, lc = r["case"], r["last_cycle"]
    f = lambda v, fmt="%.3f": "-" if v is None else fmt % v
    L = ["# Plunging NACA 0009 - results", "",
         "| item | value |", "|---|---|",
         "| airfoil | %s, AoA = %.0f deg |" % (c["airfoil"], c["AoA_deg"]),
         "| flow | air ISA sea level, Re = %.0f, U = %.3f m/s, Tu = %.1f %% |"
         % (c["Re"], c["U_inf"], c["Tu"] * 100),
         "| motion | %s; max plunge speed %.2f m/s, max effective AoA %.1f deg, k = %.3f, St = %.3f |"
         % (c["plunge"], c["max_plunge_velocity"], c["max_alpha_eff_deg"], c["reduced_frequency"],
            c["strouhal"]),
         "| model | %s, pimpleFoam (rigid mesh motion) |" % c["turbulence_model"],
         "| mesh | %d cells (%d quads), first cell %.2e m |" % (c["cells"], c["quads_2d"],
                                                               c["first_layer"]),
         "| time step | %.3e s (%d per cycle), %.2f cycles computed |"
         % (c["dt"], c["steps_per_cycle"], c["cycles_computed"]),
         "| run | %s - %s |" % (c["run"].get("status", ""), c["run"].get("reason", "")), "",
         "## Coefficients over the last cycle", "",
         "| quantity | min | max | mean |", "|---|---|---|---|",
         "| Cl | %.3f | %.3f | %.3f (rms %.3f) |" % (lc["Cl_min"], lc["Cl_max"], lc["Cl_mean"],
                                                   lc["Cl_rms"]),
         "| Cd | %.4f | %.4f | %.4f |" % (lc["Cd_min"], lc["Cd_max"], lc["Cd_mean"]),
         "| Cm (c/4, +nose-up) | %.4f | %.4f | %.4f |" % (lc["Cm_min"], lc["Cm_max"], lc["Cm_mean"]),
         "",
         "Thrust coefficient Ct = -mean(Cd) = %.4f; input power coefficient "
         "mean(-Cl dy/dt / U) = %.4f; propulsive efficiency Ct/Cp = %s."
         % (lc["thrust_coefficient_Ct"], lc["power_coefficient"],
            f(lc["propulsive_efficiency"], "%.3f")), "",
         "## Cycle-to-cycle convergence", "",
         "| cycle | Cl min | Cl max | Cl rms | mean Cd | largest Cl change vs previous cycle "
         "| vs two cycles before |", "|---|---|---|---|---|---|---|"]
    cy = r["cycles"]
    for k, s in enumerate(cy):
        c1 = "-" if k < 1 else "%.1f %%" % (100 * cycle_change(cy[k - 1], s)[0])
        c2 = "-" if k < 2 else "%.1f %%" % (100 * cycle_change(cy[k - 2], s)[0])
        L.append("| %d | %.3f | %.3f | %.3f | %.4f | %s | %s |"
                 % (s["cycle"], s["Cl_min"], s["Cl_max"], s["Cl_rms"], s["Cd_mean"], c1, c2))
    L += ["", "Convergence criterion: Cl max, Cl min and Cl rms of a cycle each within %.0f %% of "
          "their values one cycle earlier (period 1) or two cycles earlier (period 2, a pattern "
          "dynamic stall can lock into), checked from cycle %d on."
          % (100 * P.CYCLE_TOL, P.MIN_CYCLES)]
    if r.get("last_two_cycles"):
        t2 = r["last_two_cycles"]
        L += ["", "Averages over the last two cycles: Cl rms %.3f, mean Cd %.4f (Ct = %.4f), "
              "power coefficient %.4f, propulsive efficiency %s."
              % (t2["Cl_rms"], t2["Cd_mean"], -t2["Cd_mean"], t2["power_coefficient"],
                 f(t2["propulsive_efficiency"], "%.3f"))]
    L += ["", "## Four phases of the last cycle", "",
          "| phase | y [m] | dy/dt [m/s] | alpha_eff | Cl | Cd | Cm | transition upper | "
          "transition lower |", "|---|---|---|---|---|---|---|---|---|"]
    for ph, v in r["phases"].items():
        tu, tl = v["transition"]["upper"], v["transition"]["lower"]
        L.append("| %.0f deg | %+.3f | %+.2f | %+.1f deg | %+.3f | %+.4f | %+.4f | %s (%s) | %s (%s) |"
                 % (v["phase"], v["y"], v["dydt"], v["alpha_eff"], v["Cl"], v["Cd"], v["Cm"],
                    f(tu["x_transition"]), tu["method"], f(tl["x_transition"]), tl["method"]))
    L += ["", "Transition x/c from the leading edge: inner-BL intermittency > 0.5 and "
          "nut/nu > %g within %.1f mm of the wall, or the Cf minimum inside a closed laminar "
          "separation bubble. The full phase dependence is in 11_transition_vs_phase.png."
          % (P.NUT_RATIO_TR, P.TR_BAND * 1e3), "",
          "## y+ (maximum over the last cycle)", "",
          "| surface | max y+ |", "|---|---|"]
    for k, v in r.get("yplus_max", {}).items():
        L.append("| %s | %.3f |" % (k, v))
    with open(os.path.join(out, "results.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
