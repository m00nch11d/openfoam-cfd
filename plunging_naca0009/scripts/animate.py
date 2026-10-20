"""
Animation of the last plunge cycle: spanwise vorticity, the moving airfoil
and particles released continuously from 25 points at x = -0.20 m and
tracked through the time-dependent flow (periodic over the stored cycle).

Writes <resultsDir>/12_plunging_animation.gif (always) and .mp4 (when ffmpeg
is available).

usage: python3 animate.py <unsteadyCase> <resultsDir> <geom.npz> [--frames 80] [--fps 16]
"""
import argparse
import math
import os
import shutil

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402
from matplotlib import animation  # noqa: E402

import params as P  # noqa: E402
from post import X_END, X_SEED, FlowInTime, airfoil_ring, annotate, load, phase_deg  # noqa: E402


def vorticity(mesh, U):
    """Spanwise vorticity per triangle of the plotting triangulation."""
    tri = mesh.tri
    ux = mtri.LinearTriInterpolator(tri, mesh.to_points(U[:, 0]))
    uy = mtri.LinearTriInterpolator(tri, mesh.to_points(U[:, 1]))
    xc = mesh.xy[tri.triangles].mean(axis=1)
    _, duxdy = ux.gradient(xc[:, 0], xc[:, 1])
    duydx, _ = uy.gradient(xc[:, 0], xc[:, 1])
    return np.ma.filled(duydx - duxdy, 0.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("results")
    ap.add_argument("geom")
    ap.add_argument("--frames", type=int, default=80)
    ap.add_argument("--fps", type=int, default=16)
    a = ap.parse_args()
    os.makedirs(a.results, exist_ok=True)
    geom = dict(np.load(a.geom))
    mesh, snaps = load(a.case, geom)
    flow = FlowInTime(mesh, snaps)
    wz = [vorticity(mesh, s.U) for s in flow.snaps]
    ph_s = flow.ph

    def wz_at(t):
        ph = phase_deg(t)
        j = np.searchsorted(ph_s, ph) % len(ph_s)
        i = j - 1
        span = (ph_s[j] - ph_s[i]) % 360.0 or 360.0
        w = ((ph - ph_s[i]) % 360.0) / span
        return (1 - w) * wz[i] + w * wz[j]

    # particles: continuous release, advanced through one extra transit first
    seeds_y = np.linspace(-0.5, 0.5, 25)
    dt = 2.0e-4
    t_rel = 0.006
    t_last = max(s.t for s in snaps)
    t_frames = t_last - P.PERIOD + np.arange(a.frames) * P.PERIOD / a.frames
    t = t_frames[0] - 0.12
    pts = np.empty((0, 2))
    next_rel = t
    frames_pts = []
    fi = 0
    while fi < len(t_frames):
        if t >= next_rel - 1e-12:
            pts = np.vstack([pts, np.c_[np.full(len(seeds_y), X_SEED), seeds_y]])
            next_rel += t_rel
        if len(pts):
            v1 = flow(pts, t)
            v2 = flow(pts + 0.5 * dt * v1, t + 0.5 * dt)
            pts = pts + dt * v2
            pts = pts[np.isfinite(pts).all(axis=1) & (pts[:, 0] < X_END + 0.1)]
        t += dt
        while fi < len(t_frames) and t >= t_frames[fi]:
            frames_pts.append(pts.copy())
            fi += 1

    lim = 400.0                                  # 1/s, colour range of the vorticity
    fig, ax = plt.subplots(figsize=(14, 7.6))
    tri0 = mesh.tri
    pc = ax.tripcolor(mtri.Triangulation(tri0.x, tri0.y, tri0.triangles),
                      facecolors=wz_at(t_frames[0]), cmap="RdBu_r", vmin=-lim, vmax=lim,
                      shading="flat", zorder=1)
    fig.colorbar(pc, ax=ax, fraction=0.025, pad=0.01).set_label("spanwise vorticity [1/s]")
    ring = airfoil_ring(geom)
    foil = ax.fill(ring[:, 0], ring[:, 1], color="#444444", zorder=5)[0]
    dots = ax.scatter([], [], s=6, color="#111111", zorder=6)
    for hh in (-P.AMP, P.AMP):
        ax.axhline(hh, color="gray", ls=":", lw=0.7)
    ax.axvline(X_END, color="k", ls=":", lw=0.7)
    ax.set_xlim(-0.3, 1.2)
    ax.set_ylim(-0.6, 0.6)
    ax.set_aspect("equal")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    title = ax.set_title("Plunging NACA 0009", fontsize=11)
    fig.tight_layout(rect=(0, 0.07, 1, 0.97))
    annotate(fig, mesh.ncells)

    state = {"pc": pc}

    def frame(i):
        tt = t_frames[i]
        h, v = P.AMP * math.sin(P.OMEGA * tt), P.V_PLUNGE * math.cos(P.OMEGA * tt)
        state["pc"].remove()
        state["pc"] = ax.tripcolor(mtri.Triangulation(tri0.x, tri0.y + h, tri0.triangles),
                                   facecolors=wz_at(tt), cmap="RdBu_r", vmin=-lim, vmax=lim,
                                   shading="flat", zorder=1)
        foil.set_xy(ring + [0.0, h])
        dots.set_offsets(frames_pts[i] if len(frames_pts[i]) else np.empty((0, 2)))
        title.set_text("Plunging NACA 0009 - phase %3.0f deg, y = %+.3f m, dy/dt = %+.2f m/s, "
                       "alpha_eff = %+.1f deg" % (phase_deg(tt), h, v,
                                                   -math.degrees(math.atan2(v, P.U_INF))))
        return []

    ani = animation.FuncAnimation(fig, frame, frames=a.frames, interval=1000 / a.fps, blit=False)
    gif = os.path.join(a.results, "12_plunging_animation.gif")
    ani.save(gif, writer=animation.PillowWriter(fps=a.fps), dpi=70)
    print("wrote", gif)
    if shutil.which("ffmpeg"):
        mp4 = os.path.join(a.results, "12_plunging_animation.mp4")
        ani.save(mp4, writer=animation.FFMpegWriter(fps=a.fps, bitrate=3000), dpi=110)
        print("wrote", mp4)
    else:
        print("ffmpeg not found - MP4 skipped (sudo apt install ffmpeg)")
    plt.close(fig)


if __name__ == "__main__":
    main()
