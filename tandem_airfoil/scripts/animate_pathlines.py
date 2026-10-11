"""
Animated pathlines for the tandem-airfoil case (post-processing only, the
flow field is read from the latest time directory of an existing run).

Massless particles are released continuously from 26 seed points at
x = -0.20 m and advected with the steady velocity field (2nd-order
Runge-Kutta in physical time) until they are 60 cm behind the S1223
trailing edge. Every seed emits a new particle at a fixed interval, so the
animation shows streaks moving past both airfoils; the full pathlines are
drawn faintly underneath.

Writes <resultsDir>/10_pathlines_animation.gif (always) and .mp4 (when
ffmpeg is available).

usage: python3 animate_pathlines.py <caseDir> <resultsDir> <geom.npz>
           [--frames 160] [--fps 20]
"""
import argparse
import os
import shutil

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402
from matplotlib import animation  # noqa: E402

import params as P  # noqa: E402
from post import Mesh2D, annotate, draw_airfoils, latest_time, read_field  # noqa: E402

X_SEED = -0.20                    # release line [m]
X_END = P.X_TE1 + 0.60            # 60 cm behind the S1223 trailing edge
SEEDS_Y = np.r_[np.linspace(-0.15, -0.02, 7), np.linspace(-0.012, 0.03, 12),
                np.linspace(0.04, 0.16, 7)]


def trajectories(mesh, U, dt, t_max):
    """Integrate all seeds together in physical time.

    Returns positions (n_steps, n_seeds, 2); NaN once a particle has left the
    window, hit a wall (masked interpolation) or stagnated.
    """
    ux = mtri.LinearTriInterpolator(mesh.tri, mesh.to_points(U[:, 0]))
    uy = mtri.LinearTriInterpolator(mesh.tri, mesh.to_points(U[:, 1]))

    def vel(p):
        a = ux(p[:, 0], p[:, 1])
        b = uy(p[:, 0], p[:, 1])
        v = np.c_[np.ma.filled(a, np.nan), np.ma.filled(b, np.nan)]
        return v

    n = int(t_max / dt) + 1
    pos = np.full((n, len(SEEDS_Y), 2), np.nan)
    p = np.c_[np.full(len(SEEDS_Y), X_SEED), SEEDS_Y]
    pos[0] = p
    for k in range(1, n):
        v1 = vel(p)
        v2 = vel(p + 0.5 * dt * v1)
        p = p + dt * v2
        stop = ~np.isfinite(p).all(axis=1) | (p[:, 0] > X_END + 0.03) | \
            (np.linalg.norm(v2, axis=1) < 1e-3 * P.U_INF)
        p[stop] = np.nan
        pos[k] = p
        if np.isnan(p).all():
            return pos[: k + 1]
    return pos


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("results")
    ap.add_argument("geom")
    ap.add_argument("--frames", type=int, default=160)
    ap.add_argument("--fps", type=int, default=20)
    a = ap.parse_args()
    os.makedirs(a.results, exist_ok=True)

    geom = dict(np.load(a.geom))
    mesh = Mesh2D(a.case)
    t = latest_time(a.case)
    if t is None:
        raise SystemExit("no solution time directory in %s - run the case first" % a.case)
    U, _ = read_field(os.path.join(a.case, t, "U"), mesh.ncells)
    print("animating pathlines from %s/%s (%d cells)" % (a.case, t, mesh.ncells))

    # time scale: free-stream transit from the seeds to the end line
    t_transit = (X_END - X_SEED) / P.U_INF
    dt = t_transit / 1500.0
    pos = trajectories(mesh, U, dt, 3.0 * t_transit)
    n_steps = len(pos)

    # particles are released from every seed each t_release; the animation
    # covers one free-stream transit after the first wave has arrived
    t_release = t_transit / 12.0
    t_anim = np.linspace(t_transit, 2.0 * t_transit, a.frames)

    umag = mesh.to_points(np.linalg.norm(U, axis=1))
    fig, ax = plt.subplots(figsize=(15, 6.8))
    cs = ax.tricontourf(mesh.tri, np.clip(umag, 0, 1.8 * P.U_INF),
                        levels=np.linspace(0, 1.8 * P.U_INF, 37), cmap="Blues", alpha=0.85)
    fig.colorbar(cs, ax=ax, fraction=0.02, pad=0.01).set_label("|U| [m/s]")
    for s in range(len(SEEDS_Y)):
        ax.plot(pos[:, s, 0], pos[:, s, 1], color="#c0392b", lw=0.5, alpha=0.35)
    draw_airfoils(ax, geom)
    ax.axvline(X_END, color="k", ls=":", lw=0.8)
    ax.text(X_END - 0.005, 0.18, "60 cm behind S1223 TE ", fontsize=8, ha="right")
    ax.plot(np.full(len(SEEDS_Y), X_SEED), SEEDS_Y, "o", color="#c0392b", ms=3)
    ax.set_xlim(X_SEED - 0.03, X_END + 0.03)
    ax.set_ylim(-0.19, 0.2)
    ax.set_aspect("equal")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title("Pathlines of %d particle streams released at x = %.2f m (steady flow, solution "
                 "iteration %s)" % (len(SEEDS_Y), X_SEED, t))
    dots = ax.scatter([], [], s=9, color="#8b0000", zorder=6)
    clock = ax.text(0.01, 0.97, "", transform=ax.transAxes, fontsize=10, va="top",
                    bbox=dict(fc="white", ec="#999999", alpha=0.85))
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    annotate(fig, mesh.ncells)

    def frame(i):
        tt = t_anim[i]
        ages = tt - np.arange(0.0, tt + 1e-12, t_release)        # age of every release
        k = np.round(ages / dt).astype(int)
        k = k[(k >= 0) & (k < n_steps)]
        pts = pos[k].reshape(-1, 2)
        pts = pts[np.isfinite(pts).all(axis=1)]
        dots.set_offsets(pts if len(pts) else np.empty((0, 2)))
        clock.set_text("t = %6.1f ms" % ((tt - t_anim[0]) * 1e3))
        return dots, clock

    ani = animation.FuncAnimation(fig, frame, frames=a.frames, interval=1000 / a.fps, blit=True)
    gif = os.path.join(a.results, "10_pathlines_animation.gif")
    ani.save(gif, writer=animation.PillowWriter(fps=a.fps), dpi=80)
    print("wrote", gif)
    if shutil.which("ffmpeg"):
        mp4 = os.path.join(a.results, "10_pathlines_animation.mp4")
        ani.save(mp4, writer=animation.FFMpegWriter(fps=a.fps, bitrate=2400), dpi=120)
        print("wrote", mp4)
    else:
        print("ffmpeg not found - MP4 skipped (sudo apt install ffmpeg)")
    plt.close(fig)


if __name__ == "__main__":
    main()
