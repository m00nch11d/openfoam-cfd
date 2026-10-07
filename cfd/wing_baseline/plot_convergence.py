#!/usr/bin/env python3
"""
Plot residuals and force coefficients vs iteration for the wing case.

Reads (written by the solverInfo and forceCoeffs function objects):
    postProcessing/solverInfo/<t>/solverInfo.dat
    postProcessing/forceCoeffs/<t>/coefficient.dat
Restarts (several <t> folders) are stitched together, later runs win.

Usage:
    python3 plot_convergence.py            # save convergence.png and show it
    python3 plot_convergence.py --no-show  # save only (used by Allrun)
    python3 plot_convergence.py --avg 200  # average the last 200 iterations
"""
import argparse
import glob
import os
import sys

import matplotlib

# categorical palette, fixed order
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
INK, MUTED, GRID = "#1f1f1e", "#6b6a63", "#e4e3dc"


def read_dat(func, name):
    """Return (columns, rows) of a function-object .dat file, restarts merged."""
    files = sorted(glob.glob(os.path.join("postProcessing", func, "*", name)),
                   key=lambda f: float(os.path.basename(os.path.dirname(f))))
    if not files:
        return None, []
    cols, data = None, {}
    for path in files:
        with open(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    words = line[1:].split()
                    if words and words[0] == "Time":
                        cols = words
                    continue
                vals = line.split()
                if vals and cols and len(vals) == len(cols):
                    data[float(vals[0])] = vals   # later file overrides
    return cols, [data[t] for t in sorted(data)]


def column(cols, rows, name):
    i = cols.index(name)
    return [float(r[i]) for r in rows]


def style(ax, ylabel):
    ax.set_ylabel(ylabel, color=INK)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--no-show", action="store_true", help="do not open a window")
    ap.add_argument("--avg", type=int, default=100,
                    help="iterations averaged for the printed final values")
    ap.add_argument("-o", "--output", default="convergence.png")
    args = ap.parse_args()

    if args.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    rcols, rrows = read_dat("solverInfo", "solverInfo.dat")
    fcols, frows = read_dat("forceCoeffs", "coefficient.dat")
    if not rrows and not frows:
        sys.exit("No postProcessing data found - run the case first.")

    fig, axes = plt.subplots(4, 1, figsize=(9, 11), sharex=True,
                             gridspec_kw=dict(height_ratios=[1.6, 1, 1, 1]))
    fig.suptitle("S1223 wing - convergence", color=INK, fontsize=13, x=0.08,
                 ha="left")

    # residuals (initial residual of each equation, log scale)
    ax = axes[0]
    if rrows:
        it = column(rcols, rrows, "Time")
        names = [c for c in rcols if c.endswith("_initial")]
        for k, name in enumerate(names):
            ax.semilogy(it, column(rcols, rrows, name), lw=1.5,
                        color=COLORS[k % len(COLORS)],
                        label=name.replace("_initial", ""))
        ax.legend(loc="upper right", ncol=len(names), frameon=False,
                  fontsize=9, labelcolor=INK)
    style(ax, "Initial residual")

    # force coefficients: one panel each (different scales, no twin axes)
    summary = []
    for ax, (name, label) in zip(axes[1:], [("Cl", "Cl"), ("Cd", "Cd"),
                                           ("CmPitch", "Cm (pitch, nose-up +)")]):
        if frows and name in fcols:
            it = column(fcols, frows, "Time")
            y = column(fcols, frows, name)
            ax.plot(it, y, lw=1.5, color=COLORS[0])
            n = min(args.avg, len(y))
            mean = sum(y[-n:]) / n
            ax.axhline(mean, color=MUTED, lw=1, ls="--")
            ax.annotate(f"last {n} it. mean = {mean:.4f}", xy=(1, mean),
                        xycoords=("axes fraction", "data"), xytext=(-4, 4),
                        textcoords="offset points", ha="right", va="bottom",
                        color=INK, fontsize=9)
            # zoom on the converged part, ignoring the start-up transient
            tail = y[len(y) // 5:] or y
            lo, hi = min(tail), max(tail)
            pad = 0.1 * (hi - lo) + 1e-3 * max(abs(hi), 1e-3)
            ax.set_ylim(lo - pad, hi + pad)
            summary.append((name, mean, n))
        style(ax, label)
    axes[-1].set_xlabel("Iteration", color=INK)

    fig.tight_layout()
    fig.savefig(args.output, dpi=130)
    print(f"Wrote {args.output}")
    if summary:
        vals = {k: v for k, v, _ in summary}
        for k, v, n in summary:
            print(f"  {k:8s} = {v: .5f}   (mean of last {n} iterations)")
        if vals.get("Cd"):
            print(f"  L/D      = {vals['Cl'] / vals['Cd']: .2f}")
    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
