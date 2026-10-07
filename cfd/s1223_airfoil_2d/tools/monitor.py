#!/usr/bin/env python3
"""
Watch a running simpleFoam case and stop it (stopAt writeNow) when

  * DIVERGED : a residual or force coefficient becomes NaN/inf, or
               |Cl| > clDivergence, or |Cd| > 10;
  * CONVERGED: after minIterations, |Cl_i - Cl_(i-1)| / |Cl_i| < clTol for
               every one of the last clWindow iterations.

The solver itself stops at nIterations. The verdict goes to
results/run_status.txt; Allrun reads it.

Usage: python3 tools/monitor.py [case_dir]   (started by Allrun)
"""
import math
import os
import re
import sys
import time

CASE = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else
                       os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_mesh import read_parameters  # noqa: E402

P = read_parameters(os.path.join(CASE, "caseParameters"))
COEF = os.path.join(CASE, "postProcessing", "forceCoeffs", "0", "coefficient.dat")
RES = os.path.join(CASE, "postProcessing", "solverInfo", "0", "solverInfo.dat")
STATUS = os.path.join(CASE, "results", "run_status.txt")
CONTROL = os.path.join(CASE, "system", "controlDict")


def table(path):
    cols, rows = None, []
    if not os.path.exists(path):
        return cols, rows
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                w = line[1:].split()
                if w and w[0] == "Time":
                    cols = w
            elif line.strip() and cols:
                v = line.split()
                if len(v) == len(cols):
                    rows.append(v)
    return cols, rows


def stop(verdict, reason):
    os.makedirs(os.path.dirname(STATUS), exist_ok=True)
    with open(STATUS, "w") as fh:
        fh.write(f"{verdict}\n{reason}\n")
    with open(CONTROL) as fh:
        txt = fh.read()
    txt = re.sub(r"^(stopAt\s+)\w+;", r"\1writeNow;", txt, flags=re.M)
    with open(CONTROL, "w") as fh:
        fh.write(txt)
    print(f"[monitor] {verdict}: {reason}", flush=True)


def finite(x):
    try:
        return math.isfinite(float(x))
    except ValueError:
        return True          # solver names, true/false flags


def main():
    n_min, win, tol = int(P["minIterations"]), int(P["clWindow"]), P["clTol"]
    cl_max = P["clDivergence"]
    while True:
        time.sleep(2.0)
        cols, rows = table(COEF)
        rcols, rrows = table(RES)
        if rrows:
            bad = [c for c, v in zip(rcols, rrows[-1]) if not finite(v)]
            if bad:
                stop("DIVERGED", f"non-finite {', '.join(bad)} at iteration {rrows[-1][0]}")
                return
        if not rows:
            continue
        it = [float(r[0]) for r in rows]
        cl = [float(r[cols.index("Cl")]) for r in rows]
        cd = [float(r[cols.index("Cd")]) for r in rows]
        if not (math.isfinite(cl[-1]) and math.isfinite(cd[-1])) \
                or abs(cl[-1]) > cl_max or abs(cd[-1]) > 10:
            stop("DIVERGED", f"Cl = {cl[-1]}, Cd = {cd[-1]} at iteration {it[-1]:.0f}")
            return
        if it[-1] >= n_min and len(cl) > win:
            rel = [abs(cl[k] - cl[k - 1]) / max(abs(cl[k]), 1e-12)
                   for k in range(len(cl) - win, len(cl))]
            if max(rel) < tol:
                stop("CONVERGED", f"|dCl/Cl| < {tol:g} for {win} consecutive "
                                  f"iterations (iteration {it[-1]:.0f}, Cl = {cl[-1]:.5f})")
                return


if __name__ == "__main__":
    main()
