"""
Run simpleFoam (in parallel) and watch it.

 * divergence (NaN / FPE / exploding residuals or coefficients / solver crash)
   -> the run is stopped and reported, exit code 2
 * convergence: relative Cl change between two consecutive iterations
   below CL_TOL for CL_WINDOW consecutive iterations (after MIN_ITER)
   -> 'stopAt writeNow' is set in system/controlDict so the solver writes
   the current state and exits cleanly
 * otherwise the run ends at MAX_ITER

usage: python3 run_solver.py <caseDir> <nProcs> [extra mpirun args]
"""
import glob
import json
import math
import os
import re
import subprocess
import sys
import time

import params as P


def coeff_file(case, name):
    files = sorted(glob.glob(os.path.join(case, "postProcessing", name, "*", "*.dat")))
    files = [f for f in files if "coefficient" in os.path.basename(f) or "forceCoeffs" in os.path.basename(f)]
    return files[-1] if files else None


def read_coeffs(path):
    """Return dict column -> list for an OpenFOAM forceCoeffs .dat file."""
    cols, rows = None, []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                toks = line[1:].replace("\t", " ").split()
                if toks and toks[0] == "Time":
                    cols = toks
                continue
            vals = line.split()
            if vals:
                try:
                    rows.append([float(v) for v in vals])
                except ValueError:
                    pass
    if cols is None or not rows:
        return {}
    n = min(len(cols), min(len(r) for r in rows))
    return {cols[k]: [r[k] for r in rows] for k in range(n)}


def set_stop(case):
    cd = os.path.join(case, "system", "controlDict")
    s = open(cd).read()
    s = re.sub(r"stopAt\s+\w+;", "stopAt          writeNow;", s)
    with open(cd, "w") as f:
        f.write(s)
    os.utime(cd, None)


def main():
    case = os.path.abspath(sys.argv[1])
    nprocs = int(sys.argv[2])
    extra = sys.argv[3:]
    log = open(os.path.join(case, "log.simpleFoam"), "w")
    if nprocs > 1:
        cmd = ["mpirun", "-np", str(nprocs)] + extra + ["simpleFoam", "-parallel", "-case", case]
    else:
        cmd = ["simpleFoam", "-case", case]
    print("running:", " ".join(cmd), flush=True)
    proc = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT)

    status, reason = "running", ""
    t0 = time.time()
    last_report, below, it_conv = 0, 0, None
    res_re = re.compile(r"Solving for (\w+), Initial residual = ([^,]+),")
    logpath = os.path.join(case, "log.simpleFoam")
    while True:
        time.sleep(2.0)
        done = proc.poll() is not None
        # ---- divergence checks on the log
        txt = open(logpath, errors="replace").read()
        k = txt.find("Starting time loop")
        tail = txt[max(k, len(txt) - 20000):] if k >= 0 else ""
        if re.search(r"\bnan\b|sigFpe::sigHandler|Floating point exception \(core|FOAM FATAL|"
                     r"MPI_ABORT|BAD TERMINATION", tail, re.I):
            status, reason = "diverged", "NaN / floating point exception / fatal error in solver log"
        else:
            last = {}
            for m in res_re.finditer(tail):
                try:
                    last[m.group(1)] = float(m.group(2))
                except ValueError:
                    last[m.group(1)] = float("nan")
            bad = [k for k, v in last.items() if not math.isfinite(v) or v > 10.0]
            if bad:
                status, reason = "diverged", "residual of %s blew up" % ",".join(bad)
        # ---- coefficient checks
        cf = coeff_file(case, "coeffs_S1223")
        cl = read_coeffs(cf).get("Cl", []) if cf else []
        if cl and status == "running":
            if any(not math.isfinite(c) or abs(c) > 50 for c in cl[-5:]):
                status, reason = "diverged", "Cl became non-physical (%g)" % cl[-1]
            elif len(cl) > 1:
                below = 0
                for k in range(len(cl) - 1, 0, -1):
                    rel = abs(cl[k] - cl[k - 1]) / max(abs(cl[k]), 1e-6)
                    if rel < P.CL_TOL:
                        below += 1
                    else:
                        break
                if len(cl) >= P.MIN_ITER and below >= P.CL_WINDOW and it_conv is None:
                    it_conv = len(cl)
                    status = "converged"
                    reason = ("relative Cl change < %.1e for %d consecutive iterations (iteration %d)"
                              % (P.CL_TOL, P.CL_WINDOW, it_conv))
                    set_stop(case)
                    print("converged:", reason, flush=True)
        if status == "diverged":
            print("DIVERGENCE DETECTED:", reason, flush=True)
            set_stop(case)
            time.sleep(5)
            proc.kill()
            break
        if cl and len(cl) - last_report >= 25:
            last_report = len(cl)
            print("  iter %5d  Cl(S1223) = %.5f   %.0f s" % (len(cl), cl[-1], time.time() - t0), flush=True)
        if done:
            break
    proc.wait()
    end_time = int(float(re.search(r"endTime\s+([\d.eE+]+);",
                                   open(os.path.join(case, "system", "controlDict")).read()).group(1)))
    cf = coeff_file(case, "coeffs_S1223")
    n_done = len(read_coeffs(cf).get("Cl", [])) if cf else 0
    if status == "running":
        if proc.returncode != 0:
            status, reason = "diverged", "solver exited with code %d" % proc.returncode
        elif n_done >= end_time:
            status, reason = "max_iter", "reached the iteration limit (%d) without meeting the Cl criterion" % end_time
        else:
            status, reason = "aborted", "solver stopped after %d iterations" % n_done
    info = dict(status=status, reason=reason, iterations=n_done, wall_time_s=time.time() - t0, nprocs=nprocs)
    with open(os.path.join(case, "run_status.json"), "w") as f:
        json.dump(info, f, indent=2)
    print("solver finished: %s - %s" % (status, reason), flush=True)
    sys.exit(2 if status == "diverged" else (3 if status == "aborted" else 0))


if __name__ == "__main__":
    main()
