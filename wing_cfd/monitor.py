"""Watch a running simpleFoam case; stop it on divergence or CL convergence.

Stop rules (checked every iteration):
  * divergence: NaN/Inf or FATAL in the log; after iteration 100 |CL| > DIVERGE_CL,
    CD < 0 or CD > 5; an initial residual > 1 after iteration 50 -> stopAt writeNow,
    status DIVERGED
  * convergence: after MIN_ITER iterations, the relative CL change between
    consecutive iterations |CL_n - CL_n-1| / |CL_n| stays below CL_TOL for
    CL_WINDOW consecutive iterations -> stopAt writeNow, status CONVERGED
  * endTime (N_ITER_MAX) reached -> status MAX_ITER
Thresholds are set in params.py.
"""
import math
import os
import re
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import params as P  # noqa: E402

CL_TOL = P.CL_TOL
CL_WINDOW = P.CL_WINDOW
MIN_ITER = P.MIN_ITER

case = sys.argv[1] if len(sys.argv) > 1 else "case"
log = os.path.join(case, "log.simpleFoam")
coef_dir = os.path.join(case, "postProcessing/forceCoeffs")
status_file = os.path.join(case, "run_status.txt")
env = dict(os.environ)          # OpenFOAM environment set by run_all.sh


def stop(reason):
    subprocess.run(["foamDictionary", "-case", case, "-entry", "stopAt",
                    "-set", "writeNow", "system/controlDict"],
                   env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    with open(status_file, "w") as f:
        f.write(reason + "\n")
    print(reason, flush=True)


def read_coeffs():
    """forceCoeffs history over all restart directories (later wins)."""
    rows = {}
    if not os.path.isdir(coef_dir):
        return []
    for d in sorted(os.listdir(coef_dir), key=float):
        path = os.path.join(coef_dir, d, "coefficient.dat")
        if not os.path.exists(path):
            continue
        hdr = None
        with open(path) as f:
            for line in f:
                if line.startswith("#"):
                    if "Time" in line:
                        hdr = line[1:].split()
                    continue
                v = line.split()
                if hdr and len(v) == len(hdr):
                    r = dict(zip(hdr, map(float, v)))
                    rows[r["Time"]] = r
    return [rows[k] for k in sorted(rows)]


def main(pid):
    small = 0
    done_upto = 0
    first_small = None
    while True:
        alive = os.path.exists(f"/proc/{pid}")
        txt = open(log, errors="ignore").read() if os.path.exists(log) else ""
        if re.search(r"\bnan\b|\binf\b|FOAM FATAL|printStack|Floating point exception \(",
                     txt, re.I):
            stop("DIVERGED: NaN/Inf or fatal error in solver log")
            return
        res = re.findall(r"^Time = (\d+)\n(.*?)(?=^Time = |\Z)", txt[-200000:],
                         re.M | re.S)
        rows = read_coeffs()
        for r in rows[done_upto:]:
            it, cl, cd = int(r["Time"]), r["Cl"], r["Cd"]
            if not all(map(math.isfinite, (cl, cd))) or (it > 100 and (abs(cl) > P.DIVERGE_CL or cd < 0 or cd > 5)):
                stop(f"DIVERGED: unphysical coefficients at iteration {it} "
                     f"(CL={cl:.4g}, CD={cd:.4g})")
                return
            if first_small is None and done_upto > 0:
                prev = rows[done_upto - 1]["Cl"]
                if abs(cl - prev) < 0.05 * abs(cl):
                    first_small = it
            if done_upto > 0:
                prev = rows[done_upto - 1]["Cl"]
                small = small + 1 if abs(cl - prev) < CL_TOL * max(abs(cl), 1e-6) else 0
            done_upto += 1
            if it >= MIN_ITER and small >= CL_WINDOW:
                stop(f"CONVERGED: |dCL|/CL < {CL_TOL:g} for {CL_WINDOW} consecutive "
                     f"iterations at iteration {it} (CL={cl:.5f}, CD={cd:.5f}); "
                     f"literal 5% per-iteration rule first met at iteration {first_small}")
                return
        for it, block in res[-1:]:
            if int(it) > 50:
                for m in re.finditer(r"Solving for (\w+), Initial residual = ([\deE.+-]+)", block):
                    if float(m.group(2)) > 1.0:
                        stop(f"DIVERGED: initial residual of {m.group(1)} = {m.group(2)} "
                             f"at iteration {it}")
                        return
        if not alive:
            last = rows[-1] if rows else {}
            with open(status_file, "w") as f:
                f.write(f"MAX_ITER: solver finished at iteration {int(last.get('Time', 0))}; "
                        f"literal 5% rule first met at iteration {first_small}\n")
            return
        time.sleep(20)


if __name__ == "__main__":
    main(int(sys.argv[2]))
