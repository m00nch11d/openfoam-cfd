"""
Run simpleFoam (steady initialisation) or pimpleFoam (plunging airfoil) in
parallel and watch it.

 * divergence (NaN / FPE / fatal error / exploding residuals, Courant number
   or force coefficients / solver crash) -> run stopped, reported, exit code 2
 * unsteady convergence: after every completed plunge cycle the Cl
   statistics (Cl max, Cl min, Cl rms) are compared with those of the
   previous cycle; when each changes by less than CYCLE_TOL of its own value
   the run is stopped cleanly ('stopAt writeNow')
 * otherwise the run ends at the cycle limit

usage: python3 run_solver.py <caseDir> <steady|unsteady> <nProcs> [extra mpirun args]
"""
import glob
import json
import math
import os
import re
import subprocess
import sys
import time

import numpy as np

import params as P


def _start_time(path):
    try:
        return float(os.path.basename(os.path.dirname(path)))
    except ValueError:
        return -1.0


def read_coeffs(path):
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


def read_all_coeffs(case, name="coeffs"):
    """forceCoeffs history over all (re)starts; later runs override earlier ones."""
    files = glob.glob(os.path.join(case, "postProcessing", name, "*", "*.dat"))
    files = sorted([f for f in files if "coefficient" in os.path.basename(f)
                    or "forceCoeffs" in os.path.basename(f)], key=_start_time)
    rows, cols = {}, None
    for f in files:
        d = read_coeffs(f)
        if not d:
            continue
        cols = list(d)
        for k in range(len(d["Time"])):
            rows[round(d["Time"][k], 9)] = [d[c][k] for c in cols]
    if not rows:
        return {}
    ts = sorted(rows)
    return {c: np.array([rows[t][i] for t in ts]) for i, c in enumerate(cols)}


def cycle_stats(d):
    """Statistics of every complete plunge cycle."""
    t, cl, cd = d["Time"], d["Cl"], d["Cd"]
    out = []
    n_done = int(math.floor(t[-1] / P.PERIOD + 1e-6))
    for n in range(n_done):
        m = (t > n * P.PERIOD + 1e-9) & (t <= (n + 1) * P.PERIOD + 1e-9)
        if m.sum() < 0.9 * P.STEPS_PER_CYCLE:
            continue
        c, dd = cl[m], cd[m]
        out.append(dict(cycle=n + 1, Cl_max=float(c.max()), Cl_min=float(c.min()),
                        Cl_rms=float(np.sqrt(np.mean((c - c.mean()) ** 2))),
                        Cl_mean=float(c.mean()), Cd_mean=float(dd.mean()),
                        Cd_max=float(dd.max()), Cd_min=float(dd.min())))
    return out


def cycle_change(a, b):
    """Largest relative change of the Cl cycle statistics from cycle a to b,
    each relative to its own magnitude (the user criterion is on Cl). The
    change of the mean Cd is returned for information only."""
    ch = {k: abs(b[k] - a[k]) / max(abs(b[k]), 1e-6) for k in ("Cl_max", "Cl_min", "Cl_rms")}
    info = dict(ch)
    info["Cd_mean (info)"] = abs(b["Cd_mean"] - a["Cd_mean"]) / max(abs(b["Cd_mean"]), 1e-6)
    return max(ch.values()), info


def set_stop(case):
    cd = os.path.join(case, "system", "controlDict")
    s = open(cd).read()
    s = re.sub(r"stopAt\s+\w+;", "stopAt          writeNow;", s)
    with open(cd, "w") as f:
        f.write(s)
    os.utime(cd, None)


def main():
    case = os.path.abspath(sys.argv[1])
    stage = sys.argv[2]
    nprocs = int(sys.argv[3])
    extra = sys.argv[4:]
    app = "simpleFoam" if stage == "steady" else "pimpleFoam"
    logpath = os.path.join(case, "log." + app)
    log_start = os.path.getsize(logpath) if os.path.exists(logpath) else 0
    log = open(logpath, "a")
    if nprocs > 1:
        cmd = ["mpirun", "-np", str(nprocs)] + extra + [app, "-parallel", "-case", case]
    else:
        cmd = [app, "-case", case]
    print("running:", " ".join(cmd), flush=True)
    proc = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT)

    status, reason = "running", ""
    t0 = time.time()
    last_report, stats, checked = -1, [], 0
    res_re = re.compile(r"Solving for (\w+), Initial residual = ([^,]+),")
    co_re = re.compile(r"Courant Number mean: \S+ max: (\S+)")
    while True:
        time.sleep(3.0)
        done = proc.poll() is not None
        with open(logpath, errors="replace") as fh:
            fh.seek(log_start)
            txt = fh.read()
        k = txt.find("Starting time loop")
        tail = txt[max(k, len(txt) - 30000):] if k >= 0 else ""
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
            bad = [n for n, v in last.items() if not math.isfinite(v) or v > 10.0]
            cos = [float(m.group(1)) for m in co_re.finditer(tail[-5000:])]
            if bad:
                status, reason = "diverged", "residual of %s blew up" % ",".join(bad)
            elif cos and (not math.isfinite(cos[-1]) or cos[-1] > 1e4):
                status, reason = "diverged", "Courant number exploded (%g)" % cos[-1]
        d = read_all_coeffs(case)
        if d and status == "running":
            cl = d["Cl"]
            if not np.all(np.isfinite(cl[-5:])) or np.abs(cl[-5:]).max() > 50:
                status, reason = "diverged", "Cl became non-physical (%g)" % cl[-1]
            elif stage == "unsteady":
                stats = cycle_stats(d)
                if len(stats) > checked:
                    checked = len(stats)
                    s = stats[-1]
                    msg = ("  cycle %d: Cl %.3f..%.3f (rms %.3f), mean Cd %.4f"
                           % (s["cycle"], s["Cl_min"], s["Cl_max"], s["Cl_rms"], s["Cd_mean"]))
                    if len(stats) >= 2:
                        mx, info = cycle_change(stats[-2], stats[-1])
                        msg += ", change vs previous cycle: " + ", ".join(
                            "%s %.1f %%" % (k, 100 * v) for k, v in info.items())
                        if s["cycle"] >= P.MIN_CYCLES and mx < P.CYCLE_TOL:
                            status = "converged"
                            reason = ("cycle statistics changed by %.1f %% (< %.0f %%) from cycle "
                                      "%d to %d" % (100 * mx, 100 * P.CYCLE_TOL, s["cycle"] - 1,
                                                    s["cycle"]))
                            set_stop(case)
                    print(msg, flush=True)
        if status == "diverged":
            print("DIVERGENCE DETECTED:", reason, flush=True)
            set_stop(case)
            time.sleep(5)
            proc.kill()
            break
        if d:
            n = len(d["Cl"])
            step = max(P.STEPS_PER_CYCLE // 10, 1) if stage == "unsteady" else 50
            if n // step > last_report:
                last_report = n // step
                print("  t = %.4f  Cl = %+.4f  Cd = %+.5f   (%.0f s)"
                      % (d["Time"][-1], d["Cl"][-1], d["Cd"][-1], time.time() - t0), flush=True)
        if done:
            break
    proc.wait()
    d = read_all_coeffs(case)
    t_end = float(d["Time"][-1]) if d else 0.0
    end_time = float(re.search(r"endTime\s+([\d.eE+-]+);",
                               open(os.path.join(case, "system", "controlDict")).read()).group(1))
    if status == "running":
        if proc.returncode != 0:
            status, reason = "diverged", "solver exited with code %d" % proc.returncode
        elif t_end >= end_time - 0.5 * (P.DT if stage == "unsteady" else 1):
            status = "max_cycles" if stage == "unsteady" else "done"
            reason = ("reached the limit of %.0f cycles without meeting the %.0f %% cycle criterion"
                      % (end_time / P.PERIOD, 100 * P.CYCLE_TOL)) if stage == "unsteady" else \
                "steady initialisation finished"
        else:
            status, reason = "aborted", "solver stopped at t = %g" % t_end
    info = dict(stage=stage, status=status, reason=reason, end_time=t_end,
                wall_time_s=time.time() - t0, nprocs=nprocs, cycles=stats)
    with open(os.path.join(case, "run_status.json"), "w") as f:
        json.dump(info, f, indent=2)
    print("solver finished: %s - %s" % (status, reason), flush=True)
    sys.exit(2 if status == "diverged" else (3 if status == "aborted" else 0))


if __name__ == "__main__":
    main()
