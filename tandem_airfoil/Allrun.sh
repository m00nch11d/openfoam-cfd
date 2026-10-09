#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Tandem airfoil S1223 (c=0.30 m) + NACA 0009 (c=0.10 m), Re = 3e5, AoA = 0
# k-omega SST + gamma-ReTheta transition (kOmegaSSTLM), simpleFoam
#
# One command does everything (mesh -> case -> parallel run -> post):
#
#     ./Allrun.sh                 # 8 MPI ranks, max 1000 iterations
#     ./Allrun.sh -np 12 -maxiter 1500
#     ./Allrun.sh -cltol 1e-3     # Cl convergence tolerance (relative change
#                                 # per iteration, held for 50 iterations)
#
# Results (PNG figures, results.md, results.json) are written to ./results
# The OpenFOAM case is kept in ./run (open run/case.foam in ParaView).
#
# Needs: OpenFOAM (openfoam.com / ESI release, v2006 or newer), python3.
# Python packages (gmsh numpy scipy matplotlib) are installed into ./.venv
# automatically if they are missing.
#------------------------------------------------------------------------------
set -eo pipefail
cd "$(dirname "$(readlink -f "$0")")"
HERE=$PWD

NP=8
MAXIT=1000
while [ $# -gt 0 ]; do
    case "$1" in
        -np) NP="$2"; shift 2 ;;
        -maxiter) MAXIT="$2"; shift 2 ;;
        -cltol) export CL_TOL="$2"; shift 2 ;;
        -h|--help) sed -n '2,20p' "$0"; exit 0 ;;
        *) echo "unknown option $1"; exit 1 ;;
    esac
done

say() { printf '\n\033[1;34m==> %s\033[0m\n' "$*"; }
die() { printf '\n\033[1;31mERROR: %s\033[0m\n' "$*"; exit 1; }

#------------------------------------------------------------------ OpenFOAM
if ! command -v simpleFoam >/dev/null 2>&1; then
    for f in $(ls -d /usr/lib/openfoam/openfoam*/etc/bashrc /opt/openfoam*/etc/bashrc \
                     "$HOME"/OpenFOAM/OpenFOAM-v*/etc/bashrc 2>/dev/null | sort -V | tail -1); do
        say "sourcing $f"
        set +eu; source "$f"; set -e
    done
fi
command -v simpleFoam >/dev/null 2>&1 || die "OpenFOAM not found: source its etc/bashrc first"
case "${WM_PROJECT_VERSION:-}" in
    v*) ;;
    *) echo "WARNING: OpenFOAM '${WM_PROJECT_VERSION:-?}' does not look like an openfoam.com" \
            "(ESI, vYYMM) release; this case is written for v2006+." ;;
esac
say "OpenFOAM ${WM_PROJECT_VERSION:-?} ($(command -v simpleFoam))"

#------------------------------------------------------------------ python
PY=python3
if ! $PY -c "import gmsh, numpy, scipy, matplotlib" >/dev/null 2>&1; then
    if [ ! -x .venv/bin/python ]; then
        say "creating python venv (.venv) with gmsh numpy scipy matplotlib"
        python3 -m venv .venv || die "python3 -m venv failed (sudo apt install python3-venv)"
        .venv/bin/pip install -q --upgrade pip
        .venv/bin/pip install -q gmsh numpy scipy matplotlib
    fi
    PY=$HERE/.venv/bin/python
    $PY -c "import gmsh" 2>/dev/null || die "gmsh python module cannot load; try: sudo apt install libglu1-mesa libxcursor1 libxinerama1 libxft2"
fi
$PY scripts/params.py

#------------------------------------------------------------------ case
RUN=$HERE/run
rm -rf "$RUN"
mkdir -p "$RUN"

say "writing case dictionaries (np=$NP, max iterations=$MAXIT)"
$PY scripts/setup_case.py "$RUN" "$NP" "$MAXIT"

say "meshing (gmsh + structured quad layers)"
$PY scripts/mesh.py "$RUN/tandem.msh" | tee "$RUN/log.mesh"

say "converting mesh (gmshToFoam)"
(
    cd "$RUN"
    gmshToFoam tandem.msh > log.gmshToFoam 2>&1 || { tail -20 log.gmshToFoam; exit 1; }
    foamDictionary -entry entry0/frontAndBack/type -set empty constant/polyMesh/boundary >/dev/null
    for p in airfoil1 airfoil2; do
        foamDictionary -entry entry0/$p/type -set wall constant/polyMesh/boundary >/dev/null
    done
    checkMesh > log.checkMesh 2>&1 || true
    grep -E "cells:|non-orthogonality|skewness|Mesh OK|Failed" log.checkMesh
    touch case.foam
)

#------------------------------------------------------------------ solve
MPIARGS=()
if [ "$NP" -gt 1 ]; then
    say "decomposing for $NP ranks"
    (cd "$RUN" && decomposePar -force > log.decomposePar 2>&1) || die "decomposePar failed"
    if mpirun --version 2>&1 | grep -qi "open mpi"; then
        [ "$NP" -gt "$(nproc)" ] && MPIARGS+=(--oversubscribe)
        [ "$(id -u)" -eq 0 ] && MPIARGS+=(--allow-run-as-root)
    fi
fi

say "running simpleFoam (watching for divergence and Cl convergence)"
set +e
$PY scripts/run_solver.py "$RUN" "$NP" "${MPIARGS[@]}"
rc=$?
set -e
if [ $rc -eq 2 ]; then
    cat "$RUN/run_status.json"
    die "the solution DIVERGED - see $RUN/log.simpleFoam (no post-processing done)"
fi
[ $rc -eq 0 ] || die "solver driver failed (code $rc)"

if [ "$NP" -gt 1 ]; then
    (cd "$RUN" && reconstructPar -latestTime > log.reconstructPar 2>&1) || die "reconstructPar failed"
fi

#------------------------------------------------------------------ post
say "post-processing"
$PY scripts/post.py "$RUN" "$HERE/results" "$RUN/tandem_geom.npz"
say "done - results in $HERE/results"
cat "$HERE/results/results.md"
