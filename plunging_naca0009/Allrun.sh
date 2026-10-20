#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Plunging NACA 0009 (c = 0.30 m), Re = 3e5, AoA = 0,
# plunge y(t) = A sin(2 pi f t) with A = 1 chord (0.30 m), f = 3 Hz.
# k-omega SST + gamma-ReTheta transition (kOmegaSSTLM); pimpleFoam on a
# rigidly moving mesh, initialised from a steady simpleFoam solution.
#
# One command does everything (mesh -> steady init -> unsteady run -> post):
#
#     ./Allrun.sh                   # 8 MPI ranks, up to 5 plunge cycles
#     ./Allrun.sh -np 16 -maxcycles 8
#     ./Allrun.sh -cycletol 0.05    # periodic-convergence tolerance (cycle
#                                   # statistics change between two cycles)
#     ./Allrun.sh -noanim           # skip the animation
#     ./Allrun.sh -postonly         # only redo post-processing / animation
#                                   # from an existing ./run
#
# Results (PNG figures, GIF/MP4 animation, results.md, results.json) are
# written to ./results; the OpenFOAM cases are kept in ./run/steady and
# ./run/unsteady (open run/unsteady/case.foam in ParaView).
#
# Needs: OpenFOAM (openfoam.com / ESI release, v2006 or newer), python3.
# Python packages (gmsh numpy scipy matplotlib) are installed into ./.venv
# automatically if they are missing.
#------------------------------------------------------------------------------
set -eo pipefail
cd "$(dirname "$(readlink -f "$0")")"
HERE=$PWD

NP=8
MAXCYC=5
ANIM=1
POSTONLY=0
while [ $# -gt 0 ]; do
    case "$1" in
        -np) NP="$2"; shift 2 ;;
        -maxcycles) MAXCYC="$2"; shift 2 ;;
        -cycletol) export CYCLE_TOL="$2"; shift 2 ;;
        -noanim) ANIM=0; shift ;;
        -postonly) POSTONLY=1; shift ;;
        -h|--help) sed -n '2,/^[^#]/p' "$0" | sed '$d'; exit 0 ;;
        *) echo "unknown option $1"; exit 1 ;;
    esac
done

say() { printf '\n\033[1;34m==> %s\033[0m\n' "$*"; }
die() { printf '\n\033[1;31mERROR: %s\033[0m\n' "$*"; exit 1; }

#------------------------------------------------------------------ OpenFOAM
if ! command -v pimpleFoam >/dev/null 2>&1; then
    for f in $(ls -d /usr/lib/openfoam/openfoam*/etc/bashrc /opt/openfoam*/etc/bashrc \
                     "$HOME"/OpenFOAM/OpenFOAM-v*/etc/bashrc 2>/dev/null | sort -V | tail -1); do
        say "sourcing $f"
        set +eu; source "$f"; set -e
    done
fi
command -v pimpleFoam >/dev/null 2>&1 || die "OpenFOAM not found: source its etc/bashrc first"
case "${WM_PROJECT_VERSION:-}" in
    v*) ;;
    *) echo "WARNING: OpenFOAM '${WM_PROJECT_VERSION:-?}' does not look like an openfoam.com" \
            "(ESI, vYYMM) release; this case is written for v2006+." ;;
esac
say "OpenFOAM ${WM_PROJECT_VERSION:-?} ($(command -v pimpleFoam))"

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

RUN=$HERE/run
ST=$RUN/steady
UN=$RUN/unsteady

post() {
    say "post-processing"
    $PY scripts/post.py "$UN" "$HERE/results" "$RUN/plunging_geom.npz"
    if [ "$ANIM" -eq 1 ]; then
        say "animation (GIF, plus MP4 if ffmpeg is installed)"
        $PY scripts/animate.py "$UN" "$HERE/results" "$RUN/plunging_geom.npz"
    fi
    say "done - results in $HERE/results"
    cat "$HERE/results/results.md"
}

if [ "$POSTONLY" -eq 1 ]; then
    [ -d "$UN" ] || die "no existing case in $UN - run ./Allrun.sh first"
    post
    exit 0
fi

MPIARGS=()
if [ "$NP" -gt 1 ] && mpirun --version 2>&1 | grep -qi "open mpi"; then
    [ "$NP" -gt "$(nproc)" ] && MPIARGS+=(--oversubscribe)
    [ "$(id -u)" -eq 0 ] && MPIARGS+=(--allow-run-as-root)
fi

run_stage() {   # run_stage <case> <steady|unsteady>
    local case=$1 stage=$2 rc
    if [ "$NP" -gt 1 ]; then
        (cd "$case" && decomposePar -force > log.decomposePar 2>&1) || die "decomposePar failed"
    fi
    set +e
    $PY scripts/run_solver.py "$case" "$stage" "$NP" "${MPIARGS[@]}"
    rc=$?
    set -e
    if [ $rc -eq 2 ]; then
        cat "$case/run_status.json"
        die "the $stage solution DIVERGED - see $case/log.* (stopped, no post-processing)"
    fi
    [ $rc -eq 0 ] || die "solver driver failed (code $rc)"
    if [ "$NP" -gt 1 ]; then
        local opt=-latestTime
        [ "$stage" = unsteady ] && opt=""          # all snapshots of the last cycle
        (cd "$case" && reconstructPar $opt > log.reconstructPar 2>&1) || die "reconstructPar failed"
    fi
}

#------------------------------------------------------------------ mesh
rm -rf "$RUN"
mkdir -p "$ST"
say "meshing (gmsh + structured quad layer)"
$PY scripts/mesh.py "$RUN/plunging.msh" | tee "$RUN/log.mesh"

say "steady initialisation case (simpleFoam, static airfoil)"
$PY scripts/setup_case.py "$ST" steady "$NP"
(
    cd "$ST"
    gmshToFoam ../plunging.msh > log.gmshToFoam 2>&1 || { tail -20 log.gmshToFoam; exit 1; }
    foamDictionary -entry entry0/frontAndBack/type -set empty constant/polyMesh/boundary >/dev/null
    foamDictionary -entry entry0/airfoil/type -set wall constant/polyMesh/boundary >/dev/null
    checkMesh > log.checkMesh 2>&1 || true
    grep -E "cells:|non-orthogonality|skewness|Mesh OK|Failed" log.checkMesh
)
run_stage "$ST" steady

#------------------------------------------------------------------ unsteady
say "plunging case (pimpleFoam, rigid mesh motion, up to $MAXCYC cycles)"
mkdir -p "$UN/constant"
cp -r "$ST/constant/polyMesh" "$UN/constant/"
$PY scripts/setup_case.py "$UN" unsteady "$NP" "$MAXCYC"
(cd "$UN" && mapFields ../steady -consistent -sourceTime latestTime > log.mapFields 2>&1 \
    && touch case.foam) || die "mapFields failed"
run_stage "$UN" unsteady

#------------------------------------------------------------------ post
post
