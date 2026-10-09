#!/bin/bash
# =============================================================================
#  One-command pipeline for Ubuntu 24.04:
#      ./run_all.sh
#  1. installs OpenFOAM (Ubuntu package, v1912), MPI and libraries (sudo apt)
#  2. creates a Python venv (.venv) with gmsh, pyvista, numpy, scipy, matplotlib
#  3. builds the small OpenFOAM work-around library (shim/)
#  4. meshes (2D hybrid C-mesh -> spanwise extrusion), writes the case
#  5. runs simpleFoam in parallel with a monitor (divergence / convergence stop)
#  6. reconstructs and post-processes -> results/
#
#  Settings: edit params.py.   Options:
#      --no-install   skip apt / pip (already set up)
#      --clean        delete case/ and results/ first (fresh mesh + run)
#      --resume       continue the existing run from its latest time
#      --post-only    only re-run the post-processing
# =============================================================================
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")"
HERE=$PWD

INSTALL=1; CLEAN=0; RESUME=0; POST_ONLY=0
for a in "$@"; do
    case "$a" in
        --no-install) INSTALL=0 ;;
        --clean) CLEAN=1 ;;
        --resume) RESUME=1 ;;
        --post-only) POST_ONLY=1 ;;
        -h|--help) sed -n '2,19p' "$0"; exit 0 ;;
        *) echo "unknown option $a"; exit 1 ;;
    esac
done

log() { echo -e "\n=== $(date +%H:%M:%S)  $*"; }

# ------------------------------------------------------------ 1. packages
APT_PKGS="openfoam openmpi-bin build-essential python3-venv python3-dev \
libglu1-mesa libxcursor1 libxinerama1 libxft2 libxrender1 libgl1 libegl1 libosmesa6"
if [ $INSTALL = 1 ]; then
    missing=""
    for p in $APT_PKGS; do dpkg -s "$p" >/dev/null 2>&1 || missing="$missing $p"; done
    if [ -n "$missing" ]; then
        log "Installing system packages:$missing"
        SUDO=""; [ "$(id -u)" -ne 0 ] && SUDO=sudo
        $SUDO apt-get update
        $SUDO env DEBIAN_FRONTEND=noninteractive apt-get install -y $missing
    fi
fi

# ------------------------------------------------------------ 2. python venv
if [ ! -x .venv/bin/python ]; then
    log "Creating Python venv"
    python3 -m venv .venv
    INSTALL_PIP=1
else
    INSTALL_PIP=$INSTALL
fi
PY=$HERE/.venv/bin/python
if [ $INSTALL_PIP = 1 ]; then
    "$PY" -m pip install -q --upgrade pip
    "$PY" -m pip install -q -r requirements.txt
fi

# ------------------------------------------------------------ 3. OpenFOAM env
# The Ubuntu package keeps its etc/ files in /usr/share/openfoam
export WM_PROJECT_DIR=/usr/share/openfoam FOAM_ETC=/usr/share/openfoam/etc
export PYVISTA_OFF_SCREEN=true
command -v simpleFoam >/dev/null || { echo "simpleFoam not found (install openfoam)"; exit 1; }
SHIM=$HERE/shim/libdigestshim.so
if [ ! -f "$SHIM" ] || [ shim/digest_shim.c -nt "$SHIM" ]; then
    gcc -O2 -shared -fPIC -o "$SHIM" shim/digest_shim.c
fi
NP=$("$PY" -c "import params; print(params.N_PROCS)")
CORES=$(lscpu -p=Core,Socket | grep -v '^#' | sort -u | wc -l)
THREADS=$(nproc)
MPI_OPTS=""
[ "$(id -u)" -eq 0 ] && MPI_OPTS="$MPI_OPTS --allow-run-as-root"
if [ "$NP" -gt "$THREADS" ]; then
    MPI_OPTS="$MPI_OPTS --oversubscribe"
elif [ "$NP" -gt "$CORES" ]; then
    MPI_OPTS="$MPI_OPTS --use-hwthread-cpus"      # e.g. 24 ranks on 12 cores / 24 threads
fi
echo "machine: $CORES cores / $THREADS threads; running $NP MPI processes $MPI_OPTS"

if [ $POST_ONLY = 1 ]; then
    log "Post-processing"
    "$PY" postprocess.py case
    exit 0
fi

# ------------------------------------------------------------ 4. mesh + case
if [ $CLEAN = 1 ]; then
    log "Cleaning case/ and results/"
    rm -rf case results template2d.npz
fi
if [ $RESUME = 0 ]; then
    if [ -d case/processor0 ] || [ -f case/log.simpleFoam ]; then
        echo "case/ already holds a run: use --resume to continue it or --clean to start over"
        exit 1
    fi
    log "Meshing (2D template + spanwise extrusion)"
    "$PY" mesh2d.py
    "$PY" plot_mesh2d.py
    "$PY" extrude.py case
fi
"$PY" setup_case.py case
cd case
if [ $RESUME = 0 ]; then
    log "checkMesh"
    checkMesh > log.checkMesh 2>&1 || true
    grep -E "cells:|non-orthogonality Max|Mesh OK|Failed" log.checkMesh || true
    log "decomposePar ($NP subdomains)"
    decomposePar -force > log.decomposePar 2>&1
fi
rm -f run_status.txt

# ------------------------------------------------------------ 5. solve
log "simpleFoam on $NP processes (log: case/log.simpleFoam)"
mpirun $MPI_OPTS -np "$NP" -x LD_PRELOAD="$SHIM" \
    simpleFoam -parallel >> log.simpleFoam 2>&1 &
SOLVER=$!
"$PY" ../monitor.py . $SOLVER > log.monitor 2>&1 &
MON=$!
wait $SOLVER || true
wait $MON || true
echo "status: $(cat run_status.txt 2>/dev/null)"

# ------------------------------------------------------------ 6. post
log "reconstructPar"
reconstructPar -latestTime > log.reconstructPar 2>&1
cd ..
log "Post-processing -> results/"
"$PY" postprocess.py case
log "Done. Summary: results/summary.json, results/stations.csv, images in results/"
