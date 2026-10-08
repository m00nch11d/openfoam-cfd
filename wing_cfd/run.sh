#!/bin/bash
# Full pipeline: mesh -> case -> decompose -> simpleFoam (4 ranks) + monitor
# RESUME=1 ./run.sh continues from the latest decomposed time
#                -> reconstruct -> post-processing
set -e
cd "$(dirname "$0")"
export WM_PROJECT_DIR=/usr/share/openfoam FOAM_ETC=/usr/share/openfoam/etc
NP=${NP:-4}
SHIM=$PWD/shim/libdigestshim.so
[ -f "$SHIM" ] || gcc -O2 -shared -fPIC -o "$SHIM" shim/digest_shim.c

if [ ! -f case/constant/polyMesh/faces ]; then
    python3 mesh2d.py
    python3 plot_mesh2d.py
    python3 extrude.py case
fi
python3 setup_case.py case
if [ -z "$RESUME" ]; then
    (cd case && checkMesh > log.checkMesh 2>&1 || true)
    (cd case && decomposePar -force > log.decomposePar 2>&1)
fi

cd case
mpirun --allow-run-as-root -np $NP -x LD_PRELOAD=$SHIM \
    simpleFoam -parallel >> log.simpleFoam 2>&1 &
SOLVER=$!
python3 ../monitor.py . $SOLVER > log.monitor 2>&1 &
wait $SOLVER || true
wait
reconstructPar -latestTime > log.reconstructPar 2>&1
LD_PRELOAD=$SHIM simpleFoam -postProcess -latestTime \
    -func yPlus > log.yPlus 2>&1 || true
cd ..
python3 postprocess.py case
