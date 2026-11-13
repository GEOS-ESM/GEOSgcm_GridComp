#!/bin/bash
# Runs one regression case by hand (not through ctest) in a scratch directory, so
# cmpchk.py can compare it with the MAPL2 capture and the output can be kept as
# the case's checkpoints/last baseline.
#
# usage: run_case_by_hand.sh CASE RUN      e.g. turbulence-default-sa run1
# environment:
#   GEOS_BUILD  build directory of the fixture (has bin/GEOS.x and lib/)   [required]
#   DATA_ROOT   regression data tree down to GEOSturbulence_GridComp       [required]
#   WORK_DIR    scratch directory   [default: ./turbrun]
# The compiler environment (modules) must already be loaded, e.g. the one the
# build used. The run directory is $WORK_DIR/CASE-RUN; the log is run.log in it.
set -e
case=$1; run=$2
: "${GEOS_BUILD:?}" "${DATA_ROOT:?}"
here=$(cd "$(dirname "$0")/.." && pwd)
work=${WORK_DIR:-$PWD/turbrun}/$case-$run
rm -rf "$work"; mkdir -p "$work/checkpoints"
cp -r "$here/$case/$run/." "$work/"
cp -r "$DATA_ROOT/$case/$run/checkpoints/2000-04-14T22:00:00" "$work/checkpoints/"
cd "$work"
export LD_LIBRARY_PATH=$GEOS_BUILD/lib:$LD_LIBRARY_PATH OMP_NUM_THREADS=1
if mpiexec -n 6 "$GEOS_BUILD/bin/GEOS.x" mapl.yaml > run.log 2>&1; then echo RUN_OK; else echo "RUN_FAILED (see $work/run.log)"; exit 1; fi
