#!/bin/bash
rm -rf /home/charleskrop/GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/pyMoist/tests/scripts/.gt_cache_000000
export PACE_FLOAT_PRECISION=32
python -m pytest -s --disable-warnings \
    --data_path=/home/charleskrop/netcdfs \
    --backend=dace:cpu \
    --which_rank=0 \
    --which_modules=GFDL_1M \
    --grid=default \
    ..