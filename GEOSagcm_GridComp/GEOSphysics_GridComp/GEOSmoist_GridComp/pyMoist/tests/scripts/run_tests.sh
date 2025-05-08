#!/bin/bash

rm -rf ./.gt_cache_*
export PACE_FLOAT_PRECISION=32
export GT4PY_LITERAL_PRECISION=32
export PACE_TEST_N_THRESHOLD_SAMPLES=0
export GT4PY_COMPILE_OPT_LEVEL=0
export PACE_TEST_N_THRESHOLD_SAMPLES=0
python -m pytest -s -v --disable-warnings --multimodal_metric \
    --data_path=/Users/kfandric/netcdf \
    --backend=debug\
    --which_rank=0 \
    --which_modules=ComputeUwshcuInv \
    --grid=default \
   ..
