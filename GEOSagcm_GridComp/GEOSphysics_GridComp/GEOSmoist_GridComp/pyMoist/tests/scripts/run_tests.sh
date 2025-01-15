#!/bin/bash
clear
# rm -rf ./.gt_cache_*
rm -rf ./.translate-*
export GT4PY_LITERAL_PRECISION=32
export FV3_DACEMODE=Python
export PACE_FLOAT_PRECISION=32
export PACE_TEST_N_THRESHOLD_SAMPLES=0
export PACE_DACE_DEBUG=True
export GT4PY_COMPILE_OPT_LEVEL=0
# for debugging only
export OMP_NUM_THREADS=1
python -m pytest -s \
    --data_path=/Users/ckropiew/netcdfs \
    --backend=gt:cpu_kfirst \
    --grid=default \
    --multimodal_metric \
    --which_modules=icloud \
    --which_rank=0 \
    ..
