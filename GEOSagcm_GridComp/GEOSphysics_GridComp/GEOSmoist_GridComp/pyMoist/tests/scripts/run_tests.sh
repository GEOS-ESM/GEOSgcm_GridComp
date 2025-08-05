#!/bin/bash

rm -rf ./.gt_cache_*
export PACE_FLOAT_PRECISION=32
export GT4PY_LITERAL_PRECISION=32
export PACE_TEST_N_THRESHOLD_SAMPLES=0
export GT4PY_COMPILE_OPT_LEVEL=0
export FV3_DACEMODE=Python
export OPENMP_CPPFLAGS=""
export OPENMP_LDFLAGS=""
#export GT4PY_EXTRA_COMPILE_OPT_FLAGS="-fbracket-depth=512"
python -m pytest -s -v --disable-warnings --multimodal_metric \
    --data_path=/Users/kfandric/netcdf \
    --backend=dace:cpu\
    --which_rank=1 \
    --which_modules=ComputeUwshcuInv \
    --grid=default \
    ..
   
