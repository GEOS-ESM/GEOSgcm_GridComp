#!/bin/bash
rm -rf ./.gt_cache_*
export PACE_FLOAT_PRECISION=32
python -m pytest -s --disable-warnings\
    --data_path=/home/kf041796/netcdfs \
    --backend=dace:cpu \
    --which_rank=0 \
    --which_modules=Conden \
export FV3_DACEMODE=Python
python -m pytest -v -s -x\
    --data_path=../../test_data/11.5.2/Moist/TBC_C24_L72_Debug \
    --backend=dace:cpu \
    --grid=default \
    --multimodal_metric \
    ..
