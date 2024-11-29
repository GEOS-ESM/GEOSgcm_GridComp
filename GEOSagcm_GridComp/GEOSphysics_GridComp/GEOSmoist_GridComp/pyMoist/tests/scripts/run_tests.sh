#!/bin/bash
rm -rf ./.gt_cache_*
rm -rf ./.translate-*
export PACE_FLOAT_PRECISION=32
export FV3_DACEMODE=Python
python -m pytest -s \
    --data_path=/Users/ckropiew/netcdfs \
    --backend=dace:cpu \
    --grid=default \
    --multimodal_metric \
    --which_modules=warm_rain \
    --which_rank=0 \
    ..
