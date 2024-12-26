#!/bin/bash
clear
rm -rf ./.gt_cache_*
rm -rf ./.translate-*
export GT4PY_LITERAL_PRECISION=32
export FV3_DACEMODE=Python
export PACE_FLOAT_PRECISION=32
python -m pytest -s \
    --data_path=/Users/ckropiew/netcdfs \
    --backend=debug \
    --grid=default \
    --multimodal_metric \
    --which_modules=GFDL_1M_driver \
    --which_rank=0 \
    ..
