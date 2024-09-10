#!/bin/bash
rm -rf ../.gt_cache_*
export PACE_FLOAT_PRECISION=32
python -m pytest -s --disable-warnings --multimodal_metric \
    --data_path=/home/charleskrop/netcdfs \
    --backend=dace:cpu \
    --which_rank=0 \
    --which_modules=QSat \
    --grid=default \
    ..