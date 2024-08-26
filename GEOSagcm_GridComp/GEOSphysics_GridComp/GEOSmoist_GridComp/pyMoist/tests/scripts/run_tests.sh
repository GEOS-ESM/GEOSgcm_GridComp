#!/bin/bash
export PACE_FLOAT_PRECISION=32
python -m pytest -s\
    --data_path=/home/charleskrop/netcdfs \
    --backend=dace:cpu \
    --which_rank=0 \
    --which_modules=QSat \
    --grid=default \
    ..