#!/bin/bash
export PACE_FLOAT_PRECISION=32
python -m pytest -s --disable-warnings\
    --data_path=/home/charleskrop/netcdfs \
    --backend=dace:cpu \
    --which_rank=0 \
    --which_modules=find_klcl \
    --grid=default \
    ..