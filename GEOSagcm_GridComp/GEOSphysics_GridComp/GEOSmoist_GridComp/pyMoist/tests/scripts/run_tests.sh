#!/bin/bash
export PACE_FLOAT_PRECISION=32
python -m pytest -s --disable-warnings\
    --data_path=/home/kf041796/netcdfs \
    --backend=dace:cpu \
    --which_rank=0 \
    --which_modules=Conden \
    --grid=default \
    ..
