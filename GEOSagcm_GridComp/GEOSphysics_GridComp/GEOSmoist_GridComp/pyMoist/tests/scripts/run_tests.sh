#!/bin/bash
export PACE_FLOAT_PRECISION=32
python -m pytest -s\
    --data_path=/home/charleskrop/netcdfs \
    --backend=numpy \
    --which_rank=0 \
    --which_modules=evap_subl_pdf \
    --grid=default \
    ..