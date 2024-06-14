#!/bin/bash
export PACE_FLOAT_PRECISION=32
python -m pytest --pdb -v -s -x\
    --data_path=../../test_data/rad_coup_data_nc \
    --backend=numpy \
    --which_rank=0 \
    --which_modules=RadCouple \
    --grid=default \
    ..