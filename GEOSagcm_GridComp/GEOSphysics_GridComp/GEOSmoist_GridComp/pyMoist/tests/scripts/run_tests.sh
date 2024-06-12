#!/bin/bash

python -m pytest -v -s -x\
    --data_path=../../test_data/geos_11.5.2/moist \
    --backend=numpy \
    --which_rank=0 \
    --which_modules=RadCouple \
    --grid=default \
    ..