#!/bin/bash

python -m pytest -v -s \
    --data_path=../../test_data/moist \
    --backend=numpy \
    --which_rank=0 \
    --which_modules=RadCouple \
    --grid=default \
    ..
