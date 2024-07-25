#!/bin/bash
export PACE_FLOAT_PRECISION=32
python -m pytest -v -s -x\
    --pdb \
    --data_path=../../test_data/geos_11.5.2/moist \
    --backend=dace:cpu \
    --which_rank=0 \
    --which_modules=AerActivation \
    --grid=default \
    ..
