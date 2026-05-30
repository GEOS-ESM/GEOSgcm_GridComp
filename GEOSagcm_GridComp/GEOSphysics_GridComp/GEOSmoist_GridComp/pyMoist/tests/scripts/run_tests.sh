#!/bin/bash

# Usage: ./run_tests.sh [/path/to/data] [st:python:cpu:IJK|st:dace:cpu:IJK] [TranslateName]

# NDSL configuration
export NDSL_LITERAL_PRECISION=32
export GT4PY_COMPILE_OPT_LEVEL=0
export NDSL_LOGLEVEL=Debug

# pyMoist configuration
export EXP_NAME='gcm-fp'

python -m pytest -s -v --multimodal_metric \
    --data_path=$1 \
    --backend=$2\
    --which_modules=$3 \
    --grid=default \
    --threshold_overrides_file=./overrides.yml \
    --which_rank=0 \
    ../translate_tests


