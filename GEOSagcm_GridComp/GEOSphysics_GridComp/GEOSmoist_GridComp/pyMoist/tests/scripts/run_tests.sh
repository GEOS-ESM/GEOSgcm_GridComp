#!/bin/bash

# Usage: ./run_tests.sh [/path/to/data] [st:python:cpu:IJK|st:dace:cpu:IJK] [TranslateName]

# NDSL configuration
export NDSL_LITERAL_PRECISION=32
export GT4PY_COMPILE_OPT_LEVEL=0
export NDSL_LOGLEVEL=Info

# pyMoist configuration
export EXP_NAME='gcm-fp'

# UW specific
export GT4PY_EXTRA_COMPILE_OPT_FLAGS='-fconstexpr-ops-limit=1000000000'

python -m pytest -s -v --multimodal_metric \
    --data_path=$1 \
    --backend=$2\
    --which_modules=$3 \
    --grid=default \
    --no_report \
    --threshold_overrides_file=./overrides.yml \
    --which_rank=0 \
    ../translate_tests
