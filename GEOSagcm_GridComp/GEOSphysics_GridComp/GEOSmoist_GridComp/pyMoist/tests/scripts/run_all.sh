#!/bin/bash

# Usage: ./run_all.sh [/path/to/data] [st:python:cpu:IJK|st:dace:cpu:IJK]

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
    --backend=$2 \
    --grid=default \
    --threshold_overrides_file=./overrides.yml \
    --which_rank=0 \
    ../translate_tests
