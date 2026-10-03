#!/bin/bash

# Usage: ./run_tests.sh [/path/to/data] [st:python:cpu:IJK|st:dace:cpu:IJK] [TranslateName]

# NDSL configuration
export NDSL_LITERAL_PRECISION=32
export GT4PY_COMPILE_OPT_LEVEL=0
export NDSL_LOGLEVEL=Debug

# stree optimization
export NDSL_STREE_OPT=False

# verbose debugging flags
export NDSL_VERBOSE_ORCHESTRATION=False
export NDSL_VERBOSE_SCHEDULE_TREE_OPTIMIZATIONS=False

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

#    --which_savepoint=0 \
