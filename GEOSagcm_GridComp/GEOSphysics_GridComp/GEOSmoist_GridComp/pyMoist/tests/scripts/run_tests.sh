#!/bin/bash

# Usage: ./run_tests.sh [/path/to/data] [debug|dace:cpu_kfirst] [TranslateName]

export NDSL_LITERAL_PRECISION=32
export NDSL_TEST_N_THRESHOLD_SAMPLES=0
export GT4PY_COMPILE_OPT_LEVEL=0
export FV3_DACEMODE=BuildAndRun
# export OPENMP_CPPFLAGS=" "
# export OPENMP_LDFLAGS=" "

# UW specific
#export GT4PY_EXTRA_COMPILE_OPT_FLAGS='-fconstexpr-ops-limit=1000000000'

export EXP_NAME='gcm-fp'

python -m pytest -s -v --disable-warnings --multimodal_metric \
    --data_path=$1 \
    --backend=$2\
    --which_modules=$3 \
    --which_rank=0 \
    --grid=default \
    --no_report \
    --threshold_overrides_file=./overrides.yml \
    --which_rank=0 \
    ../translate_tests
