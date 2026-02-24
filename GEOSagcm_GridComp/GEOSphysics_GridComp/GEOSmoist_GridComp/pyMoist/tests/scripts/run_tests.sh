#!/bin/bash

# Usage: ./run_tests.sh [/path/to/data] [debug|dace:cpu_kfirst] [TranslateName]

export NDSL_LITERAL_PRECISION=32
export NDSL_TEST_N_THRESHOLD_SAMPLES=0
export GT4PY_COMPILE_OPT_LEVEL=0
export FV3_DACEMODE=Python

# UW specific
export GT4PY_EXTRA_COMPILE_OPT_FLAGS='-fconstexpr-ops-limit=1000000000'

python -m pytest -s -v --disable-warnings --multimodal_metric \
    -x \
    --data_path=$1 \
    --backend=$2\
    --which_modules=$3 \
    --grid=default \
    --no_report \
    --threshold_overrides_file=./overrides.yml \
    ..
