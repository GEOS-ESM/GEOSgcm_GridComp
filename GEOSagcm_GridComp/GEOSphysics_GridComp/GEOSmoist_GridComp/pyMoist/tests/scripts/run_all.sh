#!/bin/bash

# Usage: ./run_tests.sh [/path/to/data] [debug|dace:cpu_kfirst]

export NDSL_LITERAL_PRECISION=32
export NDSL_TEST_N_THRESHOLD_SAMPLES=0
export GT4PY_COMPILE_OPT_LEVEL=0
export FV3_DACEMODE=Python
export NDSL_LOGLEVEL=Critical

# UW specific
export GT4PY_EXTRA_COMPILE_OPT_FLAGS='-fconstexpr-ops-limit=1000000000'

python -m pytest -s --disable-warnings --multimodal_metric \
    -x -v \
    --data_path=$1 \
    --backend=$2\
    --grid=default \
    --no_report \
    --threshold_overrides_file=./overrides.yml \
    ..
    