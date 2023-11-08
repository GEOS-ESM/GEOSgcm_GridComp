#!/usr/bin/env bash

gfortran \
    -g -O3 \
    -fstack-arrays \
    -fopenmp \
    -foffload=nvptx-none \
    -foffload-options=nvptx-none="-O3 -lm -lgfortran -latomic -march=sm_80 -mptx=7.0" \
    gfdl_cloud_microphys.F90 driver.F90 \
    -o microphys_driver.x

#     -foffload=-msoft-stack-reserve-local=4096 \
