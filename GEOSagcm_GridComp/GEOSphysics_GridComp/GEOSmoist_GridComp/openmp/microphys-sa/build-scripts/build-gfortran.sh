#!/usr/bin/env bash

gfortran \
    -g -O3 \
    -fstack-arrays \
    -fopenmp \
    -foffload=nvptx-none \
    -foffload-options=nvptx-none="-O3 -lm -lgfortran -latomic -march=sm_80 -mptx=7.0" \
    -foffload=-msoft-stack \
    -o microphys_driver.x \
    gfdl_cloud_microphys.F90 gfdl_cloud_microphys_orig.F90 driver.F90 \
