#!/usr/bin/env bash

mpif90 \
    -g -O3 \
    -fstack-arrays \
    -fopenmp \
    -foffload=nvptx-none \
    -foffload-options=nvptx-none="-O3 -lm -lgfortran -latomic -march=sm_80 -mptx=7.0" \
    -foffload=-msoft-stack \
    -o microphys_driver.x \
    input.F90 output.F90 gfdl_cloud_microphys.F90 gfdl_cloud_microphys_cpu.F90 driver.F90
