#!/usr/bin/env bash

nvfortran \
    -g -fast -Minfo \
    -mp=gpu \
    -gpu=cc80 \
    -o microphys_driver.x \
    input.F90 gfdl_cloud_microphys.F90 gfdl_cloud_microphys_cpu.F90 driver.F90
