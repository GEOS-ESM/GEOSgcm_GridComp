#!/usr/bin/env bash

nvfortran \
    -g -fast -Minfo \
    -mp=gpu \
    -gpu=cc80 \
    -o microphys_driver.x \
    gfdl_cloud_microphys.F90 gfdl_cloud_microphys_orig.F90 driver.F90
