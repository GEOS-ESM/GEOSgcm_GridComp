#!/usr/bin/env bash

nvfortran -g -fast -Minfo -mp=gpu -gpu=cc80 gfdl_cloud_microphys.F90 driver.F90 -o microphys_driver.x
