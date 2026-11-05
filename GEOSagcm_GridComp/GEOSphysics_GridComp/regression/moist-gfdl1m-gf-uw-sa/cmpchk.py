#!/usr/bin/env python

import numpy as np
import netCDF4 as nc4

baseline_dir = "/home/pchakrab/input/moist/gfdl1m-gf-uw/C12L181/after/new-style"
# Undefined points are left out: MAPL2's (1e15) is the baseline's _FillValue, masked on read
MAPL3_UNDEF = -np.finfo(np.float32).max

for state in ["import", "internal", "export"]:
    print("\nSTATE:", state)
    baseline = f"{baseline_dir}/MOIST_{state}_after_runPhase1.nc"
    current = f"checkpoints/last/MOIST_{state}.nc"
    with nc4.Dataset(baseline) as bas, nc4.Dataset(current) as cur:
        for var in bas.variables:
            if var in cur.variables and var not in bas.dimensions:
                bas_var = bas[var][:]
                cur_var = np.ma.masked_equal(cur[var][0], MAPL3_UNDEF)
                diff = (bas_var - cur_var).filled(0.0)
                diffnorm = np.linalg.norm(diff)
                print(f" {var:<20} diff (min, max, norm): {np.min(diff)}, {np.max(diff)}, {diffnorm}")
