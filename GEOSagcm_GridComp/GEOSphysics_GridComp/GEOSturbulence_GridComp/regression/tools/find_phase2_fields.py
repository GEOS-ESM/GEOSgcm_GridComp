#!/usr/bin/env python3
"""Lists the TURBULENCE imports that differ between the capture after Run1 and
before Run2 (what Surface and Physics change between the phases): the
`phase2_fields` of run2's root.yaml and the fields of data-phase2.yaml.

Usage: find_phase2_fields.py CAPTURE_DIR      (the converted, new-style capture)
"""
import sys
import numpy as np
import netCDF4

d = sys.argv[1]
a = netCDF4.Dataset(f'{d}/TURBULENCE_import_after_runPhase1.nc')
b = netCDF4.Dataset(f'{d}/TURBULENCE_import_before_runPhase2.nc')
diff = [v for v in a.variables
        if v in b.variables and a[v].ndim >= 3 and v not in a.dimensions
        and not np.array_equal(np.ma.filled(a[v][:], 1e30), np.ma.filled(b[v][:], 1e30))]
print(' '.join(sorted(diff)))
