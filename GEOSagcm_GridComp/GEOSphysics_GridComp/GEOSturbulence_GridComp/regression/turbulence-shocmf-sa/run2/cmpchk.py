#!/usr/bin/env python

# Compares this case's MAPL3 checkpoints with the MAPL2 states captured after
# TURBULENCE's run phase PHASE: TURBULENCE's import, internal and export
# states, and the diffused tracers (DataTracers' internals, the capture's TR_*
# import fields). Prints, per variable, the min, max and L2 norm of the
# difference, and a summary of the variables that differ. Run from the run
# directory after the run.

import numpy as np
import netCDF4 as nc4

PHASE = "2"
baseline_dir = "/home/bmauer/input/turbulence/shocmf/C12L181/new-style"
# Undefined points are left out: MAPL2's (1e15) is the baseline's _FillValue, masked on read
MAPL3_UNDEF = -np.finfo(np.float32).max

differ = []

def compare(label, bas, cur, pairs):
    print(f"\n{label}")
    for bname, cname in pairs:
        bas_var = bas[bname][:]
        cur_var = np.ma.masked_equal(cur[cname][0], MAPL3_UNDEF)
        diff = (bas_var - cur_var).filled(0.0)
        diffnorm = np.linalg.norm(diff)
        print(f" {cname:<20} diff (min, max, norm): {np.min(diff)}, {np.max(diff)}, {diffnorm}")
        if diffnorm != 0.0 or np.any(np.isnan(diff)):
            differ.append(f"{label}:{cname}")

for state in ["import", "internal", "export"]:
    baseline = f"{baseline_dir}/TURBULENCE_{state}_after_runPhase{PHASE}.nc"
    current = f"checkpoints/last/TURBULENCE_{state}.nc"
    with nc4.Dataset(baseline) as bas, nc4.Dataset(current) as cur:
        pairs = [(v, v) for v in bas.variables
                 if v in cur.variables and v not in bas.dimensions and v not in cur.dimensions
                 and v not in ("lons", "lats", "corner_lons", "corner_lats")]
        compare(f"TURBULENCE {state}", bas, cur, pairs)

baseline = f"{baseline_dir}/TURBULENCE_import_after_runPhase{PHASE}.nc"
current = "checkpoints/last/DataTracers_internal.nc"
with nc4.Dataset(baseline) as bas, nc4.Dataset(current) as cur:
    pairs = [(v, v[3:].replace(".", "_")) for v in bas.variables if v.startswith("TR_")]
    compare("tracers (DataTracers internal vs TR_*)", bas, cur, pairs)

print("\nVariables that differ:", " ".join(differ) if differ else "none")
