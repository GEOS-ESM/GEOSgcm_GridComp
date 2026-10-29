#!/usr/bin/env python

import numpy as np
import netCDF4 as nc4

# There is no MAPL2 baseline for radiation with aerosol optics (see README).
# The MAPL2 radiation capture has no aerosol optics, and solar-sa and irrad-sa
# are zero-diff against it, so this compares SOLAR and IRRAD against that
# capture: everything must match except the fields that include aerosols,
# which are listed in AEROSOL. GOCART2G, DU, SS and SU have no capture at this
# step and are left to the ctest comparison against checkpoints/last.
baseline_dir = "/home/pchakrab/input/radiation/C12L181/after/new-style"
comps = {
    "SOLAR": ["import", "internal", "export"],
    "IRRAD": ["import", "internal", "export"],
}

AEROSOL = {
    "SOLAR": {
        "internal": {"FSWN", "FSCN", "FSWUN", "FSCUN", "FSWBANDN",
                     "DRUVRN", "DFUVRN", "DRPARN", "DFPARN", "DRNIRN", "DFNIRN",
                     "OSRB08RGN", "OSRB09RGN", "OSRB10RGN", "OSRB11RGN"},
        "export": {"DRUVRN", "DFUVRN", "DRPARN", "DFPARN", "DRNIRN", "DFNIRN"},
    },
    "IRRAD": {
        "internal": {"FLX", "FLXD", "FLXU", "FLC", "FLCD", "FLCU",
                     "DFDTS", "DFDTSC", "SFCEM",
                     "OLRB05RG", "OLRB06RG", "OLRB07RG", "OLRB08RG",
                     "OLRB09RG", "OLRB10RG", "OLRB11RG", "OLRB15RG",
                     "DOLRB05RGDT", "DOLRB06RGDT", "DOLRB07RGDT", "DOLRB08RGDT",
                     "DOLRB09RGDT", "DOLRB10RGDT", "DOLRB11RGDT", "DOLRB15RGDT"},
        "export": {"DFDTS", "LWS0", "SFCEM0"},
    },
}

HUGE = np.finfo(np.float32).max
SKIP = {"time", "lons", "lats", "corner_lons", "corner_lats", "anchor", "contacts", "orientation"}


def undef(x, fills):
    # NaN, MAPL_UNDEF (+huge in MAPL2, -huge in MAPL3), or either file's _FillValue
    u = ~np.isfinite(x) | (np.abs(x) == HUGE)
    for f in fills:
        u |= x == f
    return u


unexpected = []
for comp, states in comps.items():
    for state in states:
        baseline = f"{baseline_dir}/{comp}_{state}_after_runPhase1.nc"
        current = f"checkpoints/last/{comp}_{state}.nc"
        aerosol = AEROSOL[comp].get(state, set())
        print(f"\n{comp} STATE: {state}\nBASELINE: {baseline}\nCURRENT: {current}")
        with nc4.Dataset(baseline) as bas, nc4.Dataset(current) as cur:
            bas.set_auto_mask(False)
            cur.set_auto_mask(False)
            for var in bas.variables:
                if var not in cur.variables or var in bas.dimensions or var in SKIP:
                    continue
                if not np.issubdtype(bas[var].dtype, np.number):
                    continue
                b, c = bas[var][:], cur[var][0]
                fills = [getattr(v, "_FillValue") for v in (bas[var], cur[var]) if hasattr(v, "_FillValue")]
                bu, cu = undef(b, fills), undef(c, fills)
                diff = np.where(bu | cu, 0.0, np.abs(b.astype(np.float64) - c))
                differs = (bu != cu).any() or diff.max() > 0
                tag = "(aerosol)" if var in aerosol else ""
                if differs and var not in aerosol:
                    tag = "UNEXPECTED"
                    unexpected.append(f"{comp}_{state}:{var}")
                if (bu != cu).any():
                    print(f" {var:<20} UNDEF MISMATCH: {(bu & ~cu).sum()} baseline-only, {(cu & ~bu).sum()} current-only {tag}")
                print(f" {var:<20} diff (min/max/norm): {diff.min()}, {diff.max()}, {np.linalg.norm(diff)} {tag}")

print("\nUnexpected differences:", ", ".join(unexpected) if unexpected else "none")
