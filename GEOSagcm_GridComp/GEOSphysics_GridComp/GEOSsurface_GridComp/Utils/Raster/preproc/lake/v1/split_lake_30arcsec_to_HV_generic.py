#!/usr/bin/env python3
import os
import numpy as np
import xarray as xr
from datetime import datetime

"""
Split a global 30 arc-sec mask into GEOS 10°x10° H/V tiles.

This is a lightly generalized version of the HydroLAKES-TopoCat lake splitter so
that the same script can be used for lakes or reaches by changing only the
configuration block below.

Default configuration is set for lakes.

Expected global input:
  - Dimensions: lat, lon
  - Size:       21600 x 43200  (30 arc-sec global)
  - Variables:
      * one float field in [0,1]   (fraction-like field)
      * one uint8 field in {0,1}   (binary presence field)

Per-tile output:
  - Dimensions: N_lat, N_lon (1200 x 1200)
  - Coordinates: latitude(N_lat), longitude(N_lon)
  - Variables: copied using names given in CONFIG
  - Global attrs:
      N_lon_global, N_lat_global
      i_ind_offset_LL, j_ind_offset_LL
      CellSize_arc_Secs

Notes:
  - No missing values are written anywhere; NaNs are converted to zero.
  - Ocean / non-feature cells remain zero.
  - For reaches, the float field is typically an occupancy-like fraction, not
    a true area fraction.
"""

# -----------------------------------------------------------------------------
# CONFIGURATION
# Change only this block to switch between lakes and reaches.
# -----------------------------------------------------------------------------
### LAKES ####
# -----------------------------------------------------------------------------
CONFIG = {
    # Input/output locations
    "in_global": "/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/lake/lake_mask/v1/before_splitting_nc4/LakeTopoCat_Global_30arcsec.nc4",
    "out_dir": "/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/lake/lake_mask/v1/",

    # Output filename prefix
    "file_prefix": "LakeTopoCat_30arcsec",

    # Variable names in the global input file
    "frac_var_in": "lake_area_frac",
    "any_var_in":  "lake_presence_any",

    # Variable names to write in each tile file
    "frac_var_out": "lake_area_frac",
    "any_var_out":  "lake_presence_any",

    # Metadata strings
    "region_label": "LakeTopoCat",
    "source_attr": 'HydroLAKES-TopoCat v1.1 (2023); rasterized at 10" and aggregated to 30".',
    "created_by": "borescan",
}

# -----------------------------------------------------------------------------
### REACHES  ####
# -----------------------------------------------------------------------------
#CONFIG = {
    # Input/output locations
#     "in_global": "/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/lake/reach_mask/v1/before_splitting_nc4/ReachTopoCat_Global_30arcsec.nc4",
#     "out_dir": "/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/lake/reach_mask/v1/",

    # Output filename prefix
#     "file_prefix": "ReachTopoCat_30arcsec",

    # Variable names in the global input file
#     "frac_var_in": "reach_occupancy_frac",
#     "any_var_in":  "reach_presence_any",

    # Variable names to write in each tile file
#     "frac_var_out": "reach_occupancy_frac",
#     "any_var_out":  "reach_presence_any",

    # Metadata strings
#     "region_label": "ReachTopoCat",
#     "source_attr": 'HydroLAKES-TopoCat v1.1 (2023) reaches; rasterized at 10" and aggregated to 30".',
#     "created_by": "borescan",
# }

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# 30 arc-sec global grid
N_LON_GLOBAL = 43200
N_LAT_GLOBAL = 21600

# 10° x 10° tile size at 30 arc-sec
NC_10 = 1200
NR_10 = 1200


def main():
    in_global = CONFIG["in_global"]
    out_dir = CONFIG["out_dir"]
    file_prefix = CONFIG["file_prefix"]
    frac_var_in = CONFIG["frac_var_in"]
    any_var_in = CONFIG["any_var_in"]
    frac_var_out = CONFIG["frac_var_out"]
    any_var_out = CONFIG["any_var_out"]
    region_label = CONFIG["region_label"]
    source_attr = CONFIG["source_attr"]
    created_by = CONFIG["created_by"]

    os.makedirs(out_dir, exist_ok=True)

    ds = xr.open_dataset(in_global)

    # Get global 1-D coords
    lat_g = ds["lat"].astype("float32")
    lon_g = ds["lon"].astype("float32")

    # Data vars
    frac_g = ds[frac_var_in].astype("float32")
    any_g = ds[any_var_in].astype("uint8")

    # Basic sanity checks
    assert frac_g.sizes["lon"] == N_LON_GLOBAL and frac_g.sizes["lat"] == N_LAT_GLOBAL, \
        f"Unexpected grid for {frac_var_in}: {dict(frac_g.sizes)}"
    assert any_g.sizes["lon"] == N_LON_GLOBAL and any_g.sizes["lat"] == N_LAT_GLOBAL, \
        f"Unexpected grid for {any_var_in}: {dict(any_g.sizes)}"

    created = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    for jx in range(1, 18 + 1):      # V 01..18
        for ix in range(1, 36 + 1):  # H 01..36
            i0 = (ix - 1) * NC_10   # 0-based lon start
            j0 = (jx - 1) * NR_10   # 0-based lat start

            hh = f"{ix:02d}"
            vv = f"{jx:02d}"
            out = os.path.join(out_dir, f"{file_prefix}_H{hh}V{vv}.nc4")

            # Subset coords (1-D)
            lat_sub = lat_g.isel(lat=slice(j0, j0 + NR_10)).values
            lon_sub = lon_g.isel(lon=slice(i0, i0 + NC_10)).values

            # Subset data (2-D), enforce no-missing + valid ranges
            frac_sub = frac_g.isel(lat=slice(j0, j0 + NR_10), lon=slice(i0, i0 + NC_10)).values
            any_sub = any_g.isel(lat=slice(j0, j0 + NR_10), lon=slice(i0, i0 + NC_10)).values

            frac_sub = np.nan_to_num(frac_sub, nan=0.0).astype("float32")
            frac_sub = np.clip(frac_sub, 0.0, 1.0)

            # ensure binary presence is exactly 0/1 uint8
            any_sub = np.nan_to_num(any_sub, nan=0).astype("uint8")
            any_sub = np.clip(any_sub, 0, 1).astype("uint8")

            # Build soil-style dataset
            sub = xr.Dataset(
                data_vars={
                    frac_var_out: (("N_lat", "N_lon"), frac_sub),
                    any_var_out:  (("N_lat", "N_lon"), any_sub),
                },
                coords={
                    "latitude":  ("N_lat", lat_sub),
                    "longitude": ("N_lon", lon_sub),
                },
            )

            # Global attributes (soil-style)
            sub.attrs["N_lon_global"] = N_LON_GLOBAL
            sub.attrs["N_lat_global"] = N_LAT_GLOBAL
            sub.attrs["i_ind_offset_LL"] = i0 + 1   # 1-based like other products
            sub.attrs["j_ind_offset_LL"] = j0 + 1
            sub.attrs["CellSize_arc_Secs"] = 30.0
            sub.attrs["Region"] = f"{region_label} tile H{hh}V{vv}"
            sub.attrs["CreatedBy"] = created_by
            sub.attrs["Date"] = created
            sub.attrs["Source"] = source_attr

            # No fill values (data are complete everywhere after nan_to_num)
            encoding = {
                frac_var_out: {"zlib": True, "complevel": 4, "dtype": "float32", "_FillValue": None},
                any_var_out:  {"zlib": True, "complevel": 4, "dtype": "uint8",   "_FillValue": None},
                "latitude":   {"zlib": True, "complevel": 1, "dtype": "float32", "_FillValue": None},
                "longitude":  {"zlib": True, "complevel": 1, "dtype": "float32", "_FillValue": None},
            }

            sub.to_netcdf(out, format="NETCDF4", engine="netcdf4", encoding=encoding)
            print("Wrote", out)

    ds.close()


if __name__ == "__main__":
    main()

