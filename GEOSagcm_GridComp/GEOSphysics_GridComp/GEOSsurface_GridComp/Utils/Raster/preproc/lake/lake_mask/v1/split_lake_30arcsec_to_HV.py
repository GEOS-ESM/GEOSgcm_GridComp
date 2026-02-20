#!/usr/bin/env python3
import os
import numpy as np
import xarray as xr
from datetime import datetime

"""
HydroLAKES-TopoCat 30 arc-sec lake mask split into 10°x10° tiles for GEOS make_bcs.

Format intentionally matches soil-properties style:
  - Dimensions: N_lat, N_lon (1200x1200)
  - 1D coordinates: latitude(N_lat), longitude(N_lon)
  - Variables:
      lake_presence_frac(N_lat,N_lon)  float32 [0,1]
      lake_presence_any(N_lat,N_lon)   uint8 0/1
  - Global attributes include:
      N_lon_global, N_lat_global
      i_ind_offset_LL, j_ind_offset_LL
      CellSize_arc_Secs

No missing values are present anywhere.
All ocean cells are 0.
This file is designed for fast chunk-based reading inside GEOS tiling (similar to soil and snow albedo).
"""

in_global = "/discover/nobackup/borescan/brisi/lake/processing/make_bcs_preproc/LakeTopoCat_Global_30arcsec.nc4"
out_dir   = "/discover/nobackup/borescan/brisi/lake/processing/make_bcs_preproc/tiles_10deg/"
os.makedirs(out_dir, exist_ok=True)

# 30 arc-sec global grid
N_LON_GLOBAL = 43200
N_LAT_GLOBAL = 21600

# 10° x 10° tile size at 30"
nc_10 = 1200
nr_10 = 1200

ds = xr.open_dataset(in_global)

# Get global 1-D coords
lat_g = ds["lat"].astype("float32")
lon_g = ds["lon"].astype("float32")

# Data vars
frac_g = ds["lake_presence_frac"].astype("float32")
any_g  = ds["lake_presence_any"].astype("uint8")

assert frac_g.sizes["lon"] == N_LON_GLOBAL and frac_g.sizes["lat"] == N_LAT_GLOBAL

created = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

for jx in range(1, 18+1):      # V 01..18
    for ix in range(1, 36+1):  # H 01..36
        i0 = (ix - 1) * nc_10   # 0-based lon start
        j0 = (jx - 1) * nr_10   # 0-based lat start

        hh = f"{ix:02d}"
        vv = f"{jx:02d}"
        out = os.path.join(out_dir, f"LakeTopoCat_30arcsec_H{hh}V{vv}.nc4")

        # Subset coords (1-D)
        lat_sub = lat_g.isel(lat=slice(j0, j0 + nr_10)).values
        lon_sub = lon_g.isel(lon=slice(i0, i0 + nc_10)).values

        # Subset data (2-D), enforce no-missing + valid ranges
        frac_sub = frac_g.isel(lat=slice(j0, j0 + nr_10), lon=slice(i0, i0 + nc_10)).values
        any_sub  = any_g.isel (lat=slice(j0, j0 + nr_10), lon=slice(i0, i0 + nc_10)).values

        frac_sub = np.nan_to_num(frac_sub, nan=0.0).astype("float32")
        frac_sub = np.clip(frac_sub, 0.0, 1.0)

        # ensure any is exactly 0/1 uint8
        any_sub = np.nan_to_num(any_sub, nan=0).astype("uint8")
        any_sub = np.clip(any_sub, 0, 1).astype("uint8")

        # Build soil-style dataset
        sub = xr.Dataset(
            data_vars={
                "lake_presence_frac": (("N_lat", "N_lon"), frac_sub),
                "lake_presence_any":  (("N_lat", "N_lon"), any_sub),
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
        sub.attrs["Region"] = f"LakeTopoCat tile H{hh}V{vv}"
        sub.attrs["CreatedBy"] = "borescan"
        sub.attrs["Date"] = created
        sub.attrs["Source"] = "HydroLAKES-TopoCat v1.1 (2023); rasterized at 10\" and aggregated to 30\"."

        # No fill values (data are complete everywhere)
        encoding = {
            "lake_presence_frac": {"zlib": True, "complevel": 4, "dtype": "float32", "_FillValue": None},
            "lake_presence_any":  {"zlib": True, "complevel": 4, "dtype": "uint8",   "_FillValue": None},
            "latitude":           {"zlib": True, "complevel": 1, "dtype": "float32", "_FillValue": None},
            "longitude":          {"zlib": True, "complevel": 1, "dtype": "float32", "_FillValue": None},
        }

        sub.to_netcdf(out, format="NETCDF4", engine="netcdf4", encoding=encoding)
        print("Wrote", out)
