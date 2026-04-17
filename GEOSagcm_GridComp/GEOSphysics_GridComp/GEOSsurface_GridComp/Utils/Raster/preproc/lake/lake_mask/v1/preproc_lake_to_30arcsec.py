#!/usr/bin/env python3
"""
Create HydroLAKES-based lake mask for GEOS/BCS tiling from Lake-TopoCat (HydroLAKES-TopoCat v1.1 2023):

Output (the only file makebcs Fortran needs):
  LakeTopoCat_Global_30arcsec.nc4

Variables:
  - lake_area_frac(lat, lon): float32 in [0, 1]
        Fractional lake coverage of each 30 arc-sec cell, estimated by oversampling
        at 10 arc-sec and averaging 3x3 subpixels.
  - lake_presence_any(lat, lon): uint8 0/1 (optional but cheap and handy)
        1 if any 10" subpixel in the 30" cell is lake.

Method:
  1) Build 10 arc-sec binary raster in memory for a latitude band (5 degrees).
  2) Immediately aggregate 10" -> 30" in that band (3x3 mean/max).
  3) Write 30" band NetCDF.
  4) Combine all 30" bands into one global NetCDF.

Why oversample at 10"?
  Rasterization at 30" directly gives only a binary pixel mask. Oversampling at 10"
  lets us estimate fractional coverage in each 30" pixel at low cost.
"""

import os
import glob
import numpy as np
import geopandas as gpd
import rasterio.features
import xarray as xr
from affine import Affine

# ------------------------------------------------------------------
# Input lake polygon dataset
#
# Source: HydroLAKES-TopoCat v1.1 Y2023 processing)
# Files:  Lakes_pfaf_*.shp
#
# Description:
#   Global lake polygon dataset (EPSG:4326) partitioned by Pfafstetter
#   basin code. Each shapefile contains HydroLAKES polygons with attributes
#   such as lake area, lake type, and permanence.
#
# Usage in this script:
#   Only polygon geometry is used.
#   All polygons are rasterized as lake presence (value = 1).
#   No filtering by lake type, size, or permanence is applied.
# ------------------------------------------------------------
# Paths (edit as needed)
# ------------------------------------------------------------------
shapefile_dir = "/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/lake/Lake_TopoCat/v1/Lakes/"
out_dir = "/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/lake/Lake_TopoCat/v1/lake_mask_build/"
out_global = "/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/lake/lake_mask/v1/before_splitting_nc4/LakeTopoCat_Global_30arcsec.nc4"


os.makedirs(out_dir, exist_ok=True)

# ------------------------------------------------------------
# Grid definitions
# ------------------------------------------------------------
# Oversampling grid: 10 arc-sec
res10 = 10 / 3600.0
lon10 = np.arange(-180 + res10/2, 180, res10)  # 129600
lat10 = np.arange(-90  + res10/2,  90, res10)  # 64800

# Target grid: 30 arc-sec = 3 x 10 arc-sec
lon30 = lon10[1::3]         # 43200
# lat30 varies by band; global lat30 would be lat10[1::3] = 21600

# Process in latitude bands to limit memory
band_deg = 5.0

# ------------------------------------------------------------
# Shapefile list
# ------------------------------------------------------------
shp_list = sorted(glob.glob(os.path.join(shapefile_dir, "Lakes_pfaf_*.shp")))
if not shp_list:
    raise RuntimeError(f"No shapefiles found: {shapefile_dir}/Lakes_pfaf_*.shp")
print(f"Found {len(shp_list)} shapefiles")

# ------------------------------------------------------------
# Helper: 10" -> 30" aggregation (3x3 blocks)
# ------------------------------------------------------------
def agg10_to_30(binary10: np.ndarray):
    """
    binary10: uint8 array (nlat10, nlon10) with values 0/1
    Returns:
      frac30: float32 (nlat30, nlon30) = mean of 3x3 block
      any30 : uint8   (nlat30, nlon30) = max  of 3x3 block
    """
    nlat, nlon = binary10.shape
    nlat3 = (nlat // 3) * 3
    nlon3 = (nlon // 3) * 3
    binary10 = binary10[:nlat3, :nlon3]

    # reshape into 3x3 blocks
    blk = binary10.reshape(nlat3//3, 3, nlon3//3, 3)
    frac30 = blk.mean(axis=(1, 3)).astype(np.float32)
    any30  = blk.max(axis=(1, 3)).astype(np.uint8)
    return frac30, any30

# ------------------------------------------------------------
# Build 30" band files (each band is small enough to write fast)
# ------------------------------------------------------------
band_files = []

for lat_start in np.arange(-90, 90, band_deg):
    lat_end = lat_start + band_deg
    print(f"\n=== Latitude band {lat_start} to {lat_end} ===")

    # 10" lat centers in this band
    lat10_band = lat10[(lat10 >= lat_start) & (lat10 < lat_end)]
    if lat10_band.size == 0:
        continue

    out_shape10 = (lat10_band.size, lon10.size)

    # Affine transform for this band (pixel -> lon/lat)
    transform10 = (
        Affine.translation(lon10[0] - res10/2, lat10_band[0] - res10/2)
        * Affine.scale(res10, res10)
    )

    # In-memory 10" band raster (binary 0/1)
    lake10 = np.zeros(out_shape10, dtype=np.uint8)

    # Rasterize all shapefiles, but only features intersecting this lat band
    for shp in shp_list:
        gdf = gpd.read_file(shp).to_crs("EPSG:4326")
        gdf_band = gdf.cx[:, lat_start:lat_end]
        if gdf_band.empty:
            continue

        shapes = ((geom, 1) for geom in gdf_band.geometry)

        tmp = rasterio.features.rasterize(
            shapes,
            out_shape=out_shape10,
            transform=transform10,
            fill=0,
            dtype=np.uint8,
        )

        # Union across shapefiles
        lake10 = np.maximum(lake10, tmp)

    # Aggregate to 30" for this band
    frac30, any30 = agg10_to_30(lake10)

    # 30" lat centers for this band (aligned with 10" centers)
    lat30_band = lat10_band[1::3]

    # Sanity check dimensions
    if frac30.shape != (lat30_band.size, lon30.size):
        raise RuntimeError(f"Shape mismatch: frac30={frac30.shape} vs coords={(lat30_band.size, lon30.size)}")

    # Write band netcdf
    ds = xr.Dataset(
        {
            "lake_area_frac": (("lat", "lon"), frac30),
            "lake_presence_any":  (("lat", "lon"), any30),  # optional but useful
        },
        coords={"lat": lat30_band, "lon": lon30},
    )

    band_out = os.path.join(out_dir, f"LakeTopoCat_30arcsec_{lat_start}_{lat_end}.nc4")
    encoding = {
        "lake_area_frac": {"zlib": True, "complevel": 4, "dtype": "float32", "_FillValue": -9999.0},
        "lake_presence_any":  {"zlib": True, "complevel": 4, "dtype": "uint8",   "_FillValue": 255},
    }
    
    ds.to_netcdf(band_out, format="NETCDF4", engine="netcdf4", encoding=encoding)
    print("Wrote band:", band_out)
    band_files.append(band_out)

# ------------------------------------------------------------
# Combine band files into the single global 30" product
# ------------------------------------------------------------
print(f"\nCombining {len(band_files)} 30\" band files into:\n  {out_global}")

ds = xr.open_mfdataset(
    sorted(band_files),
    combine="by_coords",
    parallel=False
).sortby("lat")

# Force arrays into memory to avoid dask/lazy-write issues
frac = ds["lake_area_frac"].load().astype("float32")
anyv = ds["lake_presence_any"].load().astype("uint8")

ds_out = xr.Dataset(
    {
        "lake_area_frac": frac,
        "lake_presence_any": anyv,
    },
    coords={
        "lat": ds["lat"].values,
        "lon": ds["lon"].values,
    },
)

encoding = {
    "lake_area_frac": {
        "zlib": True,
        "complevel": 4,
        "dtype": "float32",
        "_FillValue": -9999.0,
    },
    "lake_presence_any": {
        "zlib": True,
        "complevel": 4,
        "dtype": "uint8",
        "_FillValue": 255,
    },
}

ds_out.to_netcdf(
    out_global,
    format="NETCDF4",
    engine="netcdf4",
    encoding=encoding,
)

print("Done. Global file written:")
print(" ", out_global)

