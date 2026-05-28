"""
Create HydroLAKES-TopoCat reach mask for GEOS/BCS tiling from Lake-TopoCat Reaches.

Output:
  ReachTopoCat_Global_30arcsec.nc4

Variables:
  - reach_occupancy_frac(lat, lon): float32 in [0, 1]
        Fraction of 10 arc-sec subpixels touched by any reach geometry within
        the 30 arc-sec cell. This is a line-occupancy metric, NOT true channel area.
  - reach_presence_any(lat, lon): uint8 0/1
        1 if any 10 arc-sec subpixel in the 30 arc-sec cell is touched by a reach.

Method:
  1) Build 10 arc-sec binary raster in memory for a latitude band (5 degrees).
  2) Rasterize reach lines with all_touched=True to preserve thin features.
  3) Aggregate 10" -> 30" in that band (3x3 mean/max).
  4) Write 30" band NetCDF.
  5) Combine all 30" bands into one global NetCDF.

Notes:
   This script is intentionally separate from the lake preprocessor because
   reaches are line geometries, not polygons.
   reach_occupancy_frac is useful as a density-style metric, but it should not
   be interpreted as water area fraction.
"""

import os
import glob
import numpy as np
import geopandas as gpd
import rasterio.features
import xarray as xr
from affine import Affine

# ------------------------------------------------------------------
# Input reach line dataset
#
# Source: HydroLAKES-TopoCat v1.1 (2023)
# Files:  Reaches_pfaf_*.shp
#
# Usage in this script:
#   Only line geometry is used.
#   All reaches are rasterized as binary reach presence (value = 1).
# ------------------------------------------------------------------
# Paths (edit as needed)
# ------------------------------------------------------------------
shapefile_dir = "/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/lake/Lake_TopoCat/v1/Reaches/"
out_dir       = "/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/lake/Lake_TopoCat/v1/reach_mask_build/"
out_global    = "/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/lake/reach_mask/v1/before_splitting_nc4/ReachTopoCat_Global_30arcsec.nc4"

os.makedirs(out_dir, exist_ok=True)
os.makedirs(os.path.dirname(out_global), exist_ok=True)

# ------------------------------------------------------------------
# Grid definitions
# ------------------------------------------------------------------
# Oversampling grid: 10 arc-sec
res10 = 10.0 / 3600.0
lon10 = np.arange(-180.0 + res10 / 2.0, 180.0, res10)  # 129600
lat10 = np.arange(-90.0 + res10 / 2.0,   90.0, res10)  # 64800

# Target grid: 30 arc-sec = 3 x 10 arc-sec
lon30 = lon10[1::3]  # 43200
# lat30 varies by band; global lat30 would be lat10[1::3] = 21600

# Process in latitude bands to limit memory
band_deg = 5.0

# Preserve thin line features when rasterizing
ALL_TOUCHED = True

# ------------------------------------------------------------------
# Shapefile list
# ------------------------------------------------------------------
shp_list = sorted(glob.glob(os.path.join(shapefile_dir, "Reaches_pfaf_*.shp")))
if not shp_list:
    raise RuntimeError(f"No shapefiles found: {shapefile_dir}/Reaches_pfaf_*.shp")
print(f"Found {len(shp_list)} reach shapefiles")


# ------------------------------------------------------------------
# Helper: 10" -> 30" aggregation (3x3 blocks)
# ------------------------------------------------------------------
def agg10_to_30(binary10: np.ndarray):
    """
    binary10: uint8 array (nlat10, nlon10) with values 0/1
    Returns:
      frac30: float32 (nlat30, nlon30) = mean of 3x3 block
              interpreted here as subpixel occupancy by reach lines
      any30 : uint8   (nlat30, nlon30) = max of 3x3 block
    """
    nlat, nlon = binary10.shape
    nlat3 = (nlat // 3) * 3
    nlon3 = (nlon // 3) * 3
    binary10 = binary10[:nlat3, :nlon3]

    blk = binary10.reshape(nlat3 // 3, 3, nlon3 // 3, 3)
    frac30 = blk.mean(axis=(1, 3)).astype(np.float32)
    any30 = blk.max(axis=(1, 3)).astype(np.uint8)
    return frac30, any30


# ------------------------------------------------------------------
# Build 30" band files
# ------------------------------------------------------------------
band_files = []

for lat_start in np.arange(-90.0, 90.0, band_deg):
    lat_end = lat_start + band_deg
    print(f"\n=== Latitude band {lat_start} to {lat_end} ===")

    # 10" lat centers in this band
    lat10_band = lat10[(lat10 >= lat_start) & (lat10 < lat_end)]
    if lat10_band.size == 0:
        continue

    out_shape10 = (lat10_band.size, lon10.size)

    # Affine transform for this band (pixel -> lon/lat)
    # Keep consistent with the existing lake preprocessor.
    transform10 = (
        Affine.translation(lon10[0] - res10 / 2.0, lat10_band[0] - res10 / 2.0)
        * Affine.scale(res10, res10)
    )

    # In-memory 10" band raster (binary 0/1)
    reach10 = np.zeros(out_shape10, dtype=np.uint8)

    # Rasterize all shapefiles, but only features intersecting this latitude band
    for shp in shp_list:
        gdf = gpd.read_file(shp).to_crs("EPSG:4326")
        gdf_band = gdf.cx[:, lat_start:lat_end]
        if gdf_band.empty:
            continue

        # Drop null / empty geometries
        gdf_band = gdf_band[gdf_band.geometry.notnull()]
        gdf_band = gdf_band[~gdf_band.geometry.is_empty]
        if gdf_band.empty:
            continue

        shapes = ((geom, 1) for geom in gdf_band.geometry)

        tmp = rasterio.features.rasterize(
            shapes,
            out_shape=out_shape10,
            transform=transform10,
            fill=0,
            dtype=np.uint8,
            all_touched=ALL_TOUCHED,
        )

        # Union across shapefiles
        reach10 = np.maximum(reach10, tmp)

    # Aggregate to 30"
    frac30, any30 = agg10_to_30(reach10)

    # 30" lat centers for this band (aligned with 10" centers)
    lat30_band = lat10_band[1::3]

    if frac30.shape != (lat30_band.size, lon30.size):
        raise RuntimeError(
            f"Shape mismatch: frac30={frac30.shape} vs coords={(lat30_band.size, lon30.size)}"
        )

    ds = xr.Dataset(
        {
            "reach_occupancy_frac": (("lat", "lon"), frac30),
            "reach_presence_any": (("lat", "lon"), any30),
        },
        coords={"lat": lat30_band, "lon": lon30},
    )

    ds.attrs["Source"] = "HydroLAKES-TopoCat v1.1 (2023) Reaches"
    ds.attrs["GeometryType"] = "line"
    ds.attrs["Rasterization"] = "10 arc-sec oversampling, aggregated to 30 arc-sec"
    ds.attrs["AllTouched"] = int(ALL_TOUCHED)
    ds.attrs["reach_occupancy_frac_note"] = (
        "Fraction of 10 arc-sec subpixels touched by reach lines; not true river area fraction."
    )

    band_out = os.path.join(out_dir, f"ReachTopoCat_30arcsec_{lat_start}_{lat_end}.nc4")
    encoding = {
        "reach_occupancy_frac": {
            "zlib": True,
            "complevel": 4,
            "dtype": "float32",
            "_FillValue": -9999.0,
        },
        "reach_presence_any": {
            "zlib": True,
            "complevel": 4,
            "dtype": "uint8",
            "_FillValue": 255,
        },
    }

    ds.to_netcdf(band_out, format="NETCDF4", engine="netcdf4", encoding=encoding)
    print("Wrote band:", band_out)
    band_files.append(band_out)

# ------------------------------------------------------------------
# Combine band files into the single global 30" product
# ------------------------------------------------------------------
print(f"\nCombining {len(band_files)} 30\" band files into:\n  {out_global}")

if not band_files:
    raise RuntimeError("No band files were created. Check input shapefiles and paths.")

ds = xr.open_mfdataset(
    sorted(band_files),
    combine="by_coords",
    parallel=False,
).sortby("lat")

# Force arrays into memory to avoid dask/lazy-write issues
frac = ds["reach_occupancy_frac"].load().astype("float32")
anyv = ds["reach_presence_any"].load().astype("uint8")

ds_out = xr.Dataset(
    {
        "reach_occupancy_frac": frac,
        "reach_presence_any": anyv,
    },
    coords={
        "lat": ds["lat"].values,
        "lon": ds["lon"].values,
    },
)

ds_out.attrs["Title"] = "HydroLAKES-TopoCat reaches rasterized to 30 arc-sec"
ds_out.attrs["Source"] = "HydroLAKES-TopoCat v1.1 (2023) Reaches"
ds_out.attrs["GeometryType"] = "line"
ds_out.attrs["AllTouched"] = int(ALL_TOUCHED)
ds_out.attrs["reach_occupancy_frac_note"] = (
    "Fraction of 10 arc-sec subpixels touched by reach lines; not true river area fraction."
)

encoding = {
    "reach_occupancy_frac": {
        "zlib": True,
        "complevel": 4,
        "dtype": "float32",
        "_FillValue": -9999.0,
    },
    "reach_presence_any": {
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

