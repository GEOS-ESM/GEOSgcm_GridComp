#!/usr/bin/env python3

"""
Build GPM 2.0 based peat rasters on the standard 30-arcsec global raster grid.

Supported products
------------------
conservative
    GPM 2.0 class 1 -> 1.0
    GPM 2.0 class 2 -> 0.0
    else            -> 0.0

alpha050
    GPM 2.0 class 1 -> 1.0
    GPM 2.0 class 2 -> 0.5
    else            -> 0.0

liberal
    GPM 2.0 class 1 -> 1.0
    GPM 2.0 class 2 -> 1.0
    else            -> 0.0

The default archived product is alpha050:
    peatGPA22WGS_2cl_real_30arcsec.nc4

The script reads the original GPM 2.0 GeoTIFF directly, remaps the data from
the original 0.01-degree lat/lon grid to the standard global 30-arcsec raster
grid, converts GPM classes to real values, and writes the result as NetCDF4.

The 30-arcsec target grid is constructed analytically in double precision:
    lat = -90 + dlat/2 ... 90 - dlat/2, dlat = 180 / 21600
    lon = -180 + dlon/2 ... 180 - dlon/2, dlon = 360 / 43200

No legacy PEATMAP_mask.nc4 file is required.
"""

from pathlib import Path
import numpy as np
import xarray as xr
import rasterio


# ============================================================
# USER SETTINGS
# ============================================================

GPA_TIF_FILE = Path(
    "/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/land/soil/SOIL-DATA/v2/peatGPA22WGS_2cl.tif"
)

# Default product to build/archive.
# The alpha050 product contains enough information to derive conservative or
# liberal masks later if needed.
PRODUCT_MODE = "alpha050"


OUTFILE_ALPHA050      = Path("peatGPA22WGS_2cl_real_30arcsec.nc4")
# Optional diagnostic products if PRODUCT_MODE is changed manually.
OUTFILE_CONSERVATIVE = Path("peatGPA22WGS_2cl_real_30arcsec_conservative.nc4")
OUTFILE_LIBERAL      = Path("peatGPA22WGS_2cl_real_30arcsec_liberal.nc4")

NLAT = 21600
NLON = 43200

GRID_TOL = 1.0e-10
NC_FILL  = np.float32(-9999.0)


# ============================================================
# HELPERS
# ============================================================

def arrays_match(a, b, tol=GRID_TOL):
    """Return True if two 1D coordinate arrays match within tolerance."""
    return (len(a) == len(b)) and np.allclose(a, b, rtol=0.0, atol=tol)


def get_resolution_arcsec(coord):
    """Return absolute grid spacing in arc-seconds for a 1D coordinate."""
    return abs(float(coord[1] - coord[0])) * 3600.0

def load_gpa_from_tif(tif_path):
    """
    Read GPM 2.0 (GPA22) directly from the GeoTIFF and return
    lon, lat, peatland_type(lat, lon).

    Latitude is converted to ascending order because xarray interpolation
    expects monotonic coordinates.
    """
    print(f"Reading GPM 2.0 GeoTIFF: {tif_path}")

    with rasterio.open(tif_path) as src:
        data = src.read(1)
        transform = src.transform
        width = src.width
        height = src.height

        col = np.arange(width, dtype=np.float64)
        row = np.arange(height, dtype=np.float64)

        lon = transform.c + (col + 0.5) * transform.a
        lat = transform.f + (row + 0.5) * transform.e

    if lat[1] < lat[0]:
        print("Sorting GPM 2.0 latitude to ascending order")
        lat = lat[::-1].copy()
        data = data[::-1, :]

    ptype = xr.DataArray(
        data,
        dims=("lat", "lon"),
        coords={"lat": lat.astype(np.float64), "lon": lon.astype(np.float64)},
        name="peatland_type",
    )

    return lon.astype(np.float64), lat.astype(np.float64), ptype

def build_target_grid():
    """
    Construct the standard global 30-arcsec raster cell-center grid
    in double precision.

    Latitude is ascending: south to north.
    Longitude is ascending: west to east.
    """
    dlat = np.float64(180.0) / np.float64(NLAT)
    dlon = np.float64(360.0) / np.float64(NLON)

    lat = -90.0  + dlat * (np.arange(NLAT, dtype=np.float64) + 0.5)
    lon = -180.0 + dlon * (np.arange(NLON, dtype=np.float64) + 0.5)

    return lon, lat

def build_mask(ptype, mode):
    """Build the requested peat field from GPM 2.0 peatland_type."""
    if mode == "conservative":
        field = xr.where(ptype == 1, 1.0, 0.0)
    elif mode == "alpha050":
        field = xr.zeros_like(ptype, dtype="float32")
        field = xr.where(ptype == 1, 1.0, field)
        field = xr.where(ptype == 2, 0.5, field)
    elif mode == "liberal":
        field = xr.where(ptype.isin([1, 2]), 1.0, 0.0)
    else:
        raise RuntimeError(f"Unsupported mode: {mode}")

    return field.astype("float32")


def remap_to_target(field, lon_src, lat_src, lon_tgt, lat_tgt):
    """
    Remap the GPM 2.0 field to the target PEATMAP grid using nearest neighbor.
    If the source and target grids already match exactly, remapping is skipped.
    """
    same_lon = arrays_match(lon_src, lon_tgt)
    same_lat = arrays_match(lat_src, lat_tgt)

    print(f"  Source grid: nlat={len(lat_src)} nlon={len(lon_src)}")
    print(f"  Source resolution: dlat={get_resolution_arcsec(lat_src):.6f}\" "
          f"dlon={get_resolution_arcsec(lon_src):.6f}\"")
    print(f"  Target grid: nlat={len(lat_tgt)} nlon={len(lon_tgt)}")
    print(f"  Target resolution: dlat={get_resolution_arcsec(lat_tgt):.6f}\" "
          f"dlon={get_resolution_arcsec(lon_tgt):.6f}\"")
    print(f"  Longitude arrays match: {same_lon}")
    print(f"  Latitude arrays match:  {same_lat}")

    if same_lon and same_lat:
        print("  Grids match. No remapping needed.")
        return field.transpose("lat", "lon")

    print("  Grids do not match. Remapping to standard 30-arcsec grid using nearest neighbor.")
    lon_target = xr.DataArray(lon_tgt, dims=("lon",))
    lat_target = xr.DataArray(lat_tgt, dims=("lat",))

    out = field.interp(
        lon=lon_target,
        lat=lat_target,
        method="nearest"
    ).transpose("lat", "lon")
    
    return out.fillna(0.0).astype("float32")

def write_output(outfile, lon, lat, peat_field, mode):
    """
    Write output on the standard 30-arcsec global raster grid:
      longitude(N_lon), latitude(N_lat), PEATMAP(N_lat, N_lon)
    """    
    print(f"Writing output: {outfile}")

    peat_data = peat_field.rename({"lat": "N_lat", "lon": "N_lon"})

    ds_out = xr.Dataset(
        data_vars={
            "longitude": (("N_lon",), lon.astype("float64")),
            "latitude":  (("N_lat",), lat.astype("float64")),
            "PEATMAP":   (("N_lat", "N_lon"), peat_data.data.astype("float32")),
        },
        coords={}
    )

    ds_out["longitude"].attrs["units"] = "degrees_east"
    ds_out["latitude" ].attrs["units"] = "degrees_north"
    ds_out["PEATMAP"  ].attrs["units"] = "1"

    if mode == "conservative":
        ds_out["PEATMAP"].attrs["long_name"] = "Conservative peat mask from GPM 2.0"
        ds_out["PEATMAP"].attrs["description"] = (
            "PEATMAP = 1 where GPM 2.0 peatland_type == 1 (peat dominated); "
            "PEATMAP = 0 otherwise, including GPM 2.0 class 2."
        )
        ds_out.attrs["note"] = (
            "Conservative binary peat mask: only GPM 2.0 class 1 is treated as peat; "
            "class 2 is treated as non-peat."
        )
    elif mode == "alpha050":
        ds_out["PEATMAP"].attrs["long_name"] = "Peat fraction from GPM 2.0"
        ds_out["PEATMAP"].attrs["description"] = (
            "Fractional peat cover derived from GPM 2.0: 1.0 for class 1 "
            "(peat dominated), 0.5 for class 2 (peat in soil mosaic), "
            "0.0 otherwise."
        )
        ds_out.attrs["note"] = (
            "Alpha peat mask: GPM 2.0 class 1 is assigned a value of 1 and class 2 is assigned a value of 0.5."
        )
    elif mode == "liberal":
        ds_out["PEATMAP"].attrs["long_name"] = "Peat mask from GPM 2.0 including dominant and mosaic peat"
        ds_out["PEATMAP"].attrs["description"] = (
            "PEATMAP = 1 where GPM 2.0 peatland_type is 1 or 2; "
            "PEATMAP = 0 otherwise."
        )
        ds_out.attrs["note"] = (
            "Liberal binary peat mask: GPM 2.0 classes 1 and 2 are both treated as peat."
        )
    else:
        raise RuntimeError(f"Unsupported mode: {mode}")

    ds_out.attrs["source"] = (
        "Derived from GPM 2.0 (Global Peatland Map 2.0) GeoTIFF. "
        "The original 0.01-degree lat/lon data were remapped to the standard "
        "global 30-arcsec raster grid and written as NetCDF4."
    )    
    ds_out.attrs["note"] = (
        "Alpha peat fraction: GPM 2.0 class 1 is assigned 1.0, "
        "class 2 is assigned 0.5, and all other classes are assigned 0.0. "
        "This is the default archived product."
    )    

    encoding = {
        "PEATMAP":   {"_FillValue": NC_FILL, "zlib": True, "complevel": 1},
        "longitude": {"_FillValue": None},
        "latitude":  {"_FillValue": None},
    }

    Path(outfile).unlink(missing_ok=True)
    ds_out.to_netcdf(outfile, format="NETCDF4", encoding=encoding)
    ds_out.close()


def get_requested_products(product_mode):
    """Return the ordered list of products to create."""
    if product_mode == "conservative":
        return ["conservative"]
    if product_mode == "alpha050":
        return ["alpha050"]
    if product_mode == "liberal":
        return ["liberal"]
    if product_mode == "all":
        return ["conservative", "alpha050", "liberal"]
    raise RuntimeError(f"Unsupported PRODUCT_MODE: {product_mode}")


def get_outfile(mode):
    """Return the output filename for the chosen mode."""
    if mode == "conservative":
        return OUTFILE_CONSERVATIVE
    if mode == "alpha050":
        return OUTFILE_ALPHA050
    if mode == "liberal":
        return OUTFILE_LIBERAL
    raise RuntimeError(f"Unsupported mode: {mode}")


# ============================================================
# MAIN
# ============================================================

def main():
    lon_gpa, lat_gpa, ptype = load_gpa_from_tif(GPA_TIF_FILE)

    print(f"GPM 2.0 grid: nlat={len(lat_gpa)} nlon={len(lon_gpa)}")
    print(f"GPM 2.0 resolution: dlat={get_resolution_arcsec(lat_gpa):.6f}\" "
          f"dlon={get_resolution_arcsec(lon_gpa):.6f}\"")

    lon_tgt, lat_tgt = build_target_grid()

    print(f"Target 30-arcsec grid: nlat={len(lat_tgt)} nlon={len(lon_tgt)}")
    print(f"Target resolution: dlat={get_resolution_arcsec(lat_tgt):.6f}\" "
          f"dlon={get_resolution_arcsec(lon_tgt):.6f}\"")

    for mode in get_requested_products(PRODUCT_MODE):
        print()
        print("=" * 72)
        print(f"Building product: {mode}")
        print("=" * 72)

        field = build_mask(ptype, mode)
        field_on_target = remap_to_target(field, lon_gpa, lat_gpa, lon_tgt, lat_tgt)
        write_output(get_outfile(mode), lon_tgt, lat_tgt, field_on_target, mode)

    print()
    print("Done.")

if __name__ == "__main__":
    main()

