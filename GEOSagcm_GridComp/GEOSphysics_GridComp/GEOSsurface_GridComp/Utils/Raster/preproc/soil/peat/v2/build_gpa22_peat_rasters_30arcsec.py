#!/usr/bin/env python3

"""
Build GPA22-based peat rasters on the legacy PEATMAP target grid.

Supported products
------------------
conservative
    GPA22 class 1 -> 1.0
    GPA22 class 2 -> 0.0
    else          -> 0.0

alpha050
    GPA22 class 1 -> 1.0
    GPA22 class 2 -> 0.5
    else          -> 0.0

liberal
    GPA22 class 1 -> 1.0
    GPA22 class 2 -> 1.0
    else          -> 0.0

Why this uses the legacy PEATMAP grid
-------------------------------------
The output grid is taken directly from the legacy PEATMAP file so the resulting
rasters can be used as drop-in replacements in the existing Fortran/PEATMAP
workflow, with the same dimensions, coordinate arrays, variable names, and
overall file structure expected by the current code.

This is a compatibility choice. The script does not construct a new analytic
global grid with exact 30.000000-arcsec spacing. Instead, it remaps GPA22 onto
the exact longitude/latitude arrays stored in the legacy PEATMAP product.

In practice, this means the outputs are PEATMAP-grid-compatible global rasters
with nominal ~30-arcsec resolution (21600 x 43200), rather than standalone
exact-30-arcsec products.

Because the GPA22 source grid differs from the PEATMAP target grid, the GPA22
field is remapped to the PEATMAP grid using nearest-neighbor interpolation.
"""

from pathlib import Path
import numpy as np
import xarray as xr

try:
    import rasterio
except ImportError:
    rasterio = None


# ============================================================
# USER SETTINGS
# ============================================================

# Source choice:
#   "nc"  -> read the existing GPA22 NetCDF file
#   "tif" -> read the original GeoTIFF directly
# For strict reproduction of your current files, keep this as "nc" first.
GPA_SOURCE_KIND = "nc"

GPA_NC_FILE = Path(
    "/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/land/soil/SOIL-DATA/v2/peatGPA22WGS_2cl.nc4"
)
GPA_TIF_FILE = Path(
    "/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/land/soil/SOIL-DATA/v2/peatGPA22WGS_2cl.tif"
)

# Target 30-arcsec grid used by the legacy PEATMAP product
PEATMAP_GRID_FILE = Path(
    "/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/"
    "land/soil/SOIL-DATA/PEATMAP_mask.nc4"
)

# Which products to create:
#   "conservative" -> class 1 only
#   "alpha050"     -> class 1 = 1.0, class 2 = 0.5
#   "liberal"      -> class 1 or 2 both treated as 1.0
#   "all"          -> write all three files
PRODUCT_MODE = "all"

OUTFILE_CONSERVATIVE = Path("PEATMAP_from_GPA22_like_old_conservative.nc4")
OUTFILE_ALPHA050 = Path("PEATMAP_from_GPA22_like_old_alpha0.50.nc4")
OUTFILE_LIBERAL = Path("PEATMAP_from_GPA22_like_old_liberal.nc4")

GRID_TOL = 1.0e-10
NC_FILL = np.float32(-9999.0)


# ============================================================
# HELPERS
# ============================================================

def arrays_match(a, b, tol=GRID_TOL):
    """Return True if two 1D coordinate arrays match within tolerance."""
    return (len(a) == len(b)) and np.allclose(a, b, rtol=0.0, atol=tol)


def get_resolution_arcsec(coord):
    """Return absolute grid spacing in arc-seconds for a 1D coordinate."""
    return abs(float(coord[1] - coord[0])) * 3600.0


def load_gpa_from_nc(nc_path):
    """
    Read GPA22 from the NetCDF file and return lon, lat, peatland_type(lat, lon).
    Latitude is forced to ascending order because interpolation assumes that.
    """
    print(f"Reading GPA22 NetCDF: {nc_path}")
    ds = xr.open_dataset(nc_path)

    if "lon" not in ds or "lat" not in ds or "peatland_type" not in ds:
        raise RuntimeError("Expected variables lon, lat, peatland_type in GPA22 NetCDF.")

    ptype = ds["peatland_type"]
    lon = np.asarray(ds["lon"].values, dtype=np.float64)
    lat = np.asarray(ds["lat"].values, dtype=np.float64)

    if lat[1] < lat[0]:
        print("Sorting GPA22 latitude to ascending order")
        ds = ds.sortby("lat")
        ptype = ds["peatland_type"]
        lon = np.asarray(ds["lon"].values, dtype=np.float64)
        lat = np.asarray(ds["lat"].values, dtype=np.float64)

    return ds, lon, lat, ptype


def load_gpa_from_tif(tif_path):
    """
    Read GPA22 directly from the GeoTIFF and return lon, lat, peatland_type(lat, lon).
    Latitude is converted to ascending order to match the NetCDF workflow.
    """
    if rasterio is None:
        raise RuntimeError("rasterio is required for GPA_SOURCE_KIND='tif'.")

    print(f"Reading GPA22 GeoTIFF: {tif_path}")
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
        print("Sorting GPA22 latitude to ascending order")
        lat = lat[::-1].copy()
        data = data[::-1, :]

    ptype = xr.DataArray(
        data,
        dims=("lat", "lon"),
        coords={"lat": lat, "lon": lon},
        name="peatland_type",
    )

    return None, lon.astype(np.float64), lat.astype(np.float64), ptype


def load_target_grid(peatmap_grid_file):
    """Read the target PEATMAP grid that defines the output 30-arcsec raster."""
    print(f"Reading PEATMAP target grid: {peatmap_grid_file}")
    ds = xr.open_dataset(peatmap_grid_file)

    if "longitude" not in ds or "latitude" not in ds:
        raise RuntimeError("Expected longitude and latitude in PEATMAP grid file.")

    lon = np.asarray(ds["longitude"].values, dtype=np.float64)
    lat = np.asarray(ds["latitude"].values, dtype=np.float64)

    if lat[1] < lat[0]:
        raise RuntimeError("Expected PEATMAP latitude to be ascending.")

    return ds, lon, lat


def build_mask(ptype, mode):
    """Build the requested peat field from GPA22 peatland_type."""
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
    Remap the GPA22 field to the target PEATMAP grid using nearest neighbor.
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

    print("  Grids do not match. Remapping to PEATMAP grid using nearest neighbor.")
    lon_target = xr.DataArray(lon_tgt, dims=("lon",))
    lat_target = xr.DataArray(lat_tgt, dims=("lat",))

    return field.interp(
        lon=lon_target,
        lat=lat_target,
        method="nearest"
    ).transpose("lat", "lon")


def write_output(outfile, lon, lat, peat_field, mode):
    """
    Write output in the PEATMAP-like layout:
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
    ds_out["latitude"].attrs["units"] = "degrees_north"
    ds_out["PEATMAP"].attrs["units"] = "1"

    if mode == "conservative":
        ds_out["PEATMAP"].attrs["long_name"] = "Conservative peat mask from GPA22"
        ds_out["PEATMAP"].attrs["description"] = (
            "PEATMAP = 1 where GPA22 peatland_type == 1 (peat dominated); "
            "PEATMAP = 0 otherwise, including GPA22 class 2."
        )
        ds_out.attrs["note"] = (
            "Conservative product: only GPA22 class 1 is treated as peat; "
            "class 2 is treated as non-peat."
        )
    elif mode == "alpha050":
        ds_out["PEATMAP"].attrs["long_name"] = "Peat fraction from GPA22"
        ds_out["PEATMAP"].attrs["description"] = (
            "Fractional peat cover derived from GPA22: 1.0 for class 1 "
            "(peat dominated), 0.5 for class 2 (peat in soil mosaic), "
            "0.0 otherwise."
        )
        ds_out.attrs["note"] = (
            "Alpha product with class 2 weighted by alpha = 0.5."
        )
    elif mode == "liberal":
        ds_out["PEATMAP"].attrs["long_name"] = "Peat mask from GPA22 including dominant and mosaic peat"
        ds_out["PEATMAP"].attrs["description"] = (
            "PEATMAP = 1 where GPA22 peatland_type is 1 or 2; "
            "PEATMAP = 0 otherwise."
        )
        ds_out.attrs["note"] = (
            "Liberal binary peat mask: GPA22 classes 1 and 2 are both treated as peat."
        )
    else:
        raise RuntimeError(f"Unsupported mode: {mode}")

    ds_out.attrs["source"] = "Derived from GPA22 (Global Peatland Map 2.0) on PEATMAP grid"

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
    if GPA_SOURCE_KIND == "nc":
        ds_gpa, lon_gpa, lat_gpa, ptype = load_gpa_from_nc(GPA_NC_FILE)
    elif GPA_SOURCE_KIND == "tif":
        ds_gpa, lon_gpa, lat_gpa, ptype = load_gpa_from_tif(GPA_TIF_FILE)
    else:
        raise RuntimeError(f"Unsupported GPA_SOURCE_KIND: {GPA_SOURCE_KIND}")

    print(f"GPA22 grid: nlat={len(lat_gpa)} nlon={len(lon_gpa)}")
    print(f"GPA22 resolution: dlat={get_resolution_arcsec(lat_gpa):.6f}\" "
          f"dlon={get_resolution_arcsec(lon_gpa):.6f}\"")

    ds_tgt, lon_tgt, lat_tgt = load_target_grid(PEATMAP_GRID_FILE)

    for mode in get_requested_products(PRODUCT_MODE):
        print()
        print("=" * 72)
        print(f"Building product: {mode}")
        print("=" * 72)

        field = build_mask(ptype, mode)
        field_on_target = remap_to_target(field, lon_gpa, lat_gpa, lon_tgt, lat_tgt)
        write_output(get_outfile(mode), lon_tgt, lat_tgt, field_on_target, mode)

    if ds_gpa is not None:
        ds_gpa.close()
    ds_tgt.close()

    print()
    print("Done.")


if __name__ == "__main__":
    main()

