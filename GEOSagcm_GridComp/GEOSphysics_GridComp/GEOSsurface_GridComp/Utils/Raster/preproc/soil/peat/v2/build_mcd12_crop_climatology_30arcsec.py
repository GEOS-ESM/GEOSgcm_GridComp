#!/usr/bin/env python3

from pathlib import Path
import numpy as np
from netCDF4 import Dataset


# ============================================================
# USER SETTINGS
# ============================================================

# Root directory containing one subdirectory per year, with files like:
#   ROOT/<year>/MCD12Q1.061_*.nc4
ROOT = Path("/css/modis/Collection6.1/L3/Analysis_Ready/MCD12Q1.061")

# Final global 30-arcsec climatology file.
OUTFILE = Path("/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/land/soil/SOIL-DATA/v2/MCD12Q1_LCType1_crop_climatology_30arcsec.nc4")

# Number of 30-arcsec latitude rows written at a time.
# The code reads 2x this many rows from the 15-arcsec source grid.
OUT_CHUNK_ROWS = 256

# MCD12Q1 LC_Type1 classes of interest.
STRICT_CROP_CLASS = 12
MOSAIC_CLASS = 14

# Fill value for the final floating-point climatology fields.
NC_FILL = np.float32(-9999.0)


# ============================================================
# HELPERS
# ============================================================


def discover_files(root: Path):
    """
    Find one annual MCD12Q1 file per year directory.
    """
    files = []
    for year_dir in sorted(root.iterdir()):
        if not year_dir.is_dir():
            continue
        matches = sorted(year_dir.glob("MCD12Q1.061_*.nc4"))
        if matches:
            files.append(matches[0])
    return files



def aggregate_2x2_sum(arr):
    """
    Aggregate a (2*nrows, 2*ncols) array to (nrows, ncols) using a 2x2 sum.

    This is the step that converts from the 15-arcsec source grid to a
    30-arcsec grid.
    """
    arr = np.asarray(arr)
    nlat, nlon = arr.shape

    if (nlat % 2) != 0 or (nlon % 2) != 0:
        raise RuntimeError(
            "Expected even dimensions for 2x2 aggregation, got {}".format(arr.shape)
        )

    return arr.reshape(nlat // 2, 2, nlon // 2, 2).sum(axis=(1, 3), dtype=np.uint16)



def create_output(outfile: Path, lat30, lon30, years):
    """
    Create the final 30-arcsec output file.

    The output variables are climatology fractions (0..1).
    """
    ds = Dataset(outfile, "w", format="NETCDF4")

    ds.createDimension("lat", len(lat30))
    ds.createDimension("lon", len(lon30))

    lat_var = ds.createVariable("lat", "f8", ("lat",))
    lon_var = ds.createVariable("lon", "f8", ("lon",))

    chunks = (min(OUT_CHUNK_ROWS, len(lat30)), min(len(lon30), 2048))

    crop12_var = ds.createVariable(
        "crop12_clim", "f4", ("lat", "lon"),
        fill_value=NC_FILL,
        zlib=True, complevel=1, chunksizes=chunks,
    )
    crop14_var = ds.createVariable(
        "crop14_clim", "f4", ("lat", "lon"),
        fill_value=NC_FILL,
        zlib=True, complevel=1, chunksizes=chunks,
    )
    crop1214_var = ds.createVariable(
        "crop12or14_clim", "f4", ("lat", "lon"),
        fill_value=NC_FILL,
        zlib=True, complevel=1, chunksizes=chunks,
    )

    lat_var.units = "degrees_north"
    lon_var.units = "degrees_east"

    crop12_var.units = "1"
    crop14_var.units = "1"
    crop1214_var.units = "1"

    crop12_var.long_name = "Multi-year cropland climatology from LC_Type1 == 12"
    crop14_var.long_name = "Multi-year cropland/natural vegetation mosaic climatology from LC_Type1 == 14"
    crop1214_var.long_name = "Multi-year inclusive agriculture climatology from LC_Type1 in (12, 14)"

    ds.description = "Global 30-arcsec MCD12Q1 LC_Type1 agriculture climatology"
    ds.note = (
        "Each output value is a 30-arcsec climatology fraction in the range 0..1. "
        "The code first counts annual class occurrences on the native 15-arcsec grid, "
        "then aggregates 2x2 source cells to 30 arc-seconds, and finally divides by "
        "the aggregated valid-year count."
    )
    ds.source_root = str(ROOT)
    ds.source_resolution = "15 arc-second"
    ds.output_resolution = "30 arc-second"
    ds.n_years = len(years)
    ds.source_years = ",".join(str(y) for y in years)

    lat_var[:] = lat30
    lon_var[:] = lon30

    return ds


# ============================================================
# MAIN
# ============================================================


def main():
    files = discover_files(ROOT)
    if not files:
        raise RuntimeError("No annual MCD12Q1 files found under {}".format(ROOT))

    years = [int(f.parent.name) for f in files]
    print("Found {} annual files".format(len(files)))
    print("Years: {}".format(years))

    # Read the source grid from the first annual file.
    # Assumes all years are on the same 15-arcsec global grid.
    with Dataset(files[0]) as ds0:
        ds0.set_auto_mask(False)
        lat15 = np.asarray(ds0.variables["lat"][:], dtype=np.float64)
        lon15 = np.asarray(ds0.variables["lon"][:], dtype=np.float64)
        lc = ds0.variables["LC_Type1"]
        nlat15, nlon15 = lc.shape
        fill_value = getattr(lc, "_FillValue", 255)

    print("Input grid: lat={} lon={}".format(nlat15, nlon15))
    print("Fill value:", fill_value)

    if (nlat15 % 2) != 0 or (nlon15 % 2) != 0:
        raise RuntimeError(
            "Expected even 15-arcsec dimensions for 2x2 aggregation, got lat={} lon={}"
            .format(nlat15, nlon15)
        )

    # Build the 30-arcsec coordinate arrays by averaging adjacent 15-arcsec centers.
    lat30 = 0.5 * (lat15[0::2] + lat15[1::2])
    lon30 = 0.5 * (lon15[0::2] + lon15[1::2])

    nlat30 = len(lat30)
    nlon30 = len(lon30)
    print("Output grid: lat={} lon={}".format(nlat30, nlon30))

    out = create_output(OUTFILE, lat30, lon30, years)

    crop12_out = out.variables["crop12_clim"]
    crop14_out = out.variables["crop14_clim"]
    crop1214_out = out.variables["crop12or14_clim"]

    # Process the global grid a block of 30-arcsec rows at a time.
    # Each output row corresponds to 2 source rows on the 15-arcsec grid.
    for out_row0 in range(0, nlat30, OUT_CHUNK_ROWS):
        out_row1 = min(nlat30, out_row0 + OUT_CHUNK_ROWS)
        src_row0 = 2 * out_row0
        src_row1 = 2 * out_row1
        nrows30 = out_row1 - out_row0

        print(
            "Processing output rows {}:{} from source rows {}:{} ..."
            .format(out_row0, out_row1, src_row0, src_row1)
        )

        # Accumulate 30-arcsec counts directly, rather than writing an
        # intermediate 15-arcsec count file.
        crop12_count_30 = np.zeros((nrows30, nlon30), dtype=np.uint16)
        crop14_count_30 = np.zeros((nrows30, nlon30), dtype=np.uint16)
        valid_count_30 = np.zeros((nrows30, nlon30), dtype=np.uint16)

        for f in files:
            with Dataset(f) as ds:
                ds.set_auto_mask(False)
                arr = np.asarray(
                    ds.variables["LC_Type1"][src_row0:src_row1, :],
                    dtype=np.uint8,
                )

            valid = (arr != fill_value)
            is12 = (arr == STRICT_CROP_CLASS)
            is14 = (arr == MOSAIC_CLASS)

            valid_count_30 += aggregate_2x2_sum(valid.astype(np.uint8))
            crop12_count_30 += aggregate_2x2_sum(is12.astype(np.uint8))
            crop14_count_30 += aggregate_2x2_sum(is14.astype(np.uint8))

        # Class 12 and class 14 are mutually exclusive, so the inclusive
        # agriculture count is simply their sum.
        crop1214_count_30 = crop12_count_30 + crop14_count_30

        crop12_clim = np.full((nrows30, nlon30), NC_FILL, dtype=np.float32)
        crop14_clim = np.full((nrows30, nlon30), NC_FILL, dtype=np.float32)
        crop1214_clim = np.full((nrows30, nlon30), NC_FILL, dtype=np.float32)

        mask = valid_count_30 > 0
        crop12_clim[mask] = crop12_count_30[mask] / valid_count_30[mask]
        crop14_clim[mask] = crop14_count_30[mask] / valid_count_30[mask]
        crop1214_clim[mask] = crop1214_count_30[mask] / valid_count_30[mask]

        crop12_out[out_row0:out_row1, :] = crop12_clim
        crop14_out[out_row0:out_row1, :] = crop14_clim
        crop1214_out[out_row0:out_row1, :] = crop1214_clim

        out.sync()

    out.close()
    print("Wrote {}".format(OUTFILE))


if __name__ == "__main__":
    main()

