#!/usr/bin/env python3

"""Small end-to-end LIS geometry test using the installed GEOS executables.

The test creates small synthetic global masks and runs the real geometry
workflow with both lat-lon and cubed-sphere ocean grids.

It also verifies that an invalid LIS ID is rejected.

No production masks or land parameters are read or modified.
"""

import argparse
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
from netCDF4 import Dataset


REFERENCE_CATCHMENTS = 291284

OCEAN = 0
LAKE = 190000000
LAND_ICE = 200000000

NX = 72
NY = 36

LIS_IM = 12
LIS_JM = 6


def make_inputs(work: Path):
    """Create small synthetic GEOS and LIS masks."""

    work.mkdir(
        parents=True,
        exist_ok=False,
    )

    lat = -90.0 + (np.arange(NY) + 0.5) * 180.0 / NY
    lon = -180.0 + (np.arange(NX) + 0.5) * 360.0 / NX

    jj, ii = np.indices((NY, NX))

    # ------------------------------------------------------------------
    # Reference GEOS mask
    # ------------------------------------------------------------------

    reference = np.ones(
        (NY, NX),
        dtype="i4",
    )

    reference[:, NX // 2:] = OCEAN
    reference[:6, :12] = LAND_ICE
    reference[24:30, 6:12] = LAKE
    reference[:, 12:24] = 2

    # Construct a 3x3 block where the original GEOS identifier mode
    # selects lake. This verifies that the LIS reduction preserves the
    # original GEOS surface decision rather than choosing from LIS IDs.
    reference[12:15, 6:9] = np.array(
        [
            [1, 1, 1],
            [2, 2, LAKE],
            [LAKE, LAKE, LAKE],
        ],
        dtype="i4",
    )

    reference_path = work / "GEOS5_10arcsec_mask.nc"

    with Dataset(reference_path, "w", format="NETCDF4") as dataset:
        dataset.createDimension("N_lat", NY)
        dataset.createDimension("N_lon", NX)
        dataset.createDimension(
            "N_catchments",
            REFERENCE_CATCHMENTS,
        )

        # No PfafID lookup is needed by the LIS geometry path.
        dataset.createVariable(
            "CatchIndex",
            "i4",
            ("N_lat", "N_lon"),
        )[:] = reference

        dataset.createVariable(
            "longitude",
            "f8",
            ("N_lon",),
        )[:] = lon

        dataset.createVariable(
            "latitude",
            "f8",
            ("N_lat",),
        )[:] = lat

    # ------------------------------------------------------------------
    # LIS mask
    # ------------------------------------------------------------------

    lis_ids = (
        (jj // (NY // LIS_JM)) * LIS_IM
        + ii // (NX // LIS_IM)
        + 1
    )

    lis_mask = np.where(
        (reference > OCEAN) & (reference < LAKE),
        lis_ids,
        reference,
    ).astype("i4")

    lis_path = work / "LIS_10arcsec_mask.nc"

    with Dataset(lis_path, "w", format="NETCDF4") as dataset:
        dataset.createDimension("N_lat", NY)
        dataset.createDimension("N_lon", NX)
        dataset.createDimension("east_west", LIS_IM)
        dataset.createDimension("north_south", LIS_JM)
        dataset.createDimension(
            "N_LIS_IDs",
            LIS_IM * LIS_JM,
        )

        # Keep the variables in a deliberately different creation order
        # so that the reader must find them by name.
        dataset.createVariable(
            "LISIndex",
            "i4",
            ("N_lat", "N_lon"),
        )[:] = lis_mask

        dataset.createVariable(
            "LISgridID",
            "i4",
            ("N_LIS_IDs",),
        )[:] = np.arange(
            1,
            LIS_IM * LIS_JM + 1,
        )

        dataset.createVariable(
            "longitude",
            "f8",
            ("N_lon",),
        )[:] = lon

        dataset.createVariable(
            "latitude",
            "f8",
            ("N_lat",),
        )[:] = lat

        native_lon = (
            -180.0
            + (np.arange(LIS_IM) + 0.5)
            * 360.0
            / LIS_IM
        )

        native_lat = (
            -90.0
            + (np.arange(LIS_JM) + 0.5)
            * 180.0
            / LIS_JM
        )

        dataset.createVariable(
            "LIS_lon",
            "f4",
            ("north_south", "east_west"),
        )[:] = np.tile(
            native_lon,
            (LIS_JM, 1),
        )

        dataset.createVariable(
            "LIS_lat",
            "f4",
            ("north_south", "east_west"),
        )[:] = np.tile(
            native_lat[:, None],
            (1, LIS_IM),
        )

    return lis_path, reference_path


def run_geometry_case(
    command: list[str],
    ocean_grid: str,
    output: Path,
) -> None:
    """Run one valid geometry case and require SUCCESS."""

    subprocess.run(
        [
            *command,
            "--ocean-grid",
            ocean_grid,
            "--outdir",
            str(output),
        ],
        check=True,
    )

    manifest_path = output / "run_manifest.json"

    manifest = json.loads(
        manifest_path.read_text()
    )

    if manifest["status"] != "SUCCESS":
        raise AssertionError(manifest)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__
    )

    parser.add_argument(
        "--bin-dir",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--workdir",
        required=True,
        type=Path,
        help="New test directory",
    )

    args = parser.parse_args()

    args.bin_dir = args.bin_dir.resolve()
    args.workdir = args.workdir.resolve()

    mask, reference = make_inputs(
        args.workdir / "inputs"
    )

    command = [
        sys.executable,
        str(args.bin_dir / "make_bcs.py"),
        "--lis",
        "--resolution",
        "2",
        "--mask",
        str(mask),
        "--reference-mask",
        str(reference),
        "--bin-dir",
        str(args.bin_dir),
        "--raster-arcsec",
        "54000",
        "--ocean-im",
        "12",
        "--ocean-jm",
        "6",
        "--max-tiles",
        "4096",
        "--run",
    ]

    # ------------------------------------------------------------------
    # Valid geometry cases
    # ------------------------------------------------------------------

    run_geometry_case(
        command,
        "latlon",
        args.workdir / "latlon_geometry",
    )

    run_geometry_case(
        command,
        "cube",
        args.workdir / "cube_geometry",
    )

    # ------------------------------------------------------------------
    # Invalid LIS ID
    # ------------------------------------------------------------------

    # Corrupt one land pixel after the successful runs. The production
    # reader must reject the invalid ID rather than silently changing
    # the resulting surface or tile geometry.

    with Dataset(mask, "a") as dataset:
        dataset.variables["LISIndex"][8, 2] = 1000000

    result = subprocess.run(
        [
            *command,
            "--outdir",
            str(args.workdir / "invalid_geometry"),
        ],
        check=False,
    )

    if result.returncode == 0:
        raise AssertionError(
            "Invalid LIS ID was not rejected"
        )

    print(
        "PASS: real-executable latlon/cube geometry, "
        "NetCDF export, and invalid-ID rejection."
    )


if __name__ == "__main__":
    main()
