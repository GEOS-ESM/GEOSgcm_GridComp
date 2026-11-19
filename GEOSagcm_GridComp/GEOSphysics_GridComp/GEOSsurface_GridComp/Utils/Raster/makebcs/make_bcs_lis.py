#!/usr/bin/env python3

"""Generate LIS/CF tile geometry only; leave legacy make_bcs workflows unchanged.

Run inside a compute allocation after sourcing the built GEOS environment.
By default prints a plan; --run executes. No scheduler is submitted implicitly.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys


MASK_ROOT = Path(
    "/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/shared/mask"
)

EXECUTABLES = (
    "mkCubeFVRaster.x",
    "mkLandRaster.x",
    "mkLatLonRaster.x",
    "CombineRasters.x",
    "TileFile_ASCII_to_nc4.x",
)

FULL_LONGITUDE_ARCSEC = 360 * 3600
FULL_LATITUDE_ARCSEC = 180 * 3600
MAX_RASTER_PIXELS = 2**31 - 1

DEFAULT_MAX_TILES = 8_000_000
HASH_BLOCK_SIZE = 8 * 1024 * 1024


def parser() -> argparse.ArgumentParser:
    arg_parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    arg_parser.add_argument(
        "--resolution",
        type=int,
        required=True,
        help="Uniform CF resolution, e.g. 90, 180, or 360",
    )
    arg_parser.add_argument(
        "--mask",
        type=Path,
        default=MASK_ROOT / "LIS_10arcsec_mask.nc",
    )
    arg_parser.add_argument(
        "--reference-mask",
        type=Path,
        default=MASK_ROOT / "GEOS5_10arcsec_mask.nc",
        help=(
            "Original standard GEOS mask; fixes working-raster "
            "surface classification"
        ),
    )
    arg_parser.add_argument(
        "--bin-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
    )
    arg_parser.add_argument(
        "--outdir",
        type=Path,
        required=True,
        help="New directory; existing paths are never overwritten",
    )
    arg_parser.add_argument(
        "--raster-arcsec",
        type=int,
        default=30,
    )
    arg_parser.add_argument(
        "--ocean-grid",
        choices=("latlon", "cube"),
        default="latlon",
    )
    arg_parser.add_argument(
        "--ocean-im",
        type=int,
        default=1440,
    )
    arg_parser.add_argument(
        "--ocean-jm",
        type=int,
        default=720,
    )
    arg_parser.add_argument(
        "--max-tiles",
        type=int,
        default=DEFAULT_MAX_TILES,
        help=(
            "Capacity for each CombineRasters run; "
            "raise only if it reports overflow"
        ),
    )

    action = arg_parser.add_mutually_exclusive_group()
    action.add_argument(
        "--run",
        action="store_true",
        help="Actually execute the geometry workflow",
    )
    action.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the plan (default)",
    )

    return arg_parser


def inspect_mask(path: Path) -> dict:
    from netCDF4 import Dataset

    required_dimensions = (
        "N_lon",
        "N_lat",
        "east_west",
        "north_south",
        "N_LIS_IDs",
    )

    required_variables = (
        "LISIndex",
        "LISgridID",
        "LIS_lat",
        "LIS_lon",
        "latitude",
        "longitude",
    )

    with Dataset(path) as dataset:
        sizes = {
            name: len(dataset.dimensions[name])
            for name in required_dimensions
        }

        for name in required_variables:
            if name not in dataset.variables:
                raise ValueError(f"Missing {name} in {path}")

        expected_ids = sizes["east_west"] * sizes["north_south"]

        if sizes["N_LIS_IDs"] != expected_ids:
            raise ValueError("LIS native grid and ID count differ")

    return sizes


def file_identity(path: Path, checksum: bool = False) -> dict:
    stat = path.stat()

    info = {
        "path": str(path),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }

    if checksum:
        digest = hashlib.sha256()

        with path.open("rb") as handle:
            for data in iter(
                lambda: handle.read(HASH_BLOCK_SIZE),
                b"",
            ):
                digest.update(data)

        info["sha256"] = digest.hexdigest()

    return info


def make_plan(
    args,
    native: dict,
) -> tuple[str, int, int, list[tuple[str, list[str]]]]:
    if not 1 <= args.resolution <= 9999:
        raise ValueError("CF resolution must be in 1..9999")

    if (
        args.raster_arcsec <= 0
        or FULL_LONGITUDE_ARCSEC % args.raster_arcsec
        or FULL_LATITUDE_ARCSEC % args.raster_arcsec
    ):
        raise ValueError(
            "Raster spacing must divide the full longitude and latitude spans"
        )

    nx = FULL_LONGITUDE_ARCSEC // args.raster_arcsec
    ny = FULL_LATITUDE_ARCSEC // args.raster_arcsec

    if nx * ny > MAX_RASTER_PIXELS:
        raise ValueError(
            "Current raster tile counters require <=2^31-1 pixels; "
            "use 30 arcsec or coarser"
        )

    if (
        nx % native["east_west"]
        or ny % native["north_south"]
    ):
        raise ValueError(
            "Working pixels must nest inside native LIS cells"
        )

    if (
        native["N_lon"] % nx
        or native["N_lat"] % ny
    ):
        raise ValueError(
            "Working raster must divide the input high-resolution raster"
        )

    if args.max_tiles < 1:
        raise ValueError("--max-tiles must be positive")

    # This path intentionally preserves GEOS5_10arcsec_mask.nc
    # surface semantics. The freshwater-lakes mask has separate
    # Reynolds handling and is not silently substituted here.
    if "freshwater-lakes" in args.reference_mask.name:
        raise ValueError(
            "LIS driver currently uses the standard GEOS mask, "
            "not the Reynolds freshwater-lakes variant"
        )

    cf_grid = f"CF{args.resolution:04d}x6C"

    common = [
        "-x",
        str(nx),
        "-y",
        str(ny),
    ]

    def executable(name: str) -> str:
        return str(args.bin_dir / name)

    lis_capacity = native["N_LIS_IDs"] + 3

    steps = [
        (
            "01_atmosphere",
            [
                executable("mkCubeFVRaster.x"),
                *common,
                str(args.resolution),
            ],
        ),
        (
            "02_LIS_surface",
            [
                executable("mkLandRaster.x"),
                *common,
                "-l",
                "-r",
                str(args.reference_mask),
                "-g",
                "LIS",
                "-t",
                str(lis_capacity),
                "-v",
            ],
        ),
    ]

    if args.ocean_grid == "latlon":
        if not (
            1 <= args.ocean_im <= 9999
            and 1 <= args.ocean_jm <= 9999
        ):
            raise ValueError(
                "Ocean dimensions must be in 1..9999"
            )

        if nx % args.ocean_im or ny % args.ocean_jm:
            raise ValueError(
                "Lat-lon ocean grid boundaries must align "
                "with working raster"
            )

        ocean_grid = (
            f"DE{args.ocean_im:04d}xPE{args.ocean_jm:04d}"
        )

        steps.append(
            (
                "03_ocean",
                [
                    executable("mkLatLonRaster.x"),
                    *common,
                    "-b",
                    "DE",
                    "-p",
                    "PE",
                    "-t",
                    "0",
                    str(args.ocean_im),
                    str(args.ocean_jm),
                ],
            )
        )

    else:
        ocean_grid = cf_grid

    merged_grid = f"{ocean_grid}-LIS"
    geometry_name = f"{cf_grid}_{merged_grid}"

    steps.extend(
        [
            (
                "04_surface_merge",
                [
                    executable("CombineRasters.x"),
                    "-f",
                    "0",
                    "-t",
                    str(args.max_tiles),
                    ocean_grid,
                    "LIS",
                ],
            ),
            (
                "05_CF_intersections",
                [
                    executable("CombineRasters.x"),
                    "-t",
                    str(args.max_tiles),
                    cf_grid,
                    merged_grid,
                ],
            ),
        ]
    )

    return geometry_name, nx, ny, steps


def run_step(
    label: str,
    command: list[str],
    work: Path,
    env: dict,
) -> None:
    command_text = shlex.join(command)

    print(f"{label}: {command_text}", flush=True)

    log = work / "logs" / f"{label}.log"

    with log.open("w", encoding="utf-8") as output:
        output.write(command_text + "\n")
        output.flush()

        result = subprocess.run(
            command,
            cwd=work,
            env=env,
            stdout=output,
            stderr=subprocess.STDOUT,
            check=False,
        )

    if result.returncode:
        tail = "\n".join(
            log.read_text(errors="replace").splitlines()[-25:]
        )

        raise RuntimeError(
            f"{label} failed ({result.returncode}). "
            f"Log: {log}\n{tail}"
        )


def main(argv=None) -> int:
    args = parser().parse_args(argv)

    for key in (
        "mask",
        "reference_mask",
        "bin_dir",
        "outdir",
    ):
        path = getattr(args, key)
        setattr(args, key, path.expanduser().resolve())

    manifest = {}

    try:
        native = inspect_mask(args.mask)

        if not args.reference_mask.is_file():
            raise FileNotFoundError(args.reference_mask)

        geometry_name, nx, ny, steps = make_plan(
            args,
            native,
        )

        lis_im = native["east_west"]
        lis_jm = native["north_south"]

        converter = [
            str(args.bin_dir / "TileFile_ASCII_to_nc4.x"),
            f"til/{geometry_name}.til",
            "--lis",
            str(lis_im),
            str(lis_jm),
            f"til/{geometry_name}.bounds",
        ]

        print(
            f"LIS geometry only: {geometry_name}\n"
            f"Working raster: {nx} x {ny}\n"
            f"Output: {args.outdir}"
        )

        for label, command in steps:
            print(f"{label}: {shlex.join(command)}")

        print(
            "06_validate: full working-raster and tile-weight checks; "
            "write bounds and LIS correspondence"
        )

        print(
            "07_netcdf: "
            + shlex.join(converter)
        )

        print(
            "No mkCatchParam, routing tables, Catchment README, "
            "or land-parameter packaging will run."
        )

        if not args.run:
            print(
                "Plan only. Add --run inside a compute allocation "
                "to execute."
            )
            return 0

        for executable_name in EXECUTABLES:
            executable_path = args.bin_dir / executable_name

            if (
                not executable_path.is_file()
                or not os.access(executable_path, os.X_OK)
            ):
                raise FileNotFoundError(
                    "Build/install the patched executable: "
                    f"{executable_path}"
                )

        args.outdir.mkdir(
            parents=True,
            exist_ok=False,
        )

        for directory in (
            "til",
            "rst",
            "logs",
        ):
            (args.outdir / directory).mkdir()

        # Avoid repeatedly hashing the large reference file.
        # Record its path, size, and modification time. Hash the
        # compressed LIS input itself.

        mask_id = file_identity(
            args.mask,
            checksum=True,
        )

        reference_id = file_identity(
            args.reference_mask,
        )

        manifest = {
            "status": "RUNNING",
            "started_utc": datetime.now(timezone.utc).isoformat(),
            "mask": mask_id,
            "reference_mask": reference_id,
            "native_dimensions": native,
            "resolution": args.resolution,
            "raster_shape": [ny, nx],
            "surface_policy": (
                "reference_GEOS_identifier_mode_then_LIS_ID"
            ),
            "geometry_name": geometry_name,
            "commands": steps + [
                ("07_netcdf", converter),
            ],
        }

        manifest_path = args.outdir / "run_manifest.json"

        manifest_path.write_text(
            json.dumps(manifest, indent=2) + "\n",
            encoding="utf-8",
        )

        env = os.environ.copy()
        env["MASKFILE"] = str(args.mask)
        env["MAKE_BCS_INPUT_DIR"] = str(
            args.mask.parent.parent.parent
        )

        for label, command in steps:
            run_step(
                label,
                command,
                args.outdir,
                env,
            )

        from make_bcs_lis_tools import (
            check_netcdf,
            finish_ascii_geometry,
        )

        result = finish_ascii_geometry(
            args.outdir,
            geometry_name,
            lis_im,
            lis_jm,
            args.resolution,
        )

        run_step(
            "07_netcdf",
            converter,
            args.outdir,
            env,
        )

        check_netcdf(
            args.outdir / "til" / f"{geometry_name}.nc4",
            args.outdir / "til" / f"{geometry_name}.til",
        )

        for path, original in (
            (args.mask, mask_id),
            (args.reference_mask, reference_id),
        ):
            current = file_identity(path)

            if (
                current["size"],
                current["mtime_ns"],
            ) != (
                original["size"],
                original["mtime_ns"],
            ):
                raise RuntimeError(
                    f"Input changed during geometry generation: {path}"
                )

        manifest.update(
            status="SUCCESS",
            geometry_validation=result,
            completed_utc=datetime.now(timezone.utc).isoformat(),
        )

        manifest_path.write_text(
            json.dumps(manifest, indent=2) + "\n",
            encoding="utf-8",
        )

        print(
            "Completed and validated geometry: "
            f"{args.outdir}/til/{geometry_name}.nc4",
            flush=True,
        )

        return 0

    except Exception as error:
        if manifest:
            manifest.update(
                status="FAILED",
                error=f"{type(error).__name__}: {error}",
            )

            (args.outdir / "run_manifest.json").write_text(
                json.dumps(manifest, indent=2) + "\n",
                encoding="utf-8",
            )

        print(
            f"ERROR: {error}",
            file=sys.stderr,
            flush=True,
        )

        return 1


if __name__ == "__main__":
    raise SystemExit(main())
