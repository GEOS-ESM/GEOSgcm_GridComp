#!/usr/bin/env python3

"""Geometry checks and LIS correspondence export; no land parameter generation."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Iterator

import numpy as np


# ASCII WriteTilingIR layout:
# surface, area, lon, lat, followed by four columns for each parent grid.
TILE_DTYPE = np.dtype([
    ('typ', 'i4'),
    ('area', 'f8'),
    ('lon', 'f8'),
    ('lat', 'f8'),
    ('i1', 'i4'),
    ('j1', 'i4'),
    ('frac1', 'f8'),
    ('k1', 'i4'),
    ('i2', 'i4'),
    ('j2', 'i4'),
    ('frac2', 'f8'),
    ('k2', 'i4'),
])

SURFACE_TYPES = (0, 19, 20, 100)
LAND_TYPE = 100


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def read_tiling(path: Path) -> tuple[dict, np.ndarray]:
    with Path(path).open(encoding='ascii') as handle:
        header_values = [int(value) for value in handle.readline().split()]

        require(
            len(header_values) == 4,
            f'{path}: expected current four-integer header',
        )

        ntile, nsource, nx, ny = header_values

        ngrid = int(handle.readline())

        require(
            ngrid in (1, 2),
            f'{path}: expected one or two grids',
        )

        grids = []

        for _ in range(ngrid):
            grid_name = handle.readline().strip()
            grid_im = int(handle.readline())
            grid_jm = int(handle.readline())

            grids.append((grid_name, grid_im, grid_jm))

        dtype = np.dtype(
            TILE_DTYPE.descr[:4 + 4 * ngrid]
        )

        rows = np.loadtxt(
            handle,
            dtype=dtype,
            ndmin=1,
        )

    require(
        len(rows) == ntile,
        f'{path}: tile-count mismatch',
    )

    require(
        ntile > 0 and nsource > 0 and nx > 0 and ny > 0,
        f'{path}: invalid header',
    )

    finite_fields = (
        'area',
        'lon',
        'lat',
        'frac1',
    )

    if ngrid == 2:
        finite_fields += ('frac2',)

    for name in finite_fields:
        require(
            np.isfinite(rows[name]).all(),
            f'{path}: non-finite {name}',
        )

    require(
        (rows['area'] > 0).all(),
        f'{path}: nonpositive tile area',
    )

    require(
        np.isin(rows['typ'], SURFACE_TYPES).all(),
        f'{path}: unexpected surface type',
    )

    header = {
        'ntile': ntile,
        'nsource': nsource,
        'nx': nx,
        'ny': ny,
        'grids': grids,
    }

    return header, rows


def raster_rows(
    path: Path,
    nx: int,
    ny: int,
) -> Iterator[np.ndarray]:
    """GEOS WRITERST format: little-endian int32 zonal Fortran records."""

    expected_bytes = 4 * nx

    with Path(path).open('rb') as handle:
        for j in range(ny):
            row_number = j + 1

            header = np.fromfile(
                handle,
                dtype='<i4',
                count=1,
            )

            require(
                len(header) == 1
                and int(header[0]) == expected_bytes,
                f'{path}: invalid record header at row {row_number}; '
                f'expected {expected_bytes} bytes',
            )

            values = np.fromfile(
                handle,
                dtype='<i4',
                count=nx,
            )

            trailer = np.fromfile(
                handle,
                dtype='<i4',
                count=1,
            )

            require(
                len(values) == nx
                and len(trailer) == 1
                and int(trailer[0]) == expected_bytes,
                f'{path}: short/corrupt record at row {row_number}',
            )

            yield values

        require(
            handle.read(1) == b'',
            f'{path}: trailing raster records',
        )


def check_tables(
    header: dict,
    tiles: np.ndarray,
    base_header: dict,
    base: np.ndarray,
    im: int,
    jm: int,
    nc: int,
    atol: float = 2.e-7,
) -> dict:
    native_id_count = im * jm

    require(
        header['nsource']
        == base_header['nsource']
        == native_id_count,
        'Native LIS ID bound was lost',
    )

    require(
        header['grids'][0][1:] == (nc, 6 * nc),
        'Unexpected cubed-sphere parent dimensions',
    )

    require(
        header['nx'] == base_header['nx']
        and header['ny'] == base_header['ny'],
        'Base and final working raster sizes differ',
    )

    require(
        np.isclose(
            tiles['area'].sum(),
            4 * math.pi,
            rtol=atol,
            atol=0,
        ),
        'Final tile areas do not cover the sphere',
    )

    valid_cf_indices = (
        (tiles['i1'] >= 1)
        & (tiles['i1'] <= nc)
        & (tiles['j1'] >= 1)
        & (tiles['j1'] <= 6 * nc)
    )

    require(
        valid_cf_indices.all(),
        'Invalid CF parent indices',
    )

    valid_atmospheric_weights = (
        (tiles['frac1'] > 0)
        & (tiles['frac1'] <= 1 + atol)
    )

    require(
        valid_atmospheric_weights.all(),
        'Invalid atmospheric weights',
    )

    cf_cell = (
        (tiles['j1'].astype('i8') - 1) * nc
        + tiles['i1']
        - 1
    )

    atmospheric_sums = np.bincount(
        cf_cell,
        weights=tiles['frac1'],
        minlength=6 * nc * nc,
    )

    require(
        np.allclose(
            atmospheric_sums,
            1,
            rtol=0,
            atol=atol,
        ),
        'Atmospheric fractions do not sum to one',
    )

    land = tiles['typ'] == LAND_TYPE
    lis_ids = tiles['i2'][land]

    require(
        (
            (lis_ids >= 1)
            & (lis_ids <= native_id_count)
        ).all(),
        'Land LIS IDs out of range',
    )

    require(
        (tiles['j2'][land] == 1).all(),
        'Unexpected land second-index convention',
    )

    native_area = np.bincount(
        lis_ids,
        weights=tiles['area'][land],
        minlength=native_id_count + 1,
    )

    native_fraction = np.bincount(
        lis_ids,
        weights=tiles['frac2'][land],
        minlength=native_id_count + 1,
    )

    expected_area = np.zeros(
        native_id_count + 1
    )

    base_land = base['typ'] == LAND_TYPE

    require(
        np.unique(base['i1'][base_land]).size
        == np.count_nonzero(base_land),
        'Duplicate native land IDs in base tiling',
    )

    expected_area[base['i1'][base_land]] = (
        base['area'][base_land]
    )

    require(
        np.allclose(
            native_area,
            expected_area,
            rtol=atol,
            atol=1.e-14,
        ),
        'LIS land area lost in intersections',
    )

    represented = expected_area > 0

    require(
        np.allclose(
            native_fraction[represented],
            1,
            rtol=0,
            atol=atol,
        ),
        'Fractions within a represented LIS land unit '
        'do not sum to one',
    )

    surface_area = {}

    for surface_type in SURFACE_TYPES:
        final_area = float(
            tiles['area'][
                tiles['typ'] == surface_type
            ].sum()
        )

        base_area = float(
            base['area'][
                base['typ'] == surface_type
            ].sum()
        )

        require(
            np.isclose(
                final_area,
                base_area,
                rtol=atol,
                atol=1.e-12,
            ),
            f'Surface area changed for type {surface_type}',
        )

        surface_area[str(surface_type)] = final_area

    if represented.any():
        max_lis_error = float(
            np.abs(
                native_fraction[represented] - 1
            ).max()
        )
    else:
        max_lis_error = 0.0

    return {
        'tiles': len(tiles),
        'land_tiles': int(land.sum()),
        'represented_LIS_cells': int(represented.sum()),
        'surface_area_steradians': surface_area,
        'max_atmospheric_weight_sum_error': float(
            np.abs(atmospheric_sums - 1).max()
        ),
        'max_LIS_weight_sum_error': max_lis_error,
    }


def check_raster_and_bounds(
    raster_path: Path,
    base_path: Path,
    header: dict,
    tiles: np.ndarray,
    base: np.ndarray,
    im: int,
    jm: int,
) -> tuple[np.ndarray, dict]:
    """Validate the complete working raster and compute tile bounds.

    Work on contiguous zonal runs rather than individual pixels when
    calculating extrema and areas.

    Geographic bounds use the conventional -180..180 enclosing interval.
    A seam-crossing nonland tile can therefore have a wide longitude bound.
    """

    nx = header['nx']
    ny = header['ny']
    ntile = len(tiles)

    require(
        nx % im == 0 and ny % jm == 0,
        'Working grid not nested in LIS',
    )

    dx = 360.0 / nx
    dy = 180.0 / ny

    pixels_per_lis_i = nx // im
    pixels_per_lis_j = ny // jm

    bounds = np.column_stack([
        np.full(ntile, np.inf),
        np.full(ntile, -np.inf),
        np.full(ntile, np.inf),
        np.full(ntile, -np.inf),
    ])

    area = np.zeros(ntile)
    pixels = np.zeros(
        ntile,
        dtype=np.int64,
    )

    base_rows = raster_rows(
        base_path,
        nx,
        ny,
    )

    for j, row in enumerate(
        raster_rows(raster_path, nx, ny)
    ):
        row_number = j + 1
        base_row = next(base_rows)

        require(
            (
                (row >= 1)
                & (row <= ntile)
            ).all(),
            f'Invalid final raster tile ID at row {row_number}',
        )

        require(
            (
                (base_row >= 1)
                & (base_row <= len(base))
            ).all(),
            f'Invalid base raster tile ID at row {row_number}',
        )

        require(
            np.array_equal(
                tiles['typ'][row - 1],
                base['typ'][base_row - 1],
            ),
            'Surface type changed during intersections '
            f'at row {row_number}',
        )

        starts = np.r_[
            0,
            np.flatnonzero(
                row[1:] != row[:-1]
            ) + 1,
        ]

        ends = np.r_[
            starts[1:],
            nx,
        ]

        tile_index = row[starts] - 1
        run_lengths = ends - starts

        low_lat = -90.0 + j * dy
        high_lat = -90.0 + (j + 1) * dy

        pixel_area = (
            math.radians(dx)
            * (
                math.sin(math.radians(high_lat))
                - math.sin(math.radians(low_lat))
            )
        )

        np.add.at(
            pixels,
            tile_index,
            run_lengths,
        )

        np.add.at(
            area,
            tile_index,
            run_lengths * pixel_area,
        )

        np.minimum.at(
            bounds[:, 0],
            tile_index,
            -180.0 + starts * dx,
        )

        np.maximum.at(
            bounds[:, 1],
            tile_index,
            -180.0 + ends * dx,
        )

        np.minimum.at(
            bounds[:, 2],
            tile_index,
            low_lat,
        )

        np.maximum.at(
            bounds[:, 3],
            tile_index,
            high_lat,
        )

        land_runs = (
            tiles['typ'][tile_index]
            == LAND_TYPE
        )

        lis_ids = (
            tiles['i2'][tile_index[land_runs]]
            - 1
        )

        expected_lis_j = (
            j // pixels_per_lis_j
        )

        require(
            np.all(
                lis_ids // im
                == expected_lis_j
            ),
            'LIS latitude assignment mismatch '
            f'at raster row {row_number}',
        )

        lis_i = lis_ids % im

        start_lis_i = (
            starts[land_runs]
            // pixels_per_lis_i
        )

        end_lis_i = (
            (ends[land_runs] - 1)
            // pixels_per_lis_i
        )

        require(
            np.all(lis_i == start_lis_i)
            and np.all(lis_i == end_lis_i),
            'LIS longitude assignment mismatch '
            f'at raster row {row_number}',
        )

        if row_number % 2000 == 0 or row_number == ny:
            print(
                f'Validate final raster: '
                f'{row_number:,}/{ny:,} rows',
                flush=True,
            )

    # Exhaust the base generator so that its trailing-record check
    # is also executed.
    require(
        next(base_rows, None) is None,
        'Unexpected extra base rows',
    )

    require(
        (pixels > 0).all(),
        'Output includes a tile with no raster pixels',
    )

    require(
        np.allclose(
            area,
            tiles['area'],
            rtol=2.e-7,
            atol=1.e-14,
        ),
        'ASCII tile areas do not agree with the final raster',
    )

    max_relative_area_error = float(
        np.max(
            np.abs(
                area / tiles['area'] - 1
            )
        )
    )

    stats = {
        'working_pixels_checked': int(nx) * int(ny),
        'surface_mismatches': 0,
        'LIS_assignment_mismatches': 0,
        'max_tile_relative_area_error': max_relative_area_error,
    }

    return bounds, stats


def export_mapping(
    path: Path,
    tiles: np.ndarray,
    im: int,
    jm: int,
    nc: int,
) -> None:
    columns = [
        'tile_id',
        'lis_id',
        'lis_i',
        'lis_j',
        'atm_i',
        'atm_j',
        'cf_face',
        'cf_j_on_face',
        'area_steradians',
        'fraction_atmospheric_cell',
        'fraction_LIS_GEOS_land',
        'fraction_LIS_entire_cell',
    ]

    with Path(path).open(
        'w',
        newline='',
        encoding='ascii',
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)

        land_tiles = np.flatnonzero(
            tiles['typ'] == LAND_TYPE
        )

        for tile_index in land_tiles:
            tile = tiles[tile_index]

            lis_j, lis_i = divmod(
                int(tile['i2']) - 1,
                im,
            )

            low_lat = math.radians(
                -90.0
                + lis_j * 180.0 / jm
            )

            high_lat = math.radians(
                -90.0
                + (lis_j + 1) * 180.0 / jm
            )

            lis_cell_area = (
                math.radians(360.0 / im)
                * (
                    math.sin(high_lat)
                    - math.sin(low_lat)
                )
            )

            atmospheric_j = int(
                tile['j1']
            )

            cf_face = (
                (atmospheric_j - 1) // nc
                + 1
            )

            cf_j_on_face = (
                (atmospheric_j - 1) % nc
                + 1
            )

            writer.writerow([
                int(tile_index) + 1,
                int(tile['i2']),
                lis_i + 1,
                lis_j + 1,
                int(tile['i1']),
                atmospheric_j,
                cf_face,
                cf_j_on_face,
                f"{tile['area']:.15g}",
                f"{tile['frac1']:.15g}",
                f"{tile['frac2']:.15g}",
                f"{tile['area'] / lis_cell_area:.15g}",
            ])


def finish_ascii_geometry(
    work: Path,
    name: str,
    im: int,
    jm: int,
    nc: int,
) -> dict:
    final_tiling = (
        work / 'til' / f'{name}.til'
    )

    base_tiling = (
        work / 'til' / 'LIS.til'
    )

    final_raster = (
        work / 'rst' / f'{name}.rst'
    )

    base_raster = (
        work / 'rst' / 'LIS.rst'
    )

    header, tiles = read_tiling(
        final_tiling
    )

    base_header, base = read_tiling(
        base_tiling
    )

    result = check_tables(
        header,
        tiles,
        base_header,
        base,
        im,
        jm,
        nc,
    )

    bounds, stats = check_raster_and_bounds(
        final_raster,
        base_raster,
        header,
        tiles,
        base,
        im,
        jm,
    )

    result.update(stats)

    bounds_path = (
        work / 'til' / f'{name}.bounds'
    )

    with bounds_path.open(
        'w',
        encoding='ascii',
    ) as output:
        output.write(f'{len(tiles)}\n')
        np.savetxt(
            output,
            bounds,
            fmt='%.12f',
        )

    export_mapping(
        work / f'{name}.lis_mapping.csv',
        tiles,
        im,
        jm,
        nc,
    )

    validation_path = (
        work / 'geometry_validation.json'
    )

    validation_path.write_text(
        json.dumps(
            result,
            indent=2,
        ) + '\n',
        encoding='utf-8',
    )

    return result


def check_netcdf(
    path: Path,
    ascii_path: Path,
) -> None:
    from netCDF4 import Dataset

    header, tiles = read_tiling(
        ascii_path
    )

    with Dataset(path) as dataset:
        dataset.set_auto_maskandscale(False)

        require(
            getattr(
                dataset,
                'land_discretization',
                '',
            ) == 'LIS',
            'LIS NetCDF metadata missing',
        )

        expected_lis_id = np.where(
            tiles['typ'] == LAND_TYPE,
            tiles['i2'],
            0,
        )

        expected_variables = (
            ('typ', tiles['typ']),
            ('lis_index', expected_lis_id),
        )

        for name, expected in expected_variables:
            require(
                name in dataset.variables,
                f'Converted NetCDF missing {name}',
            )

            require(
                np.array_equal(
                    dataset.variables[name][:],
                    expected,
                ),
                f'Converted NetCDF {name} mismatch',
            )

        require(
            int(dataset.N_LIS_IDs)
            == header['nsource'],
            'Converted NetCDF native count mismatch',
        )

        # MAPL compatibility slot must retain the same native
        # LIS source IDs.
        land = tiles['typ'] == LAND_TYPE

        require(
            np.array_equal(
                dataset.variables['pfaf_index'][:][land],
                expected_lis_id[land],
            ),
            'Converted NetCDF legacy source-index mismatch',
        )
