# LIS land geometry in the GEOS Raster package

This code builds GEOS tile geometry using the native LIS 0.25° land grid.

The LIS cell IDs are preserved while the LIS land area is intersected with the GEOS atmosphere and ocean grids. 
The result uses the normal GEOS `.til`, `.rst`, and NetCDF geometry formats.
This path is geometry only. It does **not** run `mkCatchParam`, create Catchment land parameters, create routing parameters, or modify the LIS runtime masks or parameters.

## Input files

The geometry workflow uses:

```text
/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/shared/mask/LIS_10arcsec_mask.nc
/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/shared/mask/GEOS5_10arcsec_mask.nc
```

`LIS_10arcsec_mask.nc` was supplied by Ronnie. It places the native LIS cell IDs on the original GEOS 10-arcsecond surface mask.

Native LIS grid:

```text
1440 x 720
0.25 degrees
1,036,800 possible native LIS IDs
```

The one-based ID convention is:

```text
LIS_ID = (row - 1) * 1440 + column
```

Surface values in the high-resolution mask are:

```text
0           ocean
1..1036800  LIS land-cell ID
190000000   lake
200000000   land ice
```

The original GEOS ocean, lake, and land-ice classifications are preserved. 
A LIS ID identifies the GEOS land portion of that 0.25° LIS cell.
The geometry workflow also reads `GEOS5_10arcsec_mask.nc` as the reference surface mask. 
This is important because the normal GEOS reduction rule is applied using the original GEOS identifiers before the selected land identifier is replaced by the LIS ID.

Ronnie also supplied:

```text
/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/shared/mask/lis_input.d01.WaterCutoff.processed.nc
```

That processed LIS file is used separately to check the LIS native grid, `LANDMASK`, and `DOMAINMASK`.
It is **not** read by the GEOS geometry driver and does not replace `LIS_10arcsec_mask.nc`.

The processed masks enable all 256,932 native LIS cells that contain GEOS land at 10-arcsecond resolution.

## Supported geometry

The current LIS path supports:

- the global regular 0.25° LIS grid
- a 30-arcsecond working raster
- a uniform cubed-sphere atmospheric grid
- a regular lat-lon ocean grid
- the same cubed-sphere grid as the ocean grid
- geometry-only output

Full real-mask runs have been completed and validated for:

```text
C90  + DE1440xPE0720 ocean
C180 + DE1440xPE0720 ocean
C360 + DE1440xPE0720 ocean
```

All three use the same 43,200 × 21,600 working raster and preserve the native LIS ID bound of 1,036,800.
The driver accepts other uniform cubed-sphere resolutions supported by the Raster utilities. 
C90, C180, and C360 are the cases validated so far.

This LIS path does not currently add support for:

- stretched cubed sphere
- EASE grids
- tripolar or freshwater/Reynolds adjustments
- Pfafstetter subdivision or routing for LIS
- Catchment land-parameter generation

### Working-raster limit

Existing Raster pixel counters are 32-bit.

The default 30-arcsecond working raster is:

```text
43,200 x 21,600 = 933,120,000 pixels
```

and fits within the current limit.

A global 10-arcsecond working raster would exceed the signed 32-bit pixel-count limit and is therefore not supported by this path.

This limit applies to the intermediate working raster, not directly to the cubed-sphere atmospheric resolution.

## How the 10″ mask is reduced to 30″

The standard GEOS NetCDF-mask path chooses the most frequent non-ocean identifier in each reduction block. If there is a tie, the lowest identifier wins.

The LIS IDs cannot simply replace the Catchment IDs before this step because that could change which surface type wins.

The LIS path therefore:

1. reads the original GEOS mask and the LIS mask;
2. applies the normal GEOS identifier-mode decision using the original GEOS mask;
3. keeps the resulting ocean, lake, land-ice, or land decision;
4. if land is selected, replaces the original GEOS land identifier with the native LIS cell ID.

The code does not average categorical IDs and does not use center sampling.

Before and during the reduction it checks that:

- the LIS and reference 10-arcsecond grids match;
- the high-resolution surface classifications agree;
- LIS land IDs fall inside the correct native 0.25° LIS cell;
- working-raster pixels nest exactly inside the LIS grid;
- invalid IDs or changed surface classes stop the run.

## Native LIS ID bound in `.til` files

The first line of a GEOS ASCII `.til` file contains:

```text
number_of_tiles  native_land_ID_bound  raster_nx  raster_ny
```

For LIS geometry, the second value must remain:

```text
1036800
```

`CombineRasters` reads this value from the second input raster as `nf2` and passes it to `WriteTiling`:

```fortran
call WriteTiling(..., N_LandUnits=nf2)
```

Without this argument, `WriteTiling` would use the legacy `SRTM_maxcat` value and the resulting LIS `.til` header would incorrectly describe the native ID range.

This changes the source-ID metadata carried in the output header. It does not change the number of intersection tiles.

## Build and tests

Build the patched GEOSldas/GEOSgcm checkout normally. The LIS changes are used by:

```text
mkLandRaster.x
CombineRasters.x
TileFile_ASCII_to_nc4.x
make_bcs.py --lis
```

Run the focused tests from this source directory:

```bash
python3 tests/test_lis_geometry.py
```

The current suite contains 22 tests. 
Some intentionally exercise error paths, so messages from deliberately failing commands can appear in the output.

A successful run ends with:

```text
Ran 22 tests
OK
```

After installing the executables, the small end-to-end smoke test can be run from a compute allocation:

```bash
python3 tests/run_lis_smoke.py \
    --bin-dir /path/to/GEOSldas/install/bin \
    --workdir /path/to/new/lis_smoke_test
```

The smoke test uses small synthetic grids and exercises:

- NetCDF reads by variable name
- LIS surface reduction
- real GEOS Raster executables
- lat-lon and cubed-sphere ocean cases
- intersections
- output validation
- NetCDF conversion
- rejection of an invalid LIS ID

It does not use the production masks.

## Running the real geometry workflow

Preview a run without executing it:

```bash
python3 /path/to/GEOSldas/install/bin/make_bcs.py --lis \
    --resolution 90 \
    --outdir /path/to/new/geometry_C90_LIS
```

Run it from an allocated compute job by adding `--run`:

```bash
python3 /path/to/GEOSldas/install/bin/make_bcs.py --lis \
    --resolution 90 \
    --outdir /path/to/new/geometry_C90_LIS \
    --run
```

Change `--resolution` for another uniform cubed-sphere resolution, for example:

```text
90
180
360
```

The driver:

- refuses to overwrite an existing output directory;
- runs in the current environment and does not submit a scheduler job;
- retains intermediate rasters and logs;
- exits nonzero if an executable or validation step fails.

The global 30-arcsecond `.rst` files are large, so real runs should use a compute allocation with sufficient disk space.

Useful options include:

```text
--mask
--reference-mask
--bin-dir
--raster-arcsec
--ocean-grid {latlon,cube}
--ocean-im
--ocean-jm
--max-tiles
```

The default ocean grid is:

```text
DE1440xPE0720
```

`--max-tiles` defaults to 8,000,000 for each `CombineRasters` operation. If the capacity is exceeded, the run fails rather than truncating the geometry.

## Main outputs

For a C90 run with the default ocean grid:

```text
til/LIS.til
rst/LIS.rst

til/CF0090x6C_DE1440xPE0720-LIS.til
til/CF0090x6C_DE1440xPE0720-LIS.nc4
rst/CF0090x6C_DE1440xPE0720-LIS.rst

CF0090x6C_DE1440xPE0720-LIS.lis_mapping.csv
geometry_validation.json
run_manifest.json
logs/
```

Atmospheric, ocean, and merged intermediate geometry files are also retained.

`run_manifest.json` records the input files, dimensions, resolution, commands, status, and run provenance.

`geometry_validation.json` records the geometry checks and summary statistics.

## IDs in the output

The `.til` files keep the existing GEOS ASCII layout.

For LIS geometry:

```text
header second integer = native LIS ID bound = 1,036,800
```

For a land row in the final two-grid `.til` file:

```text
column 9   native LIS cell ID
column 10  1
row number final GEOS tile number
```

Raster values and other internal indices are GEOS tile indices, not native LIS IDs.

A single LIS cell can produce several GEOS tiles where it intersects atmospheric cells. The native LIS IDs themselves are not compacted or renumbered.

## NetCDF LIS identifiers

The existing MAPL NetCDF geometry schema is retained.

The LIS conversion adds:

```text
lis_index
lis_i
lis_j
```

For land tiles these identify the native LIS cell. They are zero for nonland tiles.

For compatibility with the existing MAPL geometry schema, the legacy `pfaf_index` and `N_PfafCat` fields carry LIS source-ID information.

They are compatibility aliases only. They are **not** Pfafstetter hydrological catchment identifiers.

Do not use LIS geometry with unchanged Catchment or ROUTE physics as though these fields contained a Pfafstetter routing map.

## LIS mapping CSV

The mapping CSV contains one row for each final land intersection.

It records:

- final GEOS tile number
- native LIS ID
- LIS row and column
- cubed-sphere parent indices
- cubed-sphere face
- tile area in steradians
- fraction of the atmospheric cell
- fraction of the represented GEOS land in the LIS cell
- fraction of the complete 0.25° LIS rectangle

The last two fractions are intentionally different for cells that include coastlines, lakes, or land ice.

Land elevation is left undefined because this path does not generate a Catchment parameter file. LIS supplies its land parameters separately.

## Validation

Each real run checks:

- every working-raster pixel;
- surface type against the base LIS surface raster;
- native LIS geographic assignment;
- valid tile IDs;
- positive tile areas;
- total sphere area;
- agreement between raster-derived and `.til` areas;
- conservation by surface type;
- conservation by native LIS cell;
- atmospheric fractions summing to one;
- represented LIS-land fractions summing to one;
- agreement between the ASCII and NetCDF tile types and LIS IDs.

Real-mask validation completed successfully at C90, C180, and C360.

The final tile counts from the validated runs were:

```text
C90   1,492,723
C180  2,135,018
C360  3,690,940
```

All three retained:

```text
native LIS ID bound = 1,036,800
working raster      = 43,200 x 21,600
```

The C90 rerun also reproduced the previously delivered C90 `til/` geometry with no differences.

The separate 10-arcsecond LIS processed-mask check found:

```text
256,932 required LIS land-bearing cells
0 disabled required cells
0 extra enabled cells
0 missing or invalid LANDMASK/DOMAINMASK values
```

The C90 final geometry contains 256,630 land-bearing LIS cells after the 10-to-30-arcsecond surface-preserving reduction.

These checks validate the preprocessing geometry. They do not validate a coupled LIS simulation or the physical values of the LIS land parameters.
