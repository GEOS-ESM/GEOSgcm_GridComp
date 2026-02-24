## LakeTopoCat Lake Mask Preprocessing for GEOS BCS

**February 2026**


## Overview

This repository contains preprocessing scripts that convert  **HydroLAKES-TopoCat v1.1 (2023)** global lake polygons into  **30-arcsecond raster tiles** compatible with `make_bcs`.

The workflow produces:

- A native global **30″ raster** (43200 × 21600)
- 36 × 18 GEOS **10°×10° tiles** (1200 × 1200 each)
- Tile-level lake fraction and classification written directly into the tile .nc4 file during tile generation.

The output format is consistent with other GEOS static datasets  (soil properties, snow albedo, etc.).

## Workflow Summary

## Step 1 — Rasterization & Aggregation

**Script:** `preproc_lake_to_30arcsec.py`

### Method

- Rasterize HydroLAKES-TopoCat polygons at **10 arc-second resolution**
- Aggregate 3×3 blocks → **30 arc-second grid**
- Process in **5° latitude bands**
- Combine bands into global file:
LakeTopoCat_Global_30arcsec.nc4

### Why 10″ → 30″?

Oversampling at 10″ improves shoreline representation and allows estimation  of **fractional lake coverage** within each 30″ pixel.

Direct 30″ rasterization would produce only a binary mask.

## Global Output File

**Dimensions**

- 43200 (lon) × 21600 (lat)  
- Resolution: 30 arc-seconds  

### Variables

| Variable | Type | Description |
|-----------|------|------------|
| `lake_presence_frac` | float32 | Fractional lake coverage in each 30″ cell [0–1] |
| `lake_presence_any`  | uint8   | Binary lake presence (0/1) |

### Data Completeness

- Data are complete everywhere (no NaNs).
- `_FillValue` attributes exist but should not appear in valid data.


## Step 2 — Split into GEOS 10°×10° Tiles

**Script:** `split_lake_30arcsec_to_HV.py`

Produces:
LakeTopoCat_30arcsec_H##V##.nc4

### Per-Tile Properties

- 1200 × 1200 pixels
- 30 arc-second resolution
- 1D coordinates:
  - `latitude(N_lat)`
  - `longitude(N_lon)`

### Global Attributes Included

- `N_lon_global = 43200`
- `N_lat_global = 21600`
- `i_ind_offset_LL`
- `j_ind_offset_LL`
- `CellSize_arc_Secs = 30`

## Source Dataset

**HydroLAKES-TopoCat v1.1 (2023)**

- Derived from HydroLAKES v1.0
- Hydrography source: **MERIT Hydro v1.0.1 (3 arc-sec resolution)**
- Input: `Lakes_pfaf_*.shp`
- CRS: EPSG:4326

Only polygon geometry is used.  No filtering by lake size, permanence, or type is applied.


## GEOS Usage (mkCatchParam Step 01)

LakeTopoCat data are mapped directly from the native 30″ raster to GEOS tiles using the raster `tile_id` grid.

## Mapping Characteristics

- Each 30″ pixel contributes equal weight (uniform pixel-area assumption).
- Because LakeTopoCat inputs and GEOS `tile_id` raster are on the same native 30″ grid,
   **no interpolation or spatial remapping is performed**.
- Tile means are simple averages of `lake_presence_frac` over all 30″ pixels belonging to each tile.
- Processing runs **only if raster resolution = 43200 × 21600**.
- Lake variables are skipped automatically for workflows using coarser raster masks (e.g., 8640×4320).

## Tile-Level Variables Created

| Variable | Description |
|-----------|------------|
| `tile_lake_frac` | Mean lake fraction per tile [0–1] |
| `tile_is_lake_50pct` | 1 if `tile_lake_frac ≥ 0.5`, else 0 |


### Key Properties

- Native resolution: 30 arc-seconds
- Fractional representation of lakes
- Binary lake classification derived from fraction
- Fully aligned with GEOS 30″ raster geometry
- No implicit remapping anywhere in the workflow
