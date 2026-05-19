## LakeTopoCat Lake and Reach Mask Preprocessing for GEOS BCS

**February 2026**

## Overview

This repository contains preprocessing scripts that convert **HydroLAKES-TopoCat v1.1 (2023)** global **lake polygons** and **reach line features** into **30-arcsecond raster products** compatible with `make_bcs`.

The workflow produces:

- A native global **30″ raster** (43200 × 21600)
- 36 × 18 GEOS **10°×10° tiles** (1200 × 1200 each)
- Separate tiled products for **lakes** and **reaches**
- Tile-level lake and reach information that can be read later during GEOS tile aggregation

The output format is consistent with other GEOS static datasets (soil properties, snow albedo, etc.).

## Workflow Summary

## Step 1 — Rasterization & Aggregation

### Lakes

**Script:** `preproc_lake_to_30arcsec.py`

#### Method

- Rasterize HydroLAKES-TopoCat lake polygons at **10 arc-second resolution**
- Aggregate 3×3 blocks → **30 arc-second grid**
- Process in **5° latitude bands**
- Combine bands into global file:  
  `LakeTopoCat_Global_30arcsec.nc4`

#### Why 10″ → 30″?

Oversampling at 10″ improves shoreline representation and allows estimation of **fractional lake coverage** within each 30″ pixel.

Direct 30″ rasterization would produce only a binary mask.

## Global Lake Output File

**Dimensions**

- 43200 (lon) × 21600 (lat)
- Resolution: 30 arc-seconds

### Variables

| Variable | Type | Description |
|----------|------|-------------|
| `lake_area_frac` | float32 | Fractional lake coverage in each 30″ cell [0–1] |
| `lake_presence_any` | uint8 | Binary lake presence (0/1) |

### Data Completeness

- Data are complete everywhere (no NaNs).
- `_FillValue` attributes exist but should not appear in valid data.

## Step 1b — Reach Rasterization & Aggregation

### Reaches

**Script:** `preproc_reach_to_30arcsec.py`

#### Method

- Rasterize HydroLAKES-TopoCat reach line features at **10 arc-second resolution**
- Use touch-based rasterization for line geometry
- Aggregate 3×3 blocks → **30 arc-second grid**
- Process in **5° latitude bands**
- Combine bands into global file:  
  `ReachTopoCat_Global_30arcsec.nc4`

#### Why a separate reach product?

Reaches are **linear features**, not polygons, so they do not represent true area coverage in the same way as lakes.

For this reason, the reach product is stored separately from the lake product and is intended for later experimentation during GEOS tile aggregation.

## Global Reach Output File

**Dimensions**

- 43200 (lon) × 21600 (lat)
- Resolution: 30 arc-seconds

### Variables

| Variable | Type | Description |
|----------|------|-------------|
| `reach_occupancy_frac` | float32 | Fraction of 10″ subpixels within each 30″ cell touched by a reach [0–1] |
| `reach_presence_any` | uint8 | Binary reach presence (0/1) |

## Step 2 — Split into GEOS 10°×10° Tiles

**Scripts:** `split_lake_30arcsec_to_HV.py`, `split_reach_30arcsec_to_HV.py`

Produces:

- `LakeTopoCat_30arcsec_H##V##.nc4`
- `ReachTopoCat_30arcsec_H##V##.nc4`

### Per-Tile Properties

For both lake and reach tile products:

- 1200 × 1200 pixels
- 30 arc-second resolution
- 1D coordinates:
  - `latitude(N_lat)`
  - `longitude(N_lon)`

### Per-Tile Variables

**Lake tiles** (`LakeTopoCat_30arcsec_H##V##.nc4`)

| Variable | Type | Description |
|----------|------|-------------|
| `lake_area_frac` | float32 | Fractional lake coverage in each 30″ cell [0–1] |
| `lake_presence_any` | uint8 | Binary lake presence (0/1) |

**Reach tiles** (`ReachTopoCat_30arcsec_H##V##.nc4`)

| Variable | Type | Description |
|----------|------|-------------|
| `reach_occupancy_frac` | float32 | Fraction of 10″ subpixels within each 30″ cell touched by a reach [0–1] |
| `reach_presence_any` | uint8 | Binary reach presence (0/1) |

### Global Attributes Included

- `N_lon_global = 43200`
- `N_lat_global = 21600`
- `i_ind_offset_LL`
- `j_ind_offset_LL`
- `CellSize_arc_Secs = 30`

## Source Datasets

**HydroLAKES-TopoCat v1.1 (2023)**

- Derived from HydroLAKES v1.0
- Hydrography source: **MERIT Hydro v1.0.1 (3 arc-sec resolution)**
- CRS: EPSG:4326

### Lake input

- Input: `Lakes_pfaf_*.shp`
- Geometry used: polygon geometry only
- No filtering by lake size, permanence, or type is applied

### Reach input

- Input: `Reaches_pfaf_*.shp`
- Geometry used: line geometry only
- Reaches are rasterized as touch-based linear features, not polygon area

## GEOS Usage (mkCatchParam Step 01)

### Current implementation

LakeTopoCat and ReachTopoCat data are mapped directly from the native 30″ raster products to GEOS tile space using the raster `tile_id` grid.

The Fortran aggregation reads the binary any-touch fields from the 10°×10° H/V files:

- `lake_presence_any`
- `reach_presence_any`

It does not use interpolation or a separate spatial remapping step. Each 30″ raster cell is mapped directly to a GEOS tile using `tile_id(iG,jG)`.

Processing runs only when the raster resolution is `43200 × 21600`. For coarser or alternative masks, LakeTopoCat / ReachTopoCat tile aggregation is skipped to avoid an implicit remap.

### Tile-Level Variable Created

| Variable | Type | Description |
|----------|------|-------------|
| `tile_lake_type` | int | Encoded LakeTopoCat / ReachTopoCat touch type for candidate lake tiles |

### `tile_lake_type` Coding

| Value | Meaning |
|-------|---------|
| `-9999` | UNDEF / excluded tile, usually `typ==100` |
| `0` | candidate tile with no LakeTopoCat lake touch and no ReachTopoCat reach touch |
| `1` | lake touch only |
| `2` | reach touch only |
| `3` | lake + reach touch |

### Candidate Tile Types

`tile_lake_type` is defined for:

- `typ == 0`
- `typ == 19`
- `typ == 20`

and set to `-9999` for excluded tiles such as:

- `typ == 100`

For the current standard EASE file, `typ==0` is absent, so the defined candidate space is effectively `typ==19` and `typ==20`. For CF and future tile files, `typ==0` may be present and will use the same coding.

### Notes

The intermediate preprocessing products still include fractional/occupancy diagnostics:

- `lake_area_frac`
- `reach_occupancy_frac`

These are retained for QA and diagnostics, but the Fortran `tile_lake_type` product is based on the binary any-touch fields:

- `lake_presence_any`
- `reach_presence_any`
