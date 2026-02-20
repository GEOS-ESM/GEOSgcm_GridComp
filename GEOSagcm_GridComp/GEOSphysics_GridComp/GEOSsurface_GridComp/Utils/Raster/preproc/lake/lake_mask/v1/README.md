LakeTopoCat Lake Mask Preprocessing for GEOS BCS - February 2026
Overview

This directory contains preprocessing scripts to convert HydroLAKES-TopoCat v1.1 (2023) lake polygons into 30-arcsecond raster tiles compatible with make_bcs.

Two-step workflow

preproc_lake_to_30arcsec.py

split_lake_30arcsec_to_HV.py

Final output format matches GEOS static datasets (36×18 10°×10° tiles).

Source Dataset

HydroLAKES-TopoCat v1.1 (2023)

Derived from HydroLAKES v1.0

Hydrography source: MERIT Hydro v1.0.1 (3 arc-sec resolution)

Input format: polygon shapefiles Lakes_pfaf_*.shp

CRS: EPSG:4326

Only polygon geometry is used.
No filtering by lake size, permanence, or type is applied.

Step 1 — Rasterization and Aggregation

Script:

preproc_lake_to_30arcsec.py
Method

Rasterize polygons at 10 arc-second resolution

Aggregate 3×3 blocks → 30 arc-second grid

Produce global file:

LakeTopoCat_Global_30arcsec.nc4
Output Variables
Variable	Type	Description
lake_presence_frac	float32	Fractional lake coverage [0–1]
lake_presence_any	uint8	Binary lake mask (0/1)

No missing values are allowed.
All grid cells contain valid values.

Why 10" → 30"?

Oversampling at 10 arc-seconds improves shoreline representation and fractional accuracy relative to direct 30 arc-second rasterization.

Step 2 — Split into GEOS 36×18 Tiles

Script:

split_lake_30arcsec_to_HV.py

Produces:

LakeTopoCat_30arcsec_H##V##.nc4

Each file contains:

1200 × 1200 grid

30 arc-second resolution

Format consistent with GEOS static datasets (e.g., snow albedo, soil properties)

Global Grid Properties

43200 × 21600

Metadata includes:

i_ind_offset_LL

j_ind_offset_LL

CellSize_arc_Secs = 30

GEOS Usage

During tiling:

Each 30" pixel contributes area-weighted lake fraction to tiles

Tile classified as lake if:

lake_fraction ≥ 0.5

All 36×18 tiles must exist.
Missing files are not permitted.

Final Notes

Resolution: 30 arc-seconds

No missing values

lake_presence_frac ∈ [0, 1]

lake_presence_any ∈ {0, 1}
