# Legacy GEOLAND2/MODIS LAI Preprocessing

## Overview

This directory contains archived IDL and shell scripts associated with the historical preprocessing
of GEOLAND2 and MODIS Leaf Area Index (LAI) climatology datasets used in GEOS land boundary conditions.

The scripts are preserved in their original form for reference and future reproducibility efforts. 
Their execution order, dependencies, and completeness have not been independently verified.

## Source Code

The scripts were collected from:

```
land01: /l_data/model_parameters/LAI
```

Two sets of preprocessing codes are included:

- **GEOLAND2-LAI/IDL/** – Legacy scripts associated with GEOLAND2 LAI preprocessing.
- **MODIS-LAI/IDL/** – Legacy scripts associated with MODIS LAI preprocessing.

All available source files are retained without modification.

## Historical LAI Boundary Conditions

The GEOS `makebcs` workflow supports `LAIBCS='MODGEO'`, which processes GEOLAND2 10-day and 
MODIS 8-day LAI climatologies and combines them through the `merge_lai_data()` routine to 
produce the final LAI boundary conditions.

The merged MODIS/GEOLAND2 product has been used for NL3 and subsequent land BC generations.

The corresponding input datasets are located under:

```
/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/land/veg/lai_grn/v2/
```

Including:

- `GEOLAND2_10-DayClim`
- `MODIS_8-DayClim`

The relationship between the archived scripts and the exact production procedure used to
generate these datasets has not been verified.

## References

- GEOS `makebcs/mkCatchParam.F90` – LAI processing selection.
- GEOS `makebcs/mod_process_hres_data.F90` – `merge_lai_data()` implementation.
- GEOSgcm_GridComp PR #1505 – New MODIS LAI time-series processing.

## Archive Status

These scripts are provided for historical reference and preservation only. 
They have not been updated, tested, or validated as part of this archive.
No changes to the existing LAI boundary conditions or processing behavior are introduced by adding these files.
