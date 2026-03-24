# UW Shallow Water Convection physics port - from the Moist component

- This ports is aimed at the core numerics, e.g.: uwshcu.F90

- Version ported: v2.5.2 (as part of GEOS v11..5.2)

- `compute_uwshcu.py` contains the NDSL port of `UW_Run` in `GEOS_UW_InterfaceMod.F90`

- The `ComputeUwshcuInv` translate test passes on all 6 ranks with 0 ULP diff between the Fortran and DSL.

- However, we continue scientific validation with single-column runs and longer full-model runs on both CPU
and GPU.

- Last updated: 3/24/26