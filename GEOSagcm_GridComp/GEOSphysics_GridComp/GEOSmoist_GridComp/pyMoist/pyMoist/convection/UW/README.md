
- UW Shallow Water Convection physics port - from the Moist component.

- This ports is aimed at the core numerics, e.g.: uwshcu.F90

- Version ported: v2.5.2 (as part of GEOS v11..5.2)

- `compute_uwshcu.py` contains the NDSL port of `UW_Run` in `GEOS_UW_InterfaceMod.F90`

- The `ComputeUwshcuInv` translate test, which tests the contents of `compute_uwshcu.py`, fails for all 6 ranks in debug backend. However, many variables do pass the translate test and many variables that are failing have very small errors/percentages. A threshold could be applied to allow the translate test to pass under certain conditions.

- The errors that still exist likely stem from the `qtu` variable in the `buoyancy_sorting` stencil. The errors start off small, around 1e-10, then they grow and spread, resulting in translate test failure. We plan to revisit this and resolve these errors.

- `CNV_Tracers` is the worst failing variable for all 6 ranks. The errors likely originate in the `calc_pbl_fluxes` stencil, when calculating `trflx`, which is then used to update `CNV_Tracers`. We plan to revisit this and resolve these errors.

- Last updated: 2/12/26
