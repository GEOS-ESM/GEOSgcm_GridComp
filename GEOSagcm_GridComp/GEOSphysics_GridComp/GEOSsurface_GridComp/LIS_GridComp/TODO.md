# LIS_GridComp: TODO

Noted 2026-11-26. Re-check each item against the current source before fixing.

## Documents that predate the LIS merge

These documents were written before LIS was merged into `GEOSgcm_GridComp` (`3d353b60`, PR #1513) and leave LIS out entirely. `AGENTS.md` asks for documents to be kept current when a child list, build option or `components.yaml` entry changes.

Paths are relative to the `GEOSgcm_GridComp` root; `SURF/` abbreviates `GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/`.

### In `GEOSgcm_GridComp`

- [ ] `SURF/ARCHITECTURE.md`, section 2: says `CMakeLists.txt` builds Surface "with the five child directories as subcomponents". With `-DBUILD_LIS=On`, `SURF/CMakeLists.txt` appends `LIS_GridComp` as a sixth.
- [ ] `SURF/ARCHITECTURE.md`, section 2: the "Also in this directory" list omits `LIS_GridComp/`.
- [ ] `SURF/ARCHITECTURE.md`, section 9: lists one nested repository (`Topo`); `LIS_GridComp/@LISF` is a second.
- [ ] `SURF/ARCHITECTURE.md`: no link to [`LIS_GridComp/ARCHITECTURE.md`](ARCHITECTURE.md), which names Surface as its parent.
- [ ] `docs/ARCHITECTURE.md`: the Surface row of the directory table ("land, lake, landice, route, saltwater, `Shared/`, `Utils/`") omits LIS.
- [ ] `docs/ARCHITECTURE.md`: the nested-repository table lists `@ncar_topo` but not `@LISF`.
- [ ] `AGENTS.md`: the documents table does not list `LIS_GridComp/ARCHITECTURE.md`.

### In the GEOSgcm fixture (separate repository, needs its own branch)

- [ ] `docs/ARCHITECTURE.md`: the nested-repository table has Topo but not LISF (`components.yaml` entry `LISF`, sparse checkout `config/LISF.sparse`). *Unverified:* whether `BUILD_LIS` is mentioned anywhere.
- [ ] `AGENTS.md`: says "about 44 nested repositories". *Unverified:* whether that count includes LISF.

## Reusing LISF's GEOS coupling run mode (`@LISF/lis/runmodes/geos_cpl_mode`)

Read on 2026-11-26 at `@LISF` `27a6f95c`; nothing built or run. Paths below starting with `lis/` are inside `@LISF`.

### Candidate: fill exports through the geos_cpl export path

- [ ] Replace `GetDiagnostic` (MOC history entries) with LISF's export path. After `LIS_surfaceModel_run`, `LIS_surfaceModel_setexport(n)` dispatches to `NoahMP50_setgeosexport` (`lis/surfacemodels/land/noahmp.5.0/cpl_geos/noahMP50_setgeosexport.F90`). That routine calls `NoahMP50_convert_to_GEOSexport` and copies NoahMP fields into the tile arrays of `LISGEOS_export(n)`. The `RM_GEOS_COUPLING` block of `surfaceModel_setexport_noesmf` (`lis/core/LIS_surfaceModelMod.F90`) then grids them with `LIS_tile2grid`.
- [ ] In this wrapper, that would be `LISGEOS_alloc_states` once in `Initialize`, then `LIS_surfaceModel_setexport(NEST)` in `Run`, then `LISGEOS_GetExportVar(NEST, name, ptr, status)` (`LISGEOS_ExportAPIMod.F90`) for each allocated MAPL export. Its case labels accept our upper-case export names.
- [ ] Coverage: 35 of the 38 exports in `LIS_StateSpecs.rc` map by name. `USTAR` has no case in `LISGEOS_GetExportVar`. `VAI` has a case, but its arrays are never allocated (below). `SWD` has no LSM source (it would only return the import).
- [ ] Gains: removes the `MODEL_OUTPUT_LIST` and time-averaging requirements of `GetDiagnostic`. Derived fields are already computed in `NoahMP50_convert_to_GEOSexport`: total runoff, `ch_f = rhoair*ch`, `cm_f`, root-zone soil moisture, canopy water, band albedos, VAI and grid-mean ustar.
- [ ] Check differences from the RAR catalog comments: `qh_t` comes from `hfx`, `avgsurft_t` from `tsfcmean` (the catalog comment says `tg`), and `landcover` is a raw NoahMP `vegetype` passthrough.
- [ ] Check units against `LIS_StateSpecs.rc`: `totalrunoff` is `(runsf + runsb) / soiltimestep` (`NoahMP50_convert_to_GEOSexport`). The spec file declares `TOTALRUNOFF` as `kg m-2` per soil time step. *Unverified:* the units of `runsf`, `runsb` and `soiltimestep`.

### LISF bugs that block reuse (fix in `@LISF`, branch `feature/GEOS-LIS-cpl`, own PR)

- [ ] `LISGEOS_alloc_states` (`LISGEOSGridCompMod.F90`) never allocates `ustar_t`, `vai_t`, `sh2o4_t`, or the 2-D `ustar` and `vai`. `NoahMP50_setgeosexport` writes `ustar_t` and `vai_t`, and `surfaceModel_setexport_noesmf` grids `vai`/`vai_t` and `sh2o4`/`sh2o4_t`. So the first call touches unallocated arrays.
- [ ] `LISGEOS_reset_states` omits `ustar`/`ustar_t` and `vai`/`vai_t`, but sets the unallocated `sh2o4_t`. `LISGEOS_dealloc_states` only deallocates the outer arrays.
- [ ] `surfaceModel_setexport_noesmf` grids `avgsurft`, `emiss`, `qle`, `evap`, `greenness`, `landcover`, `nee`, `canopint` and `tair_t2m` twice (harmless, but redundant).
- [ ] `lis/surfacemodels/land/noahmp.5.0/cpl_geos/README.md` describes `NoahMP50_geosCplMod.F90` and `NoahMP50_geoscpl_init`, which do not exist. The fields are now members of `NoahMP50_struc`. `noahMP50_setgeosexport.F90` has a commented-out `!!call NoahMP50_struc_init(n)`.

### Configuration needed

- [ ] `lsm_setexport_noesmf` (`lis/core/LIS_lsmMod.F90`) looks up the routine by `trim(LIS_rc%lsm)//"+"//trim(LIS_rc%runmode)`. `NoahMP50_setgeosexport` is registered only under `LIS_geoscplId = "GEOS coupling"` (`lis/plugins/LIS_lsmcpl_pluginMod.F90`). So `lis.config` needs `Running mode: "GEOS coupling"`. *Unverified:* whether that run mode changes anything else in the init sequence `Initialize` calls.
- [ ] Confirm that `RM_GEOS_COUPLING` is defined and that `cpl_geos/` is in the `LIS` source list under `LISF_cmake` (`[GEOS coupling]` is `enabled: True` in `lis/make/default.cfg`; *Unverified* in the generated build inputs).

### Not reusable as is

- [ ] `lisgeosinit_coupled.F90`:
  - It calls `LIS_config_init(nx=, ny=, comm=)` instead of the `vm=`/`comm=` form used here.
  - It passes an uninitialized `coupled` to `LIS_metforcing_init`.
  - It hard-codes `LIS_rc%met_nf(:) = 17`.
  - It pulls in data assimilation, perturbation and irrigation.
  - The GEOS coupling run mode is not registered in `lis/plugins/LIS_runmode_pluginMod.F90`.
- [ ] `lisgeosrun_coupled.F90`:
  - It fills forcings from its own `LISGEOS_import` buffers, including `Emiss`, `Ch`, `Chs2`, `Cqs2`, `Q2sat`, `Forc_Hgt`, `Cosz`, `XICE` and `TMN`, which the RAR catalog does not import.
  - It reads `data_chs2`, `data_cqs2` and `data_xice`, which `LISGEOS_alloc_states` never allocates.
  - It divides precipitation by `LIS_rc%nts(n)`, treating it as an accumulated total in mm. This wrapper's `PRCP` is a rate (`kg m-2 s-1`).
  - It calls `ESMF_FieldGet` for `tmnField` twice.

### Steps this wrapper skips (investigate)

- [ ] `LIS_setDynparams(n)` before `LIS_surfaceModel_f2t`: `lisgeosrun_coupled` calls it; `Run` here does not. *Unverified:* which NoahMP 5.0 parameters then stay fixed (for example greenness or LAI climatology).
- [ ] `LIS_core_init` (updates the clock, dekadal alarms and alarm setups; `lis/core/LIS_coreMod.F90`): `lisgeosinit_coupled` calls it; `LIS_GridComp.F90` imports it but `Initialize` never calls it. *Unverified:* the effect.
- [ ] Output and restart order: `lisgeosrun_coupled` calls `LIS_surfaceModel_output` and `LIS_surfaceModel_writerestart` before the physics on every step except the first, following GEOS's write-then-run order (comment in that routine). Here, restarts are written only in `Record`, and `LIS_surfaceModel_output` is never called.
