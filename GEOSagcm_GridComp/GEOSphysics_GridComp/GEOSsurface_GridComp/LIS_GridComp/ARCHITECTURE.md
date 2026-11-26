# LIS_GridComp (LIS): architecture

**Audience:** coding agents and developers working on the GEOS wrapper around the Land Information System (LISF), and on its state specs.
**Describes:** the MAPL2 line; `GEOSgcm_GridComp` checked out on branch `develop` (the fixture's `components.yaml` pins tag `v3.1.1`), with `@LISF` on branch `feature/GEOS-LIS-cpl`.
**Verified against:** `GEOSgcm_GridComp` @ `3d353b60` plus the uncommitted state-spec update of 2026-11-25 (section 6; lost and redone 2026-11-26); `@LISF` @ `27a6f95c`. Checked 2026-11-26.
**Placement:** as in the GEOSgcm fixture (`components.yaml`): `src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/LIS_GridComp`.
Other fixtures (for example GEOSadas, GEOSldas) may place this repository elsewhere; paths inside this repository do not change.

Paths below are relative to `LIS_GridComp/` unless they start with `@LISF/` (the nested LISF clone in this directory) or `SURF/` (`GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/`).
Parent document: [Surface](../ARCHITECTURE.md). Build details: [`LISF_cmake/README.md`](LISF_cmake/README.md). Repository overview: [../../../../docs/ARCHITECTURE.md](../../../../docs/ARCHITECTURE.md).

## 1. Purpose
`LIS_GridComp` (module `LIS_GridComp`, `LIS_GridComp.F90`) is a MAPL/ESMF gridded component that runs LISF's land surface model inside GEOS. It calls LISF's core Fortran API directly (`LIS_config_init`, `LIS_domain_init`, `LIS_surfaceModel_run`, ...), the same routines that LISF's offline driver (`@LISF/lis/offline/lisdrv.F90`) and its retrospective run mode (`@LISF/lis/runmodes/retrospective/retrospective_runMod.F90`) use. It moves GEOS fields to and from LIS tile space with LIS's own `LIS_grid2tile` and `LIS_tile2grid`. No NUOPC layer is involved (module header comment). The source comments call it "the Plug".

The land model the field catalog targets is NoahMP 5.0 (`@LISF/lis/surfacemodels/land/noahmp.5.0/cpl_geos/reference/LIS_FieldCatalog_RAR.F90`). *Unverified:* which LSM a run actually uses; that is chosen in the LIS configuration file (section 5).

## 2. Place in the hierarchy and build
- **Not yet a child of anything at run time.** `SURF/GEOS_SurfaceGridComp.F90` `SetServices` adds `SALTWATER`, `LAKE`, `LANDICE`, `LAND` and `ROUTE` with `MAPL_AddChild`, but not LIS. No other `MAPL_AddChild` for `LIS_GridComp` was found in this repository. See Observation 1.
- **Built only on request.** `SURF/CMakeLists.txt` appends `LIS_GridComp` to the Surface library's `SUBCOMPONENTS` only `if (BUILD_LIS)`. `BUILD_LIS` has no `option()`; it is set on the command line (`-DBUILD_LIS=On`, `LISF_cmake/README.md`).
- **Libraries.** `CMakeLists.txt` builds the shared library `LIS_GridComp` from `LIS_GridComp.F90` alone, depending on `MAPL` and `LIS`. `LISF_cmake/CMakeLists.txt` builds `LIS` from the `@LISF` sources under CMake, using LISF's own plugin selector (`LISF_cmake/README.md`, "How the source list is derived").
- **State specs.** `CMakeLists.txt` runs `mapl_acg` on `LIS_StateSpecs.rc` with `IMPORT_SPECS EXPORT_SPECS`, generating `LIS_Import___.h` and `LIS_Export___.h`, which `SetServices` includes. No internal specs (section 4.3).
- **Nested repository.** `@LISF` is the `LISF` entry in the fixture `components.yaml`: remote `https://github.com/NASA-LIS/LISF.git`, branch `feature/GEOS-LIS-cpl`, sparse checkout `config/LISF.sparse` (excludes `docs`, `env`, `ldt`, `lvt`, `lis/utils`, `lis/configs`, `lis/testcases`). Edit it only in its own repository.
- **Files not compiled.** `LIS_FieldCatalog.F90` and `LIS_Fields.F90` are kept but not in `SRCS` (`CMakeLists.txt`, comment). See section 6.

## 3. Entry points
Registered in `SetServices` (`LIS_GridComp.F90`): `ESMF_METHOD_INITIALIZE` = `Initialize`, `ESMF_METHOD_RUN` = `Run`, `ESMF_METHOD_FINALIZE` = `Finalize`, `ESMF_METHOD_WRITERESTART` = `Record`. One run phase. One LIS nest: the module parameter `NEST = 1`.

### 3.1 SetServices
Sets the entry points, includes the generated import and export specs, adds the timers `INITIALIZE`, `RUN` and `FINALIZE`, and calls `MAPL_GenericSetServices`.

### 3.2 Initialize
1. `MAPL_GenericInitialize`.
2. Reads resource `LISF_CONFIG_FILE:` (default `lis.config`) into `LIS_rc%lis_config_file`.
3. Hands LIS the GEOS VM and its MPI communicator (`ESMF_VMGet(..., mpiCommunicator=comm)`), so LIS runs on the same PETs as its parent.
4. Calls, in this order: `LIS_config_init(cmd_args=.false., vm=vm, comm=comm)`, `LIS_domain_init`, `LIS_param_init`, `LIS_surfaceModel_init`, `LIS_surfaceModel_setup`, `LIS_metforcing_init`. This is `LIS_init_retrospective` minus data assimilation, routing, RTM, irrigation and the application model. The code comment says the run-mode dispatch `lisinit(...)` is avoided on purpose, because it would pull in those subsystems and follow the `RUNMODE:` in the config file.
5. Asserts that the LIS local domain size equals the GEOS local grid: `LIS_rc%lnc(NEST) == local_counts(1)` and `LIS_rc%lnr(NEST) == local_counts(2)`.
6. If `LIS_rc%startcode == "restart"`, calls `LIS_surfaceModel_readrestart`.
7. Sets the LIS clock from the ESMF clock (`LIS_timemgr_set`).

### 3.3 Run
1. Gets import pointers `TAIR`, `QAIR`, `SWD`, `LWD`, `U`, `V`, `PSURF`, `PRCP`, `SNOWF_F`, and export pointers `AVGSURFT` (allocated), `QH`, `QLE`, `QS`, `SNOW`.
2. Sets the LIS clock to the current GEOS time.
3. Puts each forcing into LIS's forcing state with `PutForcing`:

   | LIS forcing (`LIS_FORC_AttributesMod`) | GEOS import |
   |---|---|
   | `LIS_FORC_Tair` | `TAIR` |
   | `LIS_FORC_Qair` | `QAIR` |
   | `LIS_FORC_SWdown` | `SWD` |
   | `LIS_FORC_LWdown` | `LWD` |
   | `LIS_FORC_Wind_E` | `U` |
   | `LIS_FORC_Wind_N` | `V` |
   | `LIS_FORC_Psurf` | `PSURF` |
   | `LIS_FORC_Rainf` | `PRCP` (Observation 3) |
   | `LIS_FORC_Snowf` | `SNOWF_F` |

4. `LIS_surfaceModel_f2t(NEST)`, then `LIS_surfaceModel_run(NEST)`.
5. Fills the exports that are allocated with `GetDiagnostic`, from LIS's model-output (MOC) history entries:

   | Export | MOC index (`LIS_histDataMod`) |
   |---|---|
   | `QH` | `LIS_MOC_QH` |
   | `QLE` | `LIS_MOC_QLE` |
   | `QS` | `LIS_MOC_QS` |
   | `SNOW` | `LIS_MOC_SWE` |

   All other exports are not filled (TODO in `Run`; Observation 2).

**Helpers** (module procedures, both return `rc`):
- `PutForcing(varname, gridvar, rc)`: finds the field `varname` in `LIS_FORC_State(NEST)` and scatters the grid array into it with `LIS_grid2tile`.
- `GetDiagnostic(moc_index, moc_name, gridvar, rc)`: takes `LIS_histData(NEST)%ptr_into_lsm_list(moc_index)%dataEntryPtr` and gathers it with `LIS_tile2grid`. It fails unless the variable is enabled in LIS's `MODEL_OUTPUT_LIST` and its time-averaging option is 0 or 2. Options 1 and 3 hold a running sum, not the current value. For option 2 it reads slot 2, otherwise slot 1.

### 3.4 Finalize and Record
- `Finalize`: `LIS_finalize`, then `MAPL_GenericFinalize`.
- `Record` (write restart): `LIS_surfaceModel_writerestart(NEST)`. LIS writes its own restart files (section 4.3).

## 4. State
All fields are 2-D (`xy`), no vertical location (`N`), from `LIS_StateSpecs.rc`. The field set comes from LISF's NoahMP 5.0 reference catalog (section 6). Names are the catalog's `stateName` upper-cased; long names are its `stdName`.

### 4.1 Imports (14)
- **Used by `Run`:** `PSURF`, `TAIR`, `QAIR`, `LWD`, `SWD`, `U`, `V`, `SNOWF_F`, `PRCP`.
- **Declared, not used:** `RAINFCONV_F`, `RAINFSHCONV_F`, `RAINFLS_F`, `GRAUPELF`, `HAILF`. The catalog marks these as not mapped in NoahMP50_main.F90.

### 4.2 Exports (38)
- **Filled by `Run`:** `QH`, `QLE`, `QS`, `SNOW`.
- **Declared, not filled:** `QSB`, `TOTALRUNOFF`, `AVGSURFT`, `AVGSURFQ`, `TAIR_T2M`, `Z0`, `CH_F`, `CM_F`, `EVAP`, `ZPD`, `EMISS`, `ALBEDOVISDIR`, `ALBEDOVISDIF`, `ALBEDONIRDIR`, `ALBEDONIRDIF`, `USTARVEG`, `USTARBARE`, `USTAR`, `TAUX`, `TAUY`, `SWD`, `LAI`, `SAI`, `VAI`, `GREENNESS`, `LANDCOVER`, `CANOPINT`, `NEE`, `STC1`, `SMC1`, `ROOTMOIST`, `SNOWH`, `SNOCVR`, `RHOS`.
- `SWD` is both an import and an export (decision in section 8).

### 4.3 Internal state and restarts
There is no MAPL internal state: the catalog has no `directConn` fields. LIS keeps its model state in its own structures and writes and reads its own restarts (`Record`, `Initialize` step 6).

## 5. Configuration
- `LISF_CONFIG_FILE:` (MAPL resource, default `lis.config`): the LIS configuration file.
- In the LIS configuration file, *Inferred from the code*: the start mode (`LIS_rc%startcode`, `"restart"` triggers a restart read); a single nest whose local domain equals the GEOS local grid; and a `MODEL_OUTPUT_LIST` that enables `Qh`, `Qle`, `Qs` and `SWE` with time-averaging option 0 or 2 whenever the matching export is allocated (otherwise `GetDiagnostic` fails). *Unverified:* the configuration-file keys that set these.

## 6. Field catalog and state-spec history
- **Rules for sorting catalog fields into states** (`LIS_Fields.F90` comment; also in the `LIS_StateSpecs.rc` header): IMPORT = `adImport` and not `directConn`; EXPORT = `adExport` and not `directConn`; INTERNAL = `directConn`.
- **Old source:** `LIS_FieldCatalog.F90`, a copy of the field list and roles in `@LISF/lis/runmodes/nuopc_cpl_mode/LIS_NUOPC_Gluecode.F90` (`LIS_FieldList`, `LIS_HookupInit`), default build only. `LIS_Fields.F90` sorted it into states at run time. Both have been replaced by `LIS_StateSpecs.rc` and are no longer compiled (commit `8d008e5f`).
- **Current source:** `@LISF/lis/surfacemodels/land/noahmp.5.0/cpl_geos/reference/LIS_FieldCatalog_RAR.F90`, LISF commit `6fc91f27` (2026-11-17). It is reference only and not compiled (`@LISF/lis/surfacemodels/land/noahmp.5.0/cpl_geos/README.md`, "Files moved").
- **Spec update of 2026-11-25** (also recorded in the `LIS_StateSpecs.rc` header):
  - Field set replaced by the RAR catalog: 42 imports, 24 exports and 12 internals became 14 imports and 38 exports.
  - INTERNAL category removed (the soil fields `SMLIQFRACL1-4`, `SMFRACL1-4` and `SOILTEMPL1-4`). `INTERNAL_SPECS` was dropped from `mapl_acg` in `CMakeLists.txt`, and `#include "LIS_Internal___.h"` from `SetServices`.
  - Imports renamed in `Run`: `TAIR_F`→`TAIR`, `QAIR_F`→`QAIR`, `SWDOWN_F`→`SWD`, `LWDOWN_F`→`LWD`, `EWIND_F`→`U`, `NWIND_F`→`V`, `PSURF_F`→`PSURF`, `RAINF_F`→`PRCP`. Export `SWE`→`SNOW`.
  - Units converted to MAPL style: `kg/m2s`→`kg m-2 s-1`, `W/m2`→`W m-2`, `m/s`→`m s-1`, `N/m2`→`N m-2`, `g/m2s`→`g m-2 s-1`, `kg/kg`→`kg kg-1`, `kg/m2`→`kg m-2`, `m3/m3`→`m3 m-3`, `-`→`1`, and `mm/soil_timestep`→`kg m-2` (per LIS soil time step).
  - `SWD` kept as both import and export, as in the catalog.

## 7. Dependencies
- Fortran: `ESMF`, `MAPL`, and from LISF `LIS_coreMod`, `LIS_domainMod`, `LIS_paramsMod`, `LIS_surfaceModelMod`, `LIS_metforcingMod`, `LIS_FORC_AttributesMod`, `LIS_historyMod`, `LIS_histDataMod`, `LIS_timeMgrMod` (`use` list in `LIS_GridComp.F90`).
- CMake: `MAPL`, `LIS` (`CMakeLists.txt`). Toolchain limits (ifort stack only, HDF4 Fortran interface, linker driver, `-z notext`): `LISF_cmake/README.md` and `CMakeLists.txt` comments.

## 8. Decisions
### Call LISF's core API directly
- Pattern: the wrapper calls LIS's init, run and finalize routines itself instead of using LISF's NUOPC cap or the `lisinit` run-mode dispatch.
- Reason (code comments in `LIS_GridComp.F90` and `LIS_FieldCatalog.F90`): avoids building and linking the NUOPC cap and its dependencies, and fixes the init sequence regardless of `RUNMODE:` in the config file.

### State specs from the RAR catalog, generated by MAPL ACG
- Pattern: `LIS_StateSpecs.rc` is written by hand from the catalog and turned into specs by `mapl_acg`. It replaced building the specs at run time from `LIS_FieldCatalog`/`LIS_Fields` (commit `8d008e5f`).
- The 2026-11-25 update follows LISF's NoahMP 5.0 GEOS catalog rather than the NUOPC glue-code list.

### No MAPL internal state
- Pattern: LIS owns its state and restarts; the MAPL internal state is empty.
- Reason: the catalog has no `directConn` fields, and LIS writes its own restarts (`Record`).

### `SWD` as both import and export
- Pattern: `SWD` is in both states, as in the catalog (`adImport=.true.`, `adExport=. true.`).
- Reason: maintainer decision of 2026-11-25 to follow the catalog literally. It is recorded because the catalog's comment for `SWD` describes converting GEOS's partitioned shortwave fluxes into the bulk forcing NoahMP needs (a forcing step), and NoahMP does not change the downwelling shortwave, so the export would only return the import.

## 9. Observations (noted, not fixed) and open questions
1. **Not wired into the model.** No component adds `LIS_GridComp` as a child (section 2), so building with `BUILD_LIS` does not run LIS.
2. **Most exports are empty.** 34 of 38 exports are declared but never written in `Run`. `AVGSURFT` is even allocated (`ALLOC=.true.`) and left unset. A consumer would read uninitialized values.
3. **Snowfall is counted twice, or `Run` fails.** `PRCP` (total precipitation, rain plus snow) goes into `LIS_FORC_Rainf` and `SNOWF_F` into `LIS_FORC_Snowf` (`Run`). NoahMP 5.0 does not keep rain and snow forcings apart. `NoahMP50_f2t` (`@LISF/lis/surfacemodels/land/noahmp.5.0/NoahMP50_f2t.F90`) adds Rainf to its `prcp` and, when `LIS_FORC_Snowf%selectOpt == 1` and the value is not `LIS_rc%udef`, adds Snowf as well. `NoahMP50_main` then passes `prcp` as total precipitation (`RAINBL`), and NoahMP splits it into rain and snow itself (`IOPT_SNF`; `noahmp_driver_50.F90` sets the frozen fraction `SR = 0` and stops for `IOPT_SNF = 4`). That gives two cases, set by `Snowf:` in LIS's forcing variables list file (`get_forcingvar_attributes`, `add_forcing_fields` in `@LISF/lis/core/LIS_metforcingMod.F90`):
   - `Snowf:` selected: NoahMP receives `PRCP + SNOWF_F`, which is rain plus twice the snow.
   - `Snowf:` not selected: no Snowf field is created and its name is `NOT_SPECIFIED`, so `PutForcing` for `LIS_FORC_Snowf` fails at `ESMF_StateGet` and `Run` returns an error.

   For comparison, LISF's own GEOS run mode puts GEOS total precipitation into Rainf only and never sets Snowf (`@LISF/lis/runmodes/geos_cpl_mode/lisgeosrun_coupled.F90`). Since forcing fields start as `LIS_rc%udef` (`add_forcing_fields`) and `NoahMP50_f2t` skips undefined snowfall, leaving Snowf unset is safe. On the GEOS side, Surface imports `TPREC` (`total_precipitation`, from Moist; `SURF/GEOS_SurfaceGridComp.F90`). Catchment gets the precipitation split up instead (`PCU`, `PLS`, `SNO`, `ICE`, `FRZR`; `GEOSland_GridComp/GEOScatch_GridComp/GEOS_CatchGridComp.F90`). *Inferred:* the phase split GEOS provides is lost in NoahMP 5.0 as configured here, so the `RAINFCONV_F`, `RAINFSHCONV_F`, `RAINFLS_F`, `GRAUPELF` and `HAILF` imports (meant for `IOPT_SNF = 4`) cannot be used without changing `noahmp_driver_50.F90`.
4. **Units may not match the LIS diagnostics.** `QS`, `QSB` and `TOTALRUNOFF` are declared `kg m-2` per LIS soil time step, from the catalog's `mm/soil_timestep`. *Unverified:* the units of the `LIS_MOC_QS` entry that fills `QS`.
5. **`RHOS` (`surface_air_density`) has units `-` in the catalog**, declared here as `1`; expected `kg m-3`.
6. **The RAR catalog is not valid Fortran:** there is a stray `&.` after `adImport` for `2m_air_temperature`, and `. true.` for the `SWD` `adExport`. It is reference only, so this does not break the build.
7. **`LIS_FieldCatalog.F90` and `LIS_Fields.F90` are out of date.** They describe the old field set and no longer match `LIS_StateSpecs.rc`.
8. **Two GEOS coupling paths exist.** LISF has its own GEOS run mode, `@LISF/lis/runmodes/geos_cpl_mode` (`LISGEOSGridCompMod.F90`, `LISGEOSexport_module.F90`). It is enabled in `@LISF/lis/make/default.cfg` (`[GEOS coupling]`) and not disabled in `LISF_cmake/README.md`'s `user.cfg` list, so it is compiled into `LIS`. Its NoahMP exports go through `@LISF/lis/surfacemodels/land/noahmp.5.0/cpl_geos/noahMP50_setgeosexport.F90`, while this wrapper reads MOC history entries instead. *Unverified:* whether the wrapper should use that export path.
9. **The build README describes LISF differently from `components.yaml`.** `LISF_cmake/README.md` describes LISF as "unmodified" and `GEOSgcm_GridComp` as on branch `feature/pchakrab/integrate-lis-into-v11`. `components.yaml` pins `LISF` to branch `feature/GEOS-LIS-cpl`, which carries GEOS coupling commits, and this checkout of `GEOSgcm_GridComp` is on `develop`.
10. **Generic import names.** `U` and `V` are very generic. *Unverified:* whether they clash with other fields when the component is connected in Surface.
11. **`Record` replaces MAPL's write-restart method.** *Unverified:* whether replacing `ESMF_METHOD_WRITERESTART` skips MAPL's generic restart handling. With no internal state this probably has no effect.

### Open questions for maintainers
- Where in the Surface tree should LIS run (a child of `LAND`, a replacement for `CATCH`, or a child of Surface), and on which grid or tile space?
- Should the exports come from the `geos_cpl_mode` export path (Observation 8) instead of MOC history entries?
