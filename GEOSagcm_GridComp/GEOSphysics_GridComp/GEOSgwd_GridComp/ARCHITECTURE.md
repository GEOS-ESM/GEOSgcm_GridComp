# GEOSgwd_GridComp (GWD): architecture

**Audience:** coding agents and developers working on the gravity wave drag (GWD) gridded component (orographic, convective and background gravity-wave drag, Rayleigh friction).
**Describes:** the MAPL3 line, branch `release/MAPL-v3`.
**Verified against:** `GEOSgcm_GridComp` @ `84f0aae3` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `f50622d`; repository `GEOSgcm_GridComp`, mounted at `src/Components/@GEOSgcm_GridComp`): `src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSgwd_GridComp`.
Other fixtures (for example GEOSadas, GEOSldas) may place this repository elsewhere; paths inside this repository do not change.

Paths below are relative to `GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSgwd_GridComp/` unless they start with `GEOSagcm_GridComp/`.
**Build status on this line:** built and run. `GEOSagcm_GridComp/GEOSphysics_GridComp/CMakeLists.txt` lists `GEOSgwd_GridComp` in `alldirs` and adds it with `esma_add_subdirectories`. The Moist, Surface and Turbulence entries in the same list are commented out ("not yet ported to MAPL3"), and the Physics top-level library (`GEOS_PhysicsGridComp.F90`) is bypassed. GWD is ported to the MAPL3 API (`use MAPL, only: ...`) and has its own standalone regression test (section 3.4).

Parent: [Physics](../ARCHITECTURE.md). Siblings (not built on this line): [Moist](../GEOSmoist_GridComp/ARCHITECTURE.md), [Turbulence](../GEOSturbulence_GridComp/ARCHITECTURE.md), [Surface](../GEOSsurface_GridComp/ARCHITECTURE.md). Repository overview: [../../../docs/ARCHITECTURE.md](../../../docs/ARCHITECTURE.md).

## 1. Purpose
GWD computes the momentum and heat tendencies caused by sub-grid gravity waves that the model does not resolve. It has two independent implementations whose tendencies are added (`GEOS_GwdGridComp.F90`, Gwd_Driver):
- The **GEOS** scheme (`gw_drag.F90`, module `gw_drag`, entry `gw_intr`): an older scheme with a background (non-orographic) spectrum and orographic drag.
- The **NCAR** scheme (`ncar_gwd/`, module `gw_drag_ncar`, entry `gw_intr_ncar`): convective (Beres) and frontal/background sources, plus either an isotropic orographic scheme (`gw_oro`) or a ridge-based scheme (`gw_rdg`).

A Rayleigh friction layer (`RAYLEIGH_TAU1`) can be added in `postintr`. Which sources are active is chosen at run time (section 5). The GEOS and NCAR schemes are meant for different sources: the GEOS scheme for the extratropical background, the NCAR scheme for convective and orographic sources (`GEOS_GwdGridComp.F90`, Initialize, comments at the resource defaults).

## 2. Place in the hierarchy
- Parent is Physics. `GEOS_PhysicsGridComp.F90` still adds GWD as a leaf child `GWD` (`MAPL_AddChild(GC, NAME='GWD', SS=GwdSetServices, ...)`, SetServices) and runs it first, but that file is not compiled on this line (`GEOSphysics_GridComp/CMakeLists.txt`, "Temporarily bypass GEOSphysics top-level library"). Statements below about Physics (import sources, consumers, run order) describe that file as it is, not a running configuration. On this line the only configuration that runs GWD in this repository is the regression test (section 3.4).
- No children. The schemes are Fortran modules (`gw_drag`, `gw_drag_ncar`) called from the driver.
- `CMakeLists.txt` builds one shared library from `GEOS_GwdGridComp.F90`, `gw_drag.F90`, `machine.F90`, `gwdps.f`, `gwdc.f` and the `ncar_gwd/` sources. It runs the MAPL automatic code generator (`mapl_acg`) on `GWD_StateSpecs.rc` to generate the import, export and internal specs and the pointer declarations (`IMPORT_SPECS EXPORT_SPECS INTERNAL_SPECS GET_POINTERS DECLARE_POINTERS`). It installs `GWD_GridComp.rc` into `etc`, and adds the `regression/` subdirectory. Dependencies are `MAPL`, `ESMF::ESMF`, `NetCDF::NetCDF_Fortran`; `GEOS_Shared` was removed (comment: no GWD source uses it).

## 3. Entry points
Registered in `SetServices` (`GEOS_GwdGridComp.F90`) with `MAPL_GridCompSetEntryPoint`: `ESMF_METHOD_INITIALIZE` = `Initialize`, `ESMF_METHOD_RUN` = `Run` with `phase_name="run"`. There is one run phase. Finalize is not registered here. The shared-library entry symbol is `GWD_SetServices`, an external wrapper at the end of the file that calls the module's `SetServices`; the regression's `gwd.yaml` names it as `userRoutine: gwd_setservices_`.

### 3.1 SetServices
- Allocates a private state of derived type `GEOS_GwdGridComp` with the macros `_SET_NAMED_PRIVATE_STATE` / `_GET_NAMED_PRIVATE_STATE` under the name `"GWD_PRIVATE_STATE"` (`GEOS_GwdGridComp.F90`, SetServices).
- Reads `use_threads` with `MAPL_GridCompGetResource` (default `.false.`; `GWD_GridComp.rc` sets `.FALSE.`), sets the thread count with `MAPL_get_num_threads`, and allocates per-thread workspaces.
- Reads `NCAR_NRDG` before the generated includes (default 0) because the generated internal-spec header needs it to decide whether the ridge fields exist (`GEOS_GwdGridComp.F90`, SetServices, comment "We need to get NCAR_NRDG").
- Creates the ungridded dimension `ungrd_16` (`MAPL_UngriddedDim(16, NAME="sixteen", units="1")`) that the spec table names, then includes the generated headers `GWD_Import___.h`, `GWD_Export___.h`, `GWD_Internal___.h`. The generated headers contain `MAPL_GridCompAddSpec` calls: the generator (`src/Shared/@MAPL/apps/MAPL_GridCompSpecs_ACG.py`, `ADDSPEC = "MAPL_GridCompAddSpec"`) is run by `mapl_acg` (`src/Shared/@MAPL/cmake/mapl_acg.cmake`). The generated files are not checked in. There is no `MAPL_GenericSetServices` call on this line.

### 3.2 Initialize
Does not call `MAPL_GenericInitialize` on this line. It gets the grid, `LM` and coordinates with `MAPL_GridCompGet`, `MAPL_GridGet`, `MAPL_GridGetCoordinates`, then (`GEOS_GwdGridComp.F90`, Initialize):
1. Reads `AGCM.GRIDNAME` and `AGCM.STRETCH_FACTOR` to get a horizontal size `imsize`, and from it a tuning factor `sigma = 1 - 0.9839*exp(-0.09835*4e7*0.9/imsize/3000)`. For a cubed-sphere name, `imsize` is the digits between the third character and the `x`, times 4 when the last two characters are `CF`.
2. Reads the scheme resources and defaults (section 5).
3. Sets up the Beres convective source (`gw_common_init`, `gw_beres_init`) from the table in `BERES_FILE_NAME` (default `ExtData/g5gcm/gwd/newmfspectra40_dc25.nc`) and the `NCAR_*` resources.
4. Sets up the orographic scheme: `gw_rdg_init` if `NCAR_NRDG>0`, else `gw_oro_init`.
5. Sets up the Newtonian cooling profile with `gw_newtonian_set(LM, PREF, alpha)`; `PREF` is read from the import state with `MAPL_StateGetPointer` at this point.

### 3.3 Run
`Run` calls the internal routine `Gwd_Driver` on every call, between the timers `gwd_driver` (`MAPL_GridCompTimerStart/Stop`). The alarm test (`if ( ESMF_AlarmIsRinging( ALARM ) )`, with `ESMF_AlarmRingerOff`) is commented out, and the alarm variable is unused. The time step is read from the clock: `MAPL_ClockGet(clock, dt=DT)` in `Gwd_Driver`. State pointers come from the generated include `GWD_GetPointer___.h`, with the internal state from `MAPL_GridCompGetInternalState`. (`GEOS_GwdGridComp.F90`, Run, Gwd_Driver.)
`Gwd_Driver` does:
1. Get the vertical grid: pressure at layer midpoints (`PMID`), thickness (`PDEL`), height (`ZM`) from `pregeo` and `geopotential`.
2. **NCAR path** (`gw_intr_ncar`), if `NCAR_EFFGWORO` or `NCAR_EFFGWBKG` is non-zero. The convective source uses `HT_dc` and `(1 - CNV_FRC)^4 * HT_mi` and `WSPD_STABLE300M` when `NCAR_EFFGWBKG>0`. The orographic source needs `NCAR_EFFGWORO>0`: ridge scheme if `NCAR_NRDG>0`, isotropic otherwise.
3. **GEOS path** (`gw_intr`), if `GEOS_EFFGWORO` or `GEOS_EFFGWBKG` is non-zero, using `GEOS_PGWV` and `GEOS_BGSTRESS`.
4. Add the two sets of tendencies.
5. `postintr`: Rayleigh friction when `RAYLEIGH_TAU1>0` (with `RAYLEIGH_Z1`, `RAYLEIGH_H0`, `RAYLEIGH_HH`), and the potential- and kinetic-energy diagnostics.
6. Export the results. `DTDT` is multiplied by `PDEL` (mass-weighted, `Pa K s-1`, for Physics). `TTMGW` is the unweighted heating. The `*_EXP` aliases hold the winds and temperature after GWD (section 4).

### 3.4 Regression test
- `regression/CMakeLists.txt` reads `cases.txt` (one case, `gwd-ncar-sa`) and calls `esma_add_regression_tests(NAME gwd ... EXTRA_ENV OMP_NUM_THREADS=1)`; the data path is `GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSgwd_GridComp`.
- `regression/run_case.cmake` copies the case directory, skips the case if the regression data directory is missing, copies restarts and `newmfspectra40_dc25.nc`, runs the model on 6 processes with `run_geos`, and compares checkpoints against the stored ones only when the compiler is IntelLLVM (NaNs compare equal; coordinate variables excluded).
- Case `gwd-ncar-sa`: `cap_gridcomp.yaml` has a `Root` child (`libconfigurable_gridcomp`, config `root.yaml`) with a C12 cubed-sphere grid and 181 model levels. `Root` has children `GWD` (`libGEOSgwd_GridComp.so`, `gwd.yaml`) and `DataGwdImports` (a configurable component whose `data-gwd-imports.yaml` declares internal-state fields). `root.yaml` has 14 connections from `DataGwdImports` to the 14 GWD imports. `cap_driver.yaml` runs one 1200 s step and writes a final checkpoint.
- `gwd.yaml` sets `AGCM.GRIDNAME: PE12x72-CF`, `BERES_FILE_NAME`, `GEOS_EFFGWBKG: 0`, `NCAR_EFFGWBKG: 0.375`, `GEOS_EFFGWORO: 0`, `NCAR_EFFGWORO: 1.0`, `NCAR_NRDG: 16`, and under `mapl: misc` `activate_all_exports`, internal restart, and import/export/internal checkpoints. So the case exercises the NCAR scheme with the ridge scheme on; the GEOS scheme is off.

## 4. State
All specs come from `GWD_StateSpecs.rc` (header `schema_version: 2.0.0`, `component: GWD`; categories INTERNAL, IMPORT, EXPORT); that file is the source of truth for the list. It is a table consumed by `mapl_acg`, not YAML component config. The generated headers are not checked in. On this line the ungridded column uses the name `ungrd_16` instead of the MAPL2 `(/16/)`.

### 4.1 Imports (`GWD_StateSpecs.rc`, category IMPORT)
| Names | Restart | Source (via Physics) |
|---|---|---|
| `PLE`, `T`, `Q`, `U`, `V`, `AREA`, `PREF` | skip | `Q` from Moist; the others from Dynamics through the Agcm (`GEOS_PhysicsGridComp.F90`, SetServices, "Imports for GWD"; file not compiled on this line) |
| `PHIS`, `SGH`, `VARFLT` | skip | topography, set by Agcm in Initialize (`GEOSagcm_GridComp/GEOS_AgcmGridComp.F90`, Initialize) |
| `HT_dc`, `HT_mi` | restarted | Moist `DTDT_DC`, `DTDT_micro` |
| `CNV_FRC` | restarted | Moist |
| `WSPD_STABLE300M` | skip | Dynamics export of the same name (`GEOSagcm_GridComp/GEOSsuperdyn_GridComp/@FVdycoreCubed_GridComp/DynCore_GridCompMod.F90`, export spec). *Unverified:* the Agcm-level connection that carries it to GWD. |

`Q` is the Moist internal humidity, exported by Moist and connected to GWD (`MAPL_AddConnectivity`, `SRC_ID=MOIST`). Because GWD runs first, *Inferred:* `Q` is the value left by Dynamics and the previous Moist step. `HT_dc`, `HT_mi` and `CNV_FRC` are exports of Moist's previous run, copied by the coupler after Moist ran. *Inferred:* this is why these three imports are restarted and the others are not.

## 4.2 Exports
- Aliased copies of the state after GWD: `PLE_EXP`, `T_EXP`, `Q_EXP`, `U_EXP`, `V_EXP`, `SGH_EXP`, `PREF_EXP` (exported under the names `PLE`, `T`, `Q`, `U`, `V`, `SGH`, `PREF`). Physics connects `U`, `V`, `T` to Moist (`GEOS_PhysicsGridComp.F90`, SetServices, "Moist Imports").
- Tendencies: `DUDT`, `DVDT`, `DTDT` (mass-weighted), `TTMGW`; split by source `*_ORO`, `*_BKG`, `*_RAY`; stresses `TAUOROX/Y`, `TAUBKGX/Y`, `TAUGWX/Y`, `BKG_TAU_*`, `TAUGW*_TOT/EAST/WEST`; energy diagnostics such as `PEGWD`, `FEGW`; ridge diagnostics `RDG1_*`.
- Physics reads `DUDT`, `DVDT`, `DTDT` (as `UIG`, `VIG`, `TIG`) to build the total tendencies, and `TAUOROX`, `TAUOROY` are connected to Moist (`GEOS_PhysicsGridComp.F90`, Run, SetServices).

### 4.3 Internal state and restarts
- The ridge data `KWVRDG`, `EFFRDG`, `GBXAR`, `HWDTH`, `CLNGT`, `MXDIS`, `ANGLL`, `ANIXY` exist only if `NCAR_NRDG /= 0` (condition in `GWD_StateSpecs.rc`, category INTERNAL). All have an ungridded dimension of 16 ridges except `GBXAR`.
- They are read from the internal restart file. Restart names in the MAPL2-style `AGCM.rc`: `MODELATMGWD_INTERNAL_RESTART_FILE: gwd_internal_rst`, `MODELATMGWD_IMPORT_RESTART_FILE: gwd_import_rst` (`src/Applications/@GEOSgcm_App/AGCM.rc.tmpl`, fixture). The run script handles per-face ridge files (`gwd_internal_c${AGCM_IM}_face_${n}`) (`src/Applications/@GEOSgcm_App/gcm_run.j`, fixture).
- The private derived-type state (section 3.1) holds the scheme parameters and per-thread workspaces; it is not restarted.

## 5. Configuration that changes structure
Resources are read in Initialize with `MAPL_GridCompGetResource` (`GEOS_GwdGridComp.F90`, Initialize). Defaults depend on `LM`.

| Key | Default | Effect |
|---|---|---|
| `JASON_BKG` | `LM==72` | selects the legacy background-drag defaults |
| `GEOS_PGWV`, `GEOS_BGSTRESS`, `GEOS_EFFGWBKG`, `NCAR_EFFGWBKG`, `RAYLEIGH_TAU1` | with `JASON_BKG`: `4`, `0.9`, `0.125`, `0`, `0`. Without: `NINT(32*LM/181)`, `0`, `0`, `0.375`, `0` | which background scheme is on and with what strength |
| `JASON_ORO` | `LM==72` | selects the legacy orographic defaults |
| `GEOS_EFFGWORO`, `NCAR_EFFGWORO`, `NCAR_NRDG` | with `JASON_ORO`: `0.25`, `0`, `0`. Without: `0`, `0.25` (`1.0` if `NCAR_NRDG==16`), `0` | orographic scheme: GEOS, NCAR isotropic, or NCAR ridge |
| `NCAR_NRDG` | 0 (read in SetServices and Initialize) | 0 disables the ridge scheme; 1 to 16 enables it and adds the ridge internal state |
| `BERES_FILE_NAME`, `NCAR_DC_BERES` (default true), `NCAR_ET_FAC_DTDTM` (0.75), `NCAR_ET_FAC_WS300` (1.25), `NCAR_BKG_TAU` | | Beres convective source |
| `RAYLEIGH_TAU1`, `RAYLEIGH_Z1`, `RAYLEIGH_H0`, `RAYLEIGH_HH` | | Rayleigh friction |
| `use_threads:` in `GWD_GridComp.rc` | `.FALSE.` | OpenMP over horizontal blocks |
| `DEBUG_GWD`, `DEBUG_TQ_ERRORS` | | print diagnostics |

`NCAR_NRDG` is substituted into `AGCM.rc` by the setup script (`AGCM.rc.tmpl`, "TRB/GWD Parameterization" block). Changing it also changes the restart contents.
The names `JASON_BKG` and `JASON_ORO` are legacy tuning switches and are kept as is.

## 6. Dependencies
- Libraries: MAPL, ESMF, NetCDF Fortran, and OpenMP (not linked for NAG) (`CMakeLists.txt`). `GEOS_Shared` is no longer a dependency on this line.
- Imports from `MAPL` by name (`use MAPL, only: ...`), including `MAPL_get_current_thread` and `MAPL_find_bounds` for threading.
- Inputs: the Beres spectrum file (3.2); topography from Agcm; Moist heating and convective fraction; `WSPD_STABLE300M` from Dynamics.
- Consumed by: Physics (tendencies), Moist (`U`, `V`, `T`, `TAUOROX/Y`), and Physics's stochastic pattern code, which reads `PREF` from the GWD import (`GEOS_PhysicsGridComp.F90`, Run).

## 7. Known traps
1. **Two schemes are summed.** If both `GEOS_EFFGWBKG` and `NCAR_EFFGWBKG` (or both orographic ones) are non-zero, both act. Check all four keys when tuning.
2. **`NCAR_NRDG` changes the restart.** Turning the ridge scheme on without ridge data in `gwd_internal_rst` leaves the fields at their default values. *Unverified:* what the default is; the run script must supply the ridge files.
3. **`Initialize` parses `AGCM.GRIDNAME`.** It assumes the third character starts the digits and an `x` follows, as in `PEx...`/`CF` cubed-sphere names. Other formats make the `read` fail (`GEOS_GwdGridComp.F90`, Initialize, `imchar`). *Inferred:* the code is not tested with other grid-name forms.
4. **Run does not test an alarm on this line.** `Run` calls `Gwd_Driver` on each call, with `dt` from the clock. How often the parent calls it is not set in this repository. *Unverified:* whether any parent supercycles GWD on this line.
5. **`DTDT` is mass-weighted; `TTMGW` is not.** Use the right one.
6. **Generated headers.** `GWD_StateSpecs.rc` is the only source of the spec lists. Editing the Fortran for specs has no effect. A new ungridded dimension must also be created in `SetServices` (`ungrd_16`) for the table to name it.
7. **Run order.** GWD runs before Moist. Its `CNV_FRC`, `HT_dc`, `HT_mi` are therefore from the previous step (4.1).
8. **Per-thread workspaces.** When `use_threads` is on, any new module-level scratch array in the schemes is a race.

## 8. Decisions
### Two scheme families under one component
- Pattern: both the GEOS and the NCAR scheme are called from `Gwd_Driver` and added together; resources decide which is active (`GEOS_GwdGridComp.F90`, Gwd_Driver). [Verified]
- Likely reason: the NCAR convective and ridge sources were added alongside the older scheme without removing it, and tuning moves between them by resolution. [Inferred]
- Question for maintainers: is the long-term plan to keep both?

### State specs generated from `GWD_StateSpecs.rc`
- Pattern: imports, exports and the optional ridge internal specs are declared in a table, and CMake generates the Fortran includes (`CMakeLists.txt`, `mapl_acg`). [Verified]
- Likely reason: removes boilerplate and keeps the list of fields in one reviewable file. [Inferred]

### Ridge internal state conditional on a resource read before specs
- Pattern: `NCAR_NRDG` is read in `SetServices` so that the table's `COND` can create or omit the ridge fields (`GWD_StateSpecs.rc`, `COND`). [Verified]
- Likely reason: avoids a ~16-level ungridded set of unused fields and restart variables when the ridge scheme is off. [Inferred]

### Run first, from Dynamics' state
- Pattern: GWD is the first child in Physics `Run`, with a comment that it must use `Q` from Dynamics (`GEOS_PhysicsGridComp.F90`, Run; file not compiled on this line). [Verified]
- Likely reason: gravity-wave drag is tied to the dynamical state and not to the moist or turbulent adjustments; Moist then sees GWD-updated winds and temperature. [Inferred]

## 9. Observations (noted, not fixed) and open questions
1. **`gwdps.f` and `gwdc.f` appear unused.** They are in the CMake source list, but a search of the `.F90` files finds no call to `GWDPS` or `GWDC` (only their own definitions). *Inferred:* dead legacy code.
2. **`NCAR_NRDG` is read twice** with the same default (SetServices and in both branches of Initialize) (`GEOS_GwdGridComp.F90`).
3. **Copy-pasted description.** The `Run` and Initialize comment text refers to a "composite Moist Gridded Component" (`GEOS_GwdGridComp.F90`, Run, `DESCRIPTION`).
4. **A commented-out alarm test** (`ESMF_AlarmIsRinging`, `ESMF_AlarmRingerOff`) and an unused `alarm` variable in `Run` (3.3).
5. **A commented-out alternative default** for `RAYLEIGH_TAU1` (`172800.`) remains in Initialize.
6. **`NCAR_BKG_TAU`** starts at 6.4. It becomes `3.0 * sigma` if `NCAR_ET_FAC_DTDTM` or `NCAR_ET_FAC_WS300` is non-zero (defaults 0.75 and 1.25, so normally yes). That value is then the default of the user-settable resource `NCAR_BKG_TAU`, which overrides it (`GEOS_GwdGridComp.F90`, Initialize). *Unverified:* the physical meaning of the constant.

### Open questions for maintainers
- **Grid-name parsing.** `Initialize` derives the cubed-sphere size `imsize` by parsing the text of `AGCM.GRIDNAME` and fails on names it cannot parse (trap 3). Question for maintainers: is parsing the grid name, rather than reading the size from the grid itself, intended to remain the source of `imsize`?
