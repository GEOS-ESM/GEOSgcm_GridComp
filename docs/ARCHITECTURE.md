# GEOSgcm_GridComp Architecture

**Audience:** coding agents and developers working in the `GEOSgcm_GridComp` repository.
**Describes:** the MAPL2 line, branch `develop`.
**Verified against:** `GEOSgcm_GridComp` @ `45e1f86f` (2026-11-04).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `0f2e0fb`).
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

Paths are relative to this repository's root unless they start with `src/` or
`@` (fixture-relative). Statements without a label are Verified and carry a
`(path, routine-or-key)` citation. *Inferred:* and *Unverified:* mark the
weaker cases. Fixture-level documents (build, mepo, application) are owned by
the fixture; this file does not repeat them.

## 1. Placement

- The fixture lists this repository as `GEOSgcm_GridComp` with
  `local: ./src/Components/@GEOSgcm_GridComp` and
  `sparse: ./config/GEOSgcm_GridComp.sparse` (`components.yaml`, `GEOSgcm_GridComp`).
- The sparse file excludes the Held-Suarez test model
  (`GEOSagcm_GridComp/GEOS_AgcmSimpleGridComp.F90`, `GEOSagcm_GridComp/GEOShs_GridComp`)
  (`config/GEOSgcm_GridComp.sparse` in the fixture); those files are tracked here
  but absent from a fixture checkout.
- Parent: the fixture `GEOSgcm`, whose `src/CMakeLists.txt` and
  `src/Components/CMakeLists.txt` reach this repository
  (`esma_add_subdirectory (GEOSgcm_GridComp)`). Fixture documentation:
  [`../../../../docs/ARCHITECTURE.md`](../../../../docs/ARCHITECTURE.md) when this repository sits in the fixture.
- `GEOSgcm.x` uses this repository's root component as the MAPL root:
  `MAPL_Cap('GCM', ROOT_SetServices, ...)` with `use GEOS_GcmGridCompMod`
  (`src/Applications/@GEOSgcm_App/GEOSgcm.F90`); `ROOT_NAME: GCM`, `ROOT_CF: AGCM.rc`
  (`src/Applications/@GEOSgcm_App/CAP.rc.tmpl`).
- Children (links down): the nested repositories of section 6.

## 2. Directory map

| Directory | What it is | Module | Component document |
|---|---|---|---|
| `GEOS_GcmGridComp.F90` | Root composite "GCM": atmosphere, ocean, waves, IAU/replay, digital filter, trajectories | `GEOS_GcmGridCompMod` | section 3 here |
| `GEOSagcm_GridComp/` | Atmosphere composite (Agcm) with `GEOSsuperdyn_GridComp/` (dynamics; `GEOSdatmodyn_GridComp/` is a single-column driver per its module header; nested `@FVdycoreCubed_GridComp`) | `GEOS_AgcmGridCompMod`, `GEOS_SuperdynGridCompMod` | [`GEOSagcm_GridComp/ARCHITECTURE.md`](../GEOSagcm_GridComp/ARCHITECTURE.md) |
| `GEOSagcm_GridComp/GEOSphysics_GridComp/` | Physics composite: GWD, moist, turbulence, surface (here); chemistry, radiation (nested) | `GEOS_PhysicsGridCompMod` | section 5 and [`GEOSphysics_GridComp/ARCHITECTURE.md`](../GEOSagcm_GridComp/GEOSphysics_GridComp/ARCHITECTURE.md) |
| `GEOSogcm_GridComp/` | Ocean composite: ocean, sea ice, ocean radiation, ocean biogeochemistry; nested `@GEOS_OceanGridComp` | `GEOS_OgcmGridCompMod` | not yet documented |
| `GEOSwgcm_GridComp/` | Wave-model composite; model chosen by `wave_model:` in `WGCM.rc`; `GEOSumwm_GridComp/`, `GEOSwavewatch_GridComp/` | `GEOS_WgcmGridCompMod` | not yet documented |
| `GEOSdataatm_GridComp/` | Data atmosphere (forcing from files) used in place of the Agcm when `USE_DATA_ATM4OCN:` is true; runs a `SURFACE` child | `GEOS_DataAtmGridCompMod` | not yet documented |
| `GEOSmkiau_GridComp/` | "make IAU" (incremental analysis update) component `AIAU`; digital filter initialization `ADFI`; optional pyMLINC interface | `GEOS_mkiauGridCompMod`, `DFI_GridCompMod` | not yet documented |
| `GEOSgigatraj_GridComp/` | Trajectory component, built only with `BUILD_WITH_GIGATRAJ`; wraps nested `@GigaTraj` | `GEOS_GigatrajGridCompMod` | not yet documented |
| `.github/`, `.circleci/` | CI workflows, code owners, issue templates | | |
| `GCM_Preamble.tex` | LaTeX preamble (*Inferred:* for documents extracted from source with ProTeX; it defines `\ptx` as `protex`) | | |

Module names come from the `use ..., only: ...SetServices` lines of the parent
(`GEOS_GcmGridComp.F90`, module header; `GEOSagcm_GridComp/GEOS_AgcmGridComp.F90`,
module header). Held-Suarez: `GEOSagcm_GridComp/CMakeLists.txt` builds
`GEOS_AgcmSimpleGridComp.F90` with `GEOShs_GridComp` when that file exists,
otherwise `GEOS_AgcmGridComp.F90`.

## 3. `GEOS_GcmGridComp` (root composite)

Source: `GEOS_GcmGridComp.F90` (about 3000 lines). No separate component file.

### 3.1 Children and how they are added
`SetServices` adds children with `MAPL_AddChild` (`GEOS_GcmGridComp.F90`, `SetServices`):

| Name | SetServices from | When |
|---|---|---|
| `AGCM` | `GEOS_AgcmGridCompMod` | `USE_DATA_ATM4OCN:` false (default) |
| `DATAATM` | `GEOS_dataatmGridCompMod` | `USE_DATA_ATM4OCN:` true (takes the `AGCM` index; AIAU, ADFI, GIGATRAJ are not added) |
| `GIGATRAJ` | `GEOS_GigatrajGridCompMod` | `#ifdef HAS_GIGATRAJ` and not data atmosphere |
| `AIAU` | `GEOS_mkiauGridCompMod` | not data atmosphere |
| `ADFI` | `DFI_GridCompMod` | not data atmosphere |
| `OGCM` | `GEOS_OgcmGridCompMod` | always |
| `WGCM` | `GEOS_WgcmGridCompMod` | `USE_WAVES:` /= 0 |

When `REPLAY_MODE: Regular` and `REPLAY_HISTORY_RC:` is present, a `History_GCM` child is added
to a private parent component, not to this one (`GEOS_GcmGridComp.F90`, internal `history_setservice`).

### 3.2 Entry points
- `SetServices`: registers `Initialize` and `Run` with `MAPL_GridCompSetEntryPoint`; reads
  structural resources (3.5); adds children; borrows exports from the Agcm (`MAPL_AddExportSpec
  ... CHILD_ID = AGCM`); makes connections (3.3); creates the private state `GCM_state`
  (`ESMF_UserCompSetInternalState`); calls `MAPL_GenericSetServices`; adds timers.
- `Initialize`: creates the child grids (`MAPL_GridCreate` for AGCM, OGCM, and WGCM when waves
  are on); builds the atmosphere and ocean exchange grids from the tile file (`TILING_FILE:`,
  default `tile.data`; `MAPL_LocStreamCreate`, `MAPL_ExchangeGridSet`); creates the replay
  alarms (3.4); calls `MAPL_GenericInitialize` (initializes all children); then builds the
  atmosphere-to-ocean and ocean-to-atmosphere transforms (`MAPL_LocStreamCreateXform`, kept as
  `XFORM_A2O`, `XFORM_O2A`) and allocates the exchange exports (`AllocateExports`).
- `Run` (single phase): the replay block (3.4, `Regular` only), then each step: Gigatraj phase 1
  (if built), run `AGCM` (or `DATAATM`), Gigatraj phase 2, run `ADFI`, `RUN_OCEAN`, `RUN_WAVES`
  (if waves). The Gcm does not call `MAPL_GenericRunChildren`; it runs each child with
  `ESMF_GridCompRun` (`GEOS_GcmGridComp.F90`, `Run`).
- Finalize: the Gcm registers none itself. `MAPL_GenericSetServices` registers
  `MAPL_GenericFinalize` for any component that has not set a finalize method
  (`src/Shared/@MAPL/generic/MAPL_Generic.F90`, `register_generic_entry_points`).

### 3.3 Connections
- Agcm to AIAU: `PHIS, AK, BK, U, V, TV, PS, DELP, O3PPMV, TS, AREA`; renamed pairs `Q->QV`,
  `PPBL`, `TROPP_BLENDED`; and `QLTOT, QITOT, QRTOT, QSTOT, QGTOT`.
- AIAU to Agcm: `DUDT, DVDT, DTDT, DQVDT, DPEDT, DO3DT, DTSDT` (the analysis increments).
- Agcm to ADFI: `U_DGRID, V_DGRID, PT, PE, Q, OX`. OGCM to Agcm: `SURFSTATE` when
  `USE_CICE_Thermo:` is 2. All are `MAPL_AddConnectivity` calls in `SetServices`.
- Atmosphere-ocean exchange is not by connectivity. The OGCM imports (`TAUXW`, `PENUVR`, `LWFLX`,
  `SNOW`, ... and ice-state names depending on `USE_CICE_Thermo:` and `SEAICE_THICKNESS_EXT_DATA:`)
  are removed with `MAPL_TerminateImport` ("explicitly connected through exchange grid transforms
  Run", `SetServices`) and filled in `RUN_OCEAN` by `DO_A2O*` (`MAPL_LocStreamTransform`) from the
  Agcm's `SALTWATER` skin export state (`Initialize`: `skinname`, `expSKIN`, `impSKIN`); the
  reverse direction is `DO_O2A*`. `RUN_OCEAN` acts only when the OGCM alarm rings (`alarmOcn`,
  named `<OGCM name>_Alarm`); `BYPASS_OGCM:` skips it.
- Waves: `RUN_WAVES` copies surface fields (`U10M ... RH2M`) and ocean fields (`UW, VW, TW`) into the
  WGCM import (`DO_A2W`, `DO_O2W`), runs `WGCM`, and copies `CHARNOCK` (and `SHFX_SPRAY`, `LHFX_SPRAY`
  when `USE_SEA_SPRAY:` /= 0) back (`DO_W2A`). WGCM imports not sent to ExtData (the MAPL external
  data component) are terminated (`SetServices`).

### 3.4 IAU and replay hooks
- `SetServices` reads and writes back defaults for `ASSIMILATION_CYCLE:` (21600),
  `CORRECTOR_DURATION:`, `MKIAU_FREQUENCY:`, `REPLAY_FILE_FREQUENCY:`, `MKIAU_REFERENCE_TIME:`,
  `REPLAY_FILE_REFERENCE_TIME:` (a source comment calls the defaults "MERRA-2"); the Agcm reads them
  ([`GEOSagcm_GridComp/ARCHITECTURE.md`](../GEOSagcm_GridComp/ARCHITECTURE.md) section 3.4).
  It also defaults `CONVPAR_OPTION:` to `RAS` (relaxed Arakawa-Schubert) and `AERO_PROVIDER:` to
  `GOCART.data` for use by several children.
- `REPLAY_MODE:` is `NoReplay` (default), `Exact` or `Regular`. `DUAL_OCEAN:` asserts `Regular`.
- `Initialize` creates clock alarms in every mode: `startReplay`, `replayCycle`, `stopReplay`,
  `ReplayMKIAU`, `ReplayShutOff`, `ExactReplay09`, `RegularReplay09`, `PredictorActive`; computes
  `PREDICTOR_DURATION` and compares it with a user value. For `Regular` replay it calls
  `MAPL_AddRecord(..., MAPL_Write2Ram)` on AGCM and OGCM, and on WGCM only when waves are on, so
  their states can be restored.
- `Run`, `Regular` replay, when `startReplay` rings and `ReplayShutOff` does not:
  1. Run `AIAU` phase 2 (`Clear`: zero the IAU forcing).
  2. Turn the `PredictorActive` alarm on; save ring time and ringing state of every clock alarm;
     turn `PredictorAlarm` on; remember the current time.
  3. If `ReplayMKIAU` rings: run `AIAU` phase 1 (make the increments) and, if requested, checkpoint the
     AGCM import state (`MKIAU_CHECKPOINT_FILE:`, `MAPL_CheckpointState`).
  4. If `PREDICTOR_DURATION` is nonzero, loop until `stopReplay` rings: run ExtData, `AIAU` phase 2,
     `AGCM`, `RUN_OCEAN` (only if `USE_DATASEA:` /= 0, or as ocean phase 2 for `DUAL_OCEAN:`), advance
     the clock, `MAPL_GenericRunCouplers` for the AGCM, `AIAU` phase 1 when `ReplayMKIAU` rings.
     Waves are not run in this loop (the call is under `#if (0)`).
  5. Rewind the clock (`ESMF_DIRECTION_REVERSE`) to the remembered time, restore the alarms, call
     `MAPL_GenericRefresh` on AGCM, OGCM, WGCM (if waves), run ExtData again, then take the ordinary step.
  When `ReplayShutOff` rings, `AIAU` phase 2 is run and recording is disabled (`MAPL_DisableRecord`).
- The replay block is entered only for `Regular`; for `Exact` the Agcm reads the increment files itself.
- ExtData is not a child of the Gcm. MAPL's Cap stores its component and export state in the Gcm's
  `ExtData_state` internal state (`src/Shared/@MAPL/gridcomps/Cap/MAPL_CapGridComp.F90`,
  `initialize_gc`); the Gcm creates the empty wrapper in `SetServices`.
- With `Regular` replay and `PREDICTOR_DURATION` 0 (or MKIAU at the start), the file
  `AIAU_IMPORT_RESTART_FILE:` must exist (`Initialize`).

### 3.5 Resources that change structure
`USE_DATA_ATM4OCN:`, `USE_WAVES:`, `USE_SEA_SPRAY:`, `USE_CICE_Thermo:` (0, 1, 2; with
`CICE_N_ICE_CATEGORIES:`, `CICE_N_ICE_LAYERS:`), `USE_OCEANOBIOGEOCHEM:`, `USE_DATASEA:` (default 1),
`SEAICE_THICKNESS_EXT_DATA:`, `DUAL_OCEAN:`, `BYPASS_OGCM:`, `REPLAY_MODE:`, `TILING_FILE:`
(`GEOS_GcmGridComp.F90`, `SetServices` / `Initialize`). They come from `AGCM.rc` (template
`src/Applications/@GEOSgcm_App/AGCM.rc.tmpl`).

### 3.6 State and restarts
Private state `T_GCM_STATE` (transforms, alarms, wave-coupling states, route handles) wrapped as
`GCM_WRAP` under the name `GCM_state`. The Gcm has no `MAPL_AddInternalSpec` call, so it checkpoints
nothing of its own. Its exports are borrowed from the Agcm (about 26 `MAPL_AddExportSpec ...
CHILD_ID = AGCM` calls).

## 4. Agcm and Superdyn

Summary; detail in [`GEOSagcm_GridComp/ARCHITECTURE.md`](../GEOSagcm_GridComp/ARCHITECTURE.md).
`GEOS_AgcmGridComp` adds `SUPERDYNAMICS`, `PHYSICS` and `ORBIT` children (`SetServices`). Its `Run` applies
analysis increments or bias tendencies (run type free-run, predictor, corrector, forecast), runs
dynamics phase 1, physics, copies the physics tendencies into the dynamics imports, then runs dynamics
phase 2. `GEOS_SuperdynGridComp` selects the dynamical core with `DYCORE:` (default `FV3`) and registers
two run phases (`Run`, `RunAddIncs`).

## 5. Physics subtree (map level)

Composite: `GEOSagcm_GridComp/GEOSphysics_GridComp/GEOS_PhysicsGridComp.F90` (about 4000 lines);
component file [`GEOSphysics_GridComp/ARCHITECTURE.md`](../GEOSagcm_GridComp/GEOSphysics_GridComp/ARCHITECTURE.md).
Children are added in this order; a source comment says chemistry must precede surface so chemistry can
fill the `AERO_DP` bundle before the surface reads its restart (`GEOS_PhysicsGridComp.F90`, `SetServices`):

| Child name | Directory | Where | Component file |
|---|---|---|---|
| `GWD` | `GEOSgwd_GridComp/` | this repository | [`ARCHITECTURE.md`](../GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSgwd_GridComp/ARCHITECTURE.md) |
| `MOIST` | `GEOSmoist_GridComp/` (with `pyMoist/`) | this repository | [`ARCHITECTURE.md`](../GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/ARCHITECTURE.md) |
| `TURBULENCE` | `GEOSturbulence_GridComp/` | this repository | [`ARCHITECTURE.md`](../GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSturbulence_GridComp/ARCHITECTURE.md) |
| `CHEMISTRY` | `@GEOSchem_GridComp/` | nested repository | fixture path `src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSchem_GridComp` |
| `SURFACE` | `GEOSsurface_GridComp/` (land, lake, landice, route, saltwater, `Shared/`, `Utils/`) | this repository | [`ARCHITECTURE.md`](../GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/ARCHITECTURE.md) |
| `RADIATION` | `@GEOSradiation_GridComp/` | nested repository | fixture path `src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSradiation_GridComp` |

Run order (`GEOS_PhysicsGridComp.F90`, `Run`; children are skipped when `DOPHYSICS:` /= 1): `GWD`;
`MOIST`; stage-1 wind and temperature/humidity synchronization (`SYNCUV:`, `SYNCTQ:`) and negative-water
fill; `SURFACE` phase 1; `CHEMISTRY` phase 1; `TURBULENCE` phase 1; second fill and stage-2
synchronization; `SURFACE` phase 2; `TURBULENCE` phase 2; third fill and synchronization of the radiation
and chemistry inputs; `CHEMISTRY` phase 2; `RADIATION`. Detail in the Physics component file.

Tracer bundles: Physics exports `TRADV`, `H2ORTRI`, `MTRI`, `MCHEMTRI`, `TRANA` as `MAPL_BundleItem` and
fills them in `Initialize` from the children's friendly fields (section 7, item 5).

## 6. Nested repositories (separate git repositories in the fixture)

Each path is `src/Components/@GEOSgcm_GridComp/<path>` in the fixture (`components.yaml`). Several are
listed in `.gitignore` files here (for example `GEOSogcm_GridComp/.gitignore`).

| Path inside this repository | Repository (`components.yaml` key) |
|---|---|
| `GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSchem_GridComp` | `GEOSchem_GridComp` (and its own nested repositories) |
| `GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSradiation_GridComp` | `GEOSradiation_GridComp` (nested `RRTMGP`) |
| `GEOSagcm_GridComp/GEOSsuperdyn_GridComp/@FVdycoreCubed_GridComp` | `FVdycoreCubed_GridComp` (nested `fvdycore`) |
| `GEOSogcm_GridComp/@GEOS_OceanGridComp` | `GEOS_OceanGridComp` (nested `mom`, `mom6`, `stochastic_physics`, `mit`) |
| `GEOSogcm_GridComp/GEOSseaice_GridComp/CICE_GEOSPlug/@cice6`, `@sis2` | `cice6` (nested `icepack`), `sis2` |
| `GEOSwgcm_GridComp/GEOSwavewatch_GridComp/@ww3`, `GEOSwgcm_GridComp/GEOSumwm_GridComp/@umwm` | `ww3`, `umwm` |
| `GEOSgigatraj_GridComp/@GigaTraj` | `GigaTraj` |
| `GEOSmkiau_GridComp/pyMLINC/@geos_state_bias` | `geos_state_bias` |
| `GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/Utils/Raster/preproc/topography/@ncar_topo` | `Topo` |

Everything else in section 2 is in this repository. The ocean model plug-ins (`MOM6_GEOSPlug`,
`MOM_GEOS5PlugMod`, `MIT_GEOS5PlugMod`, `GEOSdatasea_GridComp`) live inside the nested
`GEOS_OceanGridComp`. Some children are loaded as shared objects with
`MAPL_AddChild(name, 'setservices_', parentGC, sharedObj)` (CICE6:
`GEOSogcm_GridComp/GEOSseaice_GridComp/GEOS_SeaIceGridComp.F90`, `SetServices`, `SEAICE_NAME:`).
Change nested repositories in their own repository.

`.github/CODEOWNERS`: the GCM gatekeepers own the repository; the ocean team owns `GEOSogcm_GridComp/` and
`GEOSdataatm_GridComp/`; the land, landice, surface-preprocessing and software-modernization teams own named
subdirectories; the CMake team owns every `CMakeLists.txt`.

## 7. Common MAPL2 patterns (one example each)

1. **SetServices skeleton.** Register entry points (`MAPL_GridCompSetEntryPoint`, once per run
   phase), add children and specs, make connections, call `MAPL_GenericSetServices`, add timers
   (`GEOSagcm_GridComp/GEOS_AgcmGridComp.F90`, `SetServices`; leaf:
   `GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSgwd_GridComp/GEOS_GwdGridComp.F90`, `SetServices`).
2. **`MAPL_AddChild`.** Returns an integer index kept in a module variable and used for `GCS(i)`,
   `GIM(i)`, `GEX(i)` from `MAPL_Get` (`GEOS_AgcmGridComp.F90`, `SetServices`, `Run`).
3. **`MAPL_AddConnectivity` / `MAPL_TerminateImport`.** Connect a child export to a sibling import by
   `SHORT_NAME`, or `SRC_NAME`/`DST_NAME` when names differ, with `SRC_ID`/`DST_ID`.
   `MAPL_TerminateImport` marks child imports as deliberately unconnected; where it is used the source
   comments say the parent fills them by hand (`GEOS_AgcmGridComp.F90`, `SetServices`: "We Terminate these
   IMPORTS which are manually filled"; `GEOS_GcmGridComp.F90`, `SetServices`).
4. **Borrowed exports.** `MAPL_AddExportSpec(GC, SHORT_NAME=..., CHILD_ID=...)` exposes a child's export as the
   parent's own (`GEOS_GcmGridComp.F90`, `SetServices`, `'AK'`).
5. **Friendly tracers.** `FRIENDLYTO = 'A:B:...'` on an internal or export spec tags the field as a tracer
   for the named components and also adds it to the exports (`src/Shared/@MAPL/generic/MAPL_Generic.F90`,
   `MAPL_StateAddInternalSpec`). Example tag list: `FRIENDLIES%QV = "DYNAMICS:TURBULENCE:CHEMISTRY:ANALYSIS"`
   (`GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/GEOS_BACM_1M_InterfaceMod.F90`). A parent collects
   them with `MAPL_GridCompGetFriendlies(GC, "DYNAMICS", BUNDLE)` (`GEOS_PhysicsGridComp.F90`, `Initialize`);
   `MAPL_CopyFriendliness` copies tags (`GEOS_GcmGridComp.F90`, `RUN_OCEAN`).
6. **Bundles.** A bundle export is declared with `DATATYPE = MAPL_BundleItem` and no dimensions
   (`GEOS_PhysicsGridComp.F90`, `SetServices`: `TRANA`; the Agcm's `TRADVI` also gives `DIMS = MAPL_DimsHorzVert`,
   `GEOS_AgcmGridComp.F90`, `SetServices`). It is connected like any item (`SRC_NAME='TRADV'`, Agcm
   `SetServices`). *Inferred:* bundles hold the children's fields by reference, so an in-place update through
   `ESMFL_BundleGetPointerToData` (`DO_Friendly`) changes the owning child's field.
7. **ACG (automatic code generator).** A component lists its state in `<Name>_StateSpecs.rc` (blocks
   `category: IMPORT|EXPORT|INTERNAL`). `mapl_acg(<target> <specs>.rc IMPORT_SPECS EXPORT_SPECS INTERNAL_SPECS
   GET_POINTERS DECLARE_POINTERS)` in `CMakeLists.txt` runs `MAPL_GridCompSpecs_ACG.py` and writes
   `<Name>_Import___.h`, `_Export___.h`, `_Internal___.h`, `_GetPointer___.h`, `_DeclarePointer___.h` into the
   build directory; the Fortran source `#include`s them
   (`GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSgwd_GridComp/GWD_StateSpecs.rc`, its `CMakeLists.txt`,
   `GEOS_GwdGridComp.F90`; generator `src/Shared/@MAPL/cmake/mapl_acg.cmake`). Among the components of this
   repository present in a fixture checkout, GWD is the only one using it; the rest call `MAPL_AddImportSpec` /
   `MAPL_AddExportSpec` / `MAPL_AddInternalSpec` directly. Nested repositories use it (for example
   `GEOS_OceanGridComp`, `CARMA`, `GOCART2G`). Edit the `.rc`, not the generated headers.
8. **Private state.** `ESMF_UserCompSetInternalState` with a wrapper type holding a pointer (`GEOS_GcmGridComp.F90`:
   `GCM_state`, `T_GCM_STATE`).

## 8. Build structure

Fixture-level build (mepo, presets, `parallel_build.csh`, MAPL, `esma_add_library`) is described in the
fixture documents. What is specific here:
- Each composite directory's `CMakeLists.txt` has `esma_set_this()`, a list `alldirs` of children, and two
  branches: if the composite's own source file exists, `esma_add_library(${this} SRCS <composite>.F90
  SUBCOMPONENTS ${alldirs} DEPENDENCIES ...)`; otherwise `esma_add_subdirectories(${alldirs})` (root
  `CMakeLists.txt`; `GEOSagcm_GridComp/CMakeLists.txt`; `GEOSagcm_GridComp/GEOSphysics_GridComp/CMakeLists.txt`).
  The target name is the directory name without `@` (`esma_set_this`); `OVERRIDE` renames it (sea ice:
  `GEOS_Seaice_GridComp`, `GEOSogcm_GridComp/GEOSseaice_GridComp/CMakeLists.txt`).
- A missing child directory becomes a generated stub component (`esma_add_library` in the fixture's
  `@cmake/esma_support`); `rename_<dir>` variables map a directory to a different module name
  (`GEOSagcm_GridComp/GEOSsuperdyn_GridComp/CMakeLists.txt`).
- Options, default OFF: `BUILD_WITH_GIGATRAJ` (adds `GEOSgigatraj_GridComp`, defines `HAS_GIGATRAJ` in the
  root, Agcm and Superdyn), `BUILD_WITH_PYMLINC` (defines `HAS_PYMLINC`; `GEOSmkiau_GridComp/pyMLINC.cmake`),
  `BUILD_MIT_OCEAN` (ocean and sea-ice composites) (`CMakeLists.txt`, `GEOSogcm_GridComp/CMakeLists.txt`).
- `mapl_acg` follows the `esma_add_library` of the same target.
- Tests and CI: no unit-test target for the composite components was found (`grep add_test`). The `pyMoist`
  DSL code has tests under `GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/pyMoist/tests/`
  (`numerical`, `scripts`, `translate_tests`). GitHub workflows (`.github/workflows/`): `workflow.yml` (fixture
  build tests; its matrix lists GEOSldas but passes `fixture-repo: GEOS-ESM/GEOSgcm`), `ndsl-checks.yml`,
  `enforce-labels.yml` (PR labels), `push-to-develop.yml`. CircleCI (`.circleci/config.yml`) builds the fixture
  and runs short AMIP and coupled GCM runs (job comments).

## 9. Known traps

- Order is set by the parent's `Run`, not by the order of `MAPL_AddChild`, in the Gcm, Agcm, Physics and
  Superdyn composites. Adding a child there requires adding its `ESMF_GridCompRun` call and phase. Exception:
  the Ogcm in its ordinary (not dual-ocean) mode calls `MAPL_GenericRunChildren`, which runs children in the
  order they were added (`GEOSogcm_GridComp/GEOS_OgcmGridComp.F90`, `Run`).
- A child with two run phases is run by the parent with `PHASE=1`/`PHASE=2`; *Unverified:* ESMF runs phase 1
  when no phase is given (the Physics `Run` relies on this for `GWD`, `MOIST`, `RADIATION`).
- In the existing code, imports that a parent fills by hand are terminated with `MAPL_TerminateImport` (Agcm,
  Gcm, Superdyn, Physics, Ogcm).
- Replay in the Gcm rewinds the clock and restores every clock alarm (`ESMF_ALARMLIST_ALL`), but only ring time
  and ringing state; other alarm properties (for example interval, enabled flag) are not restored
  (`GEOS_GcmGridComp.F90`, `Run`: `ESMF_AlarmGet` / `ESMF_AlarmSet` with `ringTime=`, `ringing=`).
- Module-level variables in `GEOS_GcmGridComp.F90` (`AGCM`, `OGCM`, `DO_WAVES`, `DO_OBIO`, ...) hold per-run
  configuration; *Inferred:* the component assumes one instance per process.
- Sea ice, ocean and wave components choose their implementation at `SetServices` time from resources
  (`SEAICE_NAME:`, `USE_DATASEA:`, `DUAL_OCEAN:`, `wave_model:`); structure differs per run.

## 10. Decisions

See [`decisions.md`](decisions.md).

## 11. Observations (noted, not fixed)

- `GEOS_GcmGridComp.F90`, `Run`: `PredictorIsActive` is assigned only in the `Regular` replay branch but
  tested with `ESMF_AlarmIsCreated` on every step. *Inferred:* this relies on ESMF's default initialization of
  `ESMF_Alarm` (not a defect).
- `GEOS_GcmGridComp.F90`, `Run`: the comment "fix the predictor alarm (huh?)" at the `PredictorAlarm` ringer-on.
- `GEOS_GcmGridComp.F90`, `Initialize` and `Run`: with `HAS_GIGATRAJ` defined and `USE_DATA_ATM4OCN:` true,
  `GCS(gigatraj)` is used although `SetServices` adds the `GIGATRAJ` child only in the non-data-atmosphere branch.
  *Inferred:* `gigatraj` is then unset (needs checking; the combination may be unsupported).
- `_ASSERT(..., 'needs informative message')` placeholders remain in `GEOS_GcmGridComp.F90`,
  `GEOS_AgcmGridComp.F90`, `GEOS_OgcmGridComp.F90` and `GEOS_PhysicsGridComp.F90`.
- `GEOSmkiau_GridComp/IAU_GridCompMod.F90` and `DynVec_GridComp.F90` are compiled
  (`GEOSmkiau_GridComp/CMakeLists.txt`, `srcs`) but no `use` of either was found under `src/` of the GEOSgcm
  fixture. *Inferred:* unused; other fixtures (GEOSadas) were not checked.
- The `!DESCRIPTION` of `SetServices` in `GEOS_GcmGridComp.F90` calls the component the "PhysicsGcm GC" and
  lists only AGCM and OGCM as children; the code adds more (3.1).

### Not yet documented
`GEOS_OgcmGridComp` (with sea ice, ocean radiation), `GEOS_WgcmGridComp`, `GEOS_DataAtmGridComp`,
`GEOSmkiau_GridComp` (`AIAU` with its `Run` and `Clear` phases, `ADFI`), `GEOS_GigatrajGridComp`.
Section 3 describes only how the Gcm connects to them.

### Open questions
The two maintainer questions are in [`decisions.md`](decisions.md): whether the Physics stage order is meant to
be fixed, and whether replay logic should sit outside the Gcm.
