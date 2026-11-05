# GEOSmoist_GridComp (Moist): architecture

**Audience:** coding agents and developers working on the Moist gridded component (convection, cloud microphysics, aerosol activation).
**Describes:** the MAPL2 line, branch `develop`.
**Verified against:** `GEOSgcm_GridComp` @ `45e1f86f` (2026-11-04).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `0f2e0fb`; repository `GEOSgcm_GridComp`, mounted at `src/Components/@GEOSgcm_GridComp`): `src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp`.
Other fixtures (for example GEOSadas, GEOSldas) may place this repository elsewhere; paths inside this repository do not change.

Paths below are relative to `GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/` unless they start with `GEOSagcm_GridComp/` or `src/` (the latter are relative to the GEOSgcm fixture).
Parent: [Physics](../ARCHITECTURE.md). Siblings: [Turbulence](../GEOSturbulence_GridComp/ARCHITECTURE.md), [GWD](../GEOSgwd_GridComp/ARCHITECTURE.md), [Surface](../GEOSsurface_GridComp/ARCHITECTURE.md). Repository overview: [../../../docs/ARCHITECTURE.md](../../../docs/ARCHITECTURE.md).

## 1. Purpose
Moist computes every process that changes the phase of water: deep and shallow convection, large-scale (stratiform) condensation, cloud microphysics, precipitation and cloud fraction. It also transports "friendly" tracers (chemistry and aerosol species) through convection and activates aerosol into cloud droplets and ice crystals (`GEOS_MoistGridComp.F90`, header comment, Run).
The component is one driver (`GEOS_MoistGridCompMod`) plus one *interface module* per scheme. Schemes are chosen once, by resource keys, in `SetServices`.

## 2. Place in the hierarchy
- Parent is Physics, which adds Moist as a leaf child: `MOIST = MAPL_AddChild(GC, NAME='MOIST', SS=MoistSetServices, ...)` (`GEOSagcm_GridComp/GEOSphysics_GridComp/GEOS_PhysicsGridComp.F90`, SetServices).
- Moist has **no children**. The schemes are Fortran modules called directly, not gridded components (`GEOS_MoistGridComp.F90`, `use GEOS_*_InterfaceMod`).
- `CMakeLists.txt` builds one shared library from all scheme sources, including ones that cannot be selected (NSSL, section 5.1), plus `pyMoist/pyMoist/fortran/moist_workarounds.F90`, and installs the `pyMoist` Python tree (`CMakeLists.txt`, `srcs`, `install`).

## 3. Entry points
Registered with `MAPL_GridCompSetEntryPoint` (`GEOS_MoistGridComp.F90`, SetServices). There is **one Run phase**; Physics runs Moist once per Physics step.

- **SetServices**: reads the scheme keys (5.1), calls the `*_Setup` of each selected scheme (internal specs and timers), declares imports, exports and timers, then `MAPL_GenericSetServices`.
  - Microphysics `*_Setup` is called first: a comment says binary restarts expect `Q` first in `moist_internal_rst`, so the first `MAPL_AddInternalSpec` must come from microphysics.
  - Which imports exist depends on `WSUB_OPTION` and `USE_NCLOUD_CLIM` (`select case (WSUB_OPTION)`).
- **Initialize**: `MAPL_GenericInitialize`, then each selected scheme's `*_Initialize` (scheme resources, alarms), then the radar alarm `DBZ_RunAlarm` (`GEOS_MoistGridComp.F90`, Initialize).
- **Run**: 3.1.
- **Finalize**: `GFDL_1M_Finalize`, `UW_Finalize`, `GF_Finalize` (these matter only for the Python bridge), then `MAPL_GenericFinalize` (`GEOS_MoistGridComp.F90`, FINALIZE).

### 3.1 Run (`GEOS_MoistGridComp.F90`, RUN)
The body is under `if (ESMF_AlarmIsRinging(ALARM))`, `ALARM` being the component's `RUNALARM`. Its interval is the MAPL_Generic resource `DT:`, default `RUN_DT:` (`src/Shared/@MAPL/generic/MAPL_Generic.F90`, `handle_clock_and_main_alarm`). When it rings:

1. **Prologue**: surface type `SRF_TYPE`; pressure, height, mass; lower-tropospheric stability `LTS` and estimated inversion strength `EIS`; saved winds `UMST0`, `VMST0`; convective available potential energy `CAPE`, buoyancy `BYNCY`, and the diagnosed convective fraction `CNV_FRC` (`RUN`, `MOIST_PROLOGUE`).
2. **Tracer collection**: `CNV_Tracers_Init(TR)` on the import bundle `MTR` (7).
3. **Aerosol activation**: fills `NACTL`, `NACTI` (8).
4. **Convection**: `UW_Run` and `GF_Run` (or `RAS_Run`). If `SH_MD_DP` is true: shallow, then deep; otherwise deep, then shallow (`RUN`, `if (SH_MD_DP)`).
5. **Epilogue, part 1**: `CNV_MFC`, `CNV_MFD` from the `_DC` and `_SC` exports; inactive-convection detrainment (6.3).
6. **Microphysics**: `BACM_1M_Run`, `GFDL_1M_Run`, `THOM_1M_Run` or `MGB2_2M_Run`. The microphysics `*_Run` also fills `DTDTFRIC` (section 9) and, for GFDL_1M and MGB2_2M, calls `compute_radar_diagnostics` (gated by `DBZ_RunAlarm`) (`GEOS_GFDL_1M_InterfaceMod.F90`, GFDL_1M_Run; `Process_Library.F90`, `compute_radar_diagnostics`).
7. **Epilogue, part 2**: cloud-fraction exports; rain-out of humidity above 110% of saturation (`DTDT_ER`, `DQVDT_ER`, `ER_PRCP`); sums of per-process tendencies into `DUDT`, `DVDT`, `DTDT`, `DQDT`, `DQLDT`, `DQIDT`, `DQRDT`, `DQSDT`, `DQGDT`, `DQADT`, `DPDTMST`; precipitation totals and types; `*AFMOIST` exports. The lightning export `LFR_GCC` is set to zero here; the driver imports `HEMCO_FlashRate` but does not call it (`RUN`, "Lightning Exports"; grep).

When the alarm is **not** ringing, only `UAFMOIST`, `VAFMOIST`, `TAFMOIST`, `QAFMOIST`, `THAFMOIST`, `SAFMOIST`, `RH2` are refreshed (`RUN`, `else ! Alarm ringing`). Physics reads them every step (9).

Moist **changes its imports in place**: `T`, `U`, `V` (imports) and `Q`, `QLLS`, ... (internals) are updated by each scheme in turn, and later schemes see the result (`GEOS_UW_InterfaceMod.F90`, UW_Run: `T = T + DTDT_SC*MOIST_DT`; `GEOS_GF_InterfaceMod.F90`, GF_Run).

## 4. State
- **Imports**: 58 `MAPL_AddImportSpec` calls in SetServices; only some carry `AVERAGING_INTERVAL`/`REFRESH_INTERVAL`. Couplings are made by Physics (9). Notable: bundle `MTR` (tracers), state `AERO` (aerosol), `DTDT_BL`/`DQDT_BL` (boundary-layer tendencies for GF), `KPBL_SC` (boundary-layer level for UW).
- **Internal state**, per scheme in `*_Setup`: `Q`, `QLLS`, `QLCN`, `CLLS`, `CLCN`, `QILS`, `QICN` always; precipitating species and numbers by microphysics (5.2); `NACTL`, `NACTI`; GF adds `CNV_TR`; UW adds `CUSH`.
  - `Q` is `FRIENDLYTO` `DYNAMICS:TURBULENCE:CHEMISTRY:ANALYSIS`; `QLLS`, `QLCN`, `QILS`, `QICN` to `DYNAMICS:TURBULENCE`; `CLLS`, `CLCN` to `DYNAMICS` only (`GEOS_BACM_1M_InterfaceMod.F90`, BACM_1M_Setup, `FRIENDLIES`).
- **Exports**: about 559 `MAPL_AddExportSpec` calls in the driver, including the `*_SC` and `*_DC` convection exports. Exports are also work space between driver and schemes: a scheme fetches them with `MAPL_GetPointer(EXPORT, ..., ALLOC=.TRUE.)` and the driver reads them (`UMF_DC`, `UMF_SC`, `MFD_DC`, `MFD_SC`, `CNV_FRC`, `SRF_TYPE`) (`GEOS_UW_InterfaceMod.F90`, UW_Run, "Required Exports").
- **Module-level shared state**: `Process_Library.F90` (module `GEOSmoist_Process_Library`) holds scheme tunables (`SH_MD_DP`, `USE_AEROSOL_NN`, `WSUB_OPTION`, `CNV_FRACTION_MIN/MAX/EXP`, `SIGMA_DX`, ...) and the friendly-tracer array `CNV_Tracers`; `aer_cloud.F90` holds `AeroPropsNew`; interface modules keep their own (`USE_GF2020`, `USE_PYMOIST_UW`).
- **Restarts**:
  - `moist_internal_rst`: `Q` (`MAPL_RestartRequired`) and every internal spec not marked `RestartSkip` (MAPL default `MAPL_RestartOptional`, `MAPL_Generic.F90`, internal-spec routine): cloud species, the microphysics' precipitating species and numbers, UW's `CUSH`.
  - Not restarted (`RestartSkip`): `NACTL`, `NACTI`, GF's `CNV_TR` and `DSL__GF2020_LONS/LATS` (`*_Setup`).
  - Imports without `RestartSkip` (for example `KH`, `TKE`, `WQT`, `KPBL`, `DTDT_BL`, `DQDT_BL`) go to `moist_import_rst` (`src/Applications/@GEOSgcm_App/AGCM.rc.tmpl`, `MODELATMMOIST_IMPORT_RESTART_FILE`). *Inferred:* this restores lagged imports.
  - *Unverified:* the held convective tendencies (6) are exports and are not checkpointed.

## 5. Configuration that changes structure

### 5.1 Scheme-selection keys
Read in SetServices with `MAPL_GetResource(CF, ...)` from the AGCM resource file (built from `src/Applications/@GEOSgcm_App/AGCM.rc.tmpl`), checked with `_ASSERT` (an unsupported string aborts):

| Key | Values | Default | Source |
|---|---|---|---|
| `CONVPAR_OPTION` | `GF`, `RAS`, `NONE` (deep convection) | `GF` | `GEOS_MoistGridComp.F90`, SetServices |
| `SHALLOW_OPTION` | `UW`, `NONE` (shallow convection) | `UW` | same |
| `CLDMICR_OPTION` | `BACM_1M`, `GFDL_1M`, `THOM_1M`, `MGB2_2M` | `BACM_1M` | same |

- In the template, `CONVPAR_OPTION: @CONVPAR_OPTION` is filled by the setup scripts (`GF`, or `NONE` at some very high resolutions); `CLDMICR_OPTION` and `SHALLOW_OPTION` come from prefixed template lines (`@GFDL_1M_CLDMICR_OPTION`, `@SHALLOW_UWSHALLOW_OPTION: UW`, `@SHALLOW_SHOCSHALLOW_OPTION: NONE`) (`AGCM.rc.tmpl`, `set_resolution_params.csh`). `gcm_setup` offers only `BACM_1M` (default for 72 levels), `GFDL_1M` (default otherwise), `MGB2_2M`; `THOM_1M` is accepted by Moist but not offered (`gcm_setup`).
- The choice is made by `if (adjustl(OPTION)=="...")` chains in each of SetServices, Initialize, Run and Finalize; a new scheme needs an edit in all four.
- `NSSL_2M` has an interface module and is compiled but is **not** accepted by `CLDMICR_OPTION` and not called by the driver (`GEOS_MoistGridComp.F90`, SetServices, `LCLDMICR`). *Inferred:* work in progress.
- `RAS` is accepted, but the bodies of `RAS_Initialize` and `RAS_Run` are inside `#ifdef NODISABLE`, which no CMake file defines (grep), so selecting RAS runs no convection (`GEOS_RAS_InterfaceMod.F90`). *Inferred:* RAS is disabled.

### 5.2 Interface module per scheme
Each module has `*_Setup`, `*_Initialize`, `*_Run`; some also `*_Finalize`.

| Role | Option | Interface module | Science source | State added / notes |
|---|---|---|---|---|
| Deep convection | `GF` (Grell-Freitas) | `GEOS_GF_InterfaceMod.F90` | `ConvPar_GF2020.F90`, `ConvPar_GF_GEOS5.F90`, `ConvPar_GF_Shared.F90` (also `GF2020_2M_MicrophysicsMod.F90`) | `CNV_TR`; `USE_GF2020` picks 2020 or GEOS5 code (default 1, 0 at 72 levels); alarm `GF_RunAlarm`; `USE_PYMOIST_GF2020` |
| Deep convection | `RAS` | `GEOS_RAS_InterfaceMod.F90` | `ras.F90`, `rascnvv2_v.F90`, `ddf.F90`, `RASPARAMS.F90` | disabled (5.1) |
| Shallow convection | `UW` (University of Washington) | `GEOS_UW_InterfaceMod.F90` | `uwshcu.F90` | `CUSH`; alarm `UW_RunAlarm`; `USE_PYMOIST_UW` |
| Microphysics | `BACM_1M` (one-moment) | `GEOS_BACM_1M_InterfaceMod.F90` | `cloudnew.F90`, `aer_cloud.F90` | no precipitating species in the internal state |
| Microphysics | `GFDL_1M` | `GEOS_GFDL_1M_InterfaceMod.F90` | `gfdl_cloud_microphys.F90`, `gfdl_mp.F90` | `QRAIN`, `QSNOW`, `QGRAUPEL`; `USE_PYMOIST_GFDL1M` |
| Microphysics | `THOM_1M` (Thompson) | `GEOS_THOM_1M_InterfaceMod.F90` | `module_mp_thompson.F90`, `module_mp_radar.F90` | `QRAIN`, `QSNOW`, `QGRAUPEL`, `NCPL`, `NCPI`, `NRAIN` |
| Microphysics | `MGB2_2M` (two-moment) | `GEOS_MGB2_2M_InterfaceMod.F90` | `micro_mg3_0.F90`, `micro_mg_utils.F90`, `cldwat2m_micro.F90`, `GF2020_2M_MicrophysicsMod.F90` | `MGVERSION` (default 3) sets which precipitating species are friendly to Dynamics; `NCPL`, `NCPI`, `NRAIN`, `NSNOW`, `NGRAUPEL` |
| (not selectable) | `NSSL_2M` | `GEOS_NSSL_2M_InterfaceMod.F90` | `module_mp_nssl_2mom.F90`, `mp_nssl.F90` | 5.1 |

Sources: `*_Setup` internal specs in each module; `USE_GF2020`, `GF_RunAlarm` in `GEOS_GF_InterfaceMod.F90`, GF_Initialize; `MGVERSION` in `GEOS_MGB2_2M_InterfaceMod.F90`, MGB2_2M_Setup.

### 5.3 Other keys that change structure
- `SH_MD_DP`: order of shallow versus deep (3.1 step 4); GF also reads it (`ConvPar_GF2020.F90`, plume loop). Module default `.FALSE.` (`Process_Library.F90`). **Only** `GFDL_1M_Initialize` and `NSSL_2M_Initialize` read the key (default `.TRUE.`). **Trap:** with `BACM_1M`, `THOM_1M` or `MGB2_2M` an `SH_MD_DP` line in `AGCM.rc` is ignored and the order is always deep, then shallow (grep of `GEOS_*InterfaceMod.F90`).
- `GF_DT`, `UW_DT`: supercycling intervals (6). The template sets both to `@CONV_DT`.
- `WSUB_OPTION` (0 parameterized, 1 climatology `WSUB_CLIM`, 2 `Wnet`, 3 symbolic `Wnet`): in the driver it chooses which imports exist (`SetServices`). `USE_NCLOUD_CLIM` makes `NCPL_CLIM`, `NCPI_CLIM` imports and bypasses activation. The `WSUB_OPTION` branches in code exist only in `GEOS_MGB2_2M_InterfaceMod.F90`; BACM_1M and GFDL_1M ignore it, yet the imports are declared anyway.
- `USE_AEROSOL_NN` (re-read by BACM_1M, GFDL_1M, THOM_1M, MGB2_2M and NSSL_2M Initialize), `NN_FAC_*`, `NN_MIN_*`, `NN_MAX_*`, `NCCN_OCN/LND`, `NCIN_OCN/LND`: aerosol activation (8).
- `USE_TRACER_TRANSP`, `USE_TRACER_SCAVEN` (GF), `USE_TRACER_TRANSP_UW`: convective tracer transport (7).
- `CONVECTION_TRACER` (GF, default 0): if 1, the diagnostic tracer `CNV_TR` is also friendly to `DYNAMICS`; otherwise only to Moist itself (`GEOS_GF_InterfaceMod.F90`, GF_Setup).
- `DETRAIN_INACTIVE_CNV`, `TAU_DETRAIN_CNV` (6.3). `DEBUG_MST`, `DEBUG_TQ_ERRORS`: diagnostics.
- `USE_PYMOIST_UW`, `USE_PYMOIST_GF2020`, `USE_PYMOIST_GFDL1M`: run the Python (NDSL) version through `MAPL_pybridge_gcinit/gcrun_with_internal/gcfinalize` instead of the Fortran code. The Python tree is `pyMoist/` (own README and tests).
- **Resource-label search.** `MAPL_GetResource(CF, ...)` (SetServices) goes through `MAPL_GetResource_config_scalar` without `component_name`, so only the exact label is searched. The `MAPL_GetResource(MAPL, ...)` form (Initialize and the `*_Initialize` routines) also tries `<name>_<label>`, `<type>_<label>`, `<label>` and `<name>.<label>` (`src/Shared/@MAPL/base/MAPL_Resource.F90`, `get_labels_with_prefix`; `MAPL_Generic.F90`, `MAPL_GetResourceFromMAPL_scalar`; separator `.` in `src/Shared/@MAPL/shared/Constants/InternalConstants.F90`).

## 6. Supercycling and held tendencies

### 6.1 What the code does
Deep (GF) and shallow (UW) convection may run less often than Moist.
- Each scheme creates its alarm in `*_Initialize`: `GF_RunAlarm` (interval `GF_DT`) and `UW_RunAlarm` (`UW_DT`), both defaulting to the Moist interval (`GEOS_GF_InterfaceMod.F90`, GF_Initialize; `GEOS_UW_InterfaceMod.F90`, UW_Initialize).
- In `GF_Run` and `UW_Run` the scheme call and the computation of the tendency exports (`DTDT_DC`, `DQVDT_DC`, ... / `DTDT_SC`, ...) are inside `if (alarm_is_ringing)`. The code after it runs on **every** Moist step: it adds the tendency exports (times `MOIST_DT`) to `T`, `U`, `V`, `Q`, `QLCN`, `QICN`, `CLCN` (UW also to `QLLS`, `QILS` through `QLSUB_SC`, `QLENT_SC`, ...), then fills negatives (`GF_Run`, "Add tendencies to the moist import state"; `UW_Run`, "Apply tendencies in a single fused loop").
- So between alarms the tendency is held: the exports are fetched with `ALLOC=.TRUE.` and written only at an alarm. *Inferred:* the same holds for the other exports written only inside the alarm branch (mass fluxes `UMF_*`, `MFD_*`, `CNPCPRATE`), which the driver's epilogue reads every step.
- Microphysics is **not** supercycled: it runs every Moist step.
- The schemes get their own interval (`GF_DT = DT_R8` from the alarm; `compute_uwshcu_inv(..., UW_DT, ...)`), while the held tendency is applied with `MOIST_DT`.

### 6.2 Traps
- **Alarm phase is not anchored.** Both alarms are created with `RingTime = currentTime` at Initialize; a 00Z `ringTime` is computed but not used. So the scheme runs on the first step after any start or restart and then every `GF_DT`/`UW_DT`. A restart at a time that is not a multiple of the interval shifts the phase. *Inferred:* this can break restart reproducibility when `GF_DT > RUN_DT`.
- *Unverified:* held tendencies are exports and are not in `moist_internal_rst`; since the alarm rings on the first step after a restart, they are recomputed before use.
- **Tracers are not held** (7): the scheme updates them at the alarm with the full interval; nothing re-applies them.
- *Inferred:* `GF_DT` and `UW_DT` should be multiples of `MOIST_DT`; otherwise the applied total differs from the scheme's interval.

### 6.3 Inactive-convection detrainment
If `DETRAIN_INACTIVE_CNV > 0`, wherever the layer mass flux `0.5*(CNV_MFC(L)+CNV_MFC(L+1))` is below that value the driver moves convective condensate and cloud fraction to the large-scale species at rate `(1 - MFC/DETRAIN_INACTIVE_CNV) * DT_MOIST/TAU_DETRAIN_CNV` (`TAU_DETRAIN_CNV` default 3600 s). `DETRAIN_INACTIVE_CNV` defaults to `0.0`, but to `5.0e-4` for `GFDL_1M` (`GEOS_MoistGridComp.F90`, Initialize and RUN, "Mass fluxes").

## 7. Convective transport of friendly tracers
- Physics fills Moist's import bundle `MTR` in its Initialize: `ESMF_StateGet(GIM(MOIST),'MTR',BUNDLE)` then `MAPL_GridCompGetFriendlies(GC,"MOIST",BUNDLE)` and `MAPL_GridCompGetFriendlies(GCS,"MOIST",BUNDLE)`. It also calls `MAPL_TerminateImport(... 'MTR' ..., CHILD=MOIST)`, so nothing else supplies it (`GEOSagcm_GridComp/GEOSphysics_GridComp/GEOS_PhysicsGridComp.F90`, SetServices, Initialize).
- A field is in `MTR` if its owner declared `FRIENDLYTO` containing `MOIST`; MAPL records this as the field attribute `FriendlyToMOIST` (`src/Shared/@MAPL/generic/MAPL_Generic.F90`, `'FriendlyTo'//...`). To find producers grep `FRIENDLYTO.*MOIST` (and the last column of `*StateSpecs.rc`), not `FriendlyToMOIST`. Producers in this fixture, all under `GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSchem_GridComp/`: GOCART2G (`DU2G`, `SS2G`, `SU2G`, `CA2G`, `NI2G` and `Aerosol` StateSpecs), `@QuickChem/OH_GridComp` StateSpecs, `@ACHEM/ACHEM_GridCompMod.F90`, `@MAM/MAM_GridCompMod.F90`, `@TR/TR_GridCompMod.F90` (unless `GOCART_convection` is true), `@RRG/RRGMod.F90`, GEOS-Chem (`GEOSCHEMchem_GridComp/@geos-chem/Interfaces/GEOS/geos_interface.F90`, `.../GCHP/Chem_GridCompMod.F90`). `@CARMA` declares `DYNAMICS:TURBULENCE` only (grep), so its tracers are not transported by Moist.
- Each Run the driver calls `CNV_Tracers_Init(TR)` (`Process_Library.F90`). It counts fields whose `FriendlyToMOIST` is true. The **first time** it allocates the module array `CNV_Tracers`, reads per-tracer attributes (`ScavengingFractionPerKm`, `SetofHenryLawCts`, and for GEOS-Chem washout `SetofKcScalFactors`, `IsWetDep`, `RetentionFactor`, `LiqAndGas`, `ConvFacI2G`, `OnlineCLDLIQ`, `OnlineVUD`, `UseGOCART`, `GOCARTfTempThreshold`) and stores a **pointer to the bundle data** in `CNV_Tracers(F)%Q`. Later calls only assert the count is unchanged.
- The schemes read and write `CNV_Tracers(:)%Q` directly: GF in `GF2020_INTERFACE` / `GF_GEOS5_Interface` (`Q += DT_moist*SRC_CHEM`, `DT_moist = GF_DT`, floor `mintracer`) when `USE_TRACER_TRANSP == 1` (`ConvPar_GF2020.F90`); UW in `compute_uwshcu_inv` when `USE_TRACER_TRANSP_UW == 1` (`uwshcu.F90`). With `USE_PYMOIST_*`, `CNV_Tracers_To_SOA`/`CNV_Tracers_To_AOS` copy between structure forms (`pyMoist/pyMoist/fortran/moist_workarounds.F90`).
- Tracers therefore change only at convection alarms and by the full scheme interval. If shallow runs first (`SH_MD_DP`), deep convection sees tracers already changed by UW.
- Traps:
  - **Saved pointers.** *Inferred:* anything that reallocates or reorders the bundle after the first Run breaks `CNV_Tracers(:)%Q` silently, since only the count is checked (`CNV_Tracers_Init`).
  - **Order of friendlies.** The bundle is filled in Physics Initialize (above). *Inferred:* a tracer that becomes `FRIENDLYTO MOIST` later is not seen.
  - With `CONVPAR_OPTION=NONE` deep convection does not transport tracers, but `UW` still does if `USE_TRACER_TRANSP_UW` is 1 (default) (`GEOS_UW_InterfaceMod.F90`, UW_Initialize).
  - Moist does not export a tracer tendency. Physics forms increments by differencing the bundle around the Moist call: `Initialize_IncBundle_run` and `Compute_IncBundle` in `GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSchem_GridComp/Shared/Chem_Shared/Bundle_IncrementMod.F90`, and the 2D mass-weighted `Initialize_IncMBundle_run`, `Compute_IncMBundle` in `GEOSagcm_GridComp/GEOSphysics_GridComp/MBundle_IncrementMod.F90` (called from `GEOS_PhysicsGridComp.F90`, Run).

## 8. Aerosol activation through the `AERO` state
- `AERO` is an ESMF **state** import (`DATATYPE = MAPL_StateItem`). Physics connects it from Chemistry: `MAPL_AddConnectivity(SHORT_NAME='AERO', DST_ID=MOIST, SRC_ID=CHEM)` (`GEOS_PhysicsGridComp.F90`, SetServices).
- After the tracer collection the driver gets `NACTL` (droplet number) and `NACTI` (ice number) one of three ways (`GEOS_MoistGridComp.F90`, RUN, `AERO_ACTIVATE`):
  1. `USE_NCLOUD_CLIM`: copy imports `NCPL_CLIM`, `NCPI_CLIM` (ExtData; `WSUB_ExtData.yaml`).
  2. `USE_AEROSOL_NN` (module default true): call `Aer_Activation` (`aer_actv_single_moment.F90`), scale by `NN_FAC_LIQ/ICE`, clip to `[NN_MIN, NN_MAX]`. Vertical velocity is `W`, or derived from `OMEGA` when `W` is identically zero.
  3. Otherwise constants from `NCCN_LND/OCN`, `NCIN_LND/OCN` weighted by `FRLAND`.
- `Aer_Activation` reads attributes of the `AERO` state (`number_of_aerosol_modes`, `aerosol_modes`, field names for pressure, temperature, land fraction, which it **writes**: `PLE`, `T`, `FRLAND`). For each mode it sets `aerosol_mode`, runs the state's ESMF method `aerosol_activation_properties`, and copies number, dry size, width, hygroscopicity (and optionally density, dust/soot/organic fractions) into the module array `AeroPropsNew(n)` (`aer_cloud.F90`). The method is registered by the chemistry components (grep `aerosol_activation_properties` under `GEOSchem_GridComp`).
- Activated number per mode uses updraft `w + sqrt(TKE)`; modes with hygroscopicity above 0.4 sum into `NACTL`, dust/soot modes into `NACTI`; diagnostics `NWFA`, `NIFA` (`aer_actv_single_moment.F90`, `Aer_Activation`, `GetActFrac`).
- `AeroPropsNew` is also consumed by GF (`GF_Run` passes it) and MGB2_2M (`MGB2_2M_Run`), so activation must run before convection and microphysics, as it does.
- Traps:
  - With zero modes `Aer_Activation` returns after setting only `NWFA`, `NIFA`; `NACTL`/`NACTI` keep earlier values (Observation 3).
  - `AeroPropsNew` is filled only on the `USE_AEROSOL_NN` path. *Inferred:* on the other two paths GF and MGB2 read stale or empty data.
  - The neural network `Wnet` (`Wnet.F90`, `Wnet_weights.rc`) is reached through `aer_cloud_init(use_wnet)` and is enabled only by `MGB2_2M` with `WSUB_OPTION=2`; BACM_1M and GFDL_1M call `aer_cloud_init(use_wnet=.false.)` (`GEOS_MGB2_2M_InterfaceMod.F90`, MGB2_2M_Initialize; `aer_cloud.F90`).

## 9. What Physics expects from Moist
Physics calls Moist once per step with a plain `ESMF_GridCompRun` (no phase), after GWD and before the phased Surface, Chemistry, Turbulence calls and Radiation (`GEOS_PhysicsGridComp.F90`, Run). Then:
- **Before the call**: Physics resets the increment bundles and copies the turbulence liquid and ice tendencies into Moist's `DQLDTTRB` (Observation 2). It also allocates `DTDTFRIC` of Moist's exports (`alloc=.true.`) (`Run`).
- **Lagged import.** `DTDT_BL`/`DQDT_BL` (boundary-layer tendencies for GF) are written into Moist's import state by Physics after Moist, so GF uses the previous step's values (`GEOS_PhysicsGridComp.F90`, Run, "Boundary Layer Tendencies for GF"; `GEOS_GF_InterfaceMod.F90`, GF_Run). Physics terminates the connectivity to them (`MAPL_TerminateImport`). *Inferred:* Moist restarts them from `moist_import_rst`.
- **Imports supplied by Physics** (`GEOS_PhysicsGridComp.F90`, SetServices, "Moist Imports"): `U`, `V`, `T`, `TAUOROX/Y` from GWD; turbulence moments, `KH`, `KPBL`, `KPBL_SC`, `TKE`, `ALH`, `KM`, `RI` from Turbulence; `TS`, `SNOMAS`, `FRLAND`, `FROCEAN`, `FRLANDICE`, `FRACI`, `T2M`, `Q2M`, `TA`, `QA`, `SH`, `EVAP` from Surface; `AERO` from Chemistry; `RADLW`, `RADSW` from Radiation. `PLE`, `ZLE`, `W`, `OMEGA`, `PREF` come from Dynamics through the AGCM (`GEOSagcm_GridComp/GEOS_AgcmGridComp.F90`, SetServices).
- **Exports Physics reads from Moist** (`GEX(MOIST)`): `DUDT`, `DVDT`, `DWDT`, `DTDT`, `DTDTFRIC`, `DPDTMST`, `DQDT`, `DQLDT`, `DQIDT`, `DQRDT`, `DQSDT`, `DQGDT`, and `UAFMOIST`, `VAFMOIST`, `TAFMOIST`, `THAFMOIST`, `QAFMOIST`, `SAFMOIST`. The `*AFMOIST` fields carry the post-Moist state to Surface, Turbulence and Chemistry when `SYNCUV`/`SYNCTQ` (*Inferred:* Physics switches that pass the post-Moist winds, and temperature and humidity, to the other children) are on.
- **Exports connected to siblings** (`SRC_ID=MOIST`): to Turbulence `QV`, `QLTOT`, `QITOT`, `QRTOT`, `QSTOT`, `QGTOT`, `FCLD`, `WTHV2`, `WQT_DC`, `EIS`; to Radiation `QV`, `QL`, `QI`, `QR`, `QS`, `QG`, `FCLD`, `RL`, `RI`, `RR`, `RS`, `RG`; to Surface `PCU`, `PLS`, `SNO`, `ICE`, `FRZR`, `TPREC`, `CN_PRCP`; to GWD `Q`, `CNV_FRC`, `DTDT_micro`, `DTDT_DC`; to Chemistry a long list of cloud, precipitation, mass-flux and CAPE fields.
- **Tendency convention.** `DTDT` is the sum of per-process rates times layer pressure thickness (`Pa K s-1`); it **excludes friction heating**, which comes separately as `DTDTFRIC` from the microphysics `*_Run` (`dissipative_ke_heating`, or `IKEX`/`IKEX2` in BACM_1M's `#else` branch since `KEDISS_MICRO` is undefined) and which Physics adds to the temperature budget (`GEOS_MoistGridComp.F90`, RUN, "Export Total Moist Tendencies"; `GEOS_*_InterfaceMod.F90`, `*_Run`; `GEOS_PhysicsGridComp.F90`, Run, `TICU`). The totals also omit the negative-value fills and the cloud-fraction change from 6.3 (Observation 4).
- **Dynamics and the condensate split.** In `FV_Run` (`GEOSagcm_GridComp/GEOSsuperdyn_GridComp/@FVdycoreCubed_GridComp/FV_StateMod.F90`, FVdycoreCubed @ `4deed15`), "Pull Tracers" adds `QLLS+QLCN` into total liquid, `QILS+QICN` into total ice and `CLLS+CLCN` into total cloud fraction (for `nwat` 3 or 6, chosen from the tracer count and hydrostaticity) and also copies `QLCN`, `QICN`, `CLCN` as extra tracers. "Push Tracers (Dynamic Mass Separation for Option A)" clips each advected convective tracer to its total, writes it back, sets `QLLS = total - QLCN` (likewise ice and cloud fraction) and floors at zero. Only the hydrostatic `nwat=1` path advects the six species independently. So the large-scale species are not independently conserved by Dynamics. *Inferred:* the search found no cloud-fraction-ratio re-split of the partition.
- **Negative-value fills.** Moist fills negatives itself (GF and UW: `DQVDT_FILL_DC`, `*_FILL_SC`; microphysics; `FILLQ2ZERO` calls in `*_Run`), and Physics also applies `FILLQ2ZERO` to the water species before and after the Moist call (`DQ*DT_FILL_AFMST` exports) (`GEOS_PhysicsGridComp.F90`, Run).

## 10. Dependencies
- Libraries (`CMakeLists.txt`, `esma_add_library`): `GEOS_Shared`, `GMAO_mpeu`, `MAPL`, `Chem_Shared`, `Chem_Base`, ESMF, BLAS, LAPACK, OpenMP; Python bridge `MAPL_pybridge_*` in MAPL.
- Chemistry repository: `Lightning_mod` in `GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSchem_GridComp/Shared/Chem_Shared/Lightning_mod.F90` (imported by the driver). The `AERO` activation method comes from the chemistry components.
- Shared helpers: `GEOS_UtilsMod` (saturation functions such as `GEOS_QSAT`, `GEOS_DQSAT`) and `GEOS_RadarMod`, from `GEOS_Shared`. `FILLQ2ZERO` (and `FILLQ2ZERO_OLD`) are public in `Process_Library.F90`, not in `GEOS_UtilsMod`.
- Compiler workarounds in `CMakeLists.txt`: some sources use `-O1` with GCC 15 or later; `GEOS_BACM_1M_InterfaceMod.F90` uses `-O1` with Intel Release.
- Build and test: with the GEOSgcm fixture; Python/NDSL tests in `pyMoist/tests` (see `pyMoist/README.md`).

## 11. Known traps (index)
Scheme `if` chains in four routines (5.1); first internal spec order (3); `SH_MD_DP` read by two microphysics only (5.3); unanchored alarms (6.2); cached tracer pointers (7); order inside Run is part of the contract because schemes pass data through exports and module globals (`AeroPropsNew`, `CNV_Tracers`, `SH_MD_DP`) (3.1, 8); in-place modification of imports (3.1).

## 12. Decisions
### Schemes are chosen by strings in SetServices
- Pattern: three resource strings are tested with `if` chains in SetServices, Initialize, Run and Finalize; unselected schemes declare no state (`GEOS_MoistGridComp.F90`). [Verified]
- Likely reason: keeps unused schemes' state out of restarts and History; avoids runtime polymorphism. [Inferred]
- Question for maintainers: what is the intended status of the compiled-but-inert schemes (`NSSL_2M`, `RAS`) — are they expected to be completed and used, or kept only for reference?

### Convection supercycles with held tendencies; microphysics does not
- Pattern: separate alarms for GF and UW; tendencies held in exports and applied each step (6.1). [Verified]
- Likely reason: deep convection is the most expensive scheme and is run at a longer interval at high resolution. [Inferred]
- Question for maintainers: are the unanchored alarms (6.2) intended?

### Scheme-to-scheme communication through exports (`ALLOC=.TRUE.`)
- Pattern: schemes allocate exports on demand and the driver reads them back (4). [Verified]
- Likely reason: avoids a shared derived type between independently developed schemes. [Inferred]

## 13. Observations (noted, not fixed)
1. **`WSUB_OPTION` has several defaults, one dead.** The driver reads it with default 3 (`GEOS_MoistGridComp.F90`, SetServices) and that value chooses the imports. `MGB2_2M_Setup` assigns `MGB2_WSUB_OPTION` with default 1 but never uses it; `MGB2_2M_Initialize` re-reads it with default 3 before use (`GEOS_MGB2_2M_InterfaceMod.F90`). The default 1 is dead.
2. **`DQIMST` reads the wrong label.** Physics fetches the import `'DQLDTTRB'` twice; the ice tendency should come from `'DQIDTTRB'` (`GEOS_PhysicsGridComp.F90`, Run, "Pass Turbulence QL, QI tendencies to Moist Import state"). Moist never reads either import (grep), so there is no effect today.
3. **`Aer_Activation` leaves `NACTL`/`NACTI` unassigned when there are zero modes** (`aer_actv_single_moment.F90`, `n_modes == 0`); the driver then scales and clips them.
4. **Total tendency exports omit the negative-value fills and the cloud-fraction change** from inactive-convection detrainment (which moves condensate between `QLCN`/`QLLS`, so `DQLDT`, `DQIDT` totals are unaffected, but `DQADT` is) (`GEOS_MoistGridComp.F90`, RUN).
5. **`DWDT` is declared but never computed in Moist** (`SHORT_NAME = 'DWDT'`; grep of the Moist directory), while Physics reads it as `WIM` (`GEOS_PhysicsGridComp.F90`, Run). Physics allocates it with `alloc=.true.`. *Unverified:* the initial value is zero (not checked in MAPL's `ALLOC`).
6. **Duplicate timer name.** `---CLDMACRO` is added twice (`GEOS_MoistGridComp.F90`, SetServices).
7. **Stale comments.** Initialize refers to "child Dynamics"; `pyMoist/README.md` mentions `MGB2_3M`; the module header says `CLDMACRO==MGB2_2M` while the key is `CLDMICR_OPTION`.
8. **`RAS` and `NSSL_2M` are accepted or compiled but inert** (5.1).
9. **Two implementations of `FILLQ2ZERO`.** `Process_Library.F90` (with `FILLQ2ZERO_OLD`) and a private copy in `GEOS_PhysicsGridComp.F90`; Physics does not use the Moist version.

### Open questions for maintainers
Today the GF and UW alarms start at the current time of Initialize; a 00Z ring time is computed but unused (6.2). Should the convection alarms be anchored to a fixed reference time, so that a restart at an arbitrary time continues the same supercycle phase?
