# GEOSagcm_GridComp: atmosphere composite and superdynamics

**Audience:** coding agents and developers working on `GEOS_AgcmGridComp` and `GEOS_SuperdynGridComp`.
**Describes:** the MAPL2 line, branch `develop`.
**Verified against:** `GEOSgcm_GridComp` @ `45e1f86f` (2026-11-04).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `0f2e0fb`).
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

Paths are relative to the `GEOSgcm_GridComp` repository root. Repository
overview: [`../docs/ARCHITECTURE.md`](../docs/ARCHITECTURE.md). Physics child:
[`GEOSphysics_GridComp/ARCHITECTURE.md`](GEOSphysics_GridComp/ARCHITECTURE.md).

## 1. Purpose

`GEOS_AgcmGridComp` (about 4900 lines, `GEOS_AgcmGridComp.F90`) combines
superdynamics (dynamics), physics and an orbit component into the atmosphere
model. It also applies analysis increments and bias tendencies to the
dynamics and physics tendencies (incremental analysis update, IAU) and reads
replay files. `GEOS_SuperdynGridComp` (about 850 lines,
`GEOSsuperdyn_GridComp/GEOS_SuperdynGridComp.F90`) wraps the dynamical core
and splits a dynamics step into two run phases around physics.

## 2. Place in the hierarchy

- Parent: `AGCM` child of the root Gcm (`GEOS_GcmGridComp.F90`, `SetServices`;
  [`../docs/ARCHITECTURE.md`](../docs/ARCHITECTURE.md) section 3).
- Children of the Agcm, added with `MAPL_AddChild` in `SetServices`
  (`GEOS_AgcmGridComp.F90`): `SUPERDYNAMICS` (module variable `SDYN`),
  `PHYSICS` (`PHYS`), `ORBIT` (`ORB`, `MAPL_OrbGridCompMod` from MAPL). Under
  `#ifdef SCM` the dynamics child is named `SCMDYNAMICS` and the module
  `GEOS_singcolGridCompMod` is used. With `-DHS` the physics module is
  `GEOS_hsGridCompMod` (`#define` at the top of the file); the Held-Suarez
  sources are not in a fixture checkout (see [`../docs/ARCHITECTURE.md`](../docs/ARCHITECTURE.md) section 1).
- Children of Superdyn (`GEOS_SuperdynGridComp.F90`, `SetServices`): `DYN` and
  optionally `ADV`, chosen by `DYCORE:` (default `FV3`): `FV3` adds `DYN` from
  `@FVdycoreCubed_GridComp`; `FV3+ADV` adds `DYN` and `ADV` (`AdvCore_GridCompMod`);
  `DATMO` adds `DYN` from `GEOS_DatmoDynGridCompMod`. `ADV` is -1 when absent.
- Nested repository below Superdyn: `@FVdycoreCubed_GridComp` (fixture path
  `src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSsuperdyn_GridComp/@FVdycoreCubed_GridComp`).
  `GEOSdatmodyn_GridComp/` (this repository) is a driver for a single-column
  model with profile data (its module header).

## 3. Entry points

### 3.1 `GEOS_AgcmGridComp`
- `SetServices`: registers `Initialize` and one `Run` (no extra phases); adds import
  specs `DUDT, DVDT, DTDT, DPEDT, DQVDT, DO3DT, DTSDT` (analysis increments, filled by
  the AIAU child through the Gcm); internal specs of the same names (bias tendencies,
  `FRIENDLYTO = trim(COMP_NAME)`); about 100 `MAPL_AddExportSpec` calls, most of them
  borrowed from `SDYN` or `PHYS` with `CHILD_ID`; the bundle export `TRADVI`; children;
  connections (3.3); allocates the private states `UPD_STATE` and `IAU_COEFFS`; then
  `MAPL_GenericSetServices`; then timers `INITIALIZE`, `RUN`.
- `Initialize`: `initialize_python_bridge` (from MAPL); `MAPL_GenericInitialize` (initializes
  children); `Initialize_IncBundle_init` for the tracer-increment bundle; allocates the
  physics tendency exports `DUDT ... DPEDT`; fills `PHIS`, `SGH`, `VARFLT` from topography
  files (`GEOS_TopoGet`, from `GEOS_Shared`) and copies them into the children's imports;
  marks those three exports `MAPL_InitialRestart` so that a connected import is not
  overwritten by its default (comment in `Initialize`); creates alarms
  `PredictorAlarm` (offset `PREDICTOR_OFFSET:`, default 21600 s, disabled when not in DAS mode),
  `ExactReplay` (when `REPLAY_MODE:` is `Exact`, or `Regular` with a long predictor; interval
  `REPLAY_INTERVAL:`) and `4DIAUalarm` (`4DIAU_FREQUENCY:`, default 3600 s).
- `Run` (one step), in order:
  1. Decide the run type (3.4), update the dynamics and physics imports with analysis
     increments or bias tendencies (3.4).
  2. `Initialize_IncBundle_run`; run `ORBIT`; `Pack_Chem_Groups`.
  3. Run `SDYN` phase 1 (dynamics); `Compute_IncBundle` (tracer advection increments in
     `TRADVI`); `Unpack_Chem_Groups`.
  4. Run `PHYS` phase 1 (physics).
  5. `DO_UPDATE_PHY`: copy physics exports `DUDT, DVDT, DWDT, DPEDT, DTDT` into the `SDYN`
     imports; barrier.
  6. Run `SDYN` phase 2 (`RunAddIncs`), then `MAPL_GenericRunCouplers(STATE, SDYN, ...)`.
  7. Final negative-tracer fill for the friendly advection tracers (`TRADV`), mass and
     energy diagnostics, tropopause diagnostics.
  The only `MAPL_GenericRunCouplers` call in the Agcm is the one for `SDYN` (`Run`); no call was
  found for `ORB` or `PHYS`. `MAPL_GenericRunCouplers(STATE, CHILD, ...)` runs every coupler whose
  source is `CHILD` (loop over `STATE%CCcreated(CHILD,J)`, `src/Shared/@MAPL/generic/MAPL_Generic.F90`),
  but most of the connections an agent would look at here never create a coupler in the first
  place: when `MAPL_AddConnectivity` wires a child's import to a sibling's export, the two specs
  are compared field by field (short name, `DIMS`, `VLOCATION`, `HALOWIDTH`, `PRECISION`, and, when
  both sides set them, the accumulation and coupling intervals;
  `src/Shared/@MAPL/generic/MAPL_Generic.F90`, `MAPL_WireComponent`;
  `src/Shared/@MAPL/generic/VarSpec.F90`, `MAPL_VarSpecEQ`). An exact match is wired as a direct
  connection (`MAPL_Reconnect`, comment "this a direct connection ... no additional averaging is
  needed" in `MAPL_WireComponent`): the import becomes the same field as the export, so `PHYS`
  reads what `SDYN` wrote earlier in the same step with no run-time step of its own. Only a pair
  whose specs differ gets a coupler component (`ESMF_CplCompCreate`, named `SRC_2_DST`, logged as
  "Coupler needed for ..." in `MAPL_WireComponent`), and the `MAPL_GenericRunCouplers(STATE, SDYN,
  ...)` call after `SDYN` phase 2 only runs couplers that were created this way; it has no effect
  on connections that were already direct.
  All sixteen `PHYS` imports declared directly in `GEOS_PhysicsGridComp.F90` (`U, V, W, DZ, T, S,
  ZLE, PLE, DTDTDYN, DQVDTDYN, QV_DYN_IN, T_DYN_IN, U_DYN_IN, V_DYN_IN, PLE_DYN_IN, AREA`) were
  checked against the `SDYN` export of the same name and match exactly, so these are direct
  connections: `SDYN` only promotes them unchanged from the dycore
  (`GEOSsuperdyn_GridComp/GEOS_SuperdynGridComp.F90`, `MAPL_AddExportSpec` with `CHILD_ID = DYN`),
  and the dycore declares them with the same `DIMS`/`VLOCATION`/`PRECISION`
  (`GEOSsuperdyn_GridComp/@FVdycoreCubed_GridComp/DynCore_GridCompMod.F90`). The two bundles
  connected the same way, `TRADV` (`PHYS` export to the `TRADV` import of `SDYN`'s `DYN` and `ADV`
  children) and `SATORB` (`ORB` export to the `SATORB` import of the Satsim child of `PHYS`), also
  match (same `DATATYPE = MAPL_BundleItem`, same `DIMS` where given) and are direct as well.
  *Unverified:* whether this also holds for the remaining renamed or promoted connections in 3.3
  (`TA, QA, US -> UA, VS -> VA, WSPD_STABLE300M, SPEED, PREF, TROPP_BLENDED -> TROPP, PLK, PV,
  TROPK_BLENDED -> TROPK, OMEGA, PKE, PLE -> PLEINST, LHOBS, SHOBS, DXC`): `PHYS` does not declare
  these imports itself, they are promoted up from imports declared deep in its children, and
  whether each chain ends in a matching spec was not checked. A run log line "Coupler needed for
  SUPERDYNAMICS and PHYSICS" (or "... and ORBIT") would show whether any of these need a coupler.

  The other renamed connection in 3.3, `T -> T_avg24, PLE -> PLE_avg24, ZLE -> ZLE_avg24,
  TROPP_BLENDED -> TROPP_avg24` (`GEOS_AgcmGridComp.F90`, the second `MAPL_AddConnectivity` call
  in `SetServices`), has no consumer in this checkout. The Chemistry composite's `SetServices`
  does not add `QuickChem` as a child: its `MAPL_AddChild` calls list `CHEMENV, HEMCO, PCHEM,
  ACHEM, GOCART, GOCART.data, GOCART2G, RRG, GAAS, H2O, STRATCHEM, GMICHEM, CARMA, GEOSCHEMCHEM,
  MATRIX, MAM, MAM.data, TR, DNA`, no `QuickChem`
  (`@GEOSchem_GridComp/GEOS_ChemGridComp.F90`), and its `CMakeLists.txt` `alldirs` list
  (`GEOSpchem_GridComp, GEOSCHEMchem_GridComp, H2O_GridComp, DNA_GridComp, HEMCO_GridComp, RRG`)
  does not include `@QuickChem` either (`@GEOSchem_GridComp/CMakeLists.txt`). The only
  declarations of `T_avg24`, `PLE_avg24` and `ZLE_avg24` anywhere in this repository are in
  `@GEOSchem_GridComp/@QuickChem/OH_GridComp/OH_GridCompMod.F90` (`ADD_IMPORT_24_3DC` /
  `ADD_IMPORT_24_3DE` macro calls) and `OH_StateSpecs.rc`; `TROPP_avg24` is not declared anywhere
  in the repository. Because `OH_GridComp` is never instantiated, `PHYS` never carries a
  `T_avg24` (or `PLE_avg24`, `ZLE_avg24`, `TROPP_avg24`) import up from any child, so this
  explicit `MAPL_AddConnectivity` entry is caught by the relaxation pass that runs once per
  composite right after wiring: `setup_children` calls `connect%checkReq` immediately after
  `MAPL_WireComponent` (`src/Shared/@MAPL/generic/MAPL_Generic.F90`, `setup_children`, comment
  "Relax connectivity for non-existing imports"); `checkReq` looks up the destination name in the
  destination child's import specs (`MAPL_VarSpecGetIndex`), and when it is absent there but the
  source name is present in the source child's export specs, it sets `conn%notRequired = .true.`
  (`src/Shared/@MAPL/generic/VarConn.F90`, `checkReq`). `checkUnused` then skips `notRequired`
  connections instead of logging them as unsatisfied (same file, `checkUnused`). Net effect: the
  four `*_avg24` connections are silently dropped in this checkout — no coupler is created, no
  field is wired, and no error or log line is produced for them.
- No Finalize of its own.

### 3.2 `GEOS_SuperdynGridComp`
- `SetServices`: registers `Initialize`, `Run`, and `RunAddIncs` (second run phase); reads
  `DYCORE:` and `SCM_SL:`; adds children (section 2); import specs for the summed tendencies
  `DUDT, DVDT, DWDT, DTDT`; about 50 borrowed exports from `DYN` (`MAPL_AddExportSpec`
  with `CHILD_ID = DYN`); for `FV3+ADV` the connectivity `MFX, MFY, CX, CY, PLE0, PLE1` from
  `DYN` to `ADV`; `MAPL_TerminateImport` of `DYN` imports `DUDT, DVDT, DTDT, DWDT`
  ("will be connected explicitly in the run method"); `MAPL_GenericSetServices`.
- `Initialize`: `MAPL_GenericInitialize`; allocates the `Q` export of `DYN` (needed to compute
  virtual-temperature tendencies from physics).
- `Run` (phase 1): adds the Superdyn imports `DUDT, DVDT, DTDT, DWDT` to the `DYN` imports,
  runs `DYN` phase 1, then `ADV` if present.
- `RunAddIncs` (phase 2): runs `DYN` phase 2 (adds the physics increments).
- Citations: `GEOSagcm_GridComp/GEOSsuperdyn_GridComp/GEOS_SuperdynGridComp.F90`.
- Phase numbering: the n-th `MAPL_GridCompSetEntryPoint(..., ESMF_METHOD_RUN, ...)` call gives
  phase n (`src/Shared/@MAPL/generic/MAPL_Generic.F90`, `MAPL_GridCompSetEntryPoint`:
  `MAPL_AddMethod(META%phase_run)`).

### 3.3 Connections (Agcm `SetServices`)
- `SDYN` to `PHYS`, with renames (`TROPP_BLENDED -> TROPP`, `PLE -> PLEINST`, `US -> UA`): `U, V,
  TH, T, ZLE, PS, TA, QA, PLE, W, PREF, PV, S, PLK, PKE, OMEGA, SPEED, DZ, ...`, `AREA`, `DXC`, and
  the dynamics-input diagnostics `DTDTDYN, DQVDTDYN, PLE_DYN_IN, U_DYN_IN, V_DYN_IN, QV_DYN_IN,
  T_DYN_IN`. A separate `MAPL_AddConnectivity` call renames `T, PLE, ZLE, TROPP_BLENDED` to
  `T_avg24, PLE_avg24, ZLE_avg24, TROPP_avg24`; these four have no consumer in this checkout and
  are dropped by MAPL (3.1).
- `PHYS` to `SDYN`: the bundle `TRADV`. `ORB` to `PHYS`: the bundle `SATORB`.
- `LHOBS, SHOBS` from `SDYN` to `PHYS` only when `SCM_SL:` /= 0.
- `MAPL_TerminateImport` (hand-filled imports): `SDYN` `DUDT, DVDT, DWDT, DTDT, DPEDT, DQVANA,
  DQLANA, DQIANA, DQRANA, DQSANA, DQGANA, DOXANA, PHIS, VARFLT`; `PHYS` `VARFLT, PHIS, SGH, DTSDT`
  (source comment: "We Terminate these IMPORTS which are manually filled").

### 3.4 IAU, analysis and replay (inside Agcm `Run`)
- Run type, an integer `FREERUN` (0), `PREDICTOR` (1), `CORRECTOR` (2), `FORECAST` (3):
  `IS_FCST:` /= 0 gives `FORECAST`; not in DAS mode (data assimilation system) gives `FREERUN`;
  otherwise `PredictorAlarm` ringing gives `PREDICTOR`, else `CORRECTOR`. DAS mode is on when
  `IMPORT_RESTART_FILE:` exists or `REPLAY_MODE:` is `Regular` or `Exact` (`Initialize`, `Run`).
  The `Regular` and `Exact` branches adjust the type and the last-corrector flag.
- `DO_UPDATE_ANA3D` / `DO_UPDATE_ANA2D` fill the `SDYN` (and `PHYS` for `DTSDT`) imports by type:
  `FREERUN` zero; `PREDICTOR` the bias tendency (internal state); `FORECAST` the scaled bias
  tendency; `CORRECTOR` the analysis increment (the Agcm import) times the IAU coefficient,
  with an optional pressure constraint (`CONSTRAIN_DAS:`). Bias update: `BIAS(n+1) = ALPHA*IAU(n)
  + BETA*BIAS(n)` (comment in `Run`; `ALPHA:`, `BETA:`, `ALPHAQ:`, `BETAQ:`, `ALPHAO:`, `BETAO:`).
- `get_iau_coeff`: constant `1/TAUANL`, or digital-filter weights when `IAU_DIGITAL_FILTER:` is
  `YES` (default), computed by `dfi_coeffs` and kept in `IAU_COEFFS`; uses `MKIAU_FREQUENCY:`,
  `MKIAU_RingDate:`, `MKIAU_RingTime:` (set by the Gcm, `../docs/ARCHITECTURE.md` section 3.4).
- `DO_Friendly`: adds increments to the friendly tracers in the `TRANA` bundle (water species,
  `OX`), preventing negative values (`FILL_Friendly`). `TRANA` must hold exactly two fields
  (`_ASSERT(NumFriendly==2, ...)`).
- 4D IAU: when `4DIAUalarm` rings in a corrector step, `update_ainc_` rebuilds the increments from
  the file `AINC_FILE:` (with nudging options `NUDGE_STATE:`, `SCLINC:`).
- Replay (`Exact` or `Regular` with a long predictor): `GET_REPLAY_TIME`, `MAPL_GetCurrentFile`
  (`REPLAY_FILE:`, `REPLAY_FILE09:`, or the MKIAU checkpoint template), then
  `MAPL_ESMFStateReadFromFile` into the Agcm import state. Alarms used (created by the Gcm):
  `ReplayShutOff`, `replayCycle`, `ExactReplay09`, `RegularReplay09`, `startReplay`; and
  `ExactReplay` (created in the Agcm `Initialize`). After `ReplayShutOff` rings the imports
  `DUDT ... DTSDT` are zeroed.

## 4. State

- Imports: the seven increments above (declared in `SetServices`; `DTSDT` has
  `RESTART = RST`). Internals: seven bias tendencies (`DUDT ... DTSDT`) with
  `FRIENDLYTO = trim(COMP_NAME)`, so they are also exports. Exports: mostly borrowed from children
  (`CHILD_ID`), plus diagnostics (`QTFILL`, `QVFILL`, `PEFILL`, `TROPP_*`, `DPEDT_ANA`, ...).
- Private states (`ESMF_UserCompSetInternalState`): `UPD_STATE`, `IAU_COEFFS` (digital-filter
  weights). Module types `CONNECT_ANAnBKG`/`ANAnBKG` hold analysis-grid regridders for the 4D IAU.
- Restarts: `DTSDT` is `MAPL_RestartRequired` only when `ANALYZE_TS:` /= 0, else `MAPL_RestartSkip`
  (`SetServices`, `RST`). The other internal specs pass no `RESTART`, and `MAPL_StateAddInternalSpec`
  defaults to `MAPL_RestartOptional` (`src/Shared/@MAPL/generic/MAPL_Generic.F90`,
  `MAPL_StateAddInternalSpec`: `usable_RS = MAPL_RestartOptional`; enumerator in
  `src/Shared/@MAPL/shared/Constants/InternalConstants.F90`). *Inferred:* optional means the field
  is read when present and a bootstrap start is tolerated when absent (the reader treats
  `MAPL_RestartOptional` as bootstrappable, `src/Shared/@MAPL/base/NCIO.F90`).
- The Agcm declares its state by hand (`MAPL_Add*Spec`), not with an ACG `*_StateSpecs.rc`.

## 5. Configuration that changes structure or scheme

`DYCORE:` (default `FV3`; `FV3`, `FV3+ADV`, `DATMO`), `SCM_SL:` (default 0), `REPLAY_MODE:` (default
`NoReplay`; `Exact`, `Regular`), `IS_FCST:` (default 0), `IMPORT_RESTART_FILE:` (presence selects
DAS mode), `ANALYZE_TS:` (default 0; with `REPLAY_MODE:` other than `Exact`/`Regular` it asserts
`NoReplay`), `IAU_DIGITAL_FILTER:` (default `YES`), `CONSTRAIN_DAS:` (default 1), `ANA_IS_WEIGHTED:`
(default `NO`), `AINC_FILE:`. Read from the root configuration `AGCM.rc`
(template `src/Applications/@GEOSgcm_App/AGCM.rc.tmpl`).

## 6. Dependencies

`MAPL`, `ESMF`; `GEOS_Shared` (`GEOS_TopoGetMod`, `GEOS_RemapMod`, `GEOS_UtilsMod`) and `Chem_Shared`
(`Chem_GroupMod`, `Bundle_IncrementMod`, in the nested `GEOSchem_GridComp`)
(`GEOSagcm_GridComp/CMakeLists.txt`: `DEPENDENCIES MAPL GEOS_Shared Chem_Shared ESMF::ESMF`;
`use` lines in `GEOS_AgcmGridComp.F90`). Superdyn depends on `MAPL GEOS_Shared`. `HAS_GIGATRAJ`
(from `BUILD_WITH_GIGATRAJ`) adds the exports `PL, OMEGA, TH, DTDTDYN, ZL` to Agcm and `ZL` to Superdyn.

## 7. Known traps

- Order is set by `Run`: analysis update, orbit, dynamics 1, physics 1, copy of physics
  tendencies, dynamics 2. A new child needs an explicit `ESMF_GridCompRun` here.
- Tendency units differ: `DTDT` is mass-weighted, `DUDT, DVDT, DPEDT` are not (comment in `Run`).
- Several fields (`PHIS`, `SGH`, `VARFLT`) are set in `Initialize`, not by connectivity; do not
  remove the `MAPL_InitialRestart` marking.
- The Agcm `Run` gets the replay alarms by name; renaming an alarm in the Gcm breaks it silently
  only at run time (`ESMF_ClockGetAlarm`).
- `PREDICTOR_DURATION:` and `MKIAU_FREQUENCY:` are not defaulted here; the Gcm sets them
  (comment "Initialized in GCM_GridComp"). An Agcm run without the Gcm needs them in the config.
- `Run` has `save`d locals (`REPLAY_TIME0`, `first`): *Inferred:* one Agcm instance per process.

## 8. Decisions

See [`../docs/decisions.md`](../docs/decisions.md) (explicit child running, two-phase children,
resource-selected structure).

## 9. Observations (noted, not fixed)

- Seven `_ASSERT(..., 'needs informative message')` placeholders remain in `GEOS_AgcmGridComp.F90`.
- The file comment "Held-Suarez is not in module GEOSGCM?????" at the top of `GEOS_AgcmGridComp.F90`.
- `Dry_Mass_Check` is an internal procedure of `Run` called only in `#if debug` blocks; `debug` is defined as 0 at the top of the file.
- The `SDYN` to `PHYS` connections renaming `T, PLE, ZLE, TROPP_BLENDED` to `T_avg24, PLE_avg24,
  ZLE_avg24, TROPP_avg24` (`GEOS_AgcmGridComp.F90`, second `MAPL_AddConnectivity` call in
  `SetServices`) have no consumer in this checkout: `QuickChem` is not added as a Chemistry child
  here, and `GEOS_ChemGridComp.F90` has no `enable_`-style resource for it at all (3.1). Question
  for maintainers: are these connections kept for some other fixture or configuration that does
  add `QuickChem`, or are they left over from a time when it was wired in?

Open questions: none beyond the maintainer question in
[`../docs/decisions.md`](../docs/decisions.md) on whether the order of the Physics stages is meant
to be fixed.
