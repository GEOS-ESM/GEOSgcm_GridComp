# GEOS_MoistGridComp MAPL3 porting plan

Status (2026-10-29): steps 1 and 2 done (not built; Moist is still
commented out of the MAPL3 build, `../CMakeLists.txt:5`). Step 0
(MAPL2 baseline) not done yet.

Derived from the Radiation notes
(`../@GEOSradiation_GridComp/mapl3-porting-notes.md`,
`GEOSirrad_GridComp/`, `GEOSsolar_GridComp/`, `GEOSsatsim_GridComp/`)
and GOCART's
(`../@GEOSchem_GridComp/@GOCART/ESMF/GOCART2G_GridComp/SU2G_GridComp/mapl3-porting-notes.md`).
The MAPL3 API gotchas found there are listed under "Lessons carried over"
and are not repeated in the steps.

## Reference files
- Template components: `../@GEOSradiation_GridComp/GEOSsolar_GridComp/GEOS_SolarGridComp.F90`,
  `../GEOSgwd_GridComp/GEOS_GwdGridComp.F90`
- Template state spec: `../GEOSgwd_GridComp/GWD_StateSpecs.rc`,
  `../@GEOSradiation_GridComp/GEOSsolar_GridComp/Solar_StateSpecs.rc`
- Template standalone regression: `../GEOSgwd_GridComp/regression/gwd-ncar-sa/`
  (root, cap, data-provider and component YAMLs), `../GEOSgwd_GridComp/regression/run_case.cmake`
- Service provider example: `DynCore_GridCompMod.F90:339-343` (`TRADV`
  export with `itemtype=MAPL_STATEITEM_SERVICE`) in `@FVdycoreCubed_GridComp`
- Service subscriber example: `@FVdycoreCubed_GridComp/regression/adv-dyn/data-moist.yaml`
  (`TRADV: {class: service, items: [Q]}`)
- MAPL3 docs: `src/Shared/@MAPL/docs/mapl3/api-changes.md`,
  `diffs-from-mapl2.md`, `gridcomp-template.md`

## Code overview (MAPL2, `GEOS_MoistGridComp.F90`)

`GEOS_MoistGridComp.F90` picks one deep convection scheme, one shallow
convection scheme and one cloud microphysics scheme from config, runs
them in order, and adds up their results. Each scheme lives in a
`GEOS_*_InterfaceMod.F90` with `X_Setup` / `X_Initialize` / `X_Run` /
(optional) `X_Finalize`; `X_Run` calls MAPL-free core physics
(`ConvPar_GF2020.F90`, `uwshcu.F90`, `gfdl_mp.F90`, ...).

### Module header (`:17-71`)
- Holds the three scheme choices (`CONVPAR_OPTION`, `SHALLOW_OPTION`,
  `CLDMICR_OPTION`) and a few settings. Most tuning settings are module
  variables in `Process_Library.F90` (e.g. `USE_AEROSOL_NN` `:85`,
  `SH_MD_DP` `:147`, `CNV_FRACTION_MIN` `:259`), read by the scheme
  `*_Initialize` routines, often with scheme-specific defaults (GFDL sets
  `SH_MD_DP` default `.TRUE.`, `GEOS_GFDL_1M_InterfaceMod.F90:331`).

### SetServices (`:79-5542`)
- `:126-137`: registers the Initialize, Run and Finalize entry points.
- `:146-181`: reads the scheme choices: deep GF/RAS/NONE, shallow
  UW/NONE, microphysics BACM_1M/GFDL_1M/THOM_1M/MGB2_2M. It asserts each
  one is valid. (`GEOS_NSSL_2M_InterfaceMod.F90` and `pyMoist/` exist but
  can't be selected.)
- `:189-201`: reads the aerosol settings and the fixed land/ocean CCN/IN
  numbers.
- `:206-212`: calls each chosen scheme's `*_Setup`. Microphysics goes
  first because it declares the internal state (`Q`, `QLLS/QLCN`,
  `QILS/QICN`, `CLLS/CLCN`, `NACTL/NACTI`), and binary restarts expect
  `Q` first (`:203`).
- `:214-822`: 58 imports covering the dynamics state, turbulence, surface
  fields, the `MTR` tracer bundle (`:490`) and the `AERO` state (`:526`).
- `:823-5516`: 559 exports, mostly diagnostics.
- `:5522-5537`: adds timers (`---CLDMACRO` twice, `:5526`/`:5528`) and
  calls `MAPL_GenericSetServices`.

### Initialize (`:5557-5672`)
- Runs generic init, reads the precip-type and convective-detrainment
  settings, and calls each scheme's `*_Initialize`.
- Creates a reflectivity alarm (`DBZ_RunAlarm`) on the clock that rings
  at least every 900 s (`:5649-5659`).

### Run (`:5683-6720`)
Only runs when Moist's alarm rings (`:5796`).
1. **Prologue (`:5798-6024`)**
   - Gets pointers to the state, sets the surface type, and computes
     pressure, height, layer mass and saturation humidity.
   - Computes LTS/EIS, CAPE/CIN (`BUOYANCY2`, `BUOYANCY`) and a
     convective fraction `CNV_FRC` that ramps with CAPE.
   - Saves copies of the incoming state.
2. **Tracers and aerosol activation (`:6027-6070`)**: `CNV_Tracers_Init`
   picks the convectively transported tracers out of `MTR`; gets
   `NACTL/NACTI` from a climatology, from GOCART aerosol
   (`Aer_Activation`), or from the fixed CCN numbers.
3. **Convection (`:6072-6112`)**: runs UW and GF/RAS. `SH_MD_DP` decides
   which runs first.
4. **Mass fluxes (`:6116-6165`)**: adds up the deep and shallow mass
   fluxes. Where convection is inactive, it moves convective condensate
   and cloud into the large-scale fields (`:6126-6156`).
5. **Microphysics (`:6169-6172`)**: calls the chosen `*_Run`.
6. **Epilogue (`:6183-6655`)**
   - Diagnostics: ice and liquid cloud fraction, relative humidity.
   - Rains out any vapor above 110% saturation (`:6231-6249`).
   - Adds up per-scheme tendencies into totals (`:6253-6378`). `DTDT` is
     weighted by layer pressure thickness.
   - Adds up precipitation, diagnoses precip type, and computes the
     snow:liquid ratio. The precip blocks read each other's exports by
     name and update some in place, so their order matters
     (`:6423-6536`).
   - Fills the column water paths and the after-Moist state fields
     (`*AFMOIST`) that SYNCTQ uses.
7. **When the alarm isn't ringing (`:6657-6714`)**: only refills the
   `*AFMOIST` and `RH2` exports.

### Finalize (`:6722-6750`)
Cleans up GFDL, UW and GF, then calls generic finalize.

## Starting point (survey 2026-10-28)

| File | Lines | `MAPL_GetPointer` | `MAPL_Add*Spec` | `MAPL_GetResource` | `ALLOC=.TRUE.` |
|---|---|---|---|---|---|
| `GEOS_MoistGridComp.F90` | 6752 | 229 | 618 | 28 | 45 |
| `GEOS_GFDL_1M_InterfaceMod.F90` | 1576 | 145 | 17 | 41 | 70 |
| `GEOS_GF_InterfaceMod.F90` | 1057 | 113 | 3 | 207 | 2 |
| `GEOS_UW_InterfaceMod.F90` | 781 | 96 | 1 | 55 | 51 |
| `GEOS_BACM_1M_InterfaceMod.F90` | 829 | 153 | 9 | 62 | 107 |
| `GEOS_MGB2_2M_InterfaceMod.F90` | 2710 | 236 | 17 | 102 | 141 |
| `GEOS_THOM_1M_InterfaceMod.F90` | 1112 | 121 | 15 | 19 | 64 |
| `GEOS_RAS_InterfaceMod.F90` | 347 | 7 | 0 | 38 | 0 |
| `GEOS_NSSL_2M_InterfaceMod.F90` (unused) | 1115 | 120 | 15 | 23 | 62 |
| `aer_actv_single_moment.F90` | 614 | 11 | 0 | 0 | 0 |
| `Process_Library.F90` | 7179 | 7 | 0 | 0 | 0 |

- `aer_actv_single_moment.F90` uses `ESMF_AttributeGet/Set` on the AERO
  state (14 sites, `:92-155`); `Process_Library.F90` has 23
  `ESMF_Attribute` calls.
- Driver specs: 39 `MAPL_VLocationEdge`, 27 `RESTART=`, no
  `UNGRIDDED_DIMS`, no children.
- Interface modules take a `MAPL_MetaComp` (`X_Initialize(MAPL, ...)`)
  and call `MAPL_GetObjectFromGC` / `MAPL_Get(MAPL, IM, JM, LM,
  RUNALARM=...)` themselves (e.g. `GEOS_GFDL_1M_InterfaceMod.F90:279`,
  `:571-582`).
- Default production configuration: GFDL_1M microphysics
  (`@GEOSgcm_App/gcm_setup:332`), GF deep (`gcm_setup:1045-1048`), UW
  shallow (`AGCM.rc.tmpl:109`).
- Remote branch `feature/bmauer/turbulence_mapl3_port` exists in
  `@GEOSgcm_GridComp` (Turbulence port, not inspected).

## Moist-specific MAPL3 issues

1. **Friendliness -> services.** MAPL2 `FRIENDLYTO` is stubbed in MAPL3
   (`api-changes.md` "Friendly mechanism"); the replacement is
   `MAPL_STATEITEM_SERVICE` (`diffs-from-mapl2.md` §3). Moist is on both
   sides:
   - Subscriber: microphysics internals are friendly to
     `DYNAMICS:TURBULENCE[:CHEMISTRY:ANALYSIS]`
     (`GEOS_GFDL_1M_InterfaceMod.F90:88-97`). Subscribe them to `TRADV`
     (advection) as SU did (SU2G notes step 4); turbulence once it
     provides a service.
   - Provider: Moist imports the `MTR` bundle (`:490`) and
     `CNV_Tracers_Init` selects fields whose Info has `FriendlyToMOIST`
     (`Process_Library.F90` ~`:378`). In MAPL3 Moist must export a
     convective-transport service `MTR` (as DynCore exports `TRADV`) that
     GOCART's SU/DU/SS (and GF's `CNV_TR`, `GEOS_GF_InterfaceMod.F90:68-76`)
     subscribe to. This is what GOCART's "convection service" blocker
     waits on (SU2G notes, Status).
2. **`ALLOC=.TRUE.` has no MAPL3 equivalent.** ~540 sites. MAPL3's
   `MAPL_StateGetPointer` has no `alloc` (Radiation notes), exports are
   activated only when connected/requested or by the test-only
   `activate_all_exports` (`ComponentSpec.F90:31`), and internals are
   always activated (`OuterMetaComponent/advertise_variable.F90:55-57`).
   Exports used as work arrays (`SRF_TYPE`, `CNV_FRC`, `CAPE`, `RAIN`,
   `TPREC`, `LS_PRCP`, `NWFA`, ...) and the epilogue sums that read other
   exports by name (`:6423-6536`) silently change results if those
   exports are unallocated. Handled by step 6 (decision 1).
3. **Alarms.** Moist runs on its own `RUN_DT` via `RUNALARM`
   (`:5644-5647`, `:5786-5796`); MAPL3 gives each component its own clock
   and `timestep` (`diffs-from-mapl2.md` §1, §9). `DBZ_RunAlarm` is
   created on the clock passed to Initialize (`:5655`) and is read by the
   microphysics. Both need the alarm treatment Solar/IRRAD needed (see
   Lessons).
4. **AERO state** for aerosol activation: `ESMF_Attribute*` ->
   `ESMF_Info*` (as Solar step 18), using `ESMF_InfoGetCharAlloc` /
   `ESMF_InfoGetAlloc` for allocatable outputs (SU2G step 7).
5. **Edge fields**: 39 Edge specs in the driver (`PLE`, `ZLE`, `KH`,
   `PFL_*`/`PFI_*`, mass fluxes, ...). MAPL2 indexing is mixed: most code
   is 0-based (`PLE(:,:,0:LM-1)`, `ZLE(:,:,LM)` at `:5887-5896`), but
   `TMP3Dp1(IM,JM,LM+1)` sums `UMF_DC/UMF_SC` 1-based (`:6117-6131`).
   Check each site; don't remap blindly.
6. **`MAPL_UNDEF` changed** from `1.0e15` to `-huge` (SU2G step 6). Moist
   compares imports against it (`SNOMAS /= MAPL_UNDEF` `:5841`,
   `CAPE .ne. MAPL_UNDEF` `:6007-6018`, `DETR_SC /= MAPL_UNDEF`
   `GEOS_UW_InterfaceMod.F90:587`, `v2 .eq. MAPL_UNDEF`
   `Process_Library.F90:3829`). Equality tests only work if the provider fills with
   MAPL3's value; inputs converted from MAPL2 restarts carry `1e15`.
   Grep all scheme modules and `Process_Library.F90` for inequality tests
   as well.
7. **Not checked yet**: where `Lightning_mod` (`use` at `:32`) comes
   from and whether it is MAPL3-clean; whether `pyMoist/.../moist_workarounds.F90`
   (in `CMakeLists.txt` srcs) builds against MAPL3; whether the
   "Q first" restart-order rule (`:203`) still matters with MAPL3's
   netCDF restarts.

## Steps

### 0. MAPL2 baseline (do first)
- Capture at C12L181 with `COMPONENT_TO_RECORD: MOIST`, before and after
  Run for import, internal and export, as for Solar/SU. Runs start on the
  hour; record a later step so a first-step crash doesn't lose the files
  (SU2G step 0).
- Configuration: GFDL_1M + GF + UW (production default). Moist `RUN_DT`
  = model `DT` to keep the alarm out of the first comparison.
  `USE_AEROSOL_NN: false` first, so activation uses the fixed CCN numbers
  (`:6058-6061`) and doesn't need GOCART; a second baseline with aerosol
  activation later.
- Request every Moist export in `HISTORY.rc` (MAPL2 allocates exports
  only on request; SU's first baseline had wrong totals because of this).
- Record every Moist setting used (scheme options, `SH_MD_DP`,
  `DETRAIN_INACTIVE_CNV`, `CNV_FRACTION_*`, precip-type flags, `DBZ_DT`).

### 1. Scope the schemes
Port GFDL_1M, GF and UW first. BACM_1M, THOM_1M, MGB2_2M and RAS keep
their `if` branches but `_FAIL` until ported (decision 2); NSSL_2M is
not selectable today and stays out of the build list or is removed
(decision 5).

### 2. Macro pass (no logic change)
Decided 2026-10-28: **no re-indentation, no lowercasing, no `codee
format`** for Moist (unlike Solar step 1 and SU2G steps 1-1c). Existing
indentation and keyword case stay as they are; new or rewritten code
follows the surrounding style.

Driver and the three in-scope interface modules: macros
`VERIFY_`/`RETURN_`/`__RC__` -> `_RC`/`_RETURN`/`_VERIFY`,
`RETURN_(ESMF_SUCCESS)` -> `_RETURN(_SUCCESS)`,
`MAPL_Generic.h` -> `MAPL.h`, drop `Iam`/`COMP_NAME`. Verify that only
macro spellings changed (normalize both versions and diff). Remove the
duplicated `---CLDMACRO` timer only in the functional step.

### 3. State specs
- `Moist_StateSpecs.rc` for the driver's 58 imports and 559 exports;
  `mapl_acg(... IMPORT_SPECS EXPORT_SPECS GET_POINTERS DECLARE_POINTERS)`.
- Scheme specs (decision 3, option B): the same `Moist_StateSpecs.rc`
  also holds the specs each `X_Setup` registers today. Scheme-specific
  rows get a `COND` on a scheme logical (`COND: use_gf`). ACG wraps both
  the AddSpec call and the GetPointer call of a `COND` row in
  `if (cond) then`, the GetPointer one with `else nullify(...)`
  (`@MAPL/apps/MAPL_GridCompSpecs_ACG.py:352-356,:411-420,:797-814`).
  Precedent: GWD (`GWD_StateSpecs.rc:10-14`), Ocean, DataSea, MOM6.
  - Rows every microphysics scheme declares (`Q`, `QLLS`, `QLCN`,
    `QILS`, `QICN`, `CLLS`, `CLCN`, `NACTL`, `NACTI`) are written once,
    unconditional. Unverified: that their attributes match across
    microphysics schemes. Check when the second one is ported; a field
    that differs needs one row per scheme with exclusive `COND`s.
  - `Q` stays the first internal by being the first internal row.
  - Scheme selection in one place. Today the option strings are module
    variables (`GEOS_MoistGridComp.F90:38-40`) compared 41 times
    (`adjustl(CONVPAR_OPTION)=="GF"`) in SetServices (`:206-212`),
    Initialize (`:5636-5637`), Run (`:6082-6095`) and Finalize
    (`:6742`). Set module-level logicals once, right after the options
    are read and checked (`:154-182`): `use_gfdl_1m`, `use_gf`, `use_uw`
    (one per scheme). Use them in the `COND` column and replace every
    string compare with them. They must be visible in both SetServices
    (AddSpec includes) and Run (GetPointer includes), as GWD's
    `self%NCAR_NRDG` is (`GEOS_GwdGridComp.F90:131-133,:478`). Residual
    risk: a scheme row missing its `COND` is registered (and restarted)
    unconditionally; catch it in review of the spec file.
  - Shared fields outside the driver: GF (`GEOS_GF_InterfaceMod.F90:582-589`)
    and UW (`GEOS_UW_InterfaceMod.F90:282-288`) read the microphysics
    internals and must work with any microphysics scheme. They keep
    hand-written `MAPL_StateGetPointer` calls, as GFDL_1M may
    (`GEOS_GFDL_1M_InterfaceMod.F90:596-607`). Don't include the driver's
    pointer headers there: pointer variables are named after `ALIAS` or
    `SHORT_NAME` (`MAPL_GridCompSpecs_ACG.py:570-574`), so they would
    redeclare `Q` and every other internal.
  - The driver's Run (`:5803-5811`, `:6660`) can use the generated
    pointer headers for these fields.
  - Option A (one spec file per scheme, one `mapl_acg` call each, which
    MAPL allows: `@MAPL/cmake/mapl_acg.cmake:124-139`,
    `@FVdycoreCubed_GridComp/CMakeLists.txt:121-122`) was considered and
    not chosen: shared microphysics rows would be copied into every
    microphysics file.
- Column mapping (Solar step 3, SU2G step 2): `FRIENDLYTO` dropped
  (services, step 4; `ADD2EXPORT` where it meant "export myself"),
  `DEFAULT` -> `FILL`, internals `FILL=0.0` (MAPL2 zero-filled),
  `MAPL_RestartSkip` -> `SKIP`, `MAPL_RestartOptional` -> blank,
  `AVERAGING/REFRESH_INTERVAL` dropped.
- `AERO` (`MAPL_STATEITEM_STATE`) and the `MTR`/service items stay
  hand-written `MAPL_GridCompAddSpec` calls (ACG `ITEMTYPE` only does
  F/V).

### 4. SetServices
- `MAPL_GridCompSetEntryPoint` with `phase_name='run'`; keep Initialize's
  default name (MAPL calls it as `GENERIC::INIT_USER`).
- Resources via `MAPL_GridCompGetResource(gc, "LABEL", ...)` (no colon),
  strings `character(len=:), allocatable`.
- Advection subscription for the microphysics internals (`TRADV`,
  `service_items`), as SU/DU/SS.
- Convective-transport service export `MTR`
  (`itemtype=MAPL_STATEITEM_SERVICE`, as DynCore's `TRADV`);
  `CNV_Tracers_Init` reads the service bundle instead of
  `FriendlyToMOIST` Info.
- `MAPL_GenericSetServices` removed.

### 5. Interface modules
For GFDL_1M, GF, UW: `X_Initialize(MAPL, ...)` -> take `gc`;
`MAPL_GetResource(MAPL, ...)` -> `MAPL_GridCompGetResource`;
`MAPL_GetObjectFromGC`/`MAPL_Get` -> `MAPL_GridCompGet(gc, num_levels=,
grid=)` (never `ESMF_GridCompGet(grid=)`, Solar step 18);
`MAPL_GetPointer` -> generated pointers or `MAPL_StateGetPointer`;
timers -> `MAPL_GridCompTimerStart/Stop`; `ALLOC=` per step 6.

### 6. `ALLOC=.TRUE.` replacement
Per decision 1: local work arrays copied to the export if associated;
internals with `ADD2EXPORT` (and `RESTART: SKIP`) for fields other
components always connect. Apply to the driver and the in-scope modules.
Keep the epilogue's results independent of which exports are connected.

### 7. Initialize
Drop `MAPL_GenericInitialize`. Timestep from `MAPL_ClockGet(clock,
dt=)` instead of `RUNALARM` (SU2G step 5). Recreate `DBZ_RunAlarm` the
way Solar recreated its run alarm (Solar step 18), and check it rings on
the first step (IRRAD "`irrad_lw_alarm` never rang on the first step").
Decide whether Run's `ESMF_AlarmIsRinging` branch (`:5796`) is still
needed once MAPL3 schedules Moist with its own `timestep`.

### 8. Run and Finalize
- Edge-field remaps per issue 5; `p3d => X; X(1:IM,1:JM,0:LM) => p3d`
  through a `contiguous` scratch pointer (IRRAD notes).
- `MAPL_UNDEF` comparisons per issue 6.
- `Aer_Activation`: `ESMF_Attribute*` -> `ESMF_Info*` (issue 4).
- `MAPL_MaxMin` debug prints: now a function returning `[max, min]`
  (SU2G step 9).
- Finalize: drop `MAPL_GenericFinalize`.

### 9. Build
- Uncomment `GEOSmoist_GridComp` in `../CMakeLists.txt:5`; `TYPE SHARED`;
  `mapl_acg`; external `Moist_SetServices` wrapper with a unique name
  (SU2G step 8).
- Entry points must match `I_Run`: no `intent` on `gc/import/export/clock`,
  `rc` `intent(out)` non-optional (Radiation notes; SU2G step 4).

### 10. Standalone regression `moist-sa`
After `gwd-ncar-sa`: a configurable `DataForMoist` provider restarted
from the step-0 MAPL2 import state, `activate_all_exports: true` for
Moist, a stub `AERO` (FakeGOCART) and an empty convective-transport
subscriber. Compare with `cmpchk.py`; target bit-for-bit. Then the
aerosol-activation baseline with FakeGOCART or GOCART2G.

### 11. Remaining schemes
BACM_1M, then MGB2_2M, THOM_1M, RAS, each with its own baseline.

### 12. Wiring (blocked on the Physics/AGCM ports)
Connect Moist's imports in Physics, the `TRADV` subscriptions, and the
convective-transport service to GOCART's SU/DU/SS.

## Lessons carried over (Radiation, GOCART)
- Entry-point signatures must match `I_SetServices`/`I_Run` exactly
  (Radiation notes; SU2G step 4/9).
- `MAPL_GridCompGetResource` labels have no trailing colon.
- `MAPL_GenericSetServices`/`Initialize`/`Finalize`/`RunChildren` are not
  called.
- Grid: `MAPL_GridCompGet(gc, grid=...)`, not `ESMF_GridCompGet`
  (Solar step 18); `MAPL_GridGetInterior` -> `MAPL_GridGet(interior=)`.
- Edge fields are 1-based in MAPL3; remap through a `contiguous` scratch
  pointer (IRRAD, Radiation, Solar step 18, SU2G step 6).
- `DEFAULT` -> `FILL`; MAPL2 zero-filled internals, MAPL3 does not
  (Solar step 3, SU2G step 2).
- `MAPL_UNDEF` is `-huge` in MAPL3, `1.0e15` in MAPL2 (SU2G step 6).
- `ESMF_Attribute*` -> `ESMF_Info*`; use the `...Alloc` getters for
  allocatable outputs (Solar step 18, SU2G step 7).
- Component alarms need MAPL2's main-alarm logic recreated (Solar step
  18; IRRAD first-step alarm fix).
- Unconnected child imports propagate up automatically; several
  subscribers per service provider is by design (SU2G step 4).
- MAPL2 baselines must match the standalone configuration exactly
  (`AERO_PROVIDER`, requested exports, refreshed internal restart;
  Radiation notes, Solar step 18, SU2G step 0).

## Open decisions
1. **DECIDED (2026-10-28) - `ALLOC=.TRUE.`**: by default, a local
   allocatable work array, copied to the export if it is associated.
   Fields other components always connect (`CNV_FRC`, `RAIN`, `SNO`,
   `TPREC`, ...) become internals with `ADD2EXPORT` instead (always
   active, `OuterMetaComponent/advertise_variable.F90:55-57`) and
   `RESTART: SKIP` unless they must survive a restart.
2. **DECIDED (2026-10-28) - Scheme scope**: GFDL_1M + GF + UW first;
   BACM_1M, THOM_1M, MGB2_2M and RAS `_FAIL` until ported (step 1).
3. **DECIDED (2026-10-28) - Spec files: option B.** A single
   `Moist_StateSpecs.rc` holding the driver's and every scheme's specs,
   scheme-specific rows gated with `COND` on module-level `use_*`
   logicals that also replace the 41 option-string compares, so scheme
   selection lives in one place (step 3).
4. **DECIDED (2026-10-28) - Convective-transport service name: `MTR`.** The
   advection service kept its MAPL2 bundle name: MAPL2 Physics exports the
   friendlies-to-DYNAMICS bundle as `TRADV` (`GEOS_PhysicsGridComp.F90:2425`)
   and the AGCM connects it to SUPERDYN (`GEOS_AgcmGridComp.F90:1138-1139`);
   MAPL3 DynCore exports the `TRADV` service under the same name
   (`DynCore_GridCompMod.F90:339-343`). The MAPL2 convection bundle is
   Moist's `MTR` import (`GEOS_MoistGridComp.F90:488-497`), filled by
   Physics with the friendlies to MOIST (`GEOS_PhysicsGridComp.F90:2477-2484`).
5. **NSSL_2M and pyMoist**: drop from the MAPL3 build or keep.
