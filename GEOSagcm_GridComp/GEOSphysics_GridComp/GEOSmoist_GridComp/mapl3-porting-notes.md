# GEOS_MoistGridComp MAPL3 porting plan

Status (2026-10-29): steps 1-4 done on branch
`feature/pchakrab/port-moist-to-mapl3` (off `release/MAPL-v3`; steps 1-2
in `83eb97c5`, step 3 in `befaa640`). Not built: Moist is still commented out of the MAPL3
build (`../CMakeLists.txt:5`). Step 0 (MAPL2 baseline) not done yet.

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
- `:823-5516`: 553 exports (plus 6 commented out), mostly diagnostics.
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

**Done (2026-10-29).** The BACM_1M, MGB2_2M, THOM_1M and RAS `*_Setup`
branches in SetServices (`GEOS_MoistGridComp.F90:184-205`) `_FAIL`
before their call. SetServices runs first, so this one gate also covers
their Initialize/Run/Finalize branches, which are unchanged.
- `CLDMICR_OPTION` still defaults to `"BACM_1M"` (`:166`), so a config
  without it now fails; MAPL3 configs must set `GFDL_1M`. Not changed
  (would be a logic change).
- NSSL_2M and pyMoist are still in `CMakeLists.txt` (`:13`, `:21`)
  pending decision 5.

### 2. Macro pass (no logic change)
Decided 2026-10-28: **no re-indentation, no lowercasing, no `codee
format`** for Moist (unlike Solar step 1 and SU2G steps 1-1c). Existing
indentation and keyword case stay as they are; new or rewritten code
follows the surrounding style.

Changed 2026-10-29 for the driver only: `GEOS_MoistGridComp.F90` was
formatted with `codee format --verbose --config ~/config/codee-format`
in its own commit `f9b25d72` (after step 4), so it can be reverted
(`git revert f9b25d72`) or skipped in blame (`git blame --ignore-rev
f9b25d72`). Later edits to the driver follow codee's style; the scheme
interface modules keep theirs. A second comment-only commit `ca784434`
removed the 25 decorative separator lines (`!====`, `!----`, `!!!!`);
revert or blame-ignore it the same way. A third commit ("Moist: comment cleanup in
GEOS_MoistGridComp.F90", also holding this note) is comments and blank
lines only: the 19 `! !` prefixes (`! !INTERFACE:`, ...) -> `!`, ProTeX
blocks tightened and their text re-wrapped, `ErrLog Variables` and the
file-header macro comment dropped, a stray `!EOP` removed. Verified layout/case only: comments and
cpp lines identical; code tokens identical after normalizing codee's
equivalent spellings (`endif`/`enddo` -> `end if`/`end do`, `.ne.` etc.
-> `/=` etc., `real(ESMF_KIND_R8)` -> `real(kind=ESMF_KIND_R8)`).
Without `--config`, codee uses the first `.codee-format` up the tree
(here `~/.codee-format`, which differs from `~/config/codee-format` in
`IgnoredPreprocessorDirectives`; same output for this file).

Driver and the three in-scope interface modules: macros
`VERIFY_`/`RETURN_`/`__RC__` -> `_RC`/`_RETURN`/`_VERIFY`,
`RETURN_(ESMF_SUCCESS)` -> `_RETURN(_SUCCESS)`,
`MAPL_Generic.h` -> `MAPL.h`, drop `Iam`/`COMP_NAME`. Verify that only
macro spellings changed (normalize both versions and diff). Remove the
duplicated `---CLDMACRO` timer only in the functional step.

**Done (2026-10-29).** Done with a script: `RC=STATUS )` followed by
`VERIFY_(STATUS)` on the same line (after `;`), the next line, or after
blank lines was merged into `_RC)` (~1730 sites); remaining `VERIFY_` ->
`_VERIFY`.
- Left as `_VERIFY(STATUS)`: 6 `!!$`-commented sites in the driver, the
  duplicate verify after `MAPL_TimerAdd(GC, name="--GFDL_1M")`, and the
  `compute_radar_diagnostics` call in GFDL_1M (passes `STATUS`
  positionally).
- `Iam`/`COMP_NAME` dropped everywhere; `ESMF_GridCompGet` calls that
  also get `CONFIG`/`VM` (driver Initialize, Run) keep those. GF keeps
  `COMP_NAME` in `GF_Setup` because `FRIENDLIES%CNV_TR` uses it
  (`GEOS_GF_InterfaceMod.F90:58-66`); it goes with the services rework
  (step 4).
- Check: macros expanded, comments/whitespace/`&`/`;` stripped, token
  streams of old and new diffed. Only the `Iam`/`COMP_NAME` removals
  differ.

### 3. State specs
- `Moist_StateSpecs.rc` for the driver's 58 imports and 553 exports;
  `mapl_acg(... IMPORT_SPECS EXPORT_SPECS INTERNAL_SPECS GET_POINTERS DECLARE_POINTERS)`.
- Layout (decision 3): the driver's specs in `Moist_StateSpecs.rc`;
  each scheme's fields as hand-written MAPL3 `MAPL_GridCompAddSpec`
  calls in its `X_Setup`, where MAPL2 registers them.
  - `Moist_StateSpecs.rc`: the driver's imports and exports. The driver
    has no internals, as in MAPL2.
  - Scheme fields, as each `X_Setup` registers them in MAPL2: GFDL_1M 12
    internals, including the 9 every microphysics scheme declares (`Q`,
    `QLLS`, `QLCN`, `CLLS`, `CLCN`, `QILS`, `QICN`, `NACTL`, `NACTI`;
    identical in GFDL_1M, BACM_1M, THOM_1M, MGB2_2M), and 5 exports; GF
    3 internals; UW 1 internal. Too few for ACG to save anything, and no
    `COND` needed: `X_Setup` only runs for the selected scheme.
  - `Q` stays the first internal as in MAPL2: it is the first call in
    `GFDL_1M_Setup`, and the microphysics `Setup` runs before the
    convection ones (whether order still matters is issue 7).
  - Later microphysics schemes repeat the 9 shared internals and
    `QRAIN`/`QSNOW`/`QGRAUPEL`, as MAPL2 does. `NCPL`/`NCPI` differ
    between THOM_1M and MGB2_2M (`DEFAULT` 50e6/1e3 vs 0).
  - Pointers in the scheme modules: hand-written `MAPL_StateGetPointer`,
    as in MAPL2 (step 5). Don't include the driver's generated pointer
    headers there (they declare every driver field).
  - Scheme selection: replacing the 41 option-string compares
    (`adjustl(CONVPAR_OPTION)=="GF"`; SetServices `:206-212`, Initialize
    `:5636-5637`, Run `:6082-6095`, Finalize `:6742`) with module-level
    logicals `use_gfdl_1m`, `use_gf`, `use_uw`, ... set once after the
    options are checked (`:154-182`) is optional cleanup (step 4).
  - Considered and dropped (2026-10-28/29), each step closer to MAPL2:
    - Option B, one `Moist_StateSpecs.rc` with `COND` rows per scheme:
      pointer headers are split only per state, not per scheme
      (`mapl_acg.cmake:57-69`); scheme rows away from scheme code; a
      missing `COND` registers (and restarts) a field silently; same-name
      rows need distinct `ALIAS`es (`MAPL_GridCompSpecs_ACG.py:375-403,
      :569-574`).
    - The 9 shared internals in `Moist_StateSpecs.rc`, included before
      the scheme `Setup` calls: moved them away from where MAPL2
      registers them.
    - One spec file and `mapl_acg` call per scheme
      (`GFDL_1M_StateSpecs.rc`, ...; allowed, `mapl_acg.cmake:124-139`,
      precedent `@FVdycoreCubed_GridComp/CMakeLists.txt:121-122`): only
      1-17 fields each, so inline calls are simpler and match MAPL2.
- Column mapping (Solar step 3, SU2G step 2): `FRIENDLYTO` dropped
  (services, step 4; `ADD2EXPORT` where it meant "export myself"),
  `DEFAULT` -> `FILL`, internals `FILL=0.0` (MAPL2 zero-filled),
  `MAPL_RestartSkip` -> `SKIP`, `MAPL_RestartOptional` -> blank,
  `AVERAGING/REFRESH_INTERVAL` dropped.
- `AERO` (`MAPL_STATEITEM_STATE`) and the `MTR`/service items stay
  hand-written `MAPL_GridCompAddSpec` calls (ACG `ITEMTYPE` only does
  F/V).

**Done (2026-10-29).** `Moist_StateSpecs.rc` generated by script from
the MAPL2 calls (56 imports, 553 exports) and included in SetServices
(`Moist_Import___.h`, `Moist_Export___.h`) where those calls were; one
`mapl_acg` call in `CMakeLists.txt`. The 21 scheme calls rewritten as
MAPL3 `MAPL_GridCompAddSpec` calls in place: `GFDL_1M_Setup` 12
internals (`Q` first, MAPL2 order) and 5 exports, `GF_Setup` 3
internals, `UW_Setup` 1 internal. 630 = the 630 converted calls.
- `mapl_acg` passes `INTERNAL_SPECS` although the driver has no
  internals (empty `Moist_Internal___.h`, not included): `mapl_acg`
  declares the per-state pointer outputs for all three states
  (`mapl_acg.cmake:94-98`), and ACG only writes those for the states it
  is given (`MAPL_GridCompSpecs_ACG.py:654,:680-683`).
- Still hand-written MAPL2 calls: `MTR` (bundle) and `AERO` (state)
  imports (step 4).
- `COND`: `WSUB_OPTION == 0` (`ALH`, `TAUOROX`, `TAUOROY`), `== 1`
  (`WSUB_CLIM`), `== 2` (`KM`, `RI`), `MOIST_USE_NCLOUD_CLIM`
  (`NCPL_CLIM`, `NCPI_CLIM`), replacing the `select case`/`if` around
  those calls. Both are module variables (Process_Library, driver), so
  the same expressions work in Run's GetPointer include.
- `COND: PDF_DIAG` on the 15 `PDF_*` exports, which MAPL2 registered
  under `#ifdef PDFDIAG`. ACG can't emit preprocessor conditions, so the
  driver module has `logical, parameter :: PDF_DIAG`, `.true.` only when
  `PDFDIAG` is defined. `PDFDIAG` is never defined (commented-out
  `#define` in the driver `:5`, `cloudnew.F90:4`,
  `GEOS_MGB2_2M_InterfaceMod.F90:4`; none in CMake), and only MGB2_2M
  fills these exports (`GEOS_MGB2_2M_InterfaceMod.F90:1554-1656`).
- Changes from MAPL2 beyond the column mapping:
  - `RI` is both an import (Richardson number) and an export (ice
    effective radius); the export row has `ALIAS: CLDREFFI` (the name all
    microphysics schemes use) so the generated pointers don't clash.
  - `CU2DRAINMOVE`, `CU2DSNOWMOVE`: `xy` with `VLocationCenter` in
    MAPL2; ACG rejects that (`validate_vstagger`), so `VLOC: N`.
  - Long names lose trailing spaces; short names lose MAPL2's padding
    (`'DQDT  '`).
  - 9 exports have empty MAPL2 units (`SRF_TYPE`, `CNV_FRC`, `STOCH_CNV`,
    `SIGMA_DEEP`, `SIGMA_MID`, `RAS_*`); ACG omits `units`, which MAPL3
    treats as unchecked (matches any units,
    `@MAPL/superstructure/generic/specs/UnitsAspect.F90:64-67,:98-102`).
  - `RADSW` and `NLIM_NUC` had `RC=STATUS` with no verify; the generated
    calls check `_RC`.
  - The 6 commented-out (`!!$`) export calls were dropped; the section
    comments among the export calls are `#` rows in the export table.
- Checked: ACG (Python 3.10+; system `python3` is 3.6) generates the
  headers; a script compared every generated and inlined
  `MAPL_GridCompAddSpec` with its MAPL2 call (name, units, dims,
  stagger, restart, fill, long name, condition). Only the differences
  above.
- Left for later: `RFRSHINT`/`AVRGNINT`/`DT` reads in SetServices are now
  unused (step 4 resources); `FRIENDLIES` assignments in
  `GFDL_1M_Setup`/`GF_Setup` (step 4 services); pointer headers not
  used yet (steps 5, 8). Before
  including the driver's pointer headers in Run, check generated names
  against Run's local variables.

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

**Done (2026-10-29).**
- Signature `SetServices(GC, RC)`: no intent on `GC`, `RC` `intent(out)`
  non-optional (`@MAPL/infrastructure/esmf/ESMF_Interfaces.F90:32-37`).
- Entry points: Run with `phase_name='run'`; Initialize and Finalize keep
  the default names (`GENERIC::INIT_USER`, `GENERIC::FINALIZE_USER`,
  `@MAPL/superstructure/generic/MethodPhasesMap.F90:151-157`;
  `OuterMetaComponent/finalize.F90:27-34` runs the user Finalize by that
  name).
- Resources: `MAPL_GridCompGetResource(GC, 'LABEL', var, ...)`, no colon;
  `CONVPAR_OPTION`/`SHALLOW_OPTION`/`CLDMICR_OPTION` are now
  `character(len=:), allocatable` module variables. `SIGMA_DX`/`SIGMA_EXP`
  pass `default=(X)` (a copy: `value` is `intent(inout)` and must not
  alias `default`); they had `RC=STATUS` with no verify, now `_RC`. A
  missing key without a default is an error
  (`hconfig_get_value_template.h:13-15`). `RUN_DT`, `REFRESH_INTERVAL`
  and `AVERAGING_INTERVAL` reads dropped (only used by the old spec
  arguments).
- Timers: `MAPL_TimerAdd` dropped from SetServices and the `X_Setup`s
  (MAPL3 profiler timers are created on first start,
  `OuterMetaComponent/start_timer.F90`; GWD has no TimerAdd). The
  duplicated `---CLDMACRO` timer went with them. Starts/stops: step 5.
- `AERO` import: `MAPL_GridCompAddSpec(..., itemtype=MAPL_STATEITEM_STATE)`
  as Solar.
- `MTR`: the MAPL2 bundle import is now an export with
  `itemtype=MAPL_STATEITEM_SERVICE` (as DynCore's `TRADV`). MAPL puts the
  service bundle in both the import and the export state
  (`@MAPL/superstructure/generic/specs/ServiceClassAspect.F90:166-204`),
  so Run's `ESMF_StateGet(IMPORT,'MTR',TR)` is unchanged (DynCore reads
  its `TRADV` the same way, `DynCore_GridCompMod.F90:811-817`).
  Subscribers' fields keep their own short names (no `COMP::` prefix,
  `FieldClassAspect.F90:518-528`); `CNV_Tracers(:)%CNAME/QNAME` are set
  but unused.
- `CNV_Tracers_Init` (`Process_Library.F90`): no `FriendlyToMOIST` check;
  every field of the service bundle is a subscriber's, so all are
  friendly. Per-field metadata (`ScavengingFractionPerKm`,
  `SetofHenryLawCts`) is set with `ESMF_InfoSet` by the MAPL3 SU2G/SS2G
  ports, so it is still there; `ESMF_AttributeGet` -> `ESMF_Info*` is
  step 8.
- Advection: one `TRADV` import with `itemtype=MAPL_STATEITEM_SERVICE`,
  `service_items` = a `StringVector` the scheme `Setup`s fill (one
  subscription per component; items must be the component's internals,
  SU2G notes step 4). `GFDL_1M_Setup(GC, ADVECTED, RC)` adds the 10
  internals MAPL2 made `FRIENDLYTO=DYNAMICS` (`Q`, `CLLS`, `CLCN`,
  `QLLS`, `QLCN`, `QILS`, `QICN`, `QRAIN`, `QSNOW`, `QGRAUPEL`);
  `GF_Setup(GC, ADVECTED, RC)` adds `CNV_TR` when `CONVECTION_TRACER == 1`.
  `UW_Setup(GC, RC)`. Unconnected `TRADV` imports are accepted
  standalone (SU2G notes step 4). TURBULENCE/CHEMISTRY/ANALYSIS
  friendliness (Q and the condensates): no MAPL3 service yet.
- `FRIENDLIES` types/variables removed from GFDL_1M and GF; `GF_Setup`
  no longer gets `COMP_NAME`. The unported schemes' `*_Setup(GC, CF, RC)`
  are unchanged (unreachable after their `_FAIL`); SetServices declares
  `CF` only for those calls.
- `CNV_TR` (decided 2026-10-29):
  - MAPL2: `FRIENDLYTO="MOIST"` when `CONVECTION_TRACER == 0` (default),
    `"MOISTDYNAMICS"` when `== 1`. MAPL2 splits `FRIENDLYTO` on `:` into
    `FriendlyTo<token>` (`main` checkout,
    `@MAPL/generic/MAPL_Generic.F90:6890-6916`), so `== 1` made it
    friendly to no component: neither advected nor convectively
    transported.
  - MAPL3, `== 1`: as intended (MOIST:DYNAMICS): subscribed to `TRADV`
    and convectively transported. Differs from MAPL2 only in this
    non-default setting.
  - Convective transport by Moist itself: a provider's bundle only gets
    subscribers' items (`ServiceClassAspect.F90:268-304`), and Moist
    can't subscribe to its own service, so Run adds the `CNV_TR` field
    to the `MTR` bundle before `CNV_Tracers_Init` (step 8).
  - ConvPar only touches `CNV_TR` when `CONVECTION_TRACER == 1`
    (`ConvPar_GF2020.F90:95,461,613,1253,1679`); with `== 0` it stays at
    its fill value 0, so transporting it is probably a numerical no-op
    (not checked). Step 8 may add it only when `== 1`.
- Not done (optional): the `use_*` logicals replacing the 41
  option-string compares.
- Open: several subscribers with the same item name (e.g. two
  instances) in one service bundle; `CNV_Tracers_Init` looks fields up
  by name (`ESMF_FieldBundleGet(fieldName=)`).

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
- `MTR` service bundle, as DynCore handles its `TRADV`
  (`@FVdycoreCubed_GridComp/DynCore_GridCompMod.F90:809-817`): get it from
  the EXPORT state (Moist is the provider; today's
  `ESMF_StateGet(IMPORT,'MTR',TR)` only works because MAPL adds the
  bundle to both states) and `_ASSERT(MAPL_FieldBundleSameData(...))`
  against the IMPORT copy. No copying: the fields are the subscribers'
  internals, and `CNV_Tracers_Init` takes data pointers into them
  (`ESMFL_BundleGetPointerToData`), so GF's transport updates them in
  place, as FV3 does through `PULL_Q` (`:2525-2600`). The `COMP::NAME`
  parsing in `CNV_Tracers_Init` is dead (fields keep their short names).
  No exclusion list (DynCore's `xTRADV`, `:2628-2674`) unless a scheme
  needs one.
- Add Moist's own `CNV_TR` internal to the `MTR` bundle before
  `CNV_Tracers_Init` (step 4, `CNV_TR`).
- Edge-field remaps per issue 5; `p3d => X; X(1:IM,1:JM,0:LM) => p3d`
  through a `contiguous` scratch pointer (IRRAD notes).
- `MAPL_UNDEF` comparisons per issue 6.
- `Aer_Activation`: `ESMF_Attribute*` -> `ESMF_Info*` (issue 4).
- `MAPL_MaxMin` debug prints: now a function returning `[max, min]`
  (SU2G step 9).
- Finalize: drop `MAPL_GenericFinalize`.

### 9. Build
- Uncomment `GEOSmoist_GridComp` in `../CMakeLists.txt:5`; `TYPE SHARED`;
  one `mapl_acg` per spec file (step 3); external `Moist_SetServices` wrapper with a unique name
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
3. **DECIDED (2026-10-29) - Spec layout, as in MAPL2.** (Replaces
   option B, 2026-10-28, and a per-scheme spec file version.)
   `Moist_StateSpecs.rc` for the driver's imports and exports; each
   scheme's fields, including the internals all microphysics schemes
   declare, as hand-written `MAPL_GridCompAddSpec` calls in its
   `X_Setup`, where MAPL2 registers them. No `COND` rows for schemes
   (step 3).
4. **DECIDED (2026-10-28) - Convective-transport service name: `MTR`.** The
   advection service kept its MAPL2 bundle name: MAPL2 Physics exports the
   friendlies-to-DYNAMICS bundle as `TRADV` (`GEOS_PhysicsGridComp.F90:2425`)
   and the AGCM connects it to SUPERDYN (`GEOS_AgcmGridComp.F90:1138-1139`);
   MAPL3 DynCore exports the `TRADV` service under the same name
   (`DynCore_GridCompMod.F90:339-343`). The MAPL2 convection bundle is
   Moist's `MTR` import (`GEOS_MoistGridComp.F90:488-497`), filled by
   Physics with the friendlies to MOIST (`GEOS_PhysicsGridComp.F90:2477-2484`).
5. **NSSL_2M and pyMoist**: drop from the MAPL3 build or keep.
