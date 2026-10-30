# GEOS_MoistGridComp MAPL3 porting plan

Status (2026-10-30): steps 1-8 done on branch `feature/pchakrab/port-moist-to-mapl3` (off `release/MAPL-v3`; steps 1-2 in `83eb97c5`, step 3 in `befaa640`). Not built: Moist is still commented out of the MAPL3 build (`../CMakeLists.txt:5`). Step 0 (MAPL2 baseline) not done yet.

Derived from the Radiation notes (`../@GEOSradiation_GridComp/mapl3-porting-notes.md`, `GEOSirrad_GridComp/`, `GEOSsolar_GridComp/`, `GEOSsatsim_GridComp/`) and GOCART's (`../@GEOSchem_GridComp/@GOCART/ESMF/GOCART2G_GridComp/SU2G_GridComp/mapl3-porting-notes.md`). The MAPL3 API gotchas found there are listed under "Lessons carried over" and are not repeated in the steps.

## Reference files
- Template components: `../@GEOSradiation_GridComp/GEOSsolar_GridComp/GEOS_SolarGridComp.F90`, `../GEOSgwd_GridComp/GEOS_GwdGridComp.F90`
- Template state spec: `../GEOSgwd_GridComp/GWD_StateSpecs.rc`, `../@GEOSradiation_GridComp/GEOSsolar_GridComp/Solar_StateSpecs.rc`
- Template standalone regression: `../GEOSgwd_GridComp/regression/gwd-ncar-sa/` (root, cap, data-provider and component YAMLs), `../GEOSgwd_GridComp/regression/run_case.cmake`
- Service provider example: `DynCore_GridCompMod.F90:339-343` (`TRADV` export with `itemtype=MAPL_STATEITEM_SERVICE`) in `@FVdycoreCubed_GridComp`
- Service subscriber example: `@FVdycoreCubed_GridComp/regression/adv-dyn/data-moist.yaml` (`TRADV: {class: service, items: [Q]}`)
- MAPL3 docs: `src/Shared/@MAPL/docs/mapl3/api-changes.md`, `diffs-from-mapl2.md`, `gridcomp-template.md`

## Code overview (MAPL2, `GEOS_MoistGridComp.F90`)

`GEOS_MoistGridComp.F90` picks one deep convection scheme, one shallow convection scheme and one cloud microphysics scheme from config, runs them in order, and adds up their results. Each scheme lives in a `GEOS_*_InterfaceMod.F90` with `X_Setup` / `X_Initialize` / `X_Run` / (optional) `X_Finalize`; `X_Run` calls MAPL-free core physics (`ConvPar_GF2020.F90`, `uwshcu.F90`, `gfdl_mp.F90`, ...).

### Module header (`:17-71`)
- Holds the three scheme choices (`CONVPAR_OPTION`, `SHALLOW_OPTION`, `CLDMICR_OPTION`) and a few settings. Most tuning settings are module variables in `Process_Library.F90` (e.g. `USE_AEROSOL_NN` `:85`, `SH_MD_DP` `:147`, `CNV_FRACTION_MIN` `:259`), read by the scheme `*_Initialize` routines, often with scheme-specific defaults (GFDL sets `SH_MD_DP` default `.TRUE.`, `GEOS_GFDL_1M_InterfaceMod.F90:331`).

### SetServices (`:79-5542`)
- `:126-137`: registers the Initialize, Run and Finalize entry points.
- `:146-181`: reads the scheme choices: deep GF/RAS/NONE, shallow UW/NONE, microphysics BACM_1M/GFDL_1M/THOM_1M/MGB2_2M. It asserts each one is valid. (`GEOS_NSSL_2M_InterfaceMod.F90` and `pyMoist/` exist but can't be selected.)
- `:189-201`: reads the aerosol settings and the fixed land/ocean CCN/IN numbers.
- `:206-212`: calls each chosen scheme's `*_Setup`. Microphysics goes first because it declares the internal state (`Q`, `QLLS/QLCN`, `QILS/QICN`, `CLLS/CLCN`, `NACTL/NACTI`), and binary restarts expect `Q` first (`:203`).
- `:214-822`: 58 imports covering the dynamics state, turbulence, surface fields, the `MTR` tracer bundle (`:490`) and the `AERO` state (`:526`).
- `:823-5516`: 553 exports (plus 6 commented out), mostly diagnostics.
- `:5522-5537`: adds timers (`---CLDMACRO` twice, `:5526`/`:5528`) and calls `MAPL_GenericSetServices`.

### Initialize (`:5557-5672`)
- Runs generic init, reads the precip-type and convective-detrainment settings, and calls each scheme's `*_Initialize`.
- Creates a reflectivity alarm (`DBZ_RunAlarm`) on the clock that rings at least every 900 s (`:5649-5659`).

### Run (`:5683-6720`)
Only runs when Moist's alarm rings (`:5796`).
1. **Prologue (`:5798-6024`)**
   - Gets pointers to the state, sets the surface type, and computes pressure, height, layer mass and saturation humidity.
   - Computes LTS/EIS, CAPE/CIN (`BUOYANCY2`, `BUOYANCY`) and a convective fraction `CNV_FRC` that ramps with CAPE.
   - Saves copies of the incoming state.
2. **Tracers and aerosol activation (`:6027-6070`)**: `CNV_Tracers_Init` picks the convectively transported tracers out of `MTR`; gets `NACTL/NACTI` from a climatology, from GOCART aerosol (`Aer_Activation`), or from the fixed CCN numbers.
3. **Convection (`:6072-6112`)**: runs UW and GF/RAS. `SH_MD_DP` decides which runs first.
4. **Mass fluxes (`:6116-6165`)**: adds up the deep and shallow mass fluxes. Where convection is inactive, it moves convective condensate and cloud into the large-scale fields (`:6126-6156`).
5. **Microphysics (`:6169-6172`)**: calls the chosen `*_Run`.
6. **Epilogue (`:6183-6655`)**
   - Diagnostics: ice and liquid cloud fraction, relative humidity.
   - Rains out any vapor above 110% saturation (`:6231-6249`).
   - Adds up per-scheme tendencies into totals (`:6253-6378`). `DTDT` is weighted by layer pressure thickness.
   - Adds up precipitation, diagnoses precip type, and computes the snow:liquid ratio. The precip blocks read each other's exports by name and update some in place, so their order matters (`:6423-6536`).
   - Fills the column water paths and the after-Moist state fields (`*AFMOIST`) that SYNCTQ uses.
7. **When the alarm isn't ringing (`:6657-6714`)**: only refills the `*AFMOIST` and `RH2` exports.

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

- `aer_actv_single_moment.F90` uses `ESMF_AttributeGet/Set` on the AERO state (14 sites, `:92-155`); `Process_Library.F90` has 23 `ESMF_Attribute` calls.
- Driver specs: 39 `MAPL_VLocationEdge`, 27 `RESTART=`, no `UNGRIDDED_DIMS`, no children.
- Interface modules take a `MAPL_MetaComp` (`X_Initialize(MAPL, ...)`) and call `MAPL_GetObjectFromGC` / `MAPL_Get(MAPL, IM, JM, LM, RUNALARM=...)` themselves (e.g. `GEOS_GFDL_1M_InterfaceMod.F90:279`, `:571-582`).
- Default production configuration: GFDL_1M microphysics (`@GEOSgcm_App/gcm_setup:332`), GF deep (`gcm_setup:1045-1048`), UW shallow (`AGCM.rc.tmpl:109`).
- Remote branch `feature/bmauer/turbulence_mapl3_port` exists in `@GEOSgcm_GridComp` (Turbulence port, not inspected).

## Moist-specific MAPL3 issues

1. **Friendliness -> services.** MAPL2 `FRIENDLYTO` is stubbed in MAPL3 (`api-changes.md` "Friendly mechanism"); the replacement is `MAPL_STATEITEM_SERVICE` (`diffs-from-mapl2.md` §3). Moist is on both sides:
   - Subscriber: microphysics internals are friendly to `DYNAMICS:TURBULENCE[:CHEMISTRY:ANALYSIS]` (`GEOS_GFDL_1M_InterfaceMod.F90:88-97`). Subscribe them to `TRADV` (advection) as SU did (SU2G notes step 4); turbulence once it provides a service.
   - Provider: Moist imports the `MTR` bundle (`:490`) and `CNV_Tracers_Init` selects fields whose Info has `FriendlyToMOIST` (`Process_Library.F90` ~`:378`). In MAPL3 Moist must export a convective-transport service `MTR` (as DynCore exports `TRADV`) that GOCART's SU/DU/SS (and GF's `CNV_TR`, `GEOS_GF_InterfaceMod.F90:68-76`) subscribe to. This is what GOCART's "convection service" blocker waits on (SU2G notes, Status).
2. **`ALLOC=.TRUE.` has no MAPL3 equivalent.** ~540 sites. MAPL3's `MAPL_StateGetPointer` has no `alloc` (Radiation notes), exports are activated only when connected/requested or by the test-only `activate_all_exports` (`ComponentSpec.F90:31`), and internals are always activated (`OuterMetaComponent/advertise_variable.F90:55-57`). Exports used as work arrays (`SRF_TYPE`, `CNV_FRC`, `CAPE`, `RAIN`, `TPREC`, `LS_PRCP`, `NWFA`, ...) and the epilogue sums that read other exports by name (`:6423-6536`) silently change results if those exports are unallocated. Handled by step 6 (decision 1).
3. **Alarms.** Moist runs on its own `RUN_DT` via `RUNALARM` (`:5644-5647`, `:5786-5796`); MAPL3 gives each component its own clock and `timestep` (`diffs-from-mapl2.md` §1, §9). `DBZ_RunAlarm` is created on the clock passed to Initialize (`:5655`) and is read by the microphysics. Both need the alarm treatment Solar/IRRAD needed (see Lessons).
4. **AERO state** for aerosol activation: `ESMF_Attribute*` -> `ESMF_Info*` (as Solar step 18), using `ESMF_InfoGetCharAlloc` / `ESMF_InfoGetAlloc` for allocatable outputs (SU2G step 7).
5. **Edge fields**: 39 Edge specs in the driver (`PLE`, `ZLE`, `KH`, `PFL_*`/`PFI_*`, mass fluxes, ...). MAPL2 indexing is mixed: most code is 0-based (`PLE(:,:,0:LM-1)`, `ZLE(:,:,LM)` at `:5887-5896`), but `TMP3Dp1(IM,JM,LM+1)` sums `UMF_DC/UMF_SC` 1-based (`:6117-6131`). Check each site; don't remap blindly.
6. **`MAPL_UNDEF` changed** from `1.0e15` to `-huge` (SU2G step 6). Moist compares imports against it (`SNOMAS /= MAPL_UNDEF` `:5841`, `CAPE .ne. MAPL_UNDEF` `:6007-6018`, `DETR_SC /= MAPL_UNDEF` `GEOS_UW_InterfaceMod.F90:587`, `v2 .eq. MAPL_UNDEF` `Process_Library.F90:3829`). Equality tests only work if the provider fills with MAPL3's value; inputs converted from MAPL2 restarts carry `1e15`. Grep all scheme modules and `Process_Library.F90` for inequality tests as well.
7. **Not checked yet**: where `Lightning_mod` (`use` at `:32`) comes from and whether it is MAPL3-clean; whether `pyMoist/.../moist_workarounds.F90` (in `CMakeLists.txt` srcs) builds against MAPL3; whether the "Q first" restart-order rule (`:203`) still matters with MAPL3's netCDF restarts.

## Steps

### 0. MAPL2 baseline (do first)
- Capture at C12L181 with `COMPONENT_TO_RECORD: MOIST`, before and after Run for import, internal and export, as for Solar/SU. Runs start on the hour; record a later step so a first-step crash doesn't lose the files (SU2G step 0).
- Configuration: GFDL_1M + GF + UW (production default). Moist `RUN_DT` = model `DT` to keep the alarm out of the first comparison. `USE_AEROSOL_NN: false` first, so activation uses the fixed CCN numbers (`:6058-6061`) and doesn't need GOCART; a second baseline with aerosol activation later.
- Request every Moist export in `HISTORY.rc` (MAPL2 allocates exports only on request; SU's first baseline had wrong totals because of this).
- Record every Moist setting used (scheme options, `SH_MD_DP`, `DETRAIN_INACTIVE_CNV`, `CNV_FRACTION_*`, precip-type flags, `DBZ_DT`).

### 1. Scope the schemes
Port GFDL_1M, GF and UW first. BACM_1M, THOM_1M, MGB2_2M and RAS keep their `if` branches but `_FAIL` until ported (decision 2); NSSL_2M is not selectable today and stays out of the build list or is removed (decision 5).

**Done (2026-10-29).** The BACM_1M, MGB2_2M, THOM_1M and RAS `*_Setup` branches in SetServices (`GEOS_MoistGridComp.F90:184-205`) `_FAIL` before their call. SetServices runs first, so this one gate also covers their Initialize/Run/Finalize branches, which are unchanged.
- `CLDMICR_OPTION` still defaults to `"BACM_1M"` (`:166`), so a config without it now fails; MAPL3 configs must set `GFDL_1M`. Not changed (would be a logic change).
- NSSL_2M and pyMoist are still in `CMakeLists.txt` (`:13`, `:21`) pending decision 5.

### 2. Macro pass (no logic change)
Decided 2026-10-28: **no re-indentation, no lowercasing, no `codee format`** for Moist (unlike Solar step 1 and SU2G steps 1-1c). Existing indentation and keyword case stay as they are; new or rewritten code follows the surrounding style.

Changed 2026-10-29 for the driver only: `GEOS_MoistGridComp.F90` was formatted with `codee format --verbose --config ~/config/codee-format` in its own commit `f9b25d72` (after step 4), so it can be reverted (`git revert f9b25d72`) or skipped in blame (`git blame --ignore-rev f9b25d72`). Later edits to the driver follow codee's style; the scheme interface modules keep theirs. A second comment-only commit `ca784434` removed the 25 decorative separator lines (`!====`, `!----`, `!!!!`); revert or blame-ignore it the same way. A third commit ("Moist: comment cleanup in GEOS_MoistGridComp.F90", also holding this note) is comments and blank lines only: the 19 `! !` prefixes (`! !INTERFACE:`, ...) -> `!`, ProTeX blocks tightened and their text re-wrapped, `ErrLog Variables` and the file-header macro comment dropped, a stray `!EOP` removed. Verified layout/case only: comments and cpp lines identical; code tokens identical after normalizing codee's equivalent spellings (`endif`/`enddo` -> `end if`/`end do`, `.ne.` etc. -> `/=` etc., `real(ESMF_KIND_R8)` -> `real(kind=ESMF_KIND_R8)`). Without `--config`, codee uses the first `.codee-format` up the tree (here `~/.codee-format`, which differs from `~/config/codee-format` in `IgnoredPreprocessorDirectives`; same output for this file).

Changed again 2026-10-30, for every Fortran file the port touches: commit `48d19566` (after step 6, formatting only) ran the same `codee format --config ~/config/codee-format` on the driver, `GEOS_GFDL_1M_InterfaceMod.F90`, `GEOS_GF_InterfaceMod.F90`, `GEOS_UW_InterfaceMod.F90` and `Process_Library.F90`. Revert or blame-ignore it like `f9b25d72`. All later edits to these files follow codee's style.
- UW (`#ifdef UWDIAG` inside the `compute_uwshcu_inv` argument list) and Process_Library (`#ifdef _CUDA` before `subroutine pdffrac` and `subroutine pdfcondensate`) don't parse, so codee skips them by default. They were formatted with `--on-error force`. That mode indented 16 preprocessor lines and recased `PDFDIAG` to `PDFDiAG` in some, so all preprocessor lines were put back verbatim. The two UW continuation lines after `#ifdef UWDIAG`/`#endif` were re-indented by hand.
- Checked for all five files: preprocessor lines byte-identical to the previous commit; code tokens identical after normalizing case, whitespace, continuations, semicolons, legacy relational operators, `::`, `kind=`, end-statement names, and literals split with `//` (GFDL_1M's long `WARNING_LABEL`s; characters unchanged). Only comments were re-wrapped.
- Any later re-run of codee on UW or Process_Library needs the same `--on-error force` and preprocessor-line repair.

Driver and the three in-scope interface modules: macros `VERIFY_`/`RETURN_`/`__RC__` -> `_RC`/`_RETURN`/`_VERIFY`, `RETURN_(ESMF_SUCCESS)` -> `_RETURN(_SUCCESS)`, `MAPL_Generic.h` -> `MAPL.h`, drop `Iam`/`COMP_NAME`. Verify that only macro spellings changed (normalize both versions and diff). Remove the duplicated `---CLDMACRO` timer only in the functional step.

**Done (2026-10-29).** Done with a script: `RC=STATUS )` followed by `VERIFY_(STATUS)` on the same line (after `;`), the next line, or after blank lines was merged into `_RC)` (~1730 sites); remaining `VERIFY_` -> `_VERIFY`.
- Left as `_VERIFY(STATUS)`: 6 `!!$`-commented sites in the driver, the duplicate verify after `MAPL_TimerAdd(GC, name="--GFDL_1M")`, and the `compute_radar_diagnostics` call in GFDL_1M (passes `STATUS` positionally).
- `Iam`/`COMP_NAME` dropped everywhere; `ESMF_GridCompGet` calls that also get `CONFIG`/`VM` (driver Initialize, Run) keep those. GF keeps `COMP_NAME` in `GF_Setup` because `FRIENDLIES%CNV_TR` uses it (`GEOS_GF_InterfaceMod.F90:58-66`); it goes with the services rework (step 4).
- Check: macros expanded, comments/whitespace/`&`/`;` stripped, token streams of old and new diffed. Only the `Iam`/`COMP_NAME` removals differ.

### 3. State specs
- `Moist_StateSpecs.rc` for the driver's 58 imports and 553 exports; `mapl_acg(... IMPORT_SPECS EXPORT_SPECS INTERNAL_SPECS GET_POINTERS DECLARE_POINTERS)`.
- Layout (decision 3): the driver's specs in `Moist_StateSpecs.rc`; each scheme's fields as hand-written MAPL3 `MAPL_GridCompAddSpec` calls in its `X_Setup`, where MAPL2 registers them.
  - `Moist_StateSpecs.rc`: the driver's imports and exports. The driver has no internals, as in MAPL2.
  - Scheme fields, as each `X_Setup` registers them in MAPL2: GFDL_1M 12 internals, including the 9 every microphysics scheme declares (`Q`, `QLLS`, `QLCN`, `CLLS`, `CLCN`, `QILS`, `QICN`, `NACTL`, `NACTI`; identical in GFDL_1M, BACM_1M, THOM_1M, MGB2_2M), and 5 exports; GF 3 internals; UW 1 internal. Too few for ACG to save anything, and no `COND` needed: `X_Setup` only runs for the selected scheme.
  - `Q` stays the first internal as in MAPL2: it is the first call in `GFDL_1M_Setup`, and the microphysics `Setup` runs before the convection ones (whether order still matters is issue 7).
  - Later microphysics schemes repeat the 9 shared internals and `QRAIN`/`QSNOW`/`QGRAUPEL`, as MAPL2 does. `NCPL`/`NCPI` differ between THOM_1M and MGB2_2M (`DEFAULT` 50e6/1e3 vs 0).
  - Pointers in the scheme modules: hand-written `MAPL_StateGetPointer`, as in MAPL2 (step 5). Don't include the driver's generated pointer headers there (they declare every driver field).
  - Scheme selection: replacing the 41 option-string compares (`adjustl(CONVPAR_OPTION)=="GF"`; SetServices `:206-212`, Initialize `:5636-5637`, Run `:6082-6095`, Finalize `:6742`) with module-level logicals `use_gfdl_1m`, `use_gf`, `use_uw`, ... set once after the options are checked (`:154-182`) is optional cleanup (step 4).
  - Considered and dropped (2026-10-28/29), each step closer to MAPL2:
    - Option B, one `Moist_StateSpecs.rc` with `COND` rows per scheme: pointer headers are split only per state, not per scheme (`mapl_acg.cmake:57-69`); scheme rows away from scheme code; a missing `COND` registers (and restarts) a field silently; same-name rows need distinct `ALIAS`es (`MAPL_GridCompSpecs_ACG.py:375-403, :569-574`).
    - The 9 shared internals in `Moist_StateSpecs.rc`, included before the scheme `Setup` calls: moved them away from where MAPL2 registers them.
    - One spec file and `mapl_acg` call per scheme (`GFDL_1M_StateSpecs.rc`, ...; allowed, `mapl_acg.cmake:124-139`, precedent `@FVdycoreCubed_GridComp/CMakeLists.txt:121-122`): only 1-17 fields each, so inline calls are simpler and match MAPL2.
- Column mapping (Solar step 3, SU2G step 2): `FRIENDLYTO` dropped (services, step 4; `ADD2EXPORT` where it meant "export myself"), `DEFAULT` -> `FILL`, internals `FILL=0.0` (MAPL2 zero-filled), `MAPL_RestartSkip` -> `SKIP`, `MAPL_RestartOptional` -> blank, `AVERAGING/REFRESH_INTERVAL` dropped.
- `AERO` (`MAPL_STATEITEM_STATE`) and the `MTR`/service items stay hand-written `MAPL_GridCompAddSpec` calls (ACG `ITEMTYPE` only does F/V).

**Done (2026-10-29).** `Moist_StateSpecs.rc` generated by script from the MAPL2 calls (56 imports, 553 exports) and included in SetServices (`Moist_Import___.h`, `Moist_Export___.h`) where those calls were; one `mapl_acg` call in `CMakeLists.txt`. The 21 scheme calls rewritten as MAPL3 `MAPL_GridCompAddSpec` calls in place: `GFDL_1M_Setup` 12 internals (`Q` first, MAPL2 order) and 5 exports, `GF_Setup` 3 internals, `UW_Setup` 1 internal. 630 = the 630 converted calls.
- `mapl_acg` passes `INTERNAL_SPECS` although the driver has no internals (empty `Moist_Internal___.h`, not included): `mapl_acg` declares the per-state pointer outputs for all three states (`mapl_acg.cmake:94-98`), and ACG only writes those for the states it is given (`MAPL_GridCompSpecs_ACG.py:654,:680-683`).
- Still hand-written MAPL2 calls: `MTR` (bundle) and `AERO` (state) imports (step 4).
- `COND`: `WSUB_OPTION == 0` (`ALH`, `TAUOROX`, `TAUOROY`), `== 1` (`WSUB_CLIM`), `== 2` (`KM`, `RI`), `MOIST_USE_NCLOUD_CLIM` (`NCPL_CLIM`, `NCPI_CLIM`), replacing the `select case`/`if` around those calls. Both are module variables (Process_Library, driver), so the same expressions work in Run's GetPointer include.
- `COND: PDF_DIAG` on the 15 `PDF_*` exports, which MAPL2 registered under `#ifdef PDFDIAG`. ACG can't emit preprocessor conditions, so the driver module has `logical, parameter :: PDF_DIAG`, `.true.` only when `PDFDIAG` is defined. `PDFDIAG` is never defined (commented-out `#define` in the driver `:5`, `cloudnew.F90:4`, `GEOS_MGB2_2M_InterfaceMod.F90:4`; none in CMake), and only MGB2_2M fills these exports (`GEOS_MGB2_2M_InterfaceMod.F90:1554-1656`).
- Changes from MAPL2 beyond the column mapping:
  - `RI` is both an import (Richardson number) and an export (ice effective radius); the export row has `ALIAS: CLDREFFI` (the name all microphysics schemes use) so the generated pointers don't clash.
  - `CU2DRAINMOVE`, `CU2DSNOWMOVE`: `xy` with `VLocationCenter` in MAPL2; ACG rejects that (`validate_vstagger`), so `VLOC: N`.
  - Long names lose trailing spaces; short names lose MAPL2's padding (`'DQDT  '`).
  - 9 exports have empty MAPL2 units (`SRF_TYPE`, `CNV_FRC`, `STOCH_CNV`, `SIGMA_DEEP`, `SIGMA_MID`, `RAS_*`); ACG omits `units`, which MAPL3 treats as unchecked (matches any units, `@MAPL/superstructure/generic/specs/UnitsAspect.F90:64-67,:98-102`).
  - `RADSW` and `NLIM_NUC` had `RC=STATUS` with no verify; the generated calls check `_RC`.
  - The 6 commented-out (`!!$`) export calls were dropped; the section comments among the export calls are `#` rows in the export table.
- Checked: ACG (Python 3.10+; system `python3` is 3.6) generates the headers; a script compared every generated and inlined `MAPL_GridCompAddSpec` with its MAPL2 call (name, units, dims, stagger, restart, fill, long name, condition). Only the differences above.
- Left for later: `RFRSHINT`/`AVRGNINT`/`DT` reads in SetServices are now unused (step 4 resources); `FRIENDLIES` assignments in `GFDL_1M_Setup`/`GF_Setup` (step 4 services); pointer headers not used yet (steps 5, 8). Before including the driver's pointer headers in Run, check generated names against Run's local variables.

### 4. SetServices
- `MAPL_GridCompSetEntryPoint` with `phase_name='run'`; keep Initialize's default name (MAPL calls it as `GENERIC::INIT_USER`).
- Resources via `MAPL_GridCompGetResource(gc, "LABEL", ...)` (no colon), strings `character(len=:), allocatable`.
- Advection subscription for the microphysics internals (`TRADV`, `service_items`), as SU/DU/SS.
- Convective-transport service export `MTR` (`itemtype=MAPL_STATEITEM_SERVICE`, as DynCore's `TRADV`); `CNV_Tracers_Init` reads the service bundle instead of `FriendlyToMOIST` Info.
- `MAPL_GenericSetServices` removed.

**Done (2026-10-29).**
- Signature `SetServices(GC, RC)`: no intent on `GC`, `RC` `intent(out)` non-optional (`@MAPL/infrastructure/esmf/ESMF_Interfaces.F90:32-37`).
- Entry points: Run with `phase_name='run'`; Initialize and Finalize keep the default names (`GENERIC::INIT_USER`, `GENERIC::FINALIZE_USER`, `@MAPL/superstructure/generic/MethodPhasesMap.F90:151-157`; `OuterMetaComponent/finalize.F90:27-34` runs the user Finalize by that name).
- Resources: `MAPL_GridCompGetResource(GC, 'LABEL', var, ...)`, no colon; `CONVPAR_OPTION`/`SHALLOW_OPTION`/`CLDMICR_OPTION` are now `character(len=:), allocatable` module variables. `SIGMA_DX`/`SIGMA_EXP` pass `default=(X)` (a copy: `value` is `intent(inout)` and must not alias `default`); they had `RC=STATUS` with no verify, now `_RC`. A missing key without a default is an error (`hconfig_get_value_template.h:13-15`). `RUN_DT`, `REFRESH_INTERVAL` and `AVERAGING_INTERVAL` reads dropped (only used by the old spec arguments).
- Timers: `MAPL_TimerAdd` dropped from SetServices and the `X_Setup`s (MAPL3 profiler timers are created on first start, `OuterMetaComponent/start_timer.F90`; GWD has no TimerAdd). The duplicated `---CLDMACRO` timer went with them. Starts/stops: step 5.
- `AERO` import: `MAPL_GridCompAddSpec(..., itemtype=MAPL_STATEITEM_STATE)` as Solar.
- `MTR`: the MAPL2 bundle import is now an export with `itemtype=MAPL_STATEITEM_SERVICE` (as DynCore's `TRADV`). MAPL puts the service bundle in both the import and the export state (`@MAPL/superstructure/generic/specs/ServiceClassAspect.F90:166-204`), so Run's `ESMF_StateGet(IMPORT,'MTR',TR)` is unchanged (DynCore reads its `TRADV` the same way, `DynCore_GridCompMod.F90:811-817`). Subscribers' fields keep their own short names (no `COMP::` prefix, `FieldClassAspect.F90:518-528`); `CNV_Tracers(:)%CNAME/QNAME` are set but unused.
- `CNV_Tracers_Init` (`Process_Library.F90`): no `FriendlyToMOIST` check; every field of the service bundle is a subscriber's, so all are friendly. Per-field metadata (`ScavengingFractionPerKm`, `SetofHenryLawCts`) is set with `ESMF_InfoSet` by the MAPL3 SU2G/SS2G ports, so it is still there; `ESMF_AttributeGet` -> `ESMF_Info*` is step 8.
- Advection: one `TRADV` import with `itemtype=MAPL_STATEITEM_SERVICE`, `service_items` = a `StringVector` the scheme `Setup`s fill (one subscription per component; items must be the component's internals, SU2G notes step 4). `GFDL_1M_Setup(GC, ADVECTED, RC)` adds the 10 internals MAPL2 made `FRIENDLYTO=DYNAMICS` (`Q`, `CLLS`, `CLCN`, `QLLS`, `QLCN`, `QILS`, `QICN`, `QRAIN`, `QSNOW`, `QGRAUPEL`); `GF_Setup(GC, ADVECTED, RC)` adds `CNV_TR` when `CONVECTION_TRACER == 1`. `UW_Setup(GC, RC)`. Unconnected `TRADV` imports are accepted standalone (SU2G notes step 4). TURBULENCE/CHEMISTRY/ANALYSIS friendliness (Q and the condensates): no MAPL3 service yet.
- `FRIENDLIES` types/variables removed from GFDL_1M and GF; `GF_Setup` no longer gets `COMP_NAME`. The unported schemes' `*_Setup(GC, CF, RC)` are unchanged (unreachable after their `_FAIL`); SetServices declares `CF` only for those calls.
- `CNV_TR` (decided 2026-10-29):
  - MAPL2: `FRIENDLYTO="MOIST"` when `CONVECTION_TRACER == 0` (default), `"MOISTDYNAMICS"` when `== 1`. MAPL2 splits `FRIENDLYTO` on `:` into `FriendlyTo<token>` (`main` checkout, `@MAPL/generic/MAPL_Generic.F90:6890-6916`), so `== 1` made it friendly to no component: neither advected nor convectively transported.
  - MAPL3, `== 1`: as intended (MOIST:DYNAMICS): subscribed to `TRADV` and convectively transported. Differs from MAPL2 only in this non-default setting.
  - Convective transport by Moist itself: a provider's bundle only gets subscribers' items (`ServiceClassAspect.F90:268-304`), and Moist can't subscribe to its own service, so Run adds the `CNV_TR` field to the `MTR` bundle before `CNV_Tracers_Init` (step 8).
  - ConvPar only touches `CNV_TR` when `CONVECTION_TRACER == 1` (`ConvPar_GF2020.F90:95,461,613,1253,1679`); with `== 0` it stays at its fill value 0, so transporting it is probably a numerical no-op (not checked). Step 8 may add it only when `== 1`.
- Not done (optional): the `use_*` logicals replacing the 41 option-string compares.
- Open: several subscribers with the same item name (e.g. two instances) in one service bundle; `CNV_Tracers_Init` looks fields up by name (`ESMF_FieldBundleGet(fieldName=)`).

### 5. Interface modules
For GFDL_1M, GF, UW: `X_Initialize(MAPL, ...)` -> take `gc`; `MAPL_GetResource(MAPL, ...)` -> `MAPL_GridCompGetResource`; `MAPL_GetObjectFromGC`/`MAPL_Get` -> `MAPL_GridCompGet(gc, num_levels=, grid=)` (never `ESMF_GridCompGet(grid=)`, Solar step 18); `MAPL_GetPointer` -> generated pointers or `MAPL_StateGetPointer`; timers -> `MAPL_GridCompTimerStart/Stop`; `ALLOC=` per step 6.

**Done (2026-10-29).** pyMoist not ported (see below).
- `X_Initialize(GC, CLOCK, IMPORT, EXPORT, RC)` (was `(MAPL, CF, ...)`); the driver's three calls updated. `X_Run`/`X_Finalize` already took `GC`.
- As Solar (`GEOS_SolarGridComp.F90:880-883`): `MAPL_GridCompGet(GC, num_levels=LM, grid=GRID)`, `MAPL_GridGet(GRID, IM=, JM=)`, `MAPL_GridGetCoordinates(GRID, longitudes=, latitudes=)`, `MAPL_GridCompGetInternalState(GC, INTERNAL)`. `LATS`/`LONS` are now `allocatable` (the R4 getter returns allocatables, `@MAPL/infrastructure/geom/GridGet.F90:128-131`). UW_Run's second, redundant `MAPL_Get` dropped.
- Timestep: `MAPL_ClockGet(CLOCK, dt=)` (Solar `:656`, SU2G step 5; `MAPL_Generic.F90:1296-1309`: the clock's timestep in R4 seconds) replaces the `RUNALARM` ring interval; same value when Moist's MAPL3 `timestep` equals its MAPL2 `RUN_DT`.
- 302 `MAPL_GetResource(MAPL, v, 'L:', ...)` -> `MAPL_GridCompGetResource(GC, 'L', v, ...)` (GFDL_1M 41, GF 206, UW 55), by script; checked in order against the old calls (variable, label, default). Seven self-defaults (`CCW_EVAP_EFF`, `CCI_EVAP_EFF`, `USE_AEROSOL_NN`, `SGS_W_TIMESCALE` x2, `JASON_UW`, `JASON_MFD_SC`) pass `default=(X)` (copy; `value` is `intent(inout)`). `GF_ENV_SETTING` is a `len=10` ConvPar variable (`ConvPar_GF_Shared.F90:97`): read into a local allocatable string and assigned (3 sites).
- 173 `MAPL_GetPointer(S, V, 'N', _RC)` without `ALLOC` -> `MAPL_StateGetPointer` (GFDL_1M 75, GF 53, UW 45), checked the same way. The 181 `ALLOC=.TRUE.` calls are unchanged (step 6).
- Edge fields remapped to `0:LM` right after their fetch, through a `contiguous` scratch `p3d` as IRRAD (`GEOS_IrradGridComp.F90:515-522, :648-667`). MAPL2 allocated every edge field 0-based, so remapping every one reproduces MAPL2 whatever the indexing (issue 5): GFDL_1M `ZLE`, `PLE`, `KH`; GF `PLE`, `ZLE`, `KH`, `PLE_DYN_IN`; UW `ZLE`, `PLE`, `TKE`, and (under `#ifdef UWDIAG`, `associated()`-guarded since exports may be null) `WUP_SC`, `QTUP_SC`, `THLUP_SC`, `THVUP_SC`, `UUP_SC`, `VUP_SC`.
- Timers: `MAPL_TimerOn/Off(MAPL, ...)` -> `MAPL_GridCompTimerStart/Stop(GC, ..., _RC)` (16; GF/UW had no `_RC`).
- `ESMF_GridCompGet(CONFIG=CF)` dropped where `CF` was unused (GFDL_1M, GF Run).
- Not ported (decided 2026-10-29: don't port pyMoist now): the `USE_PYMOIST_*` branches still call `MAPL_pybridge_*` with `MAPL` and `MAPL_ConfigSetAttribute(CF, ...)`, and `GFDL_1M_Finalize`, `GF_Finalize`, `UW_Finalize` (pyMoist only) still use `MAPL_GetObjectFromGC`. `MAPL`/`CF` are no longer declared in the converted routines, so these lines won't compile against MAPL3; decide at step 9 (e.g. `_FAIL` when `USE_PYMOIST_*`).
- For later steps: GF and UW create their own run alarms (`GF_RunAlarm`/`UW_RunAlarm`, period `GF_DT`/`UW_DT`, default Moist's timestep) on the clock in `X_Initialize` and look them up by name in `X_Run` (step 7). The edge exports fetched with `ALLOC` need the same remap (step 6): GFDL_1M `PFL_AN`, `PFL_LS`, `PFI_AN`, `PFI_LS`, `PFR_LS`, `PFS_LS`, `PFG_LS`; GF `UMF_DC`, `PFI_CN`, `PFL_CN`, `WQT_DC`; UW `UMF_SC`, `QTFLX_SC`, `SLFLX_SC`, `UFLX_SC`, `VFLX_SC`.

### 6. `ALLOC=.TRUE.` replacement
Per decision 1: local work arrays copied to the export if associated; internals with `ADD2EXPORT` (and `RESTART: SKIP`) for fields other components always connect. Apply to the driver and the in-scope modules. Keep the epilogue's results independent of which exports are connected.

**Survey (2026-10-29):** `mapl3-alloc-exports.md` classifies all 197 ALLOC'd exports (226 sites).
- **Group C (123), local to one routine.** Use a local allocatable named after the variable plus a `_ptr` export pointer and an `if (associated(X_ptr)) X_ptr = X` copy, as IRRAD's `CLDTT` (`GEOS_IrradGridComp.F90:2336-2337, 2401-2407, 2767`). The 22 fields that rely on the zero fill get a zero-initialized local. C3 (8) only drops `ALLOC`.
- **Group B (25), persisting across steps.** Internals with `ADD2EXPORT` (decision 1b).
- **Group A (49), shared within a step.** Internals with `ADD2EXPORT`, as group B (decision 1c). To be revisited once Tom's graph implementation is online.

**Done (2026-10-30).** All 226 `ALLOC=.TRUE.` sites in the driver, GFDL_1M, GF and UW are gone. Not built (step 9).
- Spec: the 74 A/B rows moved from `EXPORT` to a new `INTERNAL` table at the end of `Moist_StateSpecs.rc` (`RESTART: SKIP`, `FILL: 0.0`, `ADD2EXPORT: T`; `RI` keeps `ALIAS: CLDREFFI`). SetServices includes `Moist_Internal___.h` after the scheme `Setup`s, so `Q` stays the first internal.
- A/B: every fetch of those names, with or without `ALLOC`, now reads `INTERNAL` (`MAPL_StateGetPointer`), including the driver's guarded epilogue sums. The `associated()` guards there are kept; they are now always true. Edge internals are remapped to `0:LM` where the code indexes by level: `UMF_DC`, `PFI_CN`, `PFL_CN` (GF), `UMF_SC` (UW), `PFL_AN`, `PFL_LS`, `PFI_AN`, `PFI_LS` (GFDL_1M), and the six `PF*` reads in the driver's `DPDTMST` sum. The driver's whole-array uses (`CNV_MFC`, `PFL_LSAN`, `PFI_LSAN`) are not remapped.
- C (115 locals): `allocate(..., source=0.0)` where the ALLOC fetch was (zero fill as MAPL2, so C1 needs no special case), then one `X_ptr` fetch and `if (associated(X_ptr)) X_ptr = X` copy per export at the end of the block that wrote it (GF and UW: end of the alarm-ringing block; GFDL_1M: end of the non-pyMoist branch; driver: end of the ringing branch, `PTYPE` right after `DIAGNOSE_PRECIP_TYPE`). Edge locals (`WQT_DC`, `PFR_LS`, `PFS_LS`, `PFG_LS`, `QTFLX_SC`, `SLFLX_SC`, `UFLX_SC`, `VFLX_SC`) are allocated `0:LM`. No `FILL: 0.0` on the C1 export rows: the copy is unconditional, so a connected export always gets the zeroed local.
  - Exceptions: `MLCAPE`, `MLCIN`, `MUCAPE`, `MUCIN` are pointers allocated locally and deallocated after the copy, because `BUOYANCY2` takes them as pointer dummies and tests `associated()` (`Process_Library.F90:1277, 1287-1288`). `UMST0`/`VMST0` have no named variable and one use each, so they are guarded writes to the export. `CFICE`, `RAIN`, `SNO`, `TPREC` lost their `if (associated(...))` wrappers (always true under MAPL2's ALLOC); `RHCRIT3D`'s `associated()` test in GFDL_1M became unconditional.
  - C2: the epilogue re-fetches of `DTDT_ER`/`DQVDT_ER` now add the locals. `ER_PRCP` was computed through `PTR2D`; it is now a local of that name.
- C3: `ALLOC` dropped (`PREC_STRAT`, `PREC_CONV`, `*AFMOIST`).
- The driver's Run is otherwise still MAPL2 (`MAPL_GetPointer`, `MAPL_Get`, timers); step 8 converts the rest.
- Not checked: whether MAPL3 accepts an internal and an import with the same short name (`RI`; the export and import already coexisted). Seen but not changed: `DIAGNOSE_PRECIP_TYPE` receives `PLS`/`PCU` unguarded; they are plain exports, null in MAPL3 when unconnected (step 8).

### 7. Initialize
Drop `MAPL_GenericInitialize`. Timestep from `MAPL_ClockGet(clock, dt=)` instead of `RUNALARM` (SU2G step 5). Recreate `DBZ_RunAlarm` the way Solar recreated its run alarm (Solar step 18), and check it rings on the first step (IRRAD "`irrad_lw_alarm` never rang on the first step"). Decide whether Run's `ESMF_AlarmIsRinging` branch (`:5796`) is still needed once MAPL3 schedules Moist with its own `timestep`.

**Done (2026-10-30).** Not built (step 9).
- No `MAPL_GetObjectFromGC`/`MAPL_GenericInitialize`. Resources via `MAPL_GridCompGetResource` (no colon); self-defaults (`DETRAIN_INACTIVE_CNV`, `DBZ_DT`, `DBZ_VAR_INTERCP`) pass `default=(X)`.
- Timestep from Moist's clock (`ESMF_ClockGet(timeStep=)`, `MAPL_ClockGet(dt=)`). MAPL3 gives each component its own clock with its own `timestep` (`@MAPL/superstructure/generic/OuterMetaComponent/initialize_set_clock.F90:32-45`), set before `GENERIC::INIT_USER` (`@MAPL/enums/GenericPhases.F90:37-39, 85`), and advances it only on steps the component runs (`run_clock_advance.F90:31-44`). So alarms created on it ring as they did on MAPL2's heartbeat clock when Moist's `timestep` equals the heartbeat (MAPL2's `RUN_DT`).
- `DBZ_RunAlarm`: kept exactly as MAPL2 creates it (`RingTime=currTime - TINT`, `sticky=.false.`), with `TINT` now the clock's timestep. Unlike Solar/IRRAD, this was already a component-created ESMF alarm in MAPL2 (their bug was replacing MAPL2's generic main alarm), so the same call reproduces MAPL2's first-step behaviour whatever it is. Not run-tested. `GF_RunAlarm`/`UW_RunAlarm` (step 5) are likewise unchanged from MAPL2 (`ringTime=currentTime`).
- Run's `ESMF_AlarmIsRinging` branch: not needed. MAPL3 calls the user Run only when the component's run alarm rings (`OuterMetaComponent/run_user.F90:35-40`), so step 8 removes it.
- The RAS, BACM_1M, THOM_1M and MGB2_2M `*_Initialize(MAPL, ...)` calls became `_FAIL` (the old call kept in a comment): they take MAPL2's MetaComp, and SetServices already fails for those schemes (step 1).

### 8. Run and Finalize
- `MTR` service bundle, as DynCore handles its `TRADV` (`@FVdycoreCubed_GridComp/DynCore_GridCompMod.F90:809-817`): get it from the EXPORT state (Moist is the provider; today's `ESMF_StateGet(IMPORT,'MTR',TR)` only works because MAPL adds the bundle to both states) and `_ASSERT(MAPL_FieldBundleSameData(...))` against the IMPORT copy. No copying: the fields are the subscribers' internals, and `CNV_Tracers_Init` takes data pointers into them (`ESMFL_BundleGetPointerToData`), so GF's transport updates them in place, as FV3 does through `PULL_Q` (`:2525-2600`). The `COMP::NAME` parsing in `CNV_Tracers_Init` is dead (fields keep their short names). No exclusion list (DynCore's `xTRADV`, `:2628-2674`) unless a scheme needs one.
- Add Moist's own `CNV_TR` internal to the `MTR` bundle before `CNV_Tracers_Init` (step 4, `CNV_TR`).
- Edge-field remaps per issue 5; `p3d => X; X(1:IM,1:JM,0:LM) => p3d` through a `contiguous` scratch pointer (IRRAD notes).
- `MAPL_UNDEF` comparisons per issue 6.
- `Aer_Activation`: `ESMF_Attribute*` -> `ESMF_Info*` (issue 4).
- `MAPL_MaxMin` debug prints: now a function returning `[max, min]` (SU2G step 9).
- Finalize: drop `MAPL_GenericFinalize`.

**Done (2026-10-30).** Not built (step 9).
- Driver Run: `MAPL_GridCompGet(num_levels=, grid=, logger=)`, `MAPL_GridGet(IM=, JM=)`, `MAPL_GridCompGetInternalState`, `MAPL_ClockGet(dt=)`; every `MAPL_GetPointer` -> `MAPL_StateGetPointer` (`LFR_GCC`'s `NotFoundOk` dropped: it is registered); timers -> `MAPL_GridCompTimerStart/Stop`; `PMIN_CBL` via `MAPL_GridCompGetResource`.
- Alarm: the `RUNALARM` check and its not-ringing branch are gone (step 7); the ringing branch is dedented. Behaviour change only when Moist's `timestep` is longer than its parent's: MAPL2 then refilled `*AFMOIST` and `RH2` from the current state on the off steps, MAPL3 leaves them at Moist's last values. The AGCM runs Moist every heartbeat.
- Edge remaps to `0:LM`: imports `PLE`, `ZLE`, `KH`, `TKE` after their fetch; `PFI_SC`/`PFL_SC` inside their `associated()` guard in the `DPDTMST` sum (step 6 did the internals). `PREF` (1D) and the whole-array uses (`CNV_MFC`, `PFL_LSAN`, `PFI_LSAN`, `UMF_*`) need none.
- `MTR`: the export bundle, `_ASSERT(MAPL_FieldBundleSameData(...))` against the import copy, as DynCore's `TRADV`. GF's `CNV_TR` internal is added to a local bundle (`MTR`'s fields plus `CNV_TR`, created and destroyed around `CNV_Tracers_Init` each step; `ESMF_FieldBundleDestroy` leaves the fields) rather than to MAPL's `MTR` bundle. Added whenever GF is selected, as MAPL2's default `FRIENDLYTO="MOIST"` made it a convective tracer.
- `CNV_Tracers_Init` (`Process_Library.F90`): the 12 `ESMF_AttributeGet` presence checks and 11 reads -> `ESMF_InfoIsPresent`/`ESMF_InfoGet` on the field's Info (root keys, as SU2G/DU2G/SS2G set them, e.g. `SU2G_GridCompMod.F90:393-396`); `Vect_Hcts`/`KcScal` read with `values=`. `ESMFL_BundleGetPointerToData` (not in MAPL3) -> `ESMF_FieldGet(FIELD, localDe=0, farrayPtr=)`.
- `Aer_Activation` (`aer_actv_single_moment.F90`, issue 4): takes `GC` instead of the MetaComp (timers `MAPL_GridCompTimerStart/Stop`); `AERO` keys via `ESMF_InfoGetFromHost` + `ESMF_InfoGet`/`ESMF_InfoGetCharAlloc` (field names, as GOCART2G reads them, `GOCART2G_GridCompMod.F90:1427-1441`)/`ESMF_InfoGetAlloc` (`aerosol_modes`, a string list; ESMF `ESMF_Info.F90:97-104`); `aerosol_mode` set with `ESMF_InfoSet`; `MAPL_StateGetPointer`; `MAPL.h` macros. Keys match the MAPL3 provider (`@GOCART/ESMF/Shared/GOCART2G_AerosolActivationMod.F90:52-73, 183-228`).
- `compute_radar_diagnostics` (`Process_Library.F90`): its 7 `MAPL_GetPointer` -> `MAPL_StateGetPointer`.
- `PLS`/`PCU`: now local arrays copied to their exports (as group C), always computed: `DIAGNOSE_PRECIP_TYPE` takes them unguarded and `inout`, and in MAPL2 Surface always connected them. The copy is after `DIAGNOSE_PRECIP_TYPE`, which updates them.
- `MAPL_MaxMin` debug prints (`DEBUG_MST`): a module helper `report_maxmin` calls the MAPL3 function with the VM's communicator and logs through the component logger (`use pflogger`, as SU2G).
- `MAPL_UNDEF` (issue 6, decided 2026-10-30): `SNOMAS /= MAPL_UNDEF` -> `SNOMAS /= MAPL2_UNDEF`, a module parameter `1.0e15` as in SU2G/DU2G (`SU2G_GridCompMod.F90:36-37`), since SNOMAS's undefined points come from Surface or from MAPL2-converted inputs. The other `MAPL_UNDEF` uses stay: Moist fills and compares within itself with the same constant (`CAPE` from `BUOYANCY`, `Process_Library.F90:1492`; `DETR_SC` from `uwshcu.F90:1109-1110`; `v2`, `Process_Library.F90:3819, 3849`), and the export-only fills keep MAPL3's `-huge` (`@MAPL/utils/Constants/InternalConstants.F90:18-19`), which MAPL3 History treats as missing. So a MAPL2 baseline comparison (step 10) will differ at undefined points of those exports (`1e15` vs `-huge`). Check in `moist-sa` that the SNOMAS provider fills with `1e15`.
- Finalize: `MAPL_GenericFinalize` dropped.
- Left for step 9: the pyMoist branches and `GFDL_1M_Finalize`/`GF_Finalize`/`UW_Finalize` (pyMoist only) still use the MetaComp (step 5); the unported schemes' modules are MAPL2 and still in `CMakeLists.txt`.
- Not checked: whether `pflogger` needs adding to the Moist library's dependencies (SU2G uses it through MAPL); that MAPL3 tolerates `ESMF_FieldBundleCreate`/`Destroy` on fields it owns, each step.

### 9. Build
- Uncomment `GEOSmoist_GridComp` in `../CMakeLists.txt:5`; `TYPE SHARED`; one `mapl_acg` per spec file (step 3); external `Moist_SetServices` wrapper with a unique name (SU2G step 8).
- Entry points must match `I_Run`: no `intent` on `gc/import/export/clock`, `rc` `intent(out)` non-optional (Radiation notes; SU2G step 4).

### 10. Standalone regression `moist-sa`
After `gwd-ncar-sa`: a configurable `DataForMoist` provider restarted from the step-0 MAPL2 import state, `activate_all_exports: true` for Moist, a stub `AERO` (FakeGOCART) and an empty convective-transport subscriber. Compare with `cmpchk.py`; target bit-for-bit. Then the aerosol-activation baseline with FakeGOCART or GOCART2G.

### 11. Remaining schemes
BACM_1M, then MGB2_2M, THOM_1M, RAS, each with its own baseline.

### 12. Wiring (blocked on the Physics/AGCM ports)
Connect Moist's imports in Physics, the `TRADV` subscriptions, and the convective-transport service to GOCART's SU/DU/SS.

## Lessons carried over (Radiation, GOCART)
- Entry-point signatures must match `I_SetServices`/`I_Run` exactly (Radiation notes; SU2G step 4/9).
- `MAPL_GridCompGetResource` labels have no trailing colon.
- `MAPL_GenericSetServices`/`Initialize`/`Finalize`/`RunChildren` are not called.
- Grid: `MAPL_GridCompGet(gc, grid=...)`, not `ESMF_GridCompGet` (Solar step 18); `MAPL_GridGetInterior` -> `MAPL_GridGet(interior=)`.
- Edge fields are 1-based in MAPL3; remap through a `contiguous` scratch pointer (IRRAD, Radiation, Solar step 18, SU2G step 6).
- `DEFAULT` -> `FILL`; MAPL2 zero-filled internals, MAPL3 does not (Solar step 3, SU2G step 2).
- `MAPL_UNDEF` is `-huge` in MAPL3, `1.0e15` in MAPL2 (SU2G step 6).
- `ESMF_Attribute*` -> `ESMF_Info*`; use the `...Alloc` getters for allocatable outputs (Solar step 18, SU2G step 7).
- Component alarms need MAPL2's main-alarm logic recreated (Solar step 18; IRRAD first-step alarm fix).
- Unconnected child imports propagate up automatically; several subscribers per service provider is by design (SU2G step 4).
- MAPL2 baselines must match the standalone configuration exactly (`AERO_PROVIDER`, requested exports, refreshed internal restart; Radiation notes, Solar step 18, SU2G step 0).

## Open decisions
1. **DECIDED (2026-10-28) - `ALLOC=.TRUE.`**: by default, a local allocatable work array, copied to the export if it is associated. Fields other components always connect (`CNV_FRC`, `RAIN`, `SNO`, `TPREC`, ...) become internals with `ADD2EXPORT` instead (always active, `OuterMetaComponent/advertise_variable.F90:55-57`) and `RESTART: SKIP` unless they must survive a restart.

   **Amended (2026-10-29), after the survey in `mapl3-alloc-exports.md`.** The 197 ALLOC'd exports fall into three groups, A, B and C, described in that file. The local-array default above now applies to group C only. Groups A and B both become internals with `ADD2EXPORT` (1b and 1c below).

   **1b. DECIDED (2026-10-29) - Group B (25 exports) becomes internals.**
   - **Which exports:**
     - `CNPCPRATE`;
     - the GF tendencies `DUDT_DC`, `DVDT_DC`, `DTDT_DC`, `DQVDT_DC`, `DQLDT_DC`, `DQIDT_DC`, `DQADT_DC`;
     - the UW tendencies `DUDT_SC`, `DVDT_SC`, `DTDT_SC`, `DQVDT_SC`, `DQLDT_SC`, `DQIDT_SC`, `DQRDT_SC`, `DQSDT_SC`, `DQADT_SC`;
     - the UW condensate sources `QLDET_SC`, `QIDET_SC`, `QLSUB_SC`, `QISUB_SC`, `QLENT_SC`, `QIENT_SC`;
     - `SHLW_PRC3`, `SHLW_SNO3`.
   - **Why:** each is read on a later step than the one that wrote it.
     - `CNPCPRATE` in every default run: UW runs before GF and reads GF's previous value.
     - The other 24 whenever `GF_DT`/`UW_DT` (`CONV_DT`) is longer than Moist's timestep, which the production setup uses at c1080 and above.
     - 21 of the 25 feed model state.
   - **How:**
     - Move the 25 rows in `Moist_StateSpecs.rc` from `EXPORT` to `INTERNAL`, with `ADD2EXPORT: T`, `RESTART: SKIP` and `FILL: 0.0`.
     - Every `MAPL_GetPointer(EXPORT, ..., ALLOC=.TRUE.)` of these names becomes `MAPL_StateGetPointer(INTERNAL, ...)`. So do the driver's guarded epilogue reads (`GEOS_MoistGridComp.F90:878-999`, and `:1052, 1080, 1103` for `CNPCPRATE`).
     - `SKIP` matches MAPL2, which never restarted these exports. `FILL: 0.0` reproduces MAPL2's ALLOC zero fill.
     - None of the 25 is an edge field or clashes with an import name.
   - **Precedent:**
     - DynCore re-exports internals through the same ACG column, including a 3D edge field (`DynCore_StateSpecs.rc:239-252`, commit `68b4705`).
     - AgcmSimple does it with hand-written `add_to_export=.true.` calls (`GEOS_AgcmSimpleGridComp.F90:96-111`).
   - **To check later (none blocking):**
     - `RESTART: SKIP` with `ADD2EXPORT` is a first. Existing cases use `REQUIRED` or leave it unset. Check it in `moist-sa` (step 10).
     - MAPL3 writes `SKIP` fields to checkpoints anyway. See "Checkpoints" under 1c.
     - The carry-over across steps only happens if `GF_RunAlarm`/`UW_RunAlarm` still ring every `GF_DT`/`UW_DT` on Moist's MAPL3 clock (step 7).
   - **Recorded upstream issues** (kept as they are, for bit-for-bit agreement; see the report):
     - `DQLDT_SC`/`DQIDT_SC` are wrong on steps where UW doesn't run (`GEOS_UW_InterfaceMod.F90:679-680`).
     - `CNPCPRATE` may make restarted runs differ from continuous runs (unverified).

   **1c. DECIDED (2026-10-30) - Group A (49 exports) also becomes internals, as group B.**
   - **Which exports:** A1-A4 in `mapl3-alloc-exports.md` (driver inputs to the schemes, GFDL_1M macro/micro tendencies, precipitation/flux/mass-flux outputs read by the driver epilogue, and `QLTOT`, `QITOT`, `RL`, `RI`).
   - **Why:** each is shared by several routines within one step. The schemes read A1 unguarded, and the epilogue needs A2/A3 whether or not they are connected. Internals reproduce MAPL2 with no copies. The private derived-type state (the alternative in the report) was not chosen.
   - **How:** as 1b (`INTERNAL` rows with `ADD2EXPORT: T`, `RESTART: SKIP`, `FILL: 0.0`; fetches from `INTERNAL`), plus:
     - Internal `RI` needs `ALIAS: CLDREFFI`: import `RI` exists and ACG declares one pointer per short name.
     - 8 are edge fields (`PFL_AN`, `PFL_LS`, `PFI_AN`, `PFI_LS`, `PFL_CN`, `PFI_CN`, `UMF_DC`, `UMF_SC`): edge `VLOCATION` in the spec, and the `0:LM` remap after each fetch (issue 5). DynCore re-exports a 3D edge internal the same way (`DynCore_StateSpecs.rc:239-252`).
   - **Checkpoints (A and B, 74 fields): we don't want these checkpointed, but MAPL3 currently writes them.**
     - `RestartHandler%write` selects fields with `get_restart_bundle(state, is_write=.true., ...)` (`@MAPL/superstructure/generic/RestartHandler.F90:93`), which drops only incomplete fields. `SKIP` is dropped only on read (`:36-37, 204-220`). A filled internal is complete, so it is written.
     - Restarts stay correct (the read skips them); the cost is file size.
     - Fix (MAPL): in `RestartHandler%write`, also apply `SkipRestartPredicate` before `write_bundle_`. Not in `get_restart_bundle`: the in-memory checkpoint saves and restores through the same `is_write=.true.` selection (`OuterMetaComponent/write_restart.F90:132, 151`, `OuterMetaComponent/read_restart.F90:112-113`), and should keep restoring these fields so group B rewinds with the rest of the state.
     - Deferred (2026-10-30): not started. Until then Moist checkpoints carry the 74 fields.
   - **Revisit:** plan to revisit group A once Tom's graph implementation is online.
2. **DECIDED (2026-10-28) - Scheme scope**: GFDL_1M + GF + UW first; BACM_1M, THOM_1M, MGB2_2M and RAS `_FAIL` until ported (step 1).
3. **DECIDED (2026-10-29) - Spec layout, as in MAPL2.** (Replaces option B, 2026-10-28, and a per-scheme spec file version.) `Moist_StateSpecs.rc` for the driver's imports and exports; each scheme's fields, including the internals all microphysics schemes declare, as hand-written `MAPL_GridCompAddSpec` calls in its `X_Setup`, where MAPL2 registers them. No `COND` rows for schemes (step 3).
4. **DECIDED (2026-10-28) - Convective-transport service name: `MTR`.** The advection service kept its MAPL2 bundle name: MAPL2 Physics exports the friendlies-to-DYNAMICS bundle as `TRADV` (`GEOS_PhysicsGridComp.F90:2425`) and the AGCM connects it to SUPERDYN (`GEOS_AgcmGridComp.F90:1138-1139`); MAPL3 DynCore exports the `TRADV` service under the same name (`DynCore_GridCompMod.F90:339-343`). The MAPL2 convection bundle is Moist's `MTR` import (`GEOS_MoistGridComp.F90:488-497`), filled by Physics with the friendlies to MOIST (`GEOS_PhysicsGridComp.F90:2477-2484`).
5. **NSSL_2M and pyMoist**: drop from the MAPL3 build or keep.
