# GEOS_TurbulenceGridComp MAPL3 porting plan

Status (2026-11-13): step 0 done upstream (PR #1582); steps 2-7 done on branch `feature/bmauer/port-turbulence-to-mapl3` (rebased onto `release/MAPL-v3` @ `86b1a7d5`). Step 1 (MAPL2 baseline) to do; steps 8-10 not started. The driver compiles (ifx Debug, by hand against `build-debug`'s MAPL, kernels patched for issue 9 in a scratch copy only); not built by CMake, not run.

Resume here (2026-11-13): next is step 8 (build): issue 9 in the kernels (`use mapl_SatVapor_mod` in `shoc.F90`/`scm_surface.F90`, `MAPL_UNDEF` from `MAPL_Constants` in `edmf.F90`; checked to compile), `Turbulence_SetServices` wrapper, `../CMakeLists.txt`. In parallel (user): the step-1 MAPL2 capture.

Follows the Moist port (`../GEOSmoist_GridComp/mapl3-porting-notes.md`, with `mapl3-alloc-exports.md` and `mapl3-friendlies-to-services.md`) and the rules in the `GEOSgcm_GridComp` and `GEOSagcm_GridComp` `AGENTS.md`. The MAPL3 API lessons listed there ("Lessons carried over") are not repeated here.

An earlier attempt, `origin/feature/bmauer/turbulence_mapl3_port` (2026-10-28, off `release/MAPL-v3` @ `4f76b103`), is not continued: it predates the `develop` restructure of the `TR` bundle (step 0), sets `activate_all_exports`, drops the `FILL` column (it works: Moist uses `FILL: 0.0`), keeps the `FriendlyToTURBULENCE` reads instead of a service, and has no edge-field audit. Its notes are still useful for counts and for the `ESMFL_BundleGetPointerToData` replacement.

## Reference files
- Template component: `../GEOSmoist_GridComp/GEOS_MoistGridComp.F90`, `../GEOSgwd_GridComp/GEOS_GwdGridComp.F90`
- Template state spec: `../GEOSmoist_GridComp/Moist_StateSpecs.rc`, `../GEOSgwd_GridComp/GWD_StateSpecs.rc`
- Template standalone regression: `../regression/moist-gfdl1m-gf-shoc-sa/`, `../regression/run_case.cmake`
- Service provider: Moist's `MTR` (`GEOS_MoistGridComp.F90`, SetServices and Run); DynCore's `TRADV` (`@FVdycoreCubed_GridComp/DynCore_GridCompMod.F90`)
- Short names of service-bundle fields: `get_short_name` (`@FVdycoreCubed_GridComp/fv_shared.F90`)
- Two run phases registered by name: `@GOCART/ESMF/GOCART2G_GridComp/GOCART2G_GridCompMod.F90`; per-phase regression cases: `@GOCART/ESMF/GOCART2G_GridComp/regression/` (`run1`, `run2`)

## Code overview (MAPL2, `develop` @ `f31a2151`)

`GEOS_TurbulenceGridComp.F90` (7280 lines) plus `edmf.F90`, `shoc.F90`, `LockEntrain.F90`, `scm_surface.F90`. Parent: Physics (child `TURBULENCE`). No children.

- **Module header:** the `DIRECT_*` parameters list the five quantities diffused directly rather than through `TR` (`S`, `U`, `V` from the import state, `TKESHOC`, `QT2` from the internal state), in processing order (`S` first, for SHVC); `DEBUG_TRB`.
- **SetServices:** reads `DEBUG_TRB`, `USE_WAVES`, `USE_SEA_SPRAY`, `EDMF_NUMUP` (default 10), `SCM_SL`, `TURBULENCE_DO_SHOC` (label built from `COMP_NAME`); registers Initialize, Run1, Run2; 49 imports (sea-spray pair under `USE_WAVES .and. USE_SEA_SPRAY`), 159 exports (including one `<NAME>IT` per `TRI_increments::` entry, read with `GetTriIncrements`), 45 internals (`cu_scm`... under `SCM_SL /= 0`); timers; `MAPL_GenericSetServices`.
- **Initialize:** generic initialize, then fills the export bundle `TRI` with the tendency exports (`SIT`, `UIT`, `VIT`, `TKESHOCIT`, `QLLSIT`, `QILSIT`, `QVIT` as `QIT`, and the `TRI_increments` tendencies), allocating each (`alloc=.true.`) because History does not allocate bundle members.
- **Run1:** sea-spray copies; internal pointers (matrix coefficients); `REFRESH` when the generic `RUNALARM` rings (period `DT:`, default `RUN_DT:`); always `DIFFUSE` with `RUN_DT:`.
  - `REFRESH` (~2250 lines): PBL heights, EDMF (`RUN_EDMF`), SHOC (`RUN_SHOC`), Louis (`LOUIS_KS[_OPTIMIZED]`), Lock (`ENTRAIN`), Beljaars, LU decomposition (`VTRILU`, `VTRISOLVESURF`). 143 of the 165 `MAPL_GetResource` calls.
  - `DIFFUSE`: loops over the `DIRECT_*` quantities, then over the fields of `TR`. Per `TR` field: attributes `DiffuseLike` (default `Q`), `FriendlyToTURBULENCE` (default false), `WeightedTendency` (default false); scheme by name (`S`, `Q`, `QLLS`, `QILS`, `U`, `V` use the EDMF/SHOC matrices); surface value/flux/sensitivity only for `S`, `U`, `V`, `Q` (`<N>HAT` import, `<N>FLX`, `<N>DFL` exports); tendency into the export found by `GetTendency` (`QVIT` for `Q`, else `<NAME>IT`, if allocated); friendly fields updated in place. Writes `*AFDIFFUSE`.
- **Run2 / `UPDATE`:** surface increments `SDEL`, `UDEL`, `VDEL`, `QDEL` (after Surface's second phase), same direct-then-`TR` loop; dissipation `INTDIS`, `SRFDIS`, `TOPDIS`, KE diagnostics, SHVC, flux diagnostics, `SAFUPDATE`, `DPDTTRB`.
- **Module procedures:** `LOUIS_KS`, `LOUIS_KS_OPTIMIZED`, `ComputeMOScalingInterlevel`, `BELJAARS`, `VTRILU`, `VTRISOLVESURF`, `VTRISOLVE`, `GetTendency`, `GetTriIncrements`.
- Physics fills `TR` in its Initialize with the friendlies to `TURBULENCE` of Physics itself, Moist and Chemistry (`GEOS_PhysicsGridComp.F90` on `develop`, Initialize). Turbulence's own internals are not in `TR`.

## Starting point (survey 2026-11-12, `develop` @ `f31a2151`)

| Item | Count |
|---|---|
| `MAPL_AddImportSpec` / `AddExportSpec` / `AddInternalSpec` | 49 / 159 / 45 |
| `VLocationEdge` specs | 48 (below) |
| `MAPL_GetPointer` | 246 |
| `alloc=` (incl. `ALLOC=PDFALLOC`, `alloc=associated(...)`) | 45 |
| `MAPL_GetResource` | 165 (SetServices 6 incl. `ESMF_Config`, Run1 5, REFRESH 143, DIFFUSE 4, Run2 1, UPDATE 11) |
| `ESMF_Config*` | 14 (`EDMF_NUMUP`, `SCM_SL`, `TURBULENCE_DO_SHOC`, `RUN_DT` x2, `TRI_increments::`) |
| `ESMF_AttributeGet` | 6 (`DiffuseLike`, `FriendlyToTURBULENCE`, `WeightedTendency`, in DIFFUSE and UPDATE) |
| `ESMFL_BundleGetPointerToData` | 2 |
| `MAPL_MaxMin` | 9 (`DEBUG_TRB`) |
| `MAPL_TimerAdd` / `TimerOn` | 18 / 21 |
| `MAPL_UNDEF` equality tests | 7 (`ZPBL2`, `ZPBL10p`, `ZPBLRI`, `KPBL_SC`, all in REFRESH) |
| `FRIENDLYTO` | 11 internals |
| `UNGRIDDED_DIMS` | 3 (`EDMF_PLUMES_W/_THL/_QT`, `NUMUP`) |
| Bundles | `TR` (import), `TRI` (export) |
| `#ifdef _CUDA` | 3 blocks |

Edge fields:
- Imports: `PLE`, `ZLE`, `PREF` (vertical only), `WQT_DC`.
- Internals: `TKH`, `MFAW`, `SLFLXMF`, `QTFLXMF`.
- Exports: `EDMF_PLUMES_W/_THL/_QT`, `EDMF_DRY_*`/`EDMF_MOIST_*` (13), `EDMF_WQT`, `EDMF_WSL`, `EDMF_TKE`, `EDMF_MF`, `KM`, `KH`, `RI`, `DU`, `KHLS`, `KMLS`, `KHSFC`, `KHRAD`, `EKH`, `EKM`, `ALH`, `ALM`, `QTFLXTRB`, `SLFLXTRB`, `UFLXTRB`, `VFLXTRB`, `TKE`, `SHOCPRNUM`, `ISOTROPY`, `ZLES`.

`FRIENDLYTO` internals: `ZPBL`, `ZPBL_SC`, `cu_scm`, `ct_scm`, `ssurf_scm`, `qsurf_scm` (`trim(COMP_NAME)`), `TKH` (`'TURBULENCE'`), `PDF_A`, `TKESHOC`, `QT2`, `QT3` (`FRIENDLIES_SHOC`: `COMP_NAME`, or `'DYNAMICS:TURBULENCE'` when `TURBULENCE_DO_SHOC /= 0`).

## Turbulence-specific MAPL3 issues

1. **`TR` becomes a provided service.** MAPL2's `TR` import is filled by Physics with the friendlies to `TURBULENCE`. In MAPL3 Turbulence provides a service `TRDIFF` (standard name `turbulently_diffused_quantities`, the name Moist already subscribes with; Moist notes step 4), and Moist, GOCART, ... subscribe their internals. Done in steps 4 and 7; still open: short names for tracers without a known standard name (GOCART, `TRI_increments`; `diffused_short_name` returns the standard name).
   - Every field in a service bundle is a subscriber's internal, so all are updated in place: `FRIENDLY` becomes true for every `TR` field and the `FriendlyToTURBULENCE` read goes.
   - Subscribed fields carry MAPL-generated names (Moist notes, "Self-subscription"). `DIFFUSE`/`UPDATE` select by name (`Q`, `QLLS`, `QILS`, the `S`/`U`/`V`/`Q` surface terms) and `GetTendency` looks up `<NAME>IT` by name, so names must be recovered through `standard_name`, as `get_short_name` does. How tracers without a known standard name (GOCART, `TRI_increments`) are named is open.
   - `DiffuseLike`/`WeightedTendency` (`ESMF_AttributeGet`): nobody sets them for `TR` fields on `develop` (Physics only set them on `S`, `U`, `V`, now diffused directly), so the defaults always apply. Read them with `ESMF_InfoGet` (defaults kept), so a subscriber can still set them.
   - No self-subscription problem: Turbulence's own internals are not in `TR` on `develop`.
2. **`TRI` export bundle** of the tendency exports, filled in Initialize, for History only. Kept (decision 2): a `MAPL_STATEITEM_FIELDBUNDLE` export; MAPL3 does not allocate unrequested exports, so the members must be allocated (or be internals) for the tendencies to be computed, as MAPL2's `alloc=.true.` did.
3. **`FRIENDLYTO`:** to self (`COMP_NAME`, `'TURBULENCE'`) becomes `ADD2EXPORT`; `'DYNAMICS:TURBULENCE'` (SHOC on) becomes `ADD2EXPORT` plus a `TRADV` subscription of `PDF_A`, `TKESHOC`, `QT2`, `QT3`.
4. **Two run phases and the refresh alarm.** Run1 and Run2 register with `phase_name`. MAPL2 runs `REFRESH` on the generic `RUNALARM` (period `DT:`) while both phases run every `RUN_DT:`. MAPL3 calls the user Run only when the component's own alarm rings, so the refresh needs its own alarm created in Initialize, with MAPL2's first-step behaviour (IRRAD, Solar). `RUN_DT` becomes the component's clock step (`MAPL_ClockGet(dt=)`). Done in steps 4 and 5.
5. **`ALLOC=.TRUE.` (45 sites).** As Moist decision 1: locals copied to exports when connected; internals with `ADD2EXPORT` for fields shared within or across steps or always connected (`KH`, `KM`, `RI`, `KPBL_SC`, `EDMF_MF`, `*AFDIFFUSE`, `SAFUPDATE`, `INTDIS`, `TOPDIS`, `SRFDIS`, `DPDTTRB`, the tendencies). Survey needed (as `mapl3-alloc-exports.md`). Done in step 6 (42 sites; decision 5 replaces the local-and-copy plan).
6. **Edge fields (48).** Check every fetch and its uses, including edge arrays passed to `LOUIS_KS`, `BELJAARS`, `ComputeMOScalingInterlevel`, `VTRI*`, `RUN_EDMF`, `RUN_SHOC`, `ENTRAIN`; remap to `0:LM` only where the MAPL2 code is 0-based (`DIFFUSE` uses `PLE(:,:,0:LM-1)`, `ZLE(:,:,0:LM-1)`; REFRESH `KH(I,J,L-1)`). Done in step 7: every edge fetch remapped.
7. **`MAPL_UNDEF`:** the 7 equality tests compare against fills Turbulence writes itself; check that none reads a MAPL2-converted input. Checked in step 7: no change.
8. **Resources:** `TURBULENCE_DO_SHOC` is read with a `COMP_NAME` prefix; MAPL3 has no prefix lookup. `TRI_increments::` is an `ESMF_Config` table; becomes a YAML list. Decision 6; `TRI_increments` done in step 4.
9. **Sources other than the driver:** `edmf.F90`, `shoc.F90`, `scm_surface.F90` use `MAPL_SatVaporMod` (MAPL3: `mapl_SatVapor_mod`), `edmf.F90` uses `MAPL_Mod` (gone).

## Steps

### 0. Sync with `develop` (done 2026-11-12; upstream 2026-11-13)
`GEOS_TurbulenceGridComp.F90` and `edmf.F90` were first copied onto this branch from `develop` @ `f31a2151`; the other Turbulence files were already identical. The same sync then reached `release/MAPL-v3` through PR #1582 (`86b1a7d5`, merge of `develop` @ `f31a2151`), so the branch was rebased onto it and its own sync commit dropped: all Turbulence sources on `release/MAPL-v3` are byte-identical to `develop` @ `f31a2151`. Before, the `release/MAPL-v3` copy was `develop` @ `7f80bd2f` plus one release-only edit (the two `FriendlyToTURBULENCE` reads through `ESMF_InfoGet`), which the port replaces anyway (issue 1). PR #1582 also synced `GEOS_PhysicsGridComp.F90` with `develop` (not built on this line).

### 1. MAPL2 baseline (to do)
- Build MAPL2 from `develop` @ `f31a2151` (or later; then repeat step 0 with that commit).
- Capture with `COMPONENT_TO_RECORD: TURBULENCE` at C12, before and after Run1 and Run2, import, internal and export. Record a later step than the first. Request every Turbulence export in `HISTORY.rc`.
- Configurations, matching the Moist cases: L181 SHOC-MF (v12 default, as `moist-gfdl1m-gf-shoc-sa`); later L72 Louis/Lock (as `moist-bacm1m-gf-uw-sa`).
- `DT:` equal to `RUN_DT:` (refresh every step) for the first comparison.
- Record every Turbulence setting used, and the `TR` contents (Physics prints them at Initialize).

### 2. Macro pass (no logic change; done 2026-11-12)
No codee formatting (decision 1). `develop` had already converted most calls to `_RC` (#1557); this pass did the rest of `GEOS_TurbulenceGridComp.F90` by script (the other sources have no MAPL macros):
- `MAPL_Generic.h` -> `MAPL.h`; `RETURN_(ESMF_SUCCESS)` -> `_RETURN(_SUCCESS)` (10); `__STAT__` -> `_STAT` (42); `__RC__` -> `_RC` (14); `allocate(SG(0,0), stat=STATUS)` + `VERIFY_` -> `_STAT`; the two `RETURN_(ESMF_FAILURE)` for an unknown `DiffuseLike` (DIFFUSE, UPDATE) -> `_FAIL` with a message.
- 31 calls ending `RC=STATUS)` with no check (`EDMF_NUMUP`, `SCM_SL`, `TURBULENCE_DO_SHOC`, 26 `EDMF_*` reads, `SHVC_SCALING`) -> `_RC`: an error is now caught; all have defaults, so results are unchanged.
- `Iam` declarations and assignments removed (9 + 6), including the shadowing `IAm='Refresh'`/`'Diffuse'`/`'Update'` locals. `COMP_NAME` kept (labels and `FRIENDLYTO`; step 4).
- Left: `ESMF_AlarmIsRinging(ALARM, rc=status)` + `_VERIFY` (a function call inside `if`; the alarm is reworked in step 5); `ESMF_ConfigGetDim(..., rc=STATUS)` in `GetTriIncrements` (status tested: a missing table means no names; step 4); the commented-out `VERIFY_` lines.
- Seen, not changed: in the `#ifdef _CUDA` branch of REFRESH (not compiled), one-line `IF (ASSOCIATED(X)) ALLOCATE(..., _STAT)`: the status check runs whether or not the `IF` holds (same with MAPL2's `__STAT__`).
- Checked: both versions normalized for the macro spellings and diffed; only the lines above differ.

### 3. State specs (done 2026-11-12)
`Turbulence_StateSpecs.rc` (48 imports, 149 exports, 45 internals = the 242 static MAPL2 calls), generated by script from the MAPL2 calls in SetServices, in their order within each state; included in SetServices (`Turbulence_Import___.h`, `_Export___.h`, `_Internal___.h`) in place of those calls; one `mapl_acg` call in `CMakeLists.txt` (as Moist, with `GET_POINTERS DECLARE_POINTERS`; the pointer headers are not included yet). Not built (step 8).
- Column mapping: `DEFAULT` -> `FILL` (17: `SHFX_SPRAY`, `LHFX_SPRAY`, `WTHV2`, `WQT_DC` imports, `SLQT` export, 12 internals; `TKESHOC` 1e-4); the other internals `FILL: 0.0`; `MAPL_RestartSkip` -> `SKIP`; `MAPL_RestartOptional` (`SHFX_SPRAY`, `LHFX_SPRAY` imports) -> blank; `FRIENDLYTO` (11) -> `ADD2EXPORT: T` (issue 3; the `TRADV` subscription is step 4); `UNGRIDDED_DIMS = (/NUMUP/)` -> `ungrd_numup`, a `MAPL_UngriddedDim(NUMUP, name="EDMF_plume", units="1")` set in SetServices before the includes (as GWD's `ungrd_16`).
- `COND`: `DO_WAVES/=0 .and. DO_SEA_SPRAY/=0` (sea-spray imports and exports), `SCM_SL /= 0` (`SHOBS`, `LHOBS` imports; `cu_scm`, `ct_scm`, `ssurf_scm`, `qsurf_scm` internals), replacing the `if` blocks; the variables are SetServices locals read before the includes.
- `ALIAS` where a short name is in two states (pointer declarations would clash): exports `TKEDISS` (`TKEDISSx`, internal `TKEDISS` too), `SHFX_SPRAY` (`SH_SPRX`), `LHFX_SPRAY` (`LH_SPRX`), the names Run1 uses. Not checked: whether MAPL3 accepts an export and a non-exported internal with the same short name (`TKEDISS`).
- Still hand-written MAPL2 calls (step 4): `TR` import and `TRI` export (bundles), and the `TRI_increments` loop (`<NAME>IT` exports, names known only at run time).
- SetServices: the `TURBULENCE_DO_SHOC` read moved before the includes (its value only built `FRIENDLIES_SHOC`, now unused and removed; kept for the `TRADV` subscription). The 8 commented-out export calls (`SL2DIAG`, `QT2DIAG`, `SLQTDIAG`, `W3CANUTO`, `T`, `U`, `V`, `QV`) were dropped.
- Checked: ACG (`python3.12 MAPL_GridCompSpecs_ACG.py`) generates the headers; a script compared each of the 242 generated `MAPL_GridCompAddSpec` calls with its MAPL2 call (name, units, dims, stagger, restart, fill, `add_to_export`, standard name, ungridded, condition): all match after the mapping above.

### 4. SetServices (done 2026-11-13)
- Signatures as `I_SetServices`/`I_Run` (no `intent` on the ESMF arguments, `rc` `intent(out)`), for SetServices, Initialize, Run1, Run2.
- Run phases `phase_name='Run1'`/`'Run2'` (decision 4); Initialize keeps the default name.
- Resources: `MAPL_GridCompGetResource(GC, 'LABEL', ...)` for `DEBUG_TRB`, `USE_WAVES`, `USE_SEA_SPRAY`, `EDMF_NUMUP`, `SCM_SL` and `TURBULENCE_DO_SHOC` (`ESMF_ConfigGetAttribute` for the last three in MAPL2). Labels built from `COMP_NAME` are kept (issue 8).
- `TR` import -> `TRDIFF` export, `itemtype=MAPL_STATEITEM_SERVICE`, standard name `turbulently_diffused_quantities` (the name Moist subscribes with).
- `TRI`: export with `itemtype=MAPL_STATEITEM_FIELDBUNDLE` (`dims='xyz'`, centre, units `X kg m-2 s-1` as MAPL2), as SU2G's `SU_AERO_DP`.
- Tendencies (decision 2): the 7 built-in tendency exports (`SIT`, `UIT`, `VIT`, `TKESHOCIT`, `QVIT`, `QLLSIT`, `QILSIT`) moved to the end of the INTERNAL table (`SKIP`, `FILL 0.0`, `ADD2EXPORT`); the `TRI_increments` `<NAME>IT` tendencies are hand-written internals with the same attributes (MAPL2: exports, all allocated for `TRI`). MAPL2's units (`'UNITS'`, `'days s-1'`) are kept; not checked whether MAPL3 rejects `'UNITS'` at a connection.
- `TRI_increments`: a YAML string list read with `ESMF_HConfigIsDefined` + `ESMF_HConfigAsStringSeq` on the component's hconfig (as DynCore's `EXCLUDE_ADVECTION_TRACERS_LIST`), in `GetTriIncrements(GC, ...)`; built-in names and duplicates skipped as before.
- `TRADV` import (service subscription, `service_items` `PDF_A`, `TKESHOC`, `QT2`, `QT3`) when `DO_SHOC /= 0` (issue 3), as Moist's.
- Dropped: `MAPL_GetObjectFromGC`, the 18 `MAPL_TimerAdd`, `MAPL_GenericSetServices`.

### 5. Initialize (done 2026-11-13)
- No `MAPL_GenericInitialize`.
- Refresh alarm (issue 4): `turbulence_refresh_alarm` on Turbulence's clock, period `TURBULENCE_DT` (integer seconds, default the clock step; asserted positive and a multiple of the step), created with the MAPL2 `handle_clock_and_main_alarm` logic (`RUN_AT_INTERVAL_START`, `REFERENCE_DATE`, `REFERENCE_TIME`, ring time backed off by one step, `AlarmRingerOn` when it falls on the current time), the same code as IRRAD's and SOLAR's MAPL3 ports. Run1 gets it with `ESMF_ClockGetAlarm`, refreshes when it rings, and turns it off. MAPL2's label was `DT:` with the component prefix, i.e. `TURBULENCE_DT:` or `DT:`; only `TURBULENCE_DT` is read (as SOLAR reads `SOLAR_DT`). Not run-tested.
- MAPL3 advances a component's clock once per step, in its own phase (`OuterMetaComponent/run_clock_advance.F90`), so Run1 and Run2 see the same time and Run2 needs no alarm.
- `TRI`: the 7 tendency internals (`QVIT` as `QIT`, through `ESMF_NamedAlias`) and the `TRI_increments` internals are added with `MAPL_FieldBundleAdd`, taken from INTERNAL. Not run-tested: that MAPL3 accepts adding fields to an export field bundle in Initialize (SU2G does), and that History can request `TRI%NAME`.

### 6. `ALLOC` replacement (done 2026-11-13; issue 5)
Survey (2026-11-13): 42 sites (Initialize 1, REFRESH 29, DIFFUSE 4, UPDATE 8; the notes' 45 counted comments). Group B (read by DIFFUSE on steps without a refresh): `EDMF_MF`, `EDMF_ENTR`. No group A. All others are local to one routine (C1/C3/C4), including fields other components read (`KH`, `KM`, `RI`, `KPBL_SC`, `SL2`, ...): a consumer reads the export, which keeps its last refreshed value between refreshes, as in MAPL2.
- B: `EDMF_MF` (edge) and `EDMF_ENTR` moved to the INTERNAL table (`SKIP`, `FILL 0.0`, `ADD2EXPORT`); REFRESH and DIFFUSE fetch them from INTERNAL (`EDMF_MF` remapped to `0:LM`). In MAPL2, DIFFUSE read them unguarded under `EDMF_DOTRACERS` and only worked if `PDFSHAPE >= 5` or History requested `EDMF_MF`.
- C (39 sites, decision 5): `get_pointer_alloc(STATE, X, NAME, ALLOC, WORK, IM, JM[, LB, LM])` reproduces `MAPL_GetPointer(..., ALLOC=)`: `X` is the export when it is requested, else (if `ALLOC`) a zero-filled `WORK` array with MAPL2's bounds (`0:LM` for edge fields), else null. `WORK` is an `allocatable, target` local of the calling routine (`<var>_wk`), freed on return. So every `associated()` test, pointer dummy (`RI` -> `run_shoc`'s `ri_inv`, `BRUNTSHOC`, `EDMF_ENTR`), conditional `ALLOC` (`PDFALLOC`, `DO_SHOC/=0`, `associated(KETRB) .or. ...`) and side effect (`KPBL_SC` -> `ZPBL_SC`; `WQT`/`WSL` -> the in-place `QTFLXMF`/`SLFLXMF` update; `QT`/`SL` gates) behaves as in MAPL2, with no copies.
  - MAPL2's ALLOC'd buffers persisted between calls; the work arrays start from zero each call. Checked: no C site reads a previous value except `ISOTROPY` level 0 with SHOC off (never written; `max(10., ...)` gives 10 from MAPL2's persisted 10 or from 0 alike) and fields never written in a configuration (`ZCLD`, `ZSML`, `ZRADML`, `ZRADBS` without Lock, `EDMF_DQRDT`/`EDMF_DQSDT` without EDMF): zero either way.
  - No export `FILL` (AGENTS.md rule): a requested export that is never written in a configuration (the four above, the PDF moments with `PDFSHAPE < 5`, `ISOTROPY` level 0) is uninitialized in MAPL3, zero in MAPL2. A regression case must not request them in such a configuration. With SHOC off and `TKE` requested, a requested `ISOTROPY` gets level 0 from `max(10., min(2000., <uninitialized>))`, and `update_moments` reads it: request `ISOTROPY` only with SHOC on.
- Physics reads `UAFDIFFUSE`, `VAFDIFFUSE`, `SAFDIFFUSE`, `SAFUPDATE` unguarded (MAPL2 `alloc`); the Physics port must connect them (step 10).

### 7. Run1/REFRESH/DIFFUSE, Run2/UPDATE (done 2026-11-13)
- Run1/Run2: `MAPL_GridCompGet(num_levels=, grid=, logger=)`, `MAPL_GridGet(IM=, JM=)`, `MAPL_GridCompGetInternalState`; timestep `MAPL_ClockGet(dt=)` (MAPL2: `RUN_DT:`); Run2's `LATS` from `MAPL_GridGetCoordinates` (radians, as MAPL2's `MAPL_Get(LATS=)`; allocatable).
- 157 `MAPL_GetResource` -> `MAPL_GridCompGetResource(GC, 'LABEL', VAR, default=, _RC)` by script, checked against the MAPL2 calls (variable, label, default, in order). Self-defaults (`HGT_SURFACE`, `CAP_INTDIS`, `CAP_TOPDIS`) pass `default=(X)`. The 13 `SCM_*` reads without `RC` now have `_RC` (all have defaults).
- 40 `MAPL_TimerOn/Off` -> `MAPL_GridCompTimerStart/Stop(GC, ..., _RC)`, names unchanged.
- 203 `MAPL_GetPointer` without `ALLOC` -> `MAPL_StateGetPointer`; checked one-to-one with the MAPL2 fetches (state, pointer, name, order).
- Edge fields (issue 6): audit of the 48 edge fields at 54 fetch sites. MAPL2 gave every edge pointer `0:LM` (4-D: `(:,:,0:LM,1:NUMUP)`; `PREF(0:LM)`), and every bound-dependent use assumes it; the others (whole-array, explicit-shape dummies, `BELJAARS`'s `(:,:,0:)`) do not depend on bounds. So every edge fetch is remapped right after the fetch, with `zero_based_edges` (1-D, 3-D, 4-D; `X(1:,1:,0:) => X` if associated; inside `get_pointer_alloc` for the ALLOC'd ones). Sites where it matters: `PLE`, `ZLE` (Run1 and UPDATE), `TKH` (sections `TKH(:,:,1:LM)`, `KM(:,:,1:LM)`, `ISOTROPY(:,:,1:LM)` passed to `run_shoc`'s explicit-shape dummies would otherwise be one level off without an error), `KH`, `KM`, `RI`, `DU`, `KHSFC`, `TKE`, `ALH`, `ALM`, the pointer dummies of `LOUIS_KS*` (`DU_DIAG`, `ALH_DIAG`, `ALM_DIAG`), `run_shoc` (`ri_inv`, `prnum_inv`), `RUN_EDMF` (`EDMF_PLUMES_*`, `EDMF_DRY_*`/`EDMF_MOIST_*`), `EDMF_MF` in DIFFUSE, and UPDATE's `MFAW`, `QTFLXMF`, `SLFLXMF`, `UFLXTRB`, `VFLXTRB`.
- `TRDIFF` (issue 1): `GetTRDIFF` takes the bundle from EXPORT and asserts `MAPL_FieldBundleSameData` with the IMPORT copy (as Moist's `MTR`); `GetDiffusedField` replaces the per-field block of DIFFUSE and UPDATE: `ESMF_FieldBundleGet(TR, K, FIELD)` (unchanged index order), short name from `diffused_short_name` (standard name -> `Q`, `QLLS`, `QLCN`, `QILS`, `QICN`, `QRAIN`, `QSNOW`, `QGRAUPEL`, else the standard name), `DiffuseLike`/`WeightedTendency` via `ESMF_InfoGetFromHost` + `ESMF_InfoGet(default=)` (same defaults), data via `ESMF_FieldGet(localDe=0, farrayPtr=)`. `FRIENDLY = .true.` for every field (MAPL2 read `FriendlyToTURBULENCE`; Physics put only friendlies in `TR`). REFRESH's unused `TR` fetch (`NTR`) now reads `TRDIFF`. Field order: as MAPL2, `ESMF_FieldBundleGet` by index; each field is diffused independently, so order does not change results.
- `GetTendency` looks the tendency up in INTERNAL (all tendencies are internals).
- `MAPL_UNDEF` (issue 7): the 7 equality tests (`ZPBL2`, `ZPBL10p`, `ZPBLRI`, `KPBL_SC`) compare against fills REFRESH writes just before; no change. Not changed, to check in step 9: `RUN_EDMF` fills `EDMF_ENTR` with `MAPL_UNDEF` (`edmf.F90`) and DIFFUSE computes `exp(-(entdyn+EDMF_ENTX(I,J,L+1))*DZ)` (EDMF tracer transport); if an undefined level is ever read, MAPL2 gives `exp(-1e15*dz) = 0`, MAPL3 (`-huge`) overflows.
- `MAPL_MaxMin` (9 calls, `DEBUG_TRB`, in `LOUIS_KS*` and `BELJAARS`, which have no `GC`) -> `report_maxmin(label, x)`: `MAPL_MaxMin(x, comm)` on the current VM, logged through the component logger kept in a module variable (`TRB_LOGGER`, set by Run1 and Run2).
- Seen, not changed (MAPL2 behaviour): exports used unguarded without `ALLOC`, which crash if not requested: `KPBL` (REFRESH, feeds `BELJAARS` when `C_B /= 0`, the default), `SBIFRQ`/`SBITOP` (written when either is requested), `TRINVDELT`/`TRINVFRQ` (guarded by `TRINVBS` only), `SFLX`/`QFLX` with sea spray (DIFFUSE), `QFLX` in the `DPDTTRB` block (UPDATE). In the GCM, Moist and Surface connect `KPBL` and the fluxes; a standalone case must request them. `ZPBLRI` is allocated locally without a matching deallocate (leaks each refresh). `UPDATE` changes the internals `QTFLXMF`/`SLFLXMF` in place each step, so between refreshes the `-MFAW*QT` term is applied repeatedly.

### 8. Build
Uncomment `GEOSturbulence_GridComp` in `../CMakeLists.txt`; `TYPE SHARED`; `Turbulence_SetServices` wrapper; fixes in the other sources (issue 9).

### 9. Standalone regression (`../regression/`)
`DataForTurbulence` (all imports, restarted from the step-1 "before" capture), a `DataTracers` subscriber to `TRDIFF`, History collection that activates only written exports, `logging.yaml`, `cmpchk.py`. Phases: decision 3.

### 10. Wiring (blocked on the Physics/Agcm ports)
Connections with Moist and Surface (`SHAT`... `SDEL`..., `SFLX`... `SDFL`...), `TRDIFF` subscriptions from GOCART, Physics' use of `INTDIS`/`TOPDIS`, `*AFDIFFUSE`, `SAFUPDATE`.

## Decisions
1. **DECIDED (2026-11-12) - no codee format pass.** Unlike Moist: Turbulence is still changed on `develop` (EDMF, SHOC), and keeping its layout keeps later syncs mergeable. Macros only (step 2); new code follows the surrounding style.
2. **DECIDED (2026-11-12) - keep `TRI`.** An export bundle filled in Initialize with the `*IT` tendency exports, as MAPL2; closest to MAPL2, and History configurations using `TRI%NAME` keep working.
3. **Open - regression phases:** one case running Run1 then Run2 (inputs `SDEL`... supplied as data), or one case per phase as GOCART2G. Decide once the step-1 capture exists.
4. **DECIDED (2026-11-13) - run phase names `Run1`/`Run2`**, as GOCART2G (`GOCART2G_GridCompMod.F90`), whose standalone cases select a phase with `run_phases: [Run1]` (configurable component). The earlier plan said `run1`/`run2`.
5. **DECIDED (2026-11-13) - `ALLOC=` by `get_pointer_alloc`** for the local (group C) sites instead of Moist's "always a local, copy to the export" (Moist decision 1): the same results, but the export is written directly when it is requested, as in MAPL2, so `associated()` tests, pointer dummies and conditional `ALLOC` need no special cases. Group B (`EDMF_MF`, `EDMF_ENTR`) becomes internals, as in Moist.
6. **DECIDED (2026-11-13) - resource labels with `COMP_NAME` kept** (`TURBULENCE_DO_SHOC`, `TURBULENCE_LOUIS_B_KH`, ...; issue 8). MAPL3 names the user gridcomp after the child (`@MAPL/superstructure/generic/OuterMetaComponent/add_child_by_spec.F90`), so with the child named `TURBULENCE` the keys are MAPL2's. Lost: MAPL2's automatic `<COMP>_` prefix lookup for unprefixed labels (e.g. `TURBULENCE_PDFSHAPE:` overriding `PDFSHAPE:`).
