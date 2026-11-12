# GEOS_TurbulenceGridComp MAPL3 porting plan

Status (2026-11-12): step 0 done on branch `feature/bmauer/port-turbulence-to-mapl3` (off `release/MAPL-v3`). Step 1 (MAPL2 baseline) to do; steps 2-10 not started.

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

1. **`TR` becomes a provided service.** MAPL2's `TR` import is filled by Physics with the friendlies to `TURBULENCE`. In MAPL3 Turbulence provides a service `TRDIFF` (standard name `turbulently_diffused_quantities`, the name Moist already subscribes with; Moist notes step 4), and Moist, GOCART, ... subscribe their internals.
   - Every field in a service bundle is a subscriber's internal, so all are updated in place: `FRIENDLY` becomes true for every `TR` field and the `FriendlyToTURBULENCE` read goes.
   - Subscribed fields carry MAPL-generated names (Moist notes, "Self-subscription"). `DIFFUSE`/`UPDATE` select by name (`Q`, `QLLS`, `QILS`, the `S`/`U`/`V`/`Q` surface terms) and `GetTendency` looks up `<NAME>IT` by name, so names must be recovered through `standard_name`, as `get_short_name` does. How tracers without a known standard name (GOCART, `TRI_increments`) are named is open.
   - `DiffuseLike`/`WeightedTendency` (`ESMF_AttributeGet`): nobody sets them for `TR` fields on `develop` (Physics only set them on `S`, `U`, `V`, now diffused directly), so the defaults always apply. Read them with `ESMF_InfoGet` (defaults kept), so a subscriber can still set them.
   - No self-subscription problem: Turbulence's own internals are not in `TR` on `develop`.
2. **`TRI` export bundle** of the tendency exports, filled in Initialize, for History only. MAPL3 History can request the `*IT` exports directly. Open decision 2.
3. **`FRIENDLYTO`:** to self (`COMP_NAME`, `'TURBULENCE'`) becomes `ADD2EXPORT`; `'DYNAMICS:TURBULENCE'` (SHOC on) becomes `ADD2EXPORT` plus a `TRADV` subscription of `PDF_A`, `TKESHOC`, `QT2`, `QT3`.
4. **Two run phases and the refresh alarm.** Run1 and Run2 register with `phase_name`. MAPL2 runs `REFRESH` on the generic `RUNALARM` (period `DT:`) while both phases run every `RUN_DT:`. MAPL3 calls the user Run only when the component's own alarm rings, so the refresh needs its own alarm created in Initialize, with MAPL2's first-step behaviour (IRRAD, Solar). `RUN_DT` becomes the component's clock step (`MAPL_ClockGet(dt=)`).
5. **`ALLOC=.TRUE.` (45 sites).** As Moist decision 1: locals copied to exports when connected; internals with `ADD2EXPORT` for fields shared within or across steps or always connected (`KH`, `KM`, `RI`, `KPBL_SC`, `EDMF_MF`, `*AFDIFFUSE`, `SAFUPDATE`, `INTDIS`, `TOPDIS`, `SRFDIS`, `DPDTTRB`, the tendencies). Survey needed (as `mapl3-alloc-exports.md`).
6. **Edge fields (48).** Check every fetch and its uses, including edge arrays passed to `LOUIS_KS`, `BELJAARS`, `ComputeMOScalingInterlevel`, `VTRI*`, `RUN_EDMF`, `RUN_SHOC`, `ENTRAIN`; remap to `0:LM` only where the MAPL2 code is 0-based (`DIFFUSE` uses `PLE(:,:,0:LM-1)`, `ZLE(:,:,0:LM-1)`; REFRESH `KH(I,J,L-1)`).
7. **`MAPL_UNDEF`:** the 7 equality tests compare against fills Turbulence writes itself; check that none reads a MAPL2-converted input.
8. **Resources:** `TURBULENCE_DO_SHOC` is read with a `COMP_NAME` prefix; MAPL3 has no prefix lookup. `TRI_increments::` is an `ESMF_Config` table; becomes a YAML list.
9. **Sources other than the driver:** `edmf.F90`, `shoc.F90`, `scm_surface.F90` use `MAPL_SatVaporMod` (MAPL3: `mapl_SatVapor_mod`), `edmf.F90` uses `MAPL_Mod` (gone).

## Steps

### 0. Sync with `develop` (done 2026-11-12)
`GEOS_TurbulenceGridComp.F90` and `edmf.F90` copied from `develop` @ `f31a2151` (commit `9270a7e9`); the other Turbulence files were already identical. Before, the `release/MAPL-v3` copy was `develop` @ `7f80bd2f` plus one release-only edit (the two `FriendlyToTURBULENCE` reads through `ESMF_InfoGet`), which the port replaces anyway (issue 1). `GEOS_PhysicsGridComp.F90` was not synced (decided 2026-11-12): it is not built on this line.

### 1. MAPL2 baseline (to do)
- Build MAPL2 from `develop` @ `f31a2151` (or later; then repeat step 0 with that commit).
- Capture with `COMPONENT_TO_RECORD: TURBULENCE` at C12, before and after Run1 and Run2, import, internal and export. Record a later step than the first. Request every Turbulence export in `HISTORY.rc`.
- Configurations, matching the Moist cases: L181 SHOC-MF (v12 default, as `moist-gfdl1m-gf-shoc-sa`); later L72 Louis/Lock (as `moist-bacm1m-gf-uw-sa`).
- `DT:` equal to `RUN_DT:` (refresh every step) for the first comparison.
- Record every Turbulence setting used, and the `TR` contents (Physics prints them at Initialize).

### 2. Formatting and macro pass (no logic change)
Open decision 1. Macros: `MAPL_Generic.h` -> `MAPL.h`, `VERIFY_`/`RETURN_`/`__STAT__` -> `_VERIFY`/`_RETURN(_SUCCESS)`/`_STAT`, drop `Iam`/`COMP_NAME` where unused. Check by token comparison that only macro spellings changed.

### 3. State specs
`Turbulence_StateSpecs.rc` with `mapl_acg`; column mapping as Moist step 3 (`DEFAULT` -> `FILL`, internals `FILL: 0.0` where MAPL2 zero-filled, `MAPL_RestartSkip` -> `SKIP`, `MAPL_RestartOptional` -> blank, `FRIENDLYTO` per issue 3, `COND` for sea spray and `SCM_SL`). Hand-written `MAPL_GridCompAddSpec`: `TR`/`TRI`, the `NUMUP` plume exports if ACG cannot express them, the `TRI_increments` tendencies.

### 4. SetServices
Entry points (`I_SetServices`/`I_Run` signatures; `phase_name="run1"`/`"run2"`), resources, services (issue 1, 3), no timers, no `MAPL_GenericSetServices`.

### 5. Initialize
Refresh alarm (issue 4); `TRI` per decision 2; no `MAPL_GenericInitialize`.

### 6. `ALLOC` survey and replacement (issue 5)

### 7. Run1/REFRESH/DIFFUSE, Run2/UPDATE
API replacements; `TR` service handling (issue 1); edge remaps (issue 6); `MAPL_UNDEF` (issue 7); `MAPL_MaxMin`; `ESMFL_BundleGetPointerToData` -> `ESMF_FieldBundleGet` + `ESMF_FieldGet`.

### 8. Build
Uncomment `GEOSturbulence_GridComp` in `../CMakeLists.txt`; `TYPE SHARED`; `Turbulence_SetServices` wrapper; fixes in the other sources (issue 9).

### 9. Standalone regression (`../regression/`)
`DataForTurbulence` (all imports, restarted from the step-1 "before" capture), a `DataTracers` subscriber to `TRDIFF`, History collection that activates only written exports, `logging.yaml`, `cmpchk.py`. Phases: decision 3.

### 10. Wiring (blocked on the Physics/Agcm ports)
Connections with Moist and Surface (`SHAT`... `SDEL`..., `SFLX`... `SDFL`...), `TRDIFF` subscriptions from GOCART, Physics' use of `INTDIS`/`TOPDIS`, `*AFDIFFUSE`, `SAFUPDATE`.

## Open decisions
1. **codee format pass:** as Moist (formatting commit, revertible and blame-ignorable) or not (keeps later `develop` syncs mergeable; Turbulence is still changed on `develop`).
2. **`TRI` bundle:** keep (export bundle filled in Initialize, closest to MAPL2) or drop (History requests `*IT` directly).
3. **Regression phases:** one case running Run1 then Run2 (inputs `SDEL`... supplied as data), or one case per phase as GOCART2G.
