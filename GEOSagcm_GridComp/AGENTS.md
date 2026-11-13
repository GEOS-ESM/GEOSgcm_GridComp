# GEOSagcm_GridComp: guide for agents

**Describes:** the MAPL3 line, branch `release/MAPL-v3`; verified against `GEOSgcm_GridComp` @ `86b1a7d5` plus the Turbulence port branch (2026-11-13).
**Scope:** the atmosphere subtree (`GEOSagcm_GridComp`): the Agcm composite, Superdyn, and Physics with Gravity Wave Drag, Moist and Turbulence. Surface is not covered; FVdycoreCubed (`GEOSsuperdyn_GridComp/@FVdycoreCubed_GridComp`), Chemistry (`GEOSphysics_GridComp/@GEOSchem_GridComp`) and Radiation (`GEOSphysics_GridComp/@GEOSradiation_GridComp`) are separate repositories with their own `AGENTS.md`. The rules of the `GEOSgcm_GridComp` `AGENTS.md` (`../AGENTS.md`), including "Porting to MAPL3", apply here too. Paths below are relative to this directory.

> **Ported on this line: Superdyn, GWD, Moist and Turbulence.** `CMakeLists.txt` builds Superdyn (with FVdycoreCubed) and, under `GEOSphysics_GridComp`, Chemistry (GOCART only), Moist, Turbulence, GWD and Radiation, plus the regression cases. The Agcm composite (`GEOS_AgcmGridComp.F90`), the Physics composite (`GEOS_PhysicsGridComp.F90`, `MBundle_IncrementMod.F90`), Surface and `GEOSdatmodyn_GridComp` are MAPL2 code and are not compiled. The ported components therefore run only in standalone regression cases (Superdyn's dynamical core in FVdycoreCubed's). Do not port or fix the MAPL2 code unless asked.

## Rules for agents

- Read [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) (what the Agcm subtree is and how it is put together: run orders, Superdyn, the Physics children, services, contracts, the background to the traps below) before changing children, connections, exports, schemes or configuration. Update it when the structure, a contract or an invariant changes; it describes only what changes rarely, so keep status and to-dos out of it.
- Moist's port is recorded in `GEOSphysics_GridComp/GEOSmoist_GridComp/mapl3-porting-notes.md`: read it before changing Moist, and record each decision there with its reason.
- Turbulence's port is recorded in `GEOSphysics_GridComp/GEOSturbulence_GridComp/mapl3-porting-notes.md`. Same rule.
- Keep this file current: update it when status, a rule or a trap changes, and update its "verified against" commit.

## Status and to-dos

- Superdyn: only `DYCORE: FV3` is accepted; the `FV3+ADV` and data-dynamics (`DATMO`) branches are commented out, so `ADV` is never added.
- Turbulence: only the SHOC-MF configuration (`TURBULENCE_DO_SHOC: 1`, `EDMF_DOMF: 1`) is regression-tested; Louis/Lock, sea spray and the SCM surface are ported but not run. Its `TRDIFF` service maps a subscribed field's standard name back to a short name only for Moist's water species (`diffused_short_name`); other tracers keep their standard name, so a `TRI_increments` tendency for them is not found (porting notes, issue 1).
- Turbulence's regression data is only in `~/regression-data-turb` of the porter (bmauer) until copied to `/ford1` and S3; without it the cases pass without running.
- Nothing connects Superdyn and Physics yet: the Agcm is not ported, so the `TRADV` subscriptions of Moist and GOCART and the application of physics tendencies are not exercised.
- Moist schemes: GFDL_1M, BACM_1M, GF and UW are ported. MGB2_2M, THOM_1M and RAS stop with `_FAIL` in `SetServices` and `Initialize`; their interface modules are commented out of `CMakeLists.txt` (their kernels are built). NSSL_2M is not selectable.
- `CONVECTION_TRACER /= 0` is not supported (`GF_Setup` stops); the scheme owner confirmed the option is unused.
- Moist's service subscriptions (`TRADV`, `TRDIFF`, `TRCHEM`, `TRANA`) and its `MTR` provider are unconnected; the regression cases run with them unconnected.
- `FILLQ2ZERO` declares `DQDT` as an `intent(out)` pointer, so the `*_FILL_*` tendency exports are never written (same in MAPL2): GEOSgcm_GridComp#1570, `TODO` in `Process_Library.F90`.
- `UW_Run` has an OpenMP race on `ZLE0` (same in MAPL2): GEOSgcm_GridComp#1572.
- `REFL10CM_MAX` differs from the MAPL2 capture and is left out of the regression History collection (Moist porting notes, step 10).
- GWD: the `MAPL_MaxMin` debug calls are not ported (`TODO` in `GEOS_GwdGridComp.F90`).

## Known traps

Agcm and Superdyn:
- **`TRADV` changes direction.** In the MAPL2 Agcm, Physics exports a `TRADV` bundle that the Agcm connects to Superdyn. In MAPL3, Superdyn re-exports FV3's `TRADV` service and the Physics children subscribe to it: porting the Agcm means connecting Superdyn's `TRADV` to Physics and dropping Physics' MAPL2 `TRADV` export.
- **Superdyn's comments describe MAPL2.** The `Initialize` and `run` comments mention running GWD and the Topo utility; the code only runs `DYN`.
- **A `DYN` export is visible above Superdyn only if Superdyn re-exports it** (`MAPL_GridCompReexport` in `SetServices`).

Turbulence:
- **Run2 continues from Run1.** It accumulates onto Run1's `RESTART: SKIP` internals and exports (`SIT`, `SFLX`, ...), so Run2 cannot be run or tested alone; the `run2` regression case runs both phases through a test driver (`GEOSphysics_GridComp/GEOSturbulence_GridComp/regression/README.md`).
- **`ALLOC=` is emulated, not removed.** `get_pointer_alloc` points at the export when it is requested, else at a local zero-filled array; which exports are requested still changes which branches run (`associated()` tests), as in MAPL2.
- **Edge fields are remapped to `0:LM`** after every fetch (`zero_based_edges`); a new edge fetch needs the same, including pointers passed to `LOUIS_KS*`, `run_shoc` and `RUN_EDMF` pointer dummies.
- **Resource labels keep the `TURBULENCE_` prefix** built from the component name, so the child must be named `TURBULENCE` for MAPL2's keys to apply.

GWD:
- **Two schemes are summed.** If both `GEOS_EFFGWBKG` and `NCAR_EFFGWBKG` (or both orographic ones) are non-zero, both act; check all four keys when tuning.
- **`NCAR_NRDG` changes the restart.** Turning the ridge scheme on without ridge data in the internal restart leaves those fields at their defaults.
- **`Initialize` parses `AGCM.GRIDNAME`** as `PE<n>x...` (digits from the third character up to `x`); other forms make the `read` fail.
- **`Run` has no alarm:** it works on every call, with `dt` from the clock.
- **`DTDT` is mass-weighted (`Pa K s-1`); `TTMGW` is not.**
- **A new ungridded dimension** must be created in `SetServices` (as `ungrd_16`) for `GWD_StateSpecs.rc` to name it.
- **With `use_threads`, module-level scratch data is a race;** use the per-thread workspaces.

Moist:
- **A new scheme touches four places:** the `if` chains in `SetServices`, `Initialize`, `Run` and `Finalize`, plus `CMakeLists.txt`.
- **`CLDMICR_OPTION` defaults to `BACM_1M`:** a configuration without the key runs BACM_1M, not GFDL_1M.
- **Which exports are requested changes results.** Moist uses some exports as work arrays and sums others by name; MAPL3 does not allocate an export nobody requests (`GEOSphysics_GridComp/GEOSmoist_GridComp/mapl3-alloc-exports.md`). In a regression case the History collection decides; regenerate its list when the case's configuration changes.
- **Do not run Moist with more than one OpenMP thread** until #1572 is fixed: `UW_Run` races, and a Debug build traps.
- **`MTR` fields are looked up by name** in `CNV_Tracers_Init`, but service-bundle fields carry MAPL-generated names, and two subscribers may subscribe the same name; GF's own `CNV_TR` is appended to a local copy of the bundle because Moist cannot subscribe to its own service (Moist porting notes, "Self-subscription").

## Build and test

- Regression tests (label `REGRESSION`): `GEOSphysics_GridComp/regression/` (`moist-gfdl1m-gf-uw-sa`, `moist-gfdl1m-gf-shoc-sa`, `moist-bacm1m-gf-uw-sa`; its `README.md` describes them, how they were built from MAPL2 captures, `cmpchk.py`, and their run time), `GEOSphysics_GridComp/GEOSgwd_GridComp/regression/` (`gwd-ncar-sa`) and `GEOSphysics_GridComp/GEOSturbulence_GridComp/regression/` (`turbulence-shocmf-sa/run1`, `turbulence-shocmf-sa/run2`). The dynamical-core cases are in FVdycoreCubed. The general rules for regression cases are in the `GEOSgcm_GridComp` `AGENTS.md`.
- ctest sets `OMP_NUM_THREADS=1` for the Physics cases (`EXTRA_ENV` in their regression `CMakeLists.txt`); set it yourself when running a case by hand.
- `GEOSphysics_GridComp/GEOSmoist_GridComp/pyMoist/` has its own tests (`pyMoist/tests/`).
