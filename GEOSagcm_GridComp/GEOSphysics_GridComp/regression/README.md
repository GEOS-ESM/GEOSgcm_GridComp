# GEOSphysics regression tests

Standalone MAPL3 regression cases for the Physics components that have no regression directory of their own. Each case runs one component under a configurable `Root`, with every import supplied from a restarted "DataFor*" gridcomp, and compares the resulting checkpoints against a stored baseline. GWD (`../GEOSgwd_GridComp/regression/`) and Radiation (`../@GEOSradiation_GridComp/regression/`) keep their own.

| Case       | Component run | States compared by ctest       | Status |
|------------|---------------|--------------------------------|--------|
| `moist-gfdl1m-gf-uw-sa` | MOIST         | MOIST import, internal, export | matches the MAPL2 capture (cmpchk.py; REFL10CM_MAX left out); regression data in ~/regression-data only (2026-11-04), not on /ford1 |
| `moist-bacm1m-gf-uw-sa` | MOIST         | MOIST import, internal, export | matches the C12L72 MAPL2 capture (cmpchk.py, 2026-11-05); regression data in ~/regression-data only (2026-11-05), not on /ford1; reproducibility not checked with a second run |
| `moist-gfdl1m-gf-shoc-sa` | MOIST       | MOIST import, internal, export | matches the MAPL2 capture (cmpchk.py, 2026-11-05); regression data in ~/regression-data only (2026-11-05), not on /ford1; reproducibility not checked with a second run |

The cases are listed in `cases.txt` and registered by `CMakeLists.txt` through `esma_add_regression_tests` (`NAME moist`; ctest label `REGRESSION`, test name = case name). `CMakeLists.txt` is added by `add_subdirectory (regression)` in `../CMakeLists.txt`.

## Common setup

- **Geometry**: C12 cubed sphere, 181 levels (C12L181), model vertical grid defined by the `PLE` field (`root.yaml`); `moist-bacm1m-gf-uw-sa` uses 72 levels (C12L72).
- **Layout**: 6 MPI tasks (`mapl.yaml`, `run_case.cmake`), `OMP_NUM_THREADS=1`.
- **Logging**: `mapl.yaml` loads `logging.yaml` (`pflogger:`), which sets the `MAPL` logger to ERROR and logs to the console only. The run has no field dictionary, so MAPL warns once per variable spec per rank (about 5,700 warnings), and its default logging writes each one to `warnings_and_errors.log` under an MPI lock; that serialized the ranks and took about 6 s of each 17 s test (perf, 2026-11-09). With `logging.yaml` the tests take 8.5-13 s.
- **Clock**: one 20-minute step (`PT20M`) from `2000-04-14T21:20:00` (`cap_restart.yaml`); the final checkpoint is written at 21:40 (`checkpoints/last`). `moist-bacm1m-gf-uw-sa` instead takes one 450 s step (`PT7M30S`) from `2000-04-14T21:15:00` (see its section). The `start`/`stop` in `cap_driver.yaml` only bound the clock.
- **No ExtData**: `run_extdata: false`. **History only activates exports**: `run_history: true` with a `HIST` collection (`history.yaml`) whose start is after the cap stop, so it never writes a file.
- **Inputs**: the case's `DataFor*` gridcomp (`libconfigurable_gridcomp`) holds every import as an internal field, restarted from `checkpoints/2000-04-14T21:20:00/DataFor*_internal.nc`, and is connected field-by-field to the component in `root.yaml`.
- **Exports**: only those the case writes are activated, through `history.yaml` (not `activate_all_exports`: MAPL3 leaves an allocated export without a fill value uninitialized, so an export the configuration never writes would hold garbage). Regenerate the list when the configuration changes (`../GEOSmoist_GridComp/mapl3-porting-notes.md`, "Unwritten exports").

## Cases

### `moist-gfdl1m-gf-uw-sa`

MOIST alone, with GFDL_1M microphysics, GF deep convection and UW shallow convection (the L72 default; for 91 levels and up, gcm_setup uses SHOC-MF instead, see `moist-gfdl1m-gf-shoc-sa`). The other schemes `_FAIL` in the MAPL3 port.

Gridcomp hierarchy and connections (`cap_gridcomp.yaml`, `root.yaml`):

```
CAP (libMAPL.cap; ExtData and History off)
└── Root (libconfigurable_gridcomp; C12L181 geometry)
    ├── DataForMoist (libconfigurable_gridcomp)
    │     internal: 48 fields, connected directly to MOIST's imports
    ├── FakeGOCART (libFakeGOCART_GridComp)
    │     export: AERO (empty state), imported by MOIST
    ├── MOIST (libGEOSmoist_GridComp)
    │     export: MTR (service provider), imported by DataTracers
    │     import: TRADV, TRDIFF, TRCHEM, TRANA (service subscriptions, unconnected)
    └── DataTracers (libconfigurable_gridcomp)
          internal: TR1, subscribed to MOIST's MTR (items: [TR1])
```

Children of `Root`:

- `MOIST` (`moist.yaml`): internal state restarted from `checkpoints/2000-04-14T21:20:00/MOIST_internal.nc`; import, internal and export checkpoints.
- `DataForMoist` (`data-for-moist.yaml`): the 48 unconditional imports of `../GEOSmoist_GridComp/Moist_StateSpecs.rc`. The 8 conditional ones (`ALH`, `TAUOROX`, `TAUOROY`, `WSUB_CLIM`, `KM`, `RI`, `NCPL_CLIM`, `NCPI_CLIM`) are not registered with the defaults `WSUB_OPTION: 3` and `USE_NCLOUD_CLIM: false`.
- `FakeGOCART` (`fake-gocart.yaml`): an empty `AERO` state. The library is built by `../@GEOSradiation_GridComp/regression/FakeGOCART_GridComp`.
- `DataTracers` (`data-tracers.yaml`): one tracer, `TR1` (fill `1.0e-6`, no `ESMF_Info` metadata), subscribed to MOIST's convective-transport service `MTR`, so GF transports a subscriber's field as well as its own `CNV_TR`.

Settings (`moist.yaml`), from the production `AGCM.rc.tmpl`, `set_resolution_params.csh` (c12) and `gcm_setup` (L181):

- `CLDMICR_OPTION: GFDL_1M`, `CONVPAR_OPTION: GF`, `SHALLOW_OPTION: UW`.
- `GFDL_MP_KLID: 19.0`; `GF_DT` and `UW_DT` `1200`, equal to the clock step.
- GF: `DEEP: 1`, `SHALLOW: 0`, `CONGESTUS: 1`, `CLOSURE_DEEP: 0`, `CLOSURE_SHALLOW: 7`, `CLOSURE_CONGESTUS: 3`.
- `USE_AEROSOL_NN: false`: cloud droplet and ice numbers come from the fixed land/ocean CCN/IN numbers, so the empty `AERO` is enough.
- Everything else at the code defaults.

Other inputs:

- `input.nml` with an empty `&gfdl_mp_nml`, as production's `fvcore_layout.rc`. `gfdl_mp_init` stops if the file is missing.

Services:

- MOIST's four subscriptions, `TRADV` (advection), `TRDIFF` (turbulent diffusion), `TRCHEM` (chemistry) and `TRANA` (analysis), are left unconnected. Their providers change MOIST's internals between MOIST's runs, not during them, so their effect is already in the restarted internal state.
- Only MOIST's own states are compared, not the transported tracer `TR1`.

Exercised since the first run (2026-11-04), which matches the MAPL2 capture: four unconnected service subscriptions from one component, the `MTR` connection (MOIST is the provider) and an empty `AERO` with `USE_AEROSOL_NN: false`. See `../GEOSmoist_GridComp/mapl3-porting-notes.md`, step 10.

### `moist-gfdl1m-gf-shoc-sa`

MOIST alone, with GFDL_1M microphysics, GF deep convection and no UW: the stock v12 L181 setup. For `LM /= 72`, gcm_setup enables the SHOC-MF block of `AGCM.rc.tmpl` and comments out UW (`gcm_setup:341-348`). Of that block Moist reads `SHALLOW_OPTION: NONE` and `PDFSHAPE: 6` (`GEOS_GFDL_1M_InterfaceMod.F90:367`); the other keys belong to Turbulence and reach Moist only through the captured imports. So `UW_Setup` and `UW_Run` are not called (`GEOS_MoistGridComp.F90:185, 757, 787`), and the `CUSH` internal is not registered.

Everything else is copied from `moist-gfdl1m-gf-uw-sa`:
- `moist.yaml` differs only in `SHALLOW_OPTION`, `PDFSHAPE` and `UW_DT` (removed).
- `root.yaml` is unchanged, since UW adds no imports (`KPBL_SC` is a driver import).
- `history.yaml` is the UW case's list minus the 22 exports only `UW_Run` writes, confirmed by a sentinel run (2026-11-05): every listed export is written, and every written export is listed except `REFL10CM_MAX` (left out), the always-active `ADD2EXPORT` internal `CUFRC_SC` and the `MTR` bundle.
- `cmpchk.py` points at `~/input/moist/gfdl1m-gf-shoc/C12L181`.

### `moist-bacm1m-gf-uw-sa`

MOIST alone, with BACM_1M microphysics, GF deep convection and UW shallow convection at C12**L72**: the stock L72 setup, where BACM_1M is the default microphysics (`gcm_setup:356-359`) and UW is on and SHOC off (`gcm_setup:342-348`). Unlike the GFDL_1M cases (L181), this also changes GF's code path: at `LM == 72` the code defaults `USE_GF2020` to 0, so GF runs `ConvPar_GF_GEOS5.F90`, with `ZERO_DIFF_*` at 1 (`GEOS_GF_InterfaceMod.F90:139-147`), and UW defaults `JASON_UW`/`JASON_MFD_SC` to true (`GEOS_UW_InterfaceMod.F90:91, 95`). This is the first MAPL3 run of the GF_GEOS5 path. A mismatch with MAPL2 could therefore come from BACM_1M or from GF_GEOS5; GF_GEOS5's outputs are its own exports and the `*_DC` internals, so compare those first.

What MAPL2's setup changes for BACM_1M (v11.11.0):
- `set_resolution_params.csh:453-457` sets `DT = CONV_DT = CHEM_DT = 450` at every resolution: a 450 s heartbeat (`CAP.rc.tmpl:9`) and `GF_DT`/`UW_DT` 450 (`AGCM.rc.tmpl:64-65`). So `cap_driver.yaml` uses `PT7M30S`, and the case starts at `2000-04-14T21:15:00`, the second step after 21:00 (21:20 is not a step).
- `gcm_setup:1115-1123` writes `GFDL_MP_KLID` only for GFDL_1M, so `moist.yaml` has none.
- At L72, BACM_1M's ice-fall factor is 1.0 (`GEOS_BACM_1M_InterfaceMod.F90`, `BACM_1M_Initialize`; 0.25 at L181).

Otherwise copied from `moist-gfdl1m-gf-uw-sa`:
- `root.yaml`: `num_levels: 72` (73 edges). `data-for-moist.yaml` unchanged: BACM_1M and GF_GEOS5 fetch only imports the driver already declares.
- `moist.yaml`: `CLDMICR_OPTION: BACM_1M`, `GF_DT`/`UW_DT` 450, `AGCM.GRIDNAME: PE12x72-CF` (needed by `BACM_1M_Initialize` for its default `MINRHCRIT*`, as GWD), and the GF plume/closure keys, which `GF_Initialize` reads on both GF paths (`:160-165`, `:418-423`).
- `history.yaml`: 304 exports, the UW case's 281 minus the 24 only GFDL_1M writes (the 20 in its scheme module, and `NACTR`, `DBZ_1KM`, `DBZ_TOP`, `DBZ_M10C` from `compute_radar_diagnostics`, which only `GFDL_1M_Run` calls, `GEOS_GFDL_1M_InterfaceMod.F90:1561`; the first run had these four NaN, MAPL2 zeros), plus the 47 BACM_1M writes that the GFDL_1M case does not (including `PFL_SC`/`PFI_SC`). GF2020-only outputs (`ENTR_*`, `SGS_VVEL_*`, `CNV_TOPP_*`, `WQT_DC`, `LFR_GF`, `DNDCNV`, `DNICNV`, `MDNDP`) stay in: GF_GEOS5 leaves them at the zeros their locals start with (`GEOS_GF_InterfaceMod.F90:740-757`). Derived from the code and confirmed by a sentinel run (2026-11-05): every listed export is written, and every written export is listed except the `MTR` bundle.
- `cmpchk.py` points at `~/input/moist/bacm1m-gf-uw/C12L72`.

The MAPL2 capture needs a C12L72 MAPL2 experiment made with `gcm_setup` (72 levels, BACM_1M, and v14 land boundary conditions, which must be typed because the L72 default is NL3, `gcm_setup:1165-1169`; the other answers as for the L181 runs), so that the other L72 defaults (RRTMG radiation, L072 GOCART tuning, `gcm_setup:330`, `:2587`) are right, and L72 restarts for 2000-04-14 21:00. Add to its `AGCM.rc`: `COMPONENT_TO_RECORD: MOIST`, `TEST_FRAMEWORK: .true.`, `TARGET_TIME: '2000-04-14 21:15:00'`, and `USE_AEROSOL_NN: .false.`; check `HEARTBEAT_DT: 450` in its `CAP.rc` and `GF_DT`/`UW_DT: 450`; and use `~/HISTORY.bacm.rc`, which is `~/HISTORY.rc` without the five GFDL_1M-only exports `REF_DBZ`, `REF_DBZ_MAX`, `REF_DBZ_1KM`, `REF_DBZ_TOP` and `REF_DBZ_M10C` (registered only in `GFDL_1M_Setup`, `GEOS_GFDL_1M_InterfaceMod.F90:215-251`; MAPL2 History fails on them otherwise). `~/AGCM.bacm.rc` is `~/AGCM.rc` with these changes and the L72 ones (`*LM: 72`, BACM_1M on and GFDL_1M off, `GF_DT`/`UW_DT: 450`, RRTMG instead of RRTMGP; v14 needs no `AGCM.rc` change, since it keeps `NCAR_NRDG: 16`, `gcm_setup:1205`, `:2082-2084`), so it can replace that experiment's `AGCM.rc`.

## Baselines

There are two kinds of baseline: a MAPL2 reference capture, and the MAPL3 checkpoints that ctest compares against.

### MAPL2 reference captures (`~/input/moist/<case without moist- and -sa>/C12L181`: `gfdl1m-gf-uw` made 2026-11-04, `gfdl1m-gf-shoc` 2026-11-05)

- **Plan** (`../GEOSmoist_GridComp/mapl3-porting-notes.md`, step 0): C12L181, `COMPONENT_TO_RECORD: MOIST`, MOIST's import, internal and export states before and after its Run at `2000-04-14T21:20:00`, converted to the new cubed-sphere file layout with `GEOS_Restart_Utilities/convertOldCStoNewCS.py`.
- **Run configuration** must match `moist.yaml`: the same scheme options and settings, a 1200 s step and MOIST's run step equal to it, `USE_AEROSOL_NN: false`, and a `HISTORY.rc` requesting every MOIST export (MAPL2 allocates exports only when they are requested).
- **Used for**:
  - The restart files in `checkpoints/2000-04-14T21:20:00/`: `DataForMoist_internal.nc` from the "before" import state and `MOIST_internal.nc` from the "before" internal state.
  - The MAPL2-vs-MAPL3 check done by `cmpchk.py` against the "after" states (see below).

### ctest baselines

- `checkpoints/last/*.nc` in the case's regression-data directory: output of the MAPL3 port itself, generated by running the case from the restarts above once it matches the MAPL2 capture.
- `run_case.cmake` compares every `*.nc` in the baseline's `checkpoints/last` with the run's (`compare_results`, NaNs equal), excluding `lons`, `lats`, `corner_lons`, `corner_lats`. The comparison must be exact.
- The comparison runs **only for IntelLLVM (ifx) builds**. With other compilers the case runs but is not compared.

## Regression data

`esma_add_regression_tests` sets the data directory at configure time to `$LOCAL_REGRESSION_DATA_DIR/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp`. If `LOCAL_REGRESSION_DATA_DIR` is unset, it uses the same path under `<build>/regression-data/`, synced from S3 by the `moist/sync_data` fixture. Layout:

```
GEOSmoist_GridComp/
  moist-gfdl1m-gf-uw-sa/checkpoints/
    2000-04-14T21:20:00/
      DataForMoist_internal.nc    # restarts, copied into the run directory
      MOIST_internal.nc
    last/
      MOIST_import.nc             # ctest baseline
      MOIST_internal.nc
      MOIST_export.nc
  moist-gfdl1m-gf-shoc-sa/checkpoints/
    2000-04-14T21:20:00/
    last/
  moist-bacm1m-gf-uw-sa/checkpoints/
    2000-04-14T21:15:00/
    last/
```

If a case's data directory is missing, `run_case.cmake` prints "Regression data not found ... skipping" and the test **passes without running**. If it has restarts but no `checkpoints/last/`, the case runs and `compare_results` finds no baseline files, so the test passes if the run does (`@cmake/esma_support/esma_regression_run_helpers.cmake`, `compare_results`).

## Run time

Measured 2026-11-09 with `ctest -L REGRESSION` (ifx Debug, 6 PEs, `logging.yaml` in place, MAPL `develop` plus GEOS-ESM/MAPL#5500, which skips the regex in `VirtualConnectionPt%matches` for literal names, not yet merged): `moist-gfdl1m-gf-uw-sa` 7.9 s, `moist-gfdl1m-gf-shoc-sa` 7.2 s, `moist-bacm1m-gf-uw-sa` 6.5 s.

Almost all of it is MAPL, not Moist. In a perf profile of `moist-bacm1m-gf-uw-sa` (rank 3):

- MAPL's `CAP` timer is 4.7 s: Initialize 4.5 s, run 0.2 s. Moist's own `Run` is about 0.03 s.
- The other ~2 s of the test is `mpiexec`, MPI and ESMF startup, copying the restarts and the `nccmp` comparison (estimated from the difference between test time and `GEOS.x` time, not profiled).
- Of the initialize CPU time, about 33% is ESMF Info (JSON field metadata), mostly from `field_info_set_internal` / `field_info_get_internal` (and their `restart_mode` variants) in MAPL `infrastructure/field/FieldInfo.F90`. About 20% of all CPU time (~1 s) is C++ exceptions: ESMF's `has_key_json` and `update_json_pointer` (`ESMCI_Info.C`, ESMF in Baselibs 8.33.0) find a missing key by catching `json::out_of_range`, and MAPL checks for or sets missing keys on every field.
- About 15% is in `ESMF_ArrayGet` / `ESMF_DistGridGet` (callers not resolved; the profile's stacks stop in ESMF), about 10% in `StateRegistry` (extend, add_to_states, set_geometry, connections) and 6% in `StateRegistry%filter`.

History of the reductions (same cases): 17 s before `logging.yaml` (field-dictionary warnings, see Common setup), then 8.7-9.5 s with it, then the times above with the literal-name match. Further gains would have to come from ESMF (look keys up without exceptions) or MAPL (avoid probing missing keys in `FieldInfo.F90`); together these are worth at most ~1-1.5 s per test.

## Running

```
ctest -L REGRESSION    # all regression tests in the build
ctest -R moist-gfdl1m-gf-uw-sa
ctest -R moist-gfdl1m-gf-shoc-sa
ctest -R moist-bacm1m-gf-uw-sa
```

### Comparing against MAPL2 with `cmpchk.py`

After a run, each case's `cmpchk.py`, run from the run directory, compares `checkpoints/last/MOIST_<state>.nc` with the MAPL2 "after" capture (`MOIST_<state>_after_runPhase1.nc`). It prints, per variable, the min, max and L2 norm of the difference. Undefined points are left out: MAPL2's (`1.0e15`) are the capture's `_FillValue` and masked on read, MAPL3's (`-huge`) are masked explicitly. So MOIST's undefined points (e.g. `CAPE`, `DETR_SC`) are not reported as differences (porting notes, step 8, `MAPL_UNDEF`). Dimension variables are skipped, as are variables missing from the run's file.

`baseline_dir` points at the converted "after" capture under `/home/pchakrab/input/moist/`; edit it to point at your copy. The script needs `numpy` and `netCDF4` (e.g. after `module load ifx-stack`).
