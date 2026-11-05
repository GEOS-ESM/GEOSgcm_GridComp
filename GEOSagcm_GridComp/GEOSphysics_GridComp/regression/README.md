# GEOSphysics regression tests

Standalone MAPL3 regression cases for the Physics components that have no regression directory of their own. Each case runs one component under a configurable `Root`, with every import supplied from a restarted "DataFor*" gridcomp, and compares the resulting checkpoints against a stored baseline. GWD (`../GEOSgwd_GridComp/regression/`) and Radiation (`../@GEOSradiation_GridComp/regression/`) keep their own.

| Case       | Component run | States compared by ctest       | Status |
|------------|---------------|--------------------------------|--------|
| `moist-gfdl1m-gf-uw-sa` | MOIST         | MOIST import, internal, export | matches the MAPL2 capture (cmpchk.py; REFL10CM_MAX left out); regression data in ~/regression-data only (2026-11-04), not on /ford1 |

The cases are listed in `cases.txt` and registered by `CMakeLists.txt` through `esma_add_regression_tests` (`NAME moist`; ctest label `REGRESSION`, test name = case name). `CMakeLists.txt` is added by `add_subdirectory (regression)` in `../CMakeLists.txt`.

## Common setup

- **Geometry**: C12 cubed sphere, 181 levels (C12L181), model vertical grid defined by the `PLE` field (`root.yaml`).
- **Layout**: 6 MPI tasks (`mapl.yaml`, `run_case.cmake`), `OMP_NUM_THREADS=1`.
- **Clock**: one 20-minute step (`PT20M`) from `2000-04-14T21:20:00` (`cap_restart.yaml`); the final checkpoint is written at 21:40 (`checkpoints/last`). The `start`/`stop` in `cap_driver.yaml` only bound the clock.
- **No ExtData**: `run_extdata: false`. **History only activates exports**: `run_history: true` with a `HIST` collection (`history.yaml`) whose start is after the cap stop, so it never writes a file.
- **Inputs**: the case's `DataFor*` gridcomp (`libconfigurable_gridcomp`) holds every import as an internal field, restarted from `checkpoints/2000-04-14T21:20:00/DataFor*_internal.nc`, and is connected field-by-field to the component in `root.yaml`.
- **Exports**: only those the case writes are activated, through `history.yaml` (not `activate_all_exports`: MAPL3 leaves an allocated export without a fill value uninitialized, so an export the configuration never writes would hold garbage). Regenerate the list when the configuration changes (`../GEOSmoist_GridComp/mapl3-porting-notes.md`, "Unwritten exports").

## Cases

### `moist-gfdl1m-gf-uw-sa`

MOIST alone, with GFDL_1M microphysics, GF deep convection and UW shallow convection (the production default). The other schemes `_FAIL` in the MAPL3 port.

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

Untested until the first run: four unconnected service subscriptions from one component, the `MTR` connection (MOIST is the provider) and an empty `AERO` with `USE_AEROSOL_NN: false`. See `../GEOSmoist_GridComp/mapl3-porting-notes.md`, step 10.

## Baselines

There are two kinds of baseline: a MAPL2 reference capture, and the MAPL3 checkpoints that ctest compares against.

### MAPL2 reference capture (`moist-gfdl1m-gf-uw-sa`: made 2026-11-04, `~/input/moist/C12L181`)

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
```

If a case's data directory is missing, `run_case.cmake` prints "Regression data not found ... skipping" and the test **passes without running**. This is the state of `moist-gfdl1m-gf-uw-sa` today.

## Running

```
ctest -L REGRESSION    # all regression tests in the build
ctest -R moist-gfdl1m-gf-uw-sa
```

### Comparing against MAPL2 with `cmpchk.py`

After a run, `moist-gfdl1m-gf-uw-sa/cmpchk.py` from the run directory compares `checkpoints/last/MOIST_<state>.nc` with the MAPL2 "after" capture (`MOIST_<state>_after_runPhase1.nc`). It prints, per variable, the min/max/L2 norm of the absolute difference, and reports separately any cells that are undefined on one side only. A cell counts as undefined if it is NaN, `|huge|`, MAPL2's `MAPL_UNDEF` (`1.0e15`) or either file's `_FillValue`; MAPL3's `MAPL_UNDEF` is `-huge`. So MOIST's undefined points (e.g. `CAPE`, `DETR_SC`), which MAPL2 fills with `1.0e15` and MAPL3 with `-huge`, are not reported as differences (porting notes, step 8, `MAPL_UNDEF`). Grid and dimension variables are skipped, as are variables missing from either file.

`baseline_dir` is a placeholder (`/home/pchakrab/input/moist/C12L181/after`) until the capture exists; edit it to point at your copy. The script needs `numpy` and `netCDF4` (e.g. after `module load ifx-stack`).
