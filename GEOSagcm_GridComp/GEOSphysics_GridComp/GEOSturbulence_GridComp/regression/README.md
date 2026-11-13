# Turbulence regression tests

Standalone MAPL3 regression case for TURBULENCE, compared exactly against a stored MAPL3 baseline that matches a MAPL2 capture. Layout and conventions follow the Moist cases (`../../regression/README.md`) and GOCART2G's per-phase cases (`@GOCART/ESMF/GOCART2G_GridComp/regression/README.md`).

| Case | What runs | States compared by ctest | Status |
|---|---|---|---|
| `turbulence-shocmf-sa/run1` | TURBULENCE Run1 | TURBULENCE import, internal, export; DataTracers internal (the diffused tracers) | matches the MAPL2 capture (cmpchk.py, 2026-11-13); reproduced by a second run |
| `turbulence-shocmf-sa/run2` | TURBULENCE Run1, then Run2 | same | matches the MAPL2 capture (cmpchk.py, 2026-11-13); reproduced by a second run |
| `turbulence-default-sa/run1` | TURBULENCE Run1 (SHOC and EDMF off) | same | matches the MAPL2 capture (cmpchk.py, 2026-11-13, 149 variables); reproduced by a second run |
| `turbulence-default-sa/run2` | TURBULENCE Run1, then Run2 (SHOC and EDMF off) | same | matches the MAPL2 capture (cmpchk.py, 2026-11-13, 158 variables); reproduced by a second run |

The cases are listed in `cases.txt` and registered by `CMakeLists.txt` through `esma_add_regression_tests` (`NAME turbulence`; ctest label `REGRESSION`, test name = case name), added by `add_subdirectory (regression)` in `../CMakeLists.txt`. `run_case.cmake` is the Moist one.

## Common setup

- **Configuration:** `turbulence-shocmf-sa`: the stock v12 C12L181 SHOC-MF setup (`TURBULENCE_DO_SHOC: 1`, `EDMF_DOMF: 1`, `PDFSHAPE: 6`, the tuned `EDMF_*` and `TURBULENCE_SHC_LENFAC*` keys, `HGT_SURFACE: 0.0`), no waves, no SCM surface, refresh every step (no `TURBULENCE_DT`). `turbulence.yaml` holds every Turbulence key the MAPL2 run's `AGCM.rc` set; keep it in sync with that run.
- **Geometry, layout, logging:** as the Moist cases: C12L181 (vertical grid from `PLE`), 6 MPI tasks, `OMP_NUM_THREADS=1`, `logging.yaml`.
- **Clock:** one 20-minute step from `2000-04-14T22:00:00` (`cap_restart.yaml`), the fourth step of the MAPL2 run, which started at 21:00. The refresh alarm rings on the first step, so `REFRESH` runs, as in MAPL2.
- **No ExtData; History only activates exports** (`history.yaml`, a collection that never writes).
- **Inputs:** `DataForTurbulence` holds the 44 unconditional imports as internals, re-exported (`<self>` connections, as `DataForMoist`) and connected field by field to TURBULENCE (`root.yaml`).
- **Diffused tracers:** `DataTracers` holds the 47 fields of MAPL2's `TR` bundle in the capture (Physics' `QW`, Moist's 8 water species, the PChem, GOCART and TR tracers) and subscribes them to TURBULENCE's `TRDIFF` service, with the capture's long names as standard names (TURBULENCE recovers `Q`, `QLLS`, `QILS`, ... from them). Its internal state is checkpointed and compared, so the diffusion of every tracer is checked. Names: the capture's `TR_<name>` without the prefix, `.` replaced by `_` (`CA.ocphobic` -> `CA_ocphobic`).
- **Unconnected:** TURBULENCE's `TRADV` subscription (SHOC on), as Moist's.
- **Exports:** History activates the exports MAPL2 had allocated (the export capture), less the `TRI` bundle members, the fields that are internals in MAPL3, and the 7 never written with SHOC on (`DU`, `DBUOY`, `WERAD`, `ZCLD`, `ZRADBS`, `ZRADML`, `ZSML`). Which exports are requested changes what TURBULENCE's `ALLOC=` and `associated()` tests do, so the list reproduces the MAPL2 run (`../mapl3-porting-notes.md`, step 6).

## `turbulence-default-sa`

The stock v12 C12L181 default setup, with the SHOC-MF block of `AGCM.rc` commented out: `TURBULENCE_DO_SHOC` and `EDMF_DOMF` at their default 0, so Turbulence runs the Louis/Lock path. At `LM /= 72` the code defaults are `JASON_*` false, `PBLHT_OPTION` 5, `SMTH_HGT` 300 and `PDFSHAPE` 1. `turbulence.yaml` has only `HGT_SURFACE: 0.0` and `USE_WAVES: 0`. Everything else is as `turbulence-shocmf-sa` (same clock, geometry, 44 imports, 47 tracers, the same 9 `phase2_fields`; `data-for-turbulence.yaml`, `data-tracers.yaml`, `data-phase2.yaml` and the two roots differ only in comments), and so are `run1` and `run2` (above). With SHOC off there is no `TRADV` subscription.

Exports: the MAPL2 export capture lacks `BRUNTSHOC`, `EDMF_FRC`, `EDMF_MF`, `EDMF_WQT`, `EDMF_WSL` and `LOBUKHOV` (SHOC and EDMF off). `history.yaml` requests the captured exports, less the `TRI` members, the internals, and eight that must not be requested because this configuration never writes them (MAPL3 would leave them uninitialized, MAPL2 zeroed them): `DU`, `EDMF_DQRDT`, `EDMF_DQSDT`, the PDF moments `SL2`, `SL3`, `W2`, `W3` (`PDFSHAPE < 5`), and `ISOTROPY` (with SHOC off and `TKE` requested, its level 0 is read uninitialized: floating invalid trap; unrequested, it is a zero-filled work array as in MAPL2). `DBUOY`, `WERAD`, `ZCLD`, `ZRADBS`, `ZRADML` and `ZSML` are written here (Lock) and compared.

Data: converted capture in `~/input/turbulence/default/C12L181/new-style`; restart subsets as for `turbulence-shocmf-sa`, in `GEOSturbulence_GridComp/turbulence-default-sa/` of the regression data (the same layout as below).

## `run1`

Root is a configurable gridcomp running `DataForTurbulence`, `DataTracers` and TURBULENCE with `run_phases: [Run1]`. Restarts from the MAPL2 states before Run1. `history.yaml` also leaves out the 9 exports only Run2 fetches and writes (`DPDTTRB`, `INTDIS`, `QT`, `SAFUPDATE`, `SL`, `SRFDIS`, `TOPDIS`, `WQT`, `WSL`): MAPL2's capture after Run1 holds their previous-step values.

## `run2`

Run2 cannot run alone: it continues from what Run1 leaves in the `RESTART: SKIP` internals (`DKS`, `DKV`, `DKQ`, `EKV`, `FKV`, `SINC`, the tendencies) and in exports (`SFLX`, `QFLX`, `SDFL`, `QDFL`, ...), which MAPL3 cannot restart; and between the phases MAPL2's Surface (and Physics) change 9 of its imports (`EVAP`, `QDEL`, `QHAT`, `S`, `SDEL`, `SH`, `T`, `TH`, `USTAR`: every import that differs between the "after Run1" and "before Run2" captures). The configurable Root runs a child's phases back to back, so `run2`'s Root is `TurbulencePhasesDriver_GridComp` (this directory, test only):

```
CAP (libMAPL.cap; ExtData off, History on)
└── Root (libTurbulencePhasesDriver_GridComp; C12L181)
    │     imports X and X_PHASE2 for each X in phase2_fields (root.yaml)
    ├── DataForTurbulence   44 imports, restarted before Run1; X also to Root
    ├── DataPhase2          the 9 fields, restarted before Run2; to Root as X_PHASE2
    ├── DataTracers         47 tracers, subscribed to TRDIFF
    └── TURBULENCE          run by Root: Run1; X = X_PHASE2; Run2
```

Root's import `X` and TURBULENCE's import `X` are the same field (both connected to `DataForTurbulence`), so the copy changes what Run2 sees. Restarts from the MAPL2 states before Run1, plus `DataPhase2_internal.nc` from the import state before Run2; compared with the states after Run2.

## Rebuilding a case's data

`tools/` has the scripts that turn a MAPL2 capture into a case (conversion, the restart files, the data-component YAML, the History list, running by hand), and `tools/README.md` says why they are needed and how to use them. ctest does not use them.

## Regression data

`$LOCAL_REGRESSION_DATA_DIR/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSturbulence_GridComp` (set at configure time; `@cmake/esma_support/esma_regression_tests.cmake`):

```
GEOSturbulence_GridComp/turbulence-shocmf-sa/
  run1/checkpoints/
    2000-04-14T22:00:00/
      DataForTurbulence_internal.nc  # MAPL2 import state before Run1, without TR_*
      DataTracers_internal.nc        # its TR_* fields, renamed
      TURBULENCE_internal.nc         # MAPL2 internal state before Run1
    last/                            # ctest baseline: this case's MAPL3 output
      TURBULENCE_{import,internal,export}.nc
      DataTracers_internal.nc
  run2/checkpoints/
    2000-04-14T22:00:00/             # as run1, plus
      DataPhase2_internal.nc         # the 9 phase2_fields of the import state before Run2
    last/
```

If a case's data directory is missing, the test passes without running (`run_case.cmake`).

## MAPL2 capture

GEOSgcm on `develop` (Turbulence sources identical to `develop` @ `f31a2151`), C12L181 SHOC-MF, 1200 s steps from `2000-04-14T21:00:00`, with `COMPONENT_TO_RECORD: TURBULENCE`, `TEST_FRAMEWORK: .true.` and `TARGET_TIME: '2000-04-14 22:00:00'`: TURBULENCE's import, internal and export states before and after Run1 and Run2. The run's History requested only the stock collections, so the export capture holds the exports MAPL2 allocated in a normal run. Converted to the new cubed-sphere layout with `GEOS_Restart_Utilities/convertOldCStoNewCS.py`; the restarts above are subsets of the converted files.

## Comparing against MAPL2 with `cmpchk.py`

After a run, each case's `cmpchk.py`, run from the run directory, compares `checkpoints/last/TURBULENCE_<state>.nc` with the MAPL2 capture after the case's phase (`TURBULENCE_<state>_after_runPhase<N>.nc`), and `checkpoints/last/DataTracers_internal.nc` with the capture's `TR_*` imports. It prints the min, max and L2 norm of each difference and lists the variables that differ. Undefined points are left out (MAPL2's `1.0e15` is the capture's `_FillValue`; MAPL3's `-huge` is masked). `baseline_dir` points at the converted capture; edit it to point at your copy. Needs `numpy` and `netCDF4`.

## Running

```
ctest -R turbulence-shocmf-sa      # both SHOC-MF cases
ctest -R turbulence-default-sa     # both default-setup cases
ctest -R turbulence-shocmf-sa/run2
```
