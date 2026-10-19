# `adv-dyn+gocart` regression case

FV3 offline tracer advection (ADV, `AdvCore`) of GOCART2G's aerosol tracers
through the MAPL3 `TRADV` service, at C12L181, two 30-minute steps, with
GOCART2G's own physics (emissions, chemistry, deposition) running before the
dynamics in each step.

In MAPL2, tracers reached the advection core because their internals were
marked `FRIENDLYTO=DYNAMICS`, and Physics assembled them into a bundle. In
MAPL3, each tracer owner subscribes to the `TRADV` service and the provider
(ADV, or DYN) receives a bundle of the subscribers' fields. This case is the
only one in which aerosol tracers reach a provider.

It is `adv-dyn` (FVdycoreCubed) at 181 levels instead of 91, with
`DataForGOCART2G` and GOCART2G added ahead of DYN and ADV,
`ADV/TRADV -> GOCART2G/TRADV` connected, and `adv-dyn`'s `DataTracers`
(test tracers `Q001`, `Q002`) removed, so ADV advects only aerosols.

## Gridcomp hierarchy

```
Cap                      (cap_gridcomp.yaml; ExtData off)
├── Root                 (root.yaml; geometry from DYN)
│   ├── DataForGOCART2G  data component: GOCART2G's 43 imports
│   ├── GOCART2G         aerosol parent (Run1, Run2)
│   │   ├── DU           dust: subscribes DU (5 bins) to TRADV
│   │   ├── SS           sea salt: subscribes SS (5 bins) to TRADV
│   │   └── SU           sulfate: subscribes DMS, SO2, SO4, MSA to TRADV
│   ├── DataMoist        data component (Q); subscribes to TRADV
│   ├── DataPhisVarflt   data component: DYN's PHIS, VARFLT
│   ├── DataAna          data component: DYN's analysis inputs
│   ├── DYN              FV3 dycore (FVdycoreCubed), cold start
│   └── ADV              FV3 AdvCore: TRADV provider
└── History              (history.yaml)
```

GOCART2G adds DU, SS and SU at run time from `ACTIVE_INSTANCES` in
`GOCART2G_GridComp.yaml` (CA and NI are empty). The children's `TRADV`
imports propagate up to GOCART2G, and Root's `connections` link
`ADV/TRADV -> GOCART2G/TRADV`.

## Run sequence

Root (`ConfigurableGridComp`) runs its children in `root.yaml` order. Each
child runs all of its `run_phases` in order (or the single phase `run` if none
are listed) before the next child starts. GOCART2G runs each phase on its
children in `ACTIVE_INSTANCES` order (DU, SS, SU), skipping any whose name
contains "data".

Within each step:

1. `DataForGOCART2G` `run`: nothing to do (no `run` section); holds
   GOCART2G's 43 inputs from its restart.
2. `GOCART2G` `Run1`: `Run1` on DU, SS, SU (emissions).
3. `GOCART2G` `Run2`: `Run0` on DU, SS, SU; zero the total-aerosol
   diagnostics; `Run2` on DU, SS, SU (chemistry, deposition, diagnostics);
   compute the totals (AOD, PM, ...).
4. `DataMoist`, `DataPhisVarflt`, `DataAna` `run`: nothing to do (no `run`
   sections); they hold DYN's inputs.
5. `DYN` `run`: the FV3 dynamics. It advects its own `TRADV` (DataMoist's
   fields) and exports this step's `MFX`, `MFY`, `CX`, `CY`, `PLE0` and
   `PLE1` for ADV. Its `run_add_incs` phase is not called.
6. `ADV` `run`: AdvCore advects its `TRADV` bundle (14 tracers: DU's 5 bins,
   SS's 5 bins, and SU's `DMS`, `SO2`, `SO4`, `MSA`) with DYN's fluxes from
   step 5. The DU and SS bins are views of DU's and SS's internals, so
   advection updates GOCART2G's state directly (as DYN's does DataMoist's
   `Q`).

The case runs two steps (00:00 to 01:00), after which the checkpoint is
written to `checkpoints/last`.

What this means:

- **Each data component runs just before the component it feeds**, so any
  per-step input one might set later is in place when it is read.
- **Physics before dynamics.** GOCART2G updates its tracers first, and ADV
  then moves the result, so the tracers checked at the end of each step have
  been through both.
- **Advection uses the current step's winds.** ADV runs after DYN, so it uses
  the mass fluxes DYN has just computed.
- **DU's restart is essentially aerosol-free** (largest mixing ratio about
  1e-32). Dust comes from GOCART2G's emissions in the first step (up to about
  1e-6 by the end), so DU's bins carry real values through advection.

## What is tested

- **Nested subscribers.** DU, SS and SU subscribe to `TRADV` from inside
  GOCART2G (DU: `DU`; SS: `SS`; SU: `DMS`, `SO2`, `SO4`, `MSA`). Their
  unconnected imports propagate to GOCART2G's `TRADV`, and the single
  root-level connection reaches all three.
- **Providers split the tracers by name.** DataMoist (`Q`) and GOCART2G both
  connect to ADV's `TRADV`, and DataMoist also to DYN's. DYN advects only the
  cloud and water species (`field_is_cloud_water_species` in `fv_shared.F90`:
  `Q`, `QLCN`, `QLLS`, ...); ADV skips them and advects the rest. So DYN
  advects `Q` and ADV the 14 aerosol tracers, and nothing is advected
  twice.
- **Fields with an ungridded dimension.** `DU` and `SS` are single fields with
  5 bins. MAPL3 splits such fields into one 3D field per bin when it builds
  the service bundle (as MAPL2 did for friendly bundles), so AdvCore, which
  handles 3D tracers only, advects every bin. Without the split ADV fails
  with `ESMF_RC_ARG_INCOMP` reading the 4D field.
- **Advection conserves mass.** With `ADV_CORE_REPORT_TRACER_MASS: true`, ADV
  logs any tracer whose global mass changes during advection by at least
  float32 epsilon, and logs none for any of the 14 tracers.
- **Physics and advection on the same fields.** GOCART2G's Run1/Run2 and ADV
  both write to DU's, SS's and SU's internals within one step, through the
  same memory (the `TRADV` fields are views, not copies).

The comparison checks all states in `checkpoints/last`: ADV, DYN, GOCART2G,
DU, SS and SU.

## What is not tested

- **Advection on its own.** With physics running, the change in a tracer
  over a step is physics plus advection, so the baseline cannot isolate
  advection. ADV's mass report is the advection-only check. When GOCART2G
  had no run phases, every aerosol tracer changed under advection alone in
  34-56% of cells by up to about 3%; for DU that was moving values of about
  1e-32.
- **Diffusion and convection.** Turbulence and moist have no MAPL3 service
  providers yet.
- **DYN as the aerosol provider.** GOCART2G is connected to ADV only.
- **The full model's routing.** In the AGCM the connection would be
  SUPERDYN/`TRADV` -> PHYS/`TRADV`, with GOCART2G nested under CHEMISTRY and
  PHYSICS; here GOCART2G is a direct child of the root.
- **A MAPL2 baseline.** DYN cold-starts on an idealized date (as in
  `adv-dyn`), so there is no MAPL2 counterpart. The baseline in
  `checkpoints/last` is this case's own MAPL3 output; the case guards
  against changes.

## Components and inputs

- `DYN` (`dyn.yaml`, `input.nml`): FV3, cold start (`COLDSTART: true`,
  `test_case = 5`), `npz = 181`; provides `MFX`, `MFY`, `CX`, `CY`, `PLE0`,
  `PLE1` to ADV.
- `ADV` (`adv.yaml`): FV3 `AdvCore`, the `TRADV` provider.
- `DataMoist`, `DataPhisVarflt`, `DataAna`: as in `adv-dyn`. `adv-dyn`'s
  `DataTracers` is left out.
- `DataForGOCART2G` (`data-for-gocart2g.yaml`): GOCART2G's 43 imports.
- `GOCART2G` (`GOCART2G_GridComp.yaml`, `{DU,SS,SU}_instance_*.yaml`): DU, SS
  and SU active, CA and NI empty (not ported); run phases `Run1`, `Run2`.

The GOCART restarts are the MAPL2 C12L181 state (2000-04-14), used as initial
tracer values on this case's idealized date. Regression data, under
`GEOSagcm_GridComp/adv-dyn+gocart`:

- `checkpoints/1891-03-02T00:00:00/`: `DataForGOCART2G`, `DU`, `SS` and `SU`
  restarts (DYN cold-starts and has none);
- `checkpoints/last/`: the baseline.

`run_case.cmake` links `ExtData/chemistry` from
`GEOSphysics_GridComp/GEOSchem_GridComp/GOCART2G_GridComp/ExtData` (the
aerosol optics tables GOCART2G reads at initialization).

## Running

From the build directory:

```bash
ctest -R '^adv-dyn\+gocart$' --output-on-failure
```

The data directory comes from `LOCAL_REGRESSION_DATA_DIR` at configure time
(otherwise data is synced from S3). If the case's data is missing,
`run_case.cmake` prints "Regression data not found ... skipping" and CTest
still reports the test as passed.
