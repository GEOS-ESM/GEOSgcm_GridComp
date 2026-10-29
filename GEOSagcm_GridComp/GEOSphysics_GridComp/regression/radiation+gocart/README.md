# `radiation+gocart` regression case

Shortwave (SOLAR) and longwave (IRRAD) radiation with GOCART2G (DU, SS, SU)
as their aerosol provider, at C12L181, one 20-minute step from a
2000-04-14T22:00 restart.

The radiation regressions (`solar-sa`, `irrad-sa`) take `AERO` from
`FakeGOCART`, which sets `implements_aerosol_optics_method = .false.`, so
radiation skips aerosol optics entirely. The GOCART2G regressions run without
radiation. This case is the only one in which a radiation component calls
GOCART2G's aerosol optics.

## How radiation and GOCART2G are coupled

The coupling goes one way: GOCART provides aerosol optical properties, and
Solar and IR radiation use them. GOCART builds AERO as an export
- nests each child instance's own aero state inside AERO (e.g. DU, SS, SU)
- adds the fields radiation needs (PLE, RH, EXT, SSA, ASY, and optional PFMOM)
- attaches methods with ESMF_MethodAdd (radiation uses `run_aerosol_optics`)
NOTE: the rank 3 fields, e.g. SO4, is added to the bundle as is. Rank 4 fields,
e.g. DU, SS, with an ungridded dim, are split into one 3D field per bin
(DU001...DU0005). The binned fields share memory with the 4D field.

Radiation does **not** depend on GOCART. It knows noting about dust bins or
sulfate. It only known the `AERO` contract: the info keys, the method label,
and the three output fields. Radiation
- read `implements_aerosol_optics_method` from `AERO`'s info. If the flag is
  false, the run is treated as aerosol-free (`NA = 0`)
- Radiation writes its own `RH` and `PLE` into fields named by the info keys
  - `relative_humidity_for_aerosol_optics` and
  - `air_pressure_for_aerosol_optics`
- For each band, radiation
  - sets `band_for_aerosol_optics = OFFSET + band`
  - calls `ESMF_MethodExecute(AERO, label=run_aerosol_optics)`

## Gridcomp hierarchy

```
Cap                      (cap_gridcomp.yaml; ExtData and History off)
└── Root                 (root.yaml; C12 cubed sphere, vertical grid from PLE)
    ├── DataForGOCART2G  data component: GOCART2G's 43 imports
    ├── GOCART2G         aerosol parent (Run1, Run2); AERO export
    │   ├── DU           dust (5 bins)
    │   ├── SS           sea salt (5 bins)
    │   └── SU           sulfate (DMS, SO2, SO4, MSA)
    ├── DataForSolar     data component: SOLAR's 23 imports
    ├── SOLAR            shortwave radiation; AERO import from GOCART2G
    ├── DataForIrrad     data component: IRRAD's 24 imports
    └── IRRAD            longwave radiation; AERO import from GOCART2G
```

GOCART2G adds DU, SS and SU at run time from `ACTIVE_INSTANCES` in
`GOCART2G_GridComp.yaml` (CA and NI are empty). Each child's `<name>_AERO`
state goes into GOCART2G's `AERO` export, which Root connects to both SOLAR's
and IRRAD's `AERO` imports.

## Run sequence

Root (`ConfigurableGridComp`) runs its children in `root.yaml` order. Each
child runs all of its `run_phases` in order (or the single phase `run` if none
are listed) before the next child starts. GOCART2G runs each phase on its
children in `ACTIVE_INSTANCES` order (DU, SS, SU), skipping any whose name
contains "data". ExtData and History are off.

Within each step:

1. `DataForGOCART2G` `run`: nothing to do; holds GOCART2G's 43 inputs from
   its restart.
2. `GOCART2G` `Run1`: `Run1` on DU, SS, SU.
3. `GOCART2G` `Run2`: `Run0` on DU, SS, SU; zero the total-aerosol
   diagnostics; `Run2` on DU, SS, SU; compute the totals (AOD, PM, ...).
4. `DataForSolar` `run`: nothing to do; holds SOLAR's inputs.
5. `SOLAR` `run`: shortwave radiative transfer. Once per band it calls
   GOCART2G `run_aerosol_optics`, which calls `aerosol_optics` on DU, SS and
   SU.
6. `DataForIrrad` `run`: nothing to do; holds IRRAD's inputs.
7. `IRRAD` `run`: longwave radiative transfer, with the same
   `run_aerosol_optics` calls per longwave band (using its band offset).

The case runs one step (22:00 to 22:20), after which the checkpoint is
written to `checkpoints/last`.

What this means:

- **Radiation sees this step's aerosols.** GOCART2G runs first, so SOLAR's
  and IRRAD's aerosol optics use the tracers after GOCART2G's emissions,
  chemistry and deposition for the step, not the restart values.
- **The optics calls happen inside radiation's run.** SOLAR and IRRAD fill
  the `AERO` state's `PLE` and `RH` and call GOCART2G's optics method
  directly; GOCART2G's own run phases are not involved.
- **GOCART2G does not depend on radiation.** Its results are bit-identical
  whether it runs before or after SOLAR and IRRAD; only the radiation output
  depends on the order.

## What is tested

- **The `AERO` handoff.** GOCART2G's `AERO` export (a state holding its
  children's `<name>_AERO` states) connects to the `AERO` imports of both
  SOLAR and IRRAD, and both see `implements_aerosol_optics_method = .true.`.
- **GOCART2G `run_aerosol_optics`**, called once per band by SOLAR
  (shortwave bands) and by IRRAD (longwave bands):
  - the caller fills the `AERO` state's `PLE` and `RH`;
  - GOCART2G copies them to each child, calls each child's `aerosol_optics`,
    and sums the children's extinction, single-scattering albedo and
    asymmetry parameter.
- **DU2G, SS2G and SU2G `aerosol_optics`**: Mie-table lookup per band in
  `opticsBands_{DU,SS,SU}.v2.1.0.RRTMGP.nc4`, which holds both the shortwave
  and the longwave bands (IRRAD passes its band offset), with the mixing
  ratios turned into layer masses from `PLE`.
- **Both radiative-transfer passes** in each component: with aerosols
  (`include_aerosols` and the provider's optics) and without (the clean-sky
  `...NA...` fields).
- **GOCART2G Run1 and Run2** running alongside radiation (emissions,
  chemistry, deposition, diagnostics), as in `gocart+seasalt+dust+sulfate`.

The comparison checks all 17 checkpoint states in `checkpoints/last`: SOLAR
and IRRAD import/internal/export, GOCART2G import/export, and DU, SS and SU
import/internal/export.

## What is not tested

- **Advection.** DU, SS and SU subscribe to the `TRADV` service, but nothing
  here provides it (no DYN/ADV), so the subscriptions stay unconnected and
  nothing is advected. See `GEOSagcm_GridComp/regression/adv-dyn+gocart`.
- **GOCART2G `get_monochromatic_aop`** and the children's
  `monochromatic_aerosol_optics`: neither radiation component calls them.
- **Photolysis optics** (`use_photolysis_table`): radiation-band tables only.
- **A MAPL2 baseline.** The MAPL2 radiation captures used for `solar-sa` and
  `irrad-sa` had no aerosol optics, so there is none for this configuration.
  The baseline in `checkpoints/last` is this case's own MAPL3 output; the
  case guards against changes, not against MAPL2.

## How the baseline was checked

When the baseline was created, each radiation component's output was compared
with its standalone case (same inputs and restart, no aerosol optics).

SOLAR, against `solar-sa`:

- identical: all SOLAR imports, the 26 internals that do not involve aerosols
  (including the clean-sky fluxes `FSWNAN`, `FSCNAN`), and 108 of 114 exports;
- different: the 15 internals that include aerosols (`FSWN`, `FSCN`,
  `FSWUN`, `FSCUN`, the band fluxes, and the surface direct/diffuse UV, PAR
  and NIR fluxes) by up to about 10%, and the 6 exports taken directly from
  them;
- physically consistent: the surface direct beam drops in every sunlit
  column and never rises, the diffuse flux mostly rises, and net surface
  shortwave falls in almost every sunlit column (it rises in under 2% of
  cells).

IRRAD, against `irrad-sa`:

- identical: all IRRAD imports, the no-aerosol internals (`FLA`, `FLAD`,
  `FLAU`, `FLXA`, `FLXAD`, `FLXAU`, `DFDTSNA`, `DFDTSCNA`) and `TS`, and 54
  of 57 exports;
- different: the 25 internals that include aerosols (`FLX`, `FLC` and their
  up/down components, `DFDTS`, `DFDTSC`, the per-band OLR and its
  derivatives, `SFCEM`) by up to about 2%, and the exports `DFDTS`, `LWS0` and
  `SFCEM0`;
- physically consistent: the surface clear-sky downward longwave `LWS0` rises
  in 95% of cells and never falls.

GOCART2G's output (and DU's, SS's and SU's) is bit-identical to the earlier
SOLAR-only version of this case, in which radiation ran first.

## Components and inputs

- `DataForSolar` (`data-for-solar.yaml`): SOLAR's 23 imports, restarted from
  `solar-sa`'s input state; `solar.yaml` is `solar-sa`'s SOLAR configuration.
- `DataForIrrad` (`data-for-irrad.yaml`): IRRAD's 24 imports, restarted from
  `irrad-sa`'s input state; `irrad.yaml` is `irrad-sa`'s IRRAD configuration.
- `DataForGOCART2G` (`data-for-gocart2g.yaml`): GOCART2G's 43 imports (DU, SS
  and SU, GMI omitted), restarted from the MAPL2 GOCART2G import state.
- `GOCART2G` (`GOCART2G_GridComp.yaml`, `{DU,SS,SU}_instance_*.yaml`): DU, SS
  and SU active; CA and NI empty (not ported). Run phases `Run1`, `Run2`.
- `root.yaml`: the 23 SOLAR, 24 IRRAD and 43 GOCART2G input connections, plus
  `GOCART2G/AERO -> SOLAR/AERO` and `GOCART2G/AERO -> IRRAD/AERO`.

Regression data, under `GEOSagcm_GridComp/GEOSphysics_GridComp/radiation+gocart`:

- `checkpoints/2000-04-14T22:00:00/`: `DataForSolar`, `SOLAR`,
  `DataForIrrad`, `IRRAD`, `DataForGOCART2G`, `DU`, `SS` and `SU` restarts;
- `checkpoints/last/`: the baseline.

`ExtData` is not duplicated: `run_case.cmake` links `ExtData/chemistry` from
`GEOSchem_GridComp/GOCART2G_GridComp/ExtData` (aerosol optics tables, CARN
volcanic degassing) and `ExtData/g5gcm` from `GEOSradiation_GridComp/ExtData`
(RRTMGP and solar tables).

## Running

From the build directory:

```bash
ctest -R '^radiation\+gocart$' --output-on-failure
```

The data directory comes from `LOCAL_REGRESSION_DATA_DIR` at configure time
(otherwise data is synced from S3). If the case's data is missing,
`run_case.cmake` prints "Regression data not found ... skipping" and CTest
still reports the test as passed.
