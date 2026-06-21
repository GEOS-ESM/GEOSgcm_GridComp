# MGB2_2M microphysics nondeterminism forensic report

## Executive summary

This investigation mapped `CLDMICR_OPTION == "MGB2_2M"` from `GEOS_MoistGridComp` into `GEOS_MGB2_2M_InterfaceMod` and downstream microphysics (`micro_mg3_0.F90`, `micro_mg_utils.F90`, `aer_cloud.F90`). The strongest evidence for a layout-sensitive/nondeterministic failure is inside `GEOS_MGB2_2M_InterfaceMod.F90`:

* **P0-likely:** contact/soot aerosol inputs to MG are only partially initialized by column and level. For levels skipped by `kbmin`, aerosol bins with no dust/soot, and columns with `naux == 0`, the code can pass stale values or divide by zero into MG. This is consistent with same-executable nondeterminism and layout-regression failures.
* **P0-likely:** `rnsootr8` and `nsootr8` are divided by `naux` even when no soot mode is present. This can create Inf/NaN in optional/contact-nucleation inputs; later microphysics can feed large tendencies into water vapor.
* **P1-likely:** activation diagnostics/export arrays such as `SMAX_*`, `CCN*`, aerosol component diagnostics, and some masks are not fully reset for every level before use. The activation loop intentionally excludes `K=LM`, so the top/bottom edge value can remain from MAPL allocation or a prior fill.
* **P1-risk:** vapor is protected only after MG by `MAX(RAD_QV + qvlat*dt, 0)`, but this clipping does not conserve total water and can hide the source of a bad sink. If a bad `qvlatr8` is created, the state may survive while diagnostics/tendencies become inconsistent with later checks.
* **P1-risk:** convective-number coupling imports/export tendencies (`DQLDT_*`, `DQIDT_*`, `DQLDTTRB`, `DQIDTTRB`, `CNV_MFD`) and converts mass tendencies to number without explicit unit assertions. Number tendencies can alter activation and autoconversion indirectly, even though they do not directly subtract vapor.

No OpenMP directives were found in the inspected MGB2_2M interface, so the leading layout-sensitivity mechanism is uninitialized/stale local array content and edge/column-dependent optional inputs rather than an obvious shared-memory race.

## Call-tree and data-flow map

1. `GEOS_MoistGridComp` reads `CLDMICR_OPTION` and accepts `MGB2_2M` as a valid cloud microphysics option.
2. Setup dispatch calls `MGB2_2M_Setup`, and the MGB2 setup must be first internal-spec provider because binary restarts require `Q` to be first.
3. During initialize, `GEOS_MoistGridComp` skips generic `aer_cloud_init` for MGB2_2M because `MGB2_2M_Initialize` owns aerosol-cloud initialization for this option.
4. During run, after convection/shallow/turbulence have filled their exports/imports, `GEOS_MoistGridComp` calls `MGB2_2M_Run`.
5. `MGB2_2M_Run` reads internal state (`Q`, `QLLS`, `QLCN`, `QILS`, `QICN`, hydrometeor and number fields), imports meteorology/turbulence/aerosol support fields, and exports radiation, precipitation, macro/micro tendencies, activation diagnostics, number diagnostics, radar diagnostics, and `Q*TOT` fields.
6. MGB2_2M sequence:
   * derive pressure, mass, density, saturation, potential temperature;
   * `fix_up_clouds_2M` pre-cleans moisture/number/hydrometeor state;
   * compute subgrid vertical velocity (`WSUB_OPTION`);
   * activate aerosol and ice nuclei (`aerosol_activate`);
   * convert convective/shallow/turbulent condensate tendencies into droplet/ice number tendencies;
   * run macrophysics/PDF (`hystpdf`, `EVAP3`, `SUBL3`, optional `update_cld` and second `hystpdf`);
   * call MG (`mmicro_pcond` for old MG, `micro_mg_tend_interface` for MG2/3);
   * update `Q`, condensate, precipitation, number, and radiation-facing in-cloud fields;
   * diagnose precipitation/radar/cloud/radiation outputs.

Additional paths needed beyond the named subtree to complete a runtime investigation are:

* `GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/micro_mg3_0.F90` and `micro_mg_utils.F90` for MG tendency semantics and optional argument handling.
* `GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/aer_cloud.F90` for `aer_cloud_init` and `aerosol_activate` behavior.
* `GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/ConvPar_GF2020.F90`, `GEOS_RAS_InterfaceMod.F90`, and shallow/turbulence providers for `DQLDT_*`, `DQIDT_*`, `CNV_MFD`, `SHLW_*`, `DQLDTTRB`, and `DQIDTTRB` producers.
* Parent MAPL/history/restart resource files that choose `CLDMICR_OPTION`, `MGVERSION`, `WSUB_OPTION`, `USE_NCLOUD_CLIM`, `ND_CST`, `NI_CST`, `NG_CST`, `SECOND_HYSTPDF`, `DO_UPD_CLD`, and aerosol options.

## Ranked suspects

### P0: stale/uninitialized contact aerosol arrays and zero-soot divide

Evidence: the code initializes `rndstr8` and `naconr8` globally, then reinitializes them at column entry, but only fills dust bins for `K=kbmin:LM` when a mode has positive dust fraction. Soot accumulators are zeroed per level, but divided by `naux` unconditionally after the soot loop. If no soot modes are present, `rnsootr8` and `nsootr8` divide by zero. Levels outside the `kbmin:LM` loop keep prior/default values; bins not touched retain previous column values unless fully zeroed immediately before MG.

Why this can fail: optional/contact nucleation inputs (`rndstr8`, `naconr8`, `nsootr8`, `rnsootr8`, `nbincontactdust`) feed MG. Inf/NaN or stale aerosol content can produce layout-dependent ice nucleation, latent heating, and vapor sink/source tendencies. Because local allocation and column traversal vary with MPI decomposition and compiler/runtime memory state, this exactly matches nondeterministic/layout behavior.

Minimal patch proposal:

* Before each column MG call, set `rndstr8=2.0e-7`, `naconr8=0.0`, `nsootr8=0.0`, `rnsootr8=0.0`, and `nbincontactdust=0` for the full `1:LM` column.
* Replace `rnsootr8/naux` and `nsootr8/naux` with guarded division only when `naux > 0`; otherwise leave both zero.
* Do not allow `nbincontactdust` to count stale bins from earlier columns.

Long-term fix: build a small column-preparation subroutine with explicit `intent(out)` arrays and unit tests for no-dust/no-soot columns.

Diagnostic confirmation: in a debug build, assert `all(isfinite(rndstr8,naconr8,nsootr8,rnsootr8))`, `nbincontactdust >= 0`, and no zero-denominator soot averaging immediately before MG.

Regression: same initial condition with two MPI layouts and `-finit-real=snan`/bounds/FPE enabled should fail before the patch and match after the patch.

### P0/P1: edge levels in activation diagnostics not assigned every call

Evidence: activation loops only `K=1, LM-1`, while many export arrays are allocated and then later used/diagnosed as full `1:LM` arrays. `SC_ICE` is reset, but arrays such as `SMAX_LIQ`, `SMAX_ICE`, `CCN01`, `CCN02`, `CCN1`, `CDNC_NUC`, aerosol component diagnostics, `INC_NUC`, `NHET_*`, and `DUST_*` need full-field deterministic defaults before the partial loop.

Why this can fail: layout tests often expose uninitialized edges. Even if the surface/top level is rarely active, stale diagnostics can affect `INC_NUC = INC_NUC*CLLS`, `frzdepr8 = NHET_DEP/DT`, and number tendency inputs when clouds exist at the excluded level.

Minimal patch proposal: after all activation export pointers are acquired and before the activation branch, initialize full arrays to safe defaults (`0.0`, except `SC_ICE=1.0`). Then keep the existing loop overwrites.

Long-term fix: separate diagnostics from prognostic drivers; only pass arrays explicitly filled on active levels to MG.

Diagnostic confirmation: run with signaling NaNs for allocated reals and assert all activation arrays are finite before `INC_NUC` and before MG.

### P1: non-conservative post-MG vapor clipping

Evidence: after MG, vapor is updated with `RAD_QV = MAX(RAD_QV + qvlatr8*dt, 0.0)`. This prevents a negative local state but does not reduce the corresponding condensate/precipitation source or tendency diagnostics to conserve total water.

Why this can fail: if a bad sink occurs, the clip masks the immediate cause while total-water diagnostics and later physics see inconsistent moisture. A downstream negative-water-vapor crash may then appear intermittent depending on layout and ordering.

Minimal patch proposal: add an optional debug assertion/log before clipping whenever `RAD_QV + qvlatr8*dt < -eps`, including column, level, old qv, qvlatr8, condensate, number, cloud fraction, and active process options. Do not change production behavior until the source is confirmed.

Long-term fix: implement a conservative limiter that scales vapor sinks and associated condensate/latent tendencies consistently.

### P1: convective/shallow/turbulence number detrainment unit and association hazards

Evidence: `DQLDT_DC`, `DQIDT_DC`, `DQLDT_SC`, `DQIDT_SC`, `DQLDTTRB`, and `DQIDTTRB` are read opportunistically and converted to number tendencies using `AIRDEN`, fixed size assumptions, and `CFX`; `CNV_MFD` is used to diagnose convective number if associated.

Why this can fail: these fields may come from different providers and may be per-second tendencies, per-timestep increments, grid-mean or in-cloud values. A unit mismatch can overproduce number, altering autoconversion/accretion and vapor adjustment. Missing fields are tolerated by pointer association, so configurations can silently take different branches.

Minimal patch proposal: add debug-only unit/range checks and explicit comments/assertions near each conversion. Require `NotFoundOk` only where absence is intentional, and log which convective fields are associated for `MGB2_2M`.

Long-term fix: formalize a two-moment convective microphysics interface with units encoded in field metadata and a single conversion routine.

### P2: cloud partition recovery and independent clipping can violate total water

Evidence: `FQA` partitions combined cloud liquid/ice back into convective/large-scale species after MG; later `fix_up_clouds_2M`, `FILLQ2ZERO`, and radiation caps limit fields independently. Radiation fields are in-cloud values, while internal prognostics are grid-mean values, increasing the chance of diagnostic/prognostic confusion.

Why this can fail: independent clipping is usually deterministic but can amplify small layout differences into water-budget differences. It is less likely to be the primary nondeterministic trigger than uninitialized optional inputs.

Minimal patch proposal: add total-water before/after diagnostics around macro, MG, and fix-up blocks; do not alter behavior until budget deltas point to a specific limiter.

## Negative-water-vapor budget trace

Potential vapor removers/adders in MGB2_2M are:

* `hystpdf`, `EVAP3`, and `SUBL3` mutate `Q`, liquid/ice condensate, cloud fraction, and number during macrophysics.
* MG receives `qvr8`, `qcr8`, `qir8`, hydrometeor and number fields, aerosol/contact inputs, and returns `qvlatr8`, condensate tendencies, precipitation tendencies, and latent heating tendencies.
* `update_cld` and optional second `hystpdf` can further adjust `Q`/cloud condensate after MG partitioning.
* Final vapor state is assigned from `RAD_QV`, which was produced by old `Q` plus `qvlatr8 * dt`, clipped at zero.

The first diagnostic to add is a per-column water-budget check:

```fortran
water0 = Q + QLLS + QLCN + QILS + QICN + QRAIN + QSNOW + QGRAUPEL
! after macro/MG/fix-up, compare water1-water0 plus precip flux divergence
```

Print only the first failing column per rank when `min(Q) < q_floor`, `RAD_QV + qvlatr8*dt < -eps`, or `abs(total_water_residual) > tolerance`.

## Layout-regression and nondeterminism risks

* Uninitialized/stale local allocatables and export arrays are the leading risk.
* Edge-level activation omissions (`1:LM-1`) can vary with vertical clouds and memory initialization.
* Optional pointer association changes behavior silently when provider fields are absent/present.
* No explicit OpenMP race was identified in the interface; if MG internals are threaded, inspect `SAVE`/module state and any static lookup/cache arrays in `micro_mg3_0.F90` and utilities.
* Floating-point layout differences can be magnified by discontinuous thresholds (`MAX`, `MIN`, cloud fraction thresholds, autoconversion thresholds), but this is a secondary amplifier unless an upstream uninitialized or unit error seeds large differences.

## Concrete patches, small and reviewable

1. Deterministic MG column input initialization and zero-denominator guard for soot/contact aerosol inputs.
2. Full-field deterministic defaults for activation diagnostics before the partial activation loop.
3. Debug-only vapor-sink and total-water budget assertions around MG, disabled by default or controlled by existing `DEBUG_MST`.
4. Debug-only association/unit report for convective/turbulent condensate tendency fields when `CLDMICR_OPTION == "MGB2_2M"`.

## Diagnostics to add before patching

* Before MG call: finite checks for all MG inputs, especially `rndstr8`, `naconr8`, `nsootr8`, `rnsootr8`, `qvr8`, `qcr8`, `qir8`, `ncr8`, `nir8`, `cldfr8`, and `pdelr8`.
* Before clipping `RAD_QV`: print negative candidate values and all source tendencies.
* Around convective-number update: print min/max of associated `DQLDT_*`, `DQIDT_*`, `DQLDTTRB`, `DQIDTTRB`, `DNDCNV`, `DNICNV`, `CNV_MFD`, and units assumptions.
* Around macro/MG/fix-up: total water and total number conservation summaries.

## Tests to run

* Same layout repeated run: at least 5 repeats for one simulated day with identical executable and restarts; compare bitwise histories/restarts.
* Different MPI layouts: run the known layout-regression pair for at least one simulated day and compare moisture/number/history checksums.
* Restart reproducibility: stop/restart at a short interval and compare with continuous run.
* Debug/FPE/bounds/uninitialized build: compile with bounds, FPE traps, and real initialization to signaling NaN or a sentinel value.
* Conservation checks: add history or asserts for column-integrated `qv + ql + qi + qr + qs + qg` plus precipitation flux divergence across macro, MG, and cleanup.
