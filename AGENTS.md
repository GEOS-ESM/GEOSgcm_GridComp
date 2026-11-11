# GEOSgcm_GridComp: guide for agents

**Describes:** the MAPL3 line, branch `release/MAPL-v3`; verified against `GEOSgcm_GridComp` @ `b1825cc1` (2026-11-10).
**Placement:** `src/Components/@GEOSgcm_GridComp` in the GEOSgcm fixture (`components.yaml`, key `GEOSgcm_GridComp`); other fixtures may place it elsewhere. Paths below are relative to this repository.

## What this is

The component hierarchy of the GEOS general circulation model: the root composite `GEOS_GcmGridComp.F90` (atmosphere, ocean, waves, IAU/replay, digital filter) and its children. Several subdirectories are separate nested repositories (Chemistry, Radiation, FV dynamical core, Ocean, ...); edit them only in their own repository.

> **Build status on this line:** partly built; the full model is not built. On `release/MAPL-v3`, `CMakeLists.txt` adds only `GEOSagcm_GridComp`; the Ogcm, Wgcm, DataAtm and mkiau composites are commented out ("not yet ported to MAPL3") and the root library `GEOS_GcmGridComp.F90` is bypassed. The Agcm library `GEOS_AgcmGridComp.F90` and the Physics library `GEOS_PhysicsGridComp.F90` are bypassed too; Surface and Turbulence are not built. In the GEOSgcm fixture, built: Superdyn (with FVdycoreCubed), Moist, GWD, Chemistry (GOCART only), Radiation. The fixture's sparse file excludes the Held-Suarez files; a checkout that includes them builds the Held-Suarez library and the `hs` regression. The Gcm, Agcm and Physics sources are MAPL2-style code kept for a later port (`CMakeLists.txt` at each level).

## Documents

| Document | Content |
|---|---|
| `GEOSagcm_GridComp/AGENTS.md`, `GEOSagcm_GridComp/docs/ARCHITECTURE.md` | The atmosphere subtree (Agcm composite, Superdyn, Physics with GWD, Moist, Turbulence): status, traps, regression cases; map, run orders, services, contracts |
| `README.md`, `CONTRIBUTING.md` | Project basics, contribution rules |

Nested repositories (separate Git repositories, each with its own `AGENTS.md`; after `mepo clone`): `GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSchem_GridComp` (Chemistry, with GOCART and other nested repositories), `GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSradiation_GridComp` (Radiation), `GEOSagcm_GridComp/GEOSsuperdyn_GridComp/@FVdycoreCubed_GridComp` (FV3 dynamical core), `GEOSogcm_GridComp/@GEOS_OceanGridComp` (Ocean); also `cice6`, `sis2`, `ww3`, `umwm`, `GigaTraj`, `geos_state_bias` and `Topo` (`components.yaml` in the fixture).

## Rules

- Note bugs you see; do not fix them unless asked. Record a bug or an open question as a GitHub issue, not in a document; a document may link to the issue.
- Keep this file and the porting notes current: if you change a child list, run phase, connection, resource key or build option they describe, update them in the same change. Cite files and routines, not line numbers.
- Do not edit the generated `*___.h` files; edit the `*_StateSpecs.rc` table.
- Code owners are in `.github/CODEOWNERS`; `CMakeLists.txt` files belong to the CMake team.
- Error handling in code you write or touch uses the MAPL3 macros (`@MAPL/include/MAPL_ErrLog.h`): `_RC` on calls, `_STAT` on `allocate`, `_ASSERT(cond, msg)` with a message, and `_RETURN(_SUCCESS)` at every successful exit. Convert `RC=STATUS`/`VERIFY_`/`ASSERT_`/`RETURN_` in any routine you change. `_RC` expands to two statements (`rc=status);_VERIFY(status`), so it cannot end a one-line `if`; use a block `if`. Do not add decorative `!----` separator lines.

## Porting to MAPL3

- A component under port keeps `mapl3-porting-notes.md` beside its code (for example `GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/mapl3-porting-notes.md`). Read it before changing the component; record each decision there with its reason.
- Choose the option closest to the MAPL2 code. A ported component must match a MAPL2 capture bit for bit.
- Keep MAPL2 bugs as they are while porting, for bit-for-bit agreement: file an issue and leave a `TODO` that cites it (for example `FILLQ2ZERO` in `Process_Library.F90`, #1570).
- In a standalone case, activate only the exports the component writes, through a History collection that never writes a file; do not give exports fill values (`GEOSphysics_GridComp/regression/README.md`).
- MAPL3 `StateGetPointer` returns 1-based pointers for edge fields. Most MAPL2 code indexes edge fields from 0 (`PLE(:,:,0:LM)`), but not all of it (Moist sums `UMF_DC`/`UMF_SC` 1-based): check each edge-field site, and remap the bounds (through a `contiguous` pointer) only where the MAPL2 code is 0-based.
- A child inherits its parent's configuration keys (the child's value wins), so set a shared resource once, in the parent.
- Fields in a MAPL3 service bundle (`TRADV`, `MTR`) do not keep their short names; map them back through `standard_name`, as FVdycoreCubed's `get_short_name` (`fv_shared.F90`) does.
- Entry points must match `I_SetServices`/`I_Run` (`@MAPL/infrastructure/esmf/ESMF_Interfaces.F90`): no `intent` on the ESMF arguments, and `integer, intent(out) :: rc`, not `optional`. A mismatch can compile through `MAPL_GridCompSetEntryPoint` but fails once the component is added as a child.
- State specs go in `<Comp>_StateSpecs.rc`, generated by ACG (`@MAPL/apps/MAPL_GridCompSpecs_ACG.py`), which silently drops columns it does not know. MAPL2's `DEFAULT=x` must become `FILL=x` (MAPL3 does not zero-fill internals); `FRIENDLYTO=<self>` on an internal becomes `ADD2EXPORT`; `FRIENDLYTO=<other component>` becomes a service subscription. Fields whose names are known only at run time stay as `MAPL_GridCompAddSpec` calls in `SetServices`.
- `MAPL_UNDEF` is `-huge` in MAPL3 and `1.0e15` in MAPL2. Equality tests against it fail on inputs converted from MAPL2 restarts; check every comparison with `MAPL_UNDEF` (GOCART uses a local `MAPL2_UNDEF = 1.0e15`).
- MAPL3 has no `alloc=.true.`: an export nobody requests is not allocated. MAPL2 code that uses exports as work arrays, or sums other exports by name, silently changes results (`GEOSmoist_GridComp/mapl3-alloc-exports.md`).
- Each component has its own clock, and MAPL3 calls the user `Run` only when the component's run alarm rings. An alarm the component creates itself must reproduce MAPL2's first-step behaviour (IRRAD's first-step alarm fix, Solar step 18 in `@GEOSradiation_GridComp`).
- API replacements: `MAPL_MetaComp`, `MAPL_Get(MAPL, ...)` and `MAPL_GetResource(MAPL, ...)` become `MAPL_GridCompGet`, `MAPL_GridCompGetResource` (labels without the trailing colon) and `MAPL_GridCompTimerStart`/`Stop`; `MAPL_GenericSetServices`, `MAPL_GenericInitialize` and `MAPL_GenericRunChildren` are not called; `ESMF_Attribute*` becomes `ESMF_Info*`, with the `...Alloc` getters for allocatable results; private state uses `_SET_NAMED_PRIVATE_STATE`/`_GET_NAMED_PRIVATE_STATE`, which must not be mixed with `ESMF_UserCompSetInternalState`.
- A parent that read a child's export in MAPL2 declares the field as its own import and connects it with `MAPL_GridCompAddConnection(gc, src_comp=<child>, src_names=..., dst_comp="<self>")`.
- Run phases: a configurable Root (the `*-sa` cases) runs phase `"run"` by name (`@MAPL/gridcomps/configurable/ConfigurableGridComp.F90`), so register `phase_name="run"`; without `phase_name` a phase is named `GENERIC::RUN_USER`. `MAPL_GridCompAddChild` always needs an `hconfig`; for a child without configuration, pass an empty one (`ESMF_HConfigCreate(content='{}')`).
- A MAPL2 capture used as a baseline must match the standalone configuration exactly: the same providers (for example `AERO_PROVIDER: none` for Radiation), the same requested exports, and an internal restart from the same capture.

## Build and test

This repository does not build alone; build the GEOSgcm fixture (`README.md`, `parallel_build.csh`). CI: `.github/workflows/`, `.circleci/config.yml` (runs `ctest -L REGRESSION`).

Regression tests (label `REGRESSION`, one test per line of a `regression/cases.txt`):
- Physics: the Moist cases in `GEOSphysics_GridComp/regression` and `gwd-ncar-sa` in `GEOSgwd_GridComp/regression` (`GEOSagcm_GridComp/AGENTS.md`).
- `GEOShs_GridComp/regression`: `held-suarez`, only when the Held-Suarez files are in the checkout.

About the regression cases:
- Each case directory has a `logging.yaml`, loaded by `pflogger: logging.yaml` in `mapl.yaml`. Keep it in new cases: without it, MAPL's field-dictionary warnings are written under an MPI lock and add seconds to each test.
- Data is read from `$LOCAL_REGRESSION_DATA_DIR/<path>`, set at configure time (`@cmake/esma_support/esma_regression_tests.cmake`). If a case's data directory is missing, `run_case.cmake` prints "Regression data not found ... skipping" and the test passes without running.
- Results are compared exactly with `nccmp`, and only for IntelLLVM (ifx) builds; other compilers run the case without comparing.

Pull requests for this line target `release/MAPL-v3`, not `develop` (the MAPL2 line), labelled `0 diff` and `MAPL3`. This repository has no `CHANGELOG.md`.
