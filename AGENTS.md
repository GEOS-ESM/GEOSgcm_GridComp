# AGENTS.md

**Audience:** coding agents and developers working in the `GEOSgcm_GridComp` repository.
**Describes:** the MAPL3 line, branch `release/MAPL-v3`.
**Verified against:** `GEOSgcm_GridComp` @ `84f0aae3` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `f50622d`).
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

## What this is

The component hierarchy of the GEOS general circulation model: the root
composite `GEOS_GcmGridComp.F90` (atmosphere, ocean, waves, IAU/replay,
digital filter) and its children. In the GEOSgcm fixture it is at
`src/Components/@GEOSgcm_GridComp` (`components.yaml`, key
`GEOSgcm_GridComp`). Several subdirectories are separate nested repositories
(Chemistry, Radiation, FV dynamical core, Ocean, ...); edit them only in
their own repository.

> **Build status on this line:** partly built; the full model is not built.
> On `release/MAPL-v3`, `CMakeLists.txt` adds only `GEOSagcm_GridComp`; the
> Ogcm, Wgcm, DataAtm and mkiau composites are commented out ("not yet ported
> to MAPL3") and the root library `GEOS_GcmGridComp.F90` is bypassed. The
> Agcm library `GEOS_AgcmGridComp.F90` and the Physics library
> `GEOS_PhysicsGridComp.F90` are bypassed too; Moist, Surface and Turbulence
> are not built. In the GEOSgcm fixture, built: Superdyn (with
> FVdycoreCubed), GWD, Chemistry (GOCART only), Radiation. The fixture's sparse
> file excludes the Held-Suarez files; a checkout that includes them builds the
> Held-Suarez library and the `hs` regression. The Gcm, Agcm and Physics
> sources are MAPL2-style code kept for a later port. Details:
> `docs/ARCHITECTURE.md` section 2.1.

## Documents

Read `docs/ARCHITECTURE.md` first, then the file for the directory you touch.

| Document | Content |
|---|---|
| [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) | Build status, placement, directory map, the root Gcm (MAPL2 code, not compiled), Physics map, nested repositories, patterns, build |
| [`docs/decisions.md`](docs/decisions.md) | Previous architectural decisions (pattern, likely reason) |
| [`GEOSagcm_GridComp/ARCHITECTURE.md`](GEOSagcm_GridComp/ARCHITECTURE.md) | Agcm (MAPL2 code, not compiled) and Superdyn (ported, built): run order, IAU/bias/replay in the Agcm |
| [`GEOSagcm_GridComp/GEOSphysics_GridComp/ARCHITECTURE.md`](GEOSagcm_GridComp/GEOSphysics_GridComp/ARCHITECTURE.md) | Physics composite (MAPL2 code, not compiled) |
| Physics children (only GWD is built on this line): [turbulence](GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSturbulence_GridComp/ARCHITECTURE.md), [gravity wave drag](GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSgwd_GridComp/ARCHITECTURE.md), [moist](GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/ARCHITECTURE.md), [surface](GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/ARCHITECTURE.md) | Each child's component document |
| `README.md`, `CONTRIBUTING.md` | Project basics, contribution rules |

Not yet documented: Ogcm composite, Wgcm, DataAtm, mkiau, Gigatraj (see `docs/ARCHITECTURE.md` section 3).
Nested repositories (after `mepo clone`; paths in section 6 of `docs/ARCHITECTURE.md`) have their own
`AGENTS.md` or `docs/ARCHITECTURE.md`: Chemistry, GOCART, Radiation, FVdycoreCubed, Ocean.

## Rules

- Note bugs you see; do not fix them unless asked. Record them under "Observations" in the
  relevant document.
- Keep the documents current: if you change a child list, run phase, connection, resource key
  or build option described there, update the document in the same change. Cite files and
  routines, not line numbers.
- Do not edit generated files (`*_Import___.h`, `*_Export___.h`, ...); edit the
  `*_StateSpecs.rc` file.
- Code owners are in `.github/CODEOWNERS`; `CMakeLists.txt` files belong to the CMake team.
- Error handling in code you write or touch uses the MAPL3 macros (`@MAPL/include/MAPL_ErrLog.h`): `_RC` on calls, `_STAT` on `allocate`, `_ASSERT(cond, msg)` with a message, and `_RETURN(_SUCCESS)` at every successful exit. Convert `RC=STATUS`/`VERIFY_`/`ASSERT_`/`RETURN_` in any routine you change. `_RC` expands to two statements (`rc=status);_VERIFY(status`), so it cannot end a one-line `if`; use a block `if`. Do not add decorative `!----` separator lines.

## Build and test

This repository does not build alone; build the GEOSgcm fixture (`README.md`, `parallel_build.csh`);
build structure: `docs/ARCHITECTURE.md` section 8. Regression tests (label `REGRESSION`): `gwd`
(`gwd-ncar-sa`); `hs` (`held-suarez`) exists only when the Held-Suarez files are in the checkout (the fixture's
sparse file excludes them). `pyMoist` has tests (`GEOSmoist_GridComp/pyMoist/tests/`). CI:
`.github/workflows/`, `.circleci/config.yml` (runs `ctest -L REGRESSION`).
