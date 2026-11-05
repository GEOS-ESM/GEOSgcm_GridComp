# AGENTS.md

**Audience:** coding agents and developers working in the `GEOSgcm_GridComp` repository.
**Describes:** the MAPL2 line, branch `develop`.
**Verified against:** `GEOSgcm_GridComp` @ `c7d81d7c` (2026-11-05; source files identical to `45e1f86f`, the earlier base, which only differs by the added documents).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `0f2e0fb`).
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

## Documents

Read `docs/ARCHITECTURE.md` first, then the file for the directory you touch.

| Document | Content |
|---|---|
| [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) | Placement, directory map, the root Gcm in depth (children, run, connections, IAU/replay), Physics map, nested repositories, MAPL2 patterns, build |
| [`docs/decisions.md`](docs/decisions.md) | Previous architectural decisions (pattern, likely reason) |
| [`GEOSagcm_GridComp/ARCHITECTURE.md`](GEOSagcm_GridComp/ARCHITECTURE.md) | Agcm and Superdyn: run order, IAU/bias/replay in the Agcm |
| [`GEOSagcm_GridComp/GEOSphysics_GridComp/ARCHITECTURE.md`](GEOSagcm_GridComp/GEOSphysics_GridComp/ARCHITECTURE.md) | Physics composite |
| Physics children: [turbulence](GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSturbulence_GridComp/ARCHITECTURE.md), [gravity wave drag](GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSgwd_GridComp/ARCHITECTURE.md), [moist](GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/ARCHITECTURE.md), [surface](GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/ARCHITECTURE.md) | Each child's component document |
| `README.md`, `CONTRIBUTING.md` | Project basics, contribution rules |

Not yet documented: Ogcm composite, Wgcm, DataAtm, mkiau, Gigatraj (see `docs/ARCHITECTURE.md` section 3).
Nested repositories (after `mepo clone`; paths in section 6 of `docs/ARCHITECTURE.md`):
[Chemistry](GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSchem_GridComp/AGENTS.md),
[GOCART](GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSchem_GridComp/@GOCART/docs/ARCHITECTURE.md),
[Radiation](GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSradiation_GridComp/AGENTS.md),
[FVdycoreCubed](GEOSagcm_GridComp/GEOSsuperdyn_GridComp/@FVdycoreCubed_GridComp/AGENTS.md),
[Ocean](GEOSogcm_GridComp/@GEOS_OceanGridComp/AGENTS.md).

## Rules

- Note bugs you see; do not fix them unless asked. Record them under "Observations" in the
  relevant document.
- Keep the documents current: if you change a child list, run phase, connection, resource key
  or build option described there, update the document in the same change. Cite files and
  routines, not line numbers.
- Do not edit generated files (`*_Import___.h`, `*_Export___.h`, ...); edit the
  `*_StateSpecs.rc` file.
- Code owners are in `.github/CODEOWNERS`; `CMakeLists.txt` files belong to the CMake team.

## Build and test

This repository does not build alone; build the GEOSgcm fixture (`README.md`, `parallel_build.csh`);
build structure: `docs/ARCHITECTURE.md` section 8. No unit-test target for the composites was found;
`pyMoist` has tests (`GEOSmoist_GridComp/pyMoist/tests/`). CI: `.github/workflows/`, `.circleci/config.yml`.
