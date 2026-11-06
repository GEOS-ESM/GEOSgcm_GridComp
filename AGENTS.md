# AGENTS.md

**Audience:** coding agents and developers working in the `GEOSgcm_GridComp` repository.
**Describes:** the MAPL2 line, branch `develop`.
**Verified against:** `GEOSgcm_GridComp` @ `45e1f86f` (2026-11-04).
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
| [`.../GEOSturbulence_GridComp/ARCHITECTURE.md`](GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSturbulence_GridComp/ARCHITECTURE.md) | Turbulence |
| [`.../GEOSgwd_GridComp/ARCHITECTURE.md`](GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSgwd_GridComp/ARCHITECTURE.md) | Gravity wave drag (the ACG example) |
| [`.../GEOSmoist_GridComp/ARCHITECTURE.md`](GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/ARCHITECTURE.md) | Moist processes |
| [`.../GEOSsurface_GridComp/ARCHITECTURE.md`](GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/ARCHITECTURE.md) | Surface (land, lake, landice, route, saltwater) |
| `README.md`, `CONTRIBUTING.md` | Project basics, contribution rules |

Not yet documented (only their connection to the Gcm is in `docs/ARCHITECTURE.md`
section 3): `GEOSogcm_GridComp` (ocean composite, sea ice), `GEOSwgcm_GridComp`,
`GEOSdataatm_GridComp`, `GEOSmkiau_GridComp` (`AIAU`, `ADFI`), `GEOSgigatraj_GridComp`.
Documents for the nested repositories live in those repositories (fixture paths in
section 6 of `docs/ARCHITECTURE.md`).

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

This repository does not build alone; build the GEOSgcm fixture (its `README.md`,
`parallel_build.csh`). Build structure: `docs/ARCHITECTURE.md` section 8. No unit-test target
for the composite components was found. Tests exist for `pyMoist`
(`GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/pyMoist/tests/`: `numerical`,
`scripts`, `translate_tests`). CI: `.github/workflows/` (`workflow.yml`, `ndsl-checks.yml`,
`enforce-labels.yml`, `push-to-develop.yml`) and `.circleci/config.yml`.
