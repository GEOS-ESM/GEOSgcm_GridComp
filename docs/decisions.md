# GEOSgcm_GridComp: previous architectural decisions

**Audience:** coding agents and developers working in the `GEOSgcm_GridComp` repository.
**Describes:** the MAPL3 line, branch `release/MAPL-v3`.
**Verified against:** `GEOSgcm_GridComp` @ `84f0aae3` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `f50622d`).
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

Only what is easy to see in the code now. Entries grow over time. Context:
[`ARCHITECTURE.md`](ARCHITECTURE.md). The entries below describe the MAPL2-style
composites (Gcm, Agcm, Physics), which are not compiled on this line; the last
entry is specific to this line.

### Parents run their children explicitly
- Pattern: the Gcm, Agcm, Physics and Superdyn composites call `ESMF_GridCompRun` on each
  child themselves, choosing order and `PHASE` (`GEOS_GcmGridComp.F90`, `Run`;
  `GEOSagcm_GridComp/GEOS_AgcmGridComp.F90`, `Run`;
  `GEOSagcm_GridComp/GEOSphysics_GridComp/GEOS_PhysicsGridComp.F90`, `Run`;
  `GEOSagcm_GridComp/GEOSsuperdyn_GridComp/GEOS_SuperdynGridComp.F90`, `Run`). Coupler runs
  (`MAPL_GenericRunCouplers`) are made explicitly in the Physics `Run` after most children, in
  the Agcm only for the dynamics child, in the Gcm only inside the replay loop (for the Agcm),
  and not in Superdyn. The Ogcm uses `MAPL_GenericRunChildren` in its ordinary mode and explicit
  calls only for dual-ocean (`GEOSogcm_GridComp/GEOS_OgcmGridComp.F90`, `Run`).  [Verified]
- Likely reason: the parent must insert its own work between children (analysis increments
  between dynamics phases, synchronization of winds and temperature between physics stages).  [Inferred]
- Question for maintainers: is the order of the Physics stages (surface 1, chemistry 1,
  turbulence 1, surface 2, ...) meant to be fixed, or should it be data-driven?

### Two-phase children for coupled stages
- Pattern: a child that must exchange data with a sibling mid-step registers two Run phases
  by calling `MAPL_GridCompSetEntryPoint(..., ESMF_METHOD_RUN, ...)` twice (Superdyn
  `Run`/`RunAddIncs`, AIAU `Run`/`Clear`, Gigatraj `GetInitVars`/`Run`); the parent passes
  `PHASE=n` (`GEOSagcm_GridComp/GEOSsuperdyn_GridComp/GEOS_SuperdynGridComp.F90`,
  `SetServices`; `GEOSmkiau_GridComp/GEOS_mkiauGridComp.F90`, `SetServices`).  [Verified]
- Likely reason: dynamics needs the physics tendencies between its two halves; AIAU needs a
  "make increments" and a "zero increments" action.  [Inferred]

### Atmosphere-ocean coupling by exchange-grid transforms in the Gcm
- Pattern: the Gcm builds exchange grids from the tile file and moves skin-layer fields with
  `MAPL_LocStreamTransform` in `RUN_OCEAN`; the OGCM imports that are filled this way are
  removed from automatic connectivity with `MAPL_TerminateImport` (`GEOS_GcmGridComp.F90`,
  `SetServices`, `Initialize`, `Run`).  [Verified]
- Likely reason: atmosphere and ocean decompose the tile space differently, so a plain
  field copy between components is not enough.  [Inferred]

### Replay by rewinding the clock in the Gcm
- Pattern: "Regular" replay runs a predictor loop of AGCM and ocean steps, rewinds the clock
  with `ESMF_DIRECTION_REVERSE`, restores saved alarms and refreshes children from recorded
  states (`MAPL_AddRecord`, `MAPL_GenericRefresh`) (`GEOS_GcmGridComp.F90`, `Initialize`,
  `Run`).  [Verified]
- Likely reason: the corrector step needs the model state at the start of the predictor
  window and the increments computed from it.  [Inferred]
- Question for maintainers: what is the reason replay logic sits in the Gcm?

### Resource-selected structure at `SetServices` time
- Pattern: which children exist is decided by run-time resources read in `SetServices`:
  `USE_DATA_ATM4OCN:`, `USE_WAVES:`, `DYCORE:`, `SEAICE_NAME:`, `USE_OCEANOBIOGEOCHEM:`,
  `wave_model:` in `WGCM.rc` (`GEOS_GcmGridComp.F90`; `GEOS_SuperdynGridComp.F90`;
  `GEOSogcm_GridComp/GEOSseaice_GridComp/GEOS_SeaIceGridComp.F90`;
  `GEOSwgcm_GridComp/GEOS_WgcmGridComp.F90`, each `SetServices`).  [Verified]
- Likely reason: one executable serves AMIP, coupled, ocean-only and wave configurations.  [Inferred]

### Build: composite builds its own library and its children, or only the children
- Pattern: each composite's `CMakeLists.txt` tests for its own source file; if present it builds a
  library with `SUBCOMPONENTS`, otherwise it only adds the child directories
  (`CMakeLists.txt`; `GEOSagcm_GridComp/CMakeLists.txt`;
  `GEOSagcm_GridComp/GEOSphysics_GridComp/CMakeLists.txt`). Missing child directories become
  generated stub components (`esma_add_library`, fixture `@cmake/esma_support`).  [Verified]
- On this line the root and Physics `CMakeLists.txt` no longer have the library branch (commented out), and the
  Agcm directory has a Held-Suarez library branch, taken only when its files are present (see the last entry).
- Likely reason: the same tree builds in sparse checkouts and in fixtures that omit parts.
  The CI workflow's matrix names GEOSldas as well as GEOSgcm, but passes
  `fixture-repo: GEOS-ESM/GEOSgcm` (`.github/workflows/workflow.yml`), so this is not
  evidence of a GEOSldas build.  [Inferred]

### State specs are moving to the ACG
- Pattern: GWD declares its state in `GWD_StateSpecs.rc` and generates the Fortran include
  files with `mapl_acg`; on this line Held-Suarez does the same (`GEOShs_GridComp/HS_StateSpecs.rc`,
  `GEOShs_GridComp/CMakeLists.txt`), while Superdyn declares no specs of its own and
  the other components of this repository (Agcm, Physics, Moist, Turbulence, Surface, Gcm, ...) call `MAPL_AddImportSpec` /
  `MAPL_AddExportSpec` / `MAPL_AddInternalSpec` by hand
  (`GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSgwd_GridComp/CMakeLists.txt`, `mapl_acg`). Nested
  repositories already use `mapl_acg`: `GEOS_OceanGridComp` (`GEOS_Ocean_StateSpecs.rc`),
  `MOM6_GEOSPlug`, `MOM_GEOS5PlugMod`, `GEOSdatasea_GridComp`, `CARMA`, `GOCART2G` and `CA2G`
  (their `CMakeLists.txt`). The repository's issue template "Type 1: ACG Migration" names GWD as the
  reference implementation (`.github/ISSUE_TEMPLATE/type1_acg_migration.md`).  [Verified]
- Likely reason: a machine-readable state description keeps declarations, pointer fetches and
  documentation in one file, which is why the migration is under way.  [Inferred]

### Port in steps: bypass composites that use non-ported components (MAPL3 line)
- Pattern: composites whose children are not ported are not compiled; their library blocks are commented out
  with the comment "Temporarily bypass ... top-level library: ... USE-s non-ported components. Restore the block
  below when all SUBCOMPONENTS are ported to MAPL3", the unported child directories are commented out of `alldirs`
  ("not yet ported to MAPL3 - restore when ported"), and a checkout that includes the Held-Suarez files
  (the fixture's sparse file excludes them) builds a Held-Suarez composite
  (`GEOS_AgcmSimpleGridComp.F90`) in place of the Agcm; children are configured by per-component YAML files
  (`superdyn.yaml`, `held-suarez.yaml`, `dyn.yaml`) (`CMakeLists.txt`; `GEOSagcm_GridComp/CMakeLists.txt`;
  `GEOSagcm_GridComp/GEOSphysics_GridComp/CMakeLists.txt`; `GEOSagcm_GridComp/GEOSsuperdyn_GridComp/CMakeLists.txt`;
  `GEOSagcm_GridComp/GEOS_AgcmSimpleGridComp.F90`, `SetServices`).  [Verified]
- Likely reason: let ported components (dynamics, GWD, and Held-Suarez where its files are present) be built
  and tested by regression tests while the rest is ported.  [Inferred]
