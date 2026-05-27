---
name: "Type 1: ACG Migration"
about: "Migrate a GridComp from manual MAPL_Add*Spec calls to mapl_acg"
title: "[ACG] ComponentName: migrate to mapl_acg"
labels: ["acg-migration", "mapl3-readiness"]
---

## Component

**Name:** <!-- e.g. MoistGridComp -->
**Primary file:** <!-- e.g. GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/GEOS_MoistGridComp.F90 -->
**Sub-repo:** <!-- GEOSgcm_GridComp | GEOSradiation_GridComp | GEOS_OceanGridComp | FVdycoreCubed_GridComp -->
**Effort:** <!-- S | M | L | XL -->
**MAPL depth:** <!-- leaf | depth-1 | depth-2 | depth-3 | root -->

## Current state

- Manual `MAPL_AddImportSpec` calls: <!-- count -->
- Manual `MAPL_AddExportSpec` calls: <!-- count -->
- Manual `MAPL_AddInternalSpec` calls: <!-- count -->
- Existing `*_StateSpecs.rc` file: <!-- Yes / No -->
- Existing ACG includes (`*___.h`): <!-- Yes / No / Partial -->

## Reference implementation

**GwdGridComp** is the gold standard. Study these files before starting:

- `GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSgwd_GridComp/GWD_StateSpecs.rc`
- `GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSgwd_GridComp/CMakeLists.txt`
- `GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSgwd_GridComp/GEOS_GwdGridComp.F90`

## Checklist

### Preparation
- [ ] Read `GWD_StateSpecs.rc` to understand `schema_version: 2.0.0` format
- [ ] Identify all `MAPL_AddImportSpec` / `MAPL_AddExportSpec` / `MAPL_AddInternalSpec` calls
- [ ] Check whether a `*_StateSpecs.rc` file already exists (partial migration)

### StateSpecs.rc
- [ ] Create (or complete) `<PREFIX>_StateSpecs.rc` with `schema_version: 2.0.0`
- [ ] Populate `IMPORT:` section — one entry per `MAPL_AddImportSpec` call
- [ ] Populate `EXPORT:` section — one entry per `MAPL_AddExportSpec` call
- [ ] Populate `INTERNAL:` section — one entry per `MAPL_AddInternalSpec` call (omit if none)
- [ ] Verify field count matches: `grep -c MAPL_Add.*Spec` before == entries after

### CMakeLists.txt
- [ ] Add `mapl_acg` invocation:
  ```cmake
  mapl_acg(${this} <PREFIX>_StateSpecs.rc
           IMPORT_SPECS EXPORT_SPECS INTERNAL_SPECS
           GET_POINTERS DECLARE_POINTERS)
  ```
- [ ] Remove any old manual code-generation targets if present

### GridComp F90
- [ ] Add `#include "<PREFIX>_Import___.h"`
- [ ] Add `#include "<PREFIX>_Export___.h"`
- [ ] Add `#include "<PREFIX>_Internal___.h"` (if INTERNAL section exists)
- [ ] Remove all `MAPL_AddImportSpec(...)` calls
- [ ] Remove all `MAPL_AddExportSpec(...)` calls
- [ ] Remove all `MAPL_AddInternalSpec(...)` calls
- [ ] Verify none remain: `grep -c MAPL_Add.*Spec <file>` == 0

### Build & test
- [ ] CMake configure succeeds
- [ ] Build succeeds
- [ ] CI passes

### PR
- [ ] PR title: `[ACG] ComponentName: migrate to mapl_acg`
- [ ] PR targets `main`
- [ ] PR links to this issue and references GwdGridComp as template

## Notes

<!-- Component-specific notes: tricky fields, conditional specs, CHILD_ID usage, etc. -->

## Links

- Parent epic: #<!-- meta-epic issue number -->
- Related Type 2 audit: #<!-- audit issue number -->
