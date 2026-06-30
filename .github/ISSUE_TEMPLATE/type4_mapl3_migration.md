---
name: "Type 4: MAPL3 Migration"
about: "Migrate a GridComp to MAPL3 APIs (created LATE, after all Type 1/2/3 are closed)"
title: "[MAPL3] ComponentName: migrate to MAPL3 APIs"
labels: ["mapl3-migration", "mapl3-readiness"]
---

## Component

**Name:** <!-- e.g. MoistGridComp -->
**Primary file:** <!-- e.g. GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/GEOS_MoistGridComp.F90 -->
**Sub-repo:** <!-- GEOSgcm_GridComp | GEOSradiation_GridComp | GEOS_OceanGridComp | FVdycoreCubed_GridComp -->
**Effort:** <!-- S | M | L | XL — re-estimate at creation time -->
**MAPL depth:** <!-- leaf | depth-1 | depth-2 | depth-3 | root -->

## Pre-conditions (must all be closed before this issue is opened)

- [ ] Type 1 ACG migration: #<!-- N -->
- [ ] Type 2 pre-migration audit: #<!-- N -->
- [ ] Type 3 fixes: <!-- list all #NNN, or "none spawned" -->
- [ ] Children's Type 4 issues (if this is a parent component):
  - [ ] <!-- child component Type 4 #NNN -->
  - <!-- repeat for each child -->

## MAPL3 API changes affecting this component

<!-- 
At issue-creation time, list the specific MAPL3 API changes that apply.
Common categories:

- MAPL_Cap / MAPL_InternalState changes
- MAPL_GetPointer signature changes
- SetServices / Initialize / Run / Finalize interface changes
- Clock/alarm API changes
- New required fields or removed fields
- Any component-specific MAPL2 calls that have MAPL3 equivalents

Link to MAPL3 migration guide: <!-- URL when available -->
-->

## Checklist

### Pre-work
- [ ] Confirm all pre-conditions above are closed
- [ ] Read MAPL3 migration guide for the relevant API surface
- [ ] Identify all MAPL calls in this component that have MAPL3 equivalents:
      `grep -n 'MAPL_' <file> | grep -v '! '`

### Migration
- [ ] Update `use` statements to MAPL3 module names (if changed)
- [ ] Update `SetServices` signature/registration if required
- [ ] Update `Initialize` / `Run` / `Finalize` interface signatures if required
- [ ] Update all `MAPL_GetPointer` / `MAPL_Get` calls to MAPL3 equivalents
- [ ] Update any `MAPL_AddChild` / `MAPL_SetEntryPoint` calls if API changed
- [ ] Update `MAPL_Cap` / `MAPL_InternalState` usage if changed
- [ ] Replace any remaining MAPL2-only calls with MAPL3 equivalents

### Validation
- [ ] Build succeeds against MAPL3
- [ ] Unit/regression tests pass
- [ ] Parent component (if any) is notified that this child is ready
- [ ] No MAPL2 API calls remain: `grep -c 'MAPL2_\|mapl2_' <file>` == 0

### PR
- [ ] PR title: `[MAPL3] ComponentName: migrate to MAPL3 APIs`
- [ ] PR targets the MAPL3 migration branch (not `main`)
- [ ] PR description lists all API changes made
- [ ] PR description links to this issue and all pre-condition issues

## Sub-issues

<!-- If this migration turns out to be complex, open sub-issues and list them here. -->
<!-- Sub-issue title: [MAPL3] ComponentName: <specific aspect> -->

## Notes

<!-- Any component-specific complexity, known tricky interactions,
     or deferred decisions. -->

## Links

- Parent epic: #<!-- meta-epic issue number -->
- Type 1 ACG: #<!-- N -->
- Type 2 audit: #<!-- N -->
- Type 3 fixes: <!-- list -->
- Children's Type 4s: <!-- list -->
- Parent's Type 4 (blocked by this): #<!-- N if known -->
