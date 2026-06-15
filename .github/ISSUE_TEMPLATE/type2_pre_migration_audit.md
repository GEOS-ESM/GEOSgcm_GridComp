---
name: "Type 2: Pre-Migration Audit"
about: "Read a GridComp and identify anything that would complicate the MAPL3 migration"
title: "[Audit] ComponentName: pre-migration audit"
labels: ["pre-migration-audit", "mapl3-readiness"]
---

## Component

**Name:** <!-- e.g. MoistGridComp -->
**Primary file:** <!-- e.g. GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/GEOS_MoistGridComp.F90 -->
**Sub-repo:** <!-- GEOSgcm_GridComp | GEOSradiation_GridComp | GEOS_OceanGridComp | FVdycoreCubed_GridComp -->
**Effort:** <!-- S | M | L | XL -->
**MAPL depth:** <!-- leaf | depth-1 | depth-2 | depth-3 | root -->

## Purpose

This is a **read-and-report** issue, not a fix issue. The auditor reads the
component and identifies anything that would complicate or obscure the MAPL3
migration: code smells, VERIFY_ patterns, repeated blocks, lurking bugs,
tight coupling, anything relevant.

**Output:** a written-up findings list posted as a comment on this issue,
plus zero or more new `[Fix]` issues spawned from those findings.

This issue closes when the audit writeup is complete — **not** when fixes are done.

## Audit checklist

### Error-handling patterns
- [ ] Count `VERIFY_(` calls — are they using the old macro form?
- [ ] Any bare `STATUS` checks that should use `_RC`?
- [ ] Any error paths that silently swallow failures?

### Code structure
- [ ] Repeated code blocks that could be extracted into a helper subroutine
- [ ] Very long subroutines (>200 lines) that should be split
- [ ] Dead code / commented-out blocks that should be removed

### MAPL usage
- [ ] Use of deprecated MAPL2 APIs that will break in MAPL3
- [ ] Unusual `MAPL_AddChild` patterns
- [ ] Pointer acquisition patterns that may need updating

### State/coupling
- [ ] Over-broad imports (fields imported but never used)
- [ ] Restart fields that appear inconsistently handled
- [ ] Missing or wrong units/DIMS/VLOCATION on specs

### General quality
- [ ] `IMPLICIT NONE` present?
- [ ] Any obviously wrong default values on spec fields

## Findings

<!-- Auditor fills this in. One sub-section per finding. -->
<!-- Severity: BLOCKER (must fix before MAPL3) / SHOULD / NICE -->

### Finding 1: <!-- short title -->
**Severity:** <!-- BLOCKER | SHOULD | NICE -->
**Location:** `<!-- file:line -->`
**Description:** <!-- What is the problem? -->
**Proposed fix:** <!-- What should be done? -->
**Spawns Type 3 issue:** <!-- Yes → #NNN | No -->

## Summary

**Total findings:** <!-- N -->
**BLOCKERs:** <!-- N -->
**Type 3 issues opened:** <!-- list issue numbers -->
**Audit complete:** <!-- date -->
**Audited by:** <!-- @handle -->

## Links

- Parent epic: #<!-- meta-epic issue number -->
- Related Type 1 ACG: #<!-- ACG issue number -->
- Spawned Type 3 fixes: <!-- list #NNN -->
