---
name: "Type 3: Pre-Migration Fix"
about: "Fix a specific finding from a pre-migration audit (Type 2)"
title: "[Fix] ComponentName: <specific description>"
labels: ["pre-migration-fix", "mapl3-readiness"]
---

## Component

**Name:** <!-- e.g. SurfaceGridComp -->
**Primary file:** <!-- e.g. path/to/GEOS_SurfaceGridComp.F90 -->
**Sub-repo:** <!-- GEOSgcm_GridComp | GEOSradiation_GridComp | GEOS_OceanGridComp | FVdycoreCubed_GridComp -->
**Fix category:** <!-- _rc-cleanup | refactor | helper-function | bug | cleanup | other -->
**Effort:** <!-- S | M | L | XL -->

## Parent audit

**Spawned by:** #<!-- Type 2 audit issue number -->
**Finding title:** <!-- copy the finding title from the audit -->
**Finding severity:** <!-- BLOCKER | SHOULD | NICE -->

## Problem description

<!-- 
Describe the specific problem identified in the audit. Be concrete:
- What is the current code doing / not doing?
- Why is this a problem for the MAPL3 migration or for correctness?
- Include file + line number reference(s) if helpful.
-->

## Proposed solution

<!--
Describe what the fix looks like. Enough detail that anyone picking up
the issue can start without needing to re-read the audit. Examples:

For _rc-cleanup:
  "Replace all VERIFY_(STATUS) calls with _RC. There are 316 occurrences.
   Use: sed -i 's/VERIFY_(STATUS)/_RC/g' ... as a starting point, then
   verify by inspection that no STATUS variable remains."

For refactor:
  "Lines 450-510 and 620-680 are near-identical restart read blocks.
   Extract into a private subroutine read_restart_tile(GC, IMPORT, CLOCK, RC)."

For helper-function:
  "The MAPL_GetPointer / bounds-check / copy pattern appears 14 times.
   Introduce get_field_pointer(GC, name, ptr, RC) helper."

For bug:
  "Line 233: the VERIFY_ check fires after the pointer is dereferenced,
   not before. Reorder to check before dereference."
-->

## Acceptance criteria

- [ ] <!-- specific, testable criterion 1 -->
- [ ] <!-- specific, testable criterion 2 -->
- [ ] CI passes on `main`
- [ ] No new `VERIFY_(` calls introduced (for `_rc-cleanup` issues)

## Notes

<!-- Any caveats, dependencies on other Type 3 fixes, or component-specific context. -->

## Links

- Parent Type 2 audit: #<!-- N -->
- Parent epic: #<!-- meta-epic issue number -->
- Related fixes (if ordering matters): #<!-- N -->
