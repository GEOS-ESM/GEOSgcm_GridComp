# GEOS_TurbulenceGridComp MAPL3 porting plan

Status: in progress, session 1. This component is substantially bigger and
structurally different from the already-ported `GEOSirrad_GridComp`/
`GEOSsolar_GridComp` (7302 lines vs ~4500/~8000 lines but with two Run
phases, six dynamically-sized `ESMF_FieldBundle` import/export specs, and
a per-field "friendly" bundle-item mechanism that IRRAD/Solar don't have).
Treat those two components' own `mapl3-porting-notes.md` files as the
baseline reference for every generic MAPL2->MAPL3 API mapping (resource
reads, pointer fetches, HConfig, private state, Edge 0-based remap, ACG
schema/columns) - this file only records what's *specific* to Turbulence
or *new* beyond what IRRAD/Solar already discovered.

There is a real, already-configured MAPL3 build at
`GEOSgcm/build-debug` (ifx compiler) with `GEOSirrad_GridComp` already
built (`build-debug/lib/libGEOSirrad_GridComp.a` exists) - unlike the
IRRAD/Solar sessions, this port CAN and SHOULD be incrementally
build-verified rather than trusted by inspection alone.

## Structural summary (this component, as of session start)

- Module `GEOS_TurbulenceGridCompMod`, single public symbol `SetServices`.
- No custom `Initialize`/`Finalize` - relies entirely on
  `MAPL_GenericSetServices`'s default versions (comment at top: "uses the
  default Initialize and Finalize methods of GEOS_Generic"). Since
  `MAPL_GenericSetServices` is stubbed/gone in MAPL3 and the framework
  drives lifecycle automatically (see IRRAD notes' `api-changes.md`
  citation), **this component likely needs NO custom Initialize at all**
  (unlike IRRAD/Solar, which both needed one for alarm/orbit setup this
  component doesn't have). Just delete the `MAPL_GenericSetServices` call
  and the `MAPL_TimerAdd` calls (see Timers section below).
- **Two Run phases**, not one: `RUN1` (line 2513-5564, contains nested
  `REFRESH` at 2845-5110 and `DIFFUSE` at 5121-5562) and `RUN2` (line
  5578-6443, contains nested `UPDATE` at 5683-6441). Registered today via
  `MAPL_GridCompSetEntryPoint(GC, ESMF_METHOD_RUN, Run1, RC=STATUS)` /
  `...Run2, RC=STATUS)` with **no `phase_name`** (relies on MAPL2's
  implicit call-order-based phase numbering).
- Module-level (non-contained) physics helpers after `end module`... no,
  actually *before* `end module`, as siblings of `SetServices`/`RUN1`/
  `RUN2`: `LOUIS_KS`, `LOUIS_KS_OPTIMIZED`, `ComputeMOScalingInterlevel`,
  `BELJAARS`, `VTRILU`, `VTRISOLVESURF`, `VTRISOLVE` (lines 6460-7298) -
  these are plain numerical kernels with explicit dummy-argument
  signatures, no MAPL/ESMF API usage at all found in them (not yet
  fully verified line-by-line, but a grep for `MAPL_`/`ESMF_` across
  lines 6460-7298 should be done before assuming - **TODO: verify these
  5 module-level routines are 100% MAPL/ESMF-API-free** before treating
  them as needing zero changes).
- **No manual internal-state / private-state pattern found**
  (`ESMF_UserCompSetInternalState`/`GetInternalState`: 0 matches,
  `MAPL_VarSpec`: 0 matches) - simpler than IRRAD/Solar in this respect,
  no `_SET_NAMED_PRIVATE_STATE` conversion needed.
- `MAPL_MetaComp` used 3x (all in `SetServices`, via
  `type(MAPL_MetaComp), pointer :: MAPL` + `MAPL_GetObjectFromGC(GC, MAPL,
  RC=STATUS)`, feeding the `MAPL_GetResource(MAPL, ...)` calls) - convert
  per IRRAD/Solar precedent: drop the `MAPL_MetaComp`/`MAPL_GetObjectFromGC`
  entirely, convert every `MAPL_GetResource(MAPL, ...)` to
  `MAPL_GridCompGetResource(GC, "LABEL", var, default=..., _RC)` (no
  trailing colon on the label - the old labels have a trailing `:`,
  e.g. `"DEBUG_TRB:"` -> `"DEBUG_TRB"`).
- `ESMF_ConfigGetAttribute` used 5x, all in `SetServices`, reading from a
  `type(ESMF_Config) :: CF` obtained via
  `ESMF_GridCompGet(GC, CONFIG=CF, NAME=COMP_NAME, RC=STATUS)`. Convert
  each to `MAPL_GridCompGetResource(GC, "LABEL", var, default=..., _RC)`
  too (same treatment as the `MAPL_GetResource` calls - simpler than
  going through `ESMF_HConfig` directly, matching IRRAD's stated
  preference "PREFER MAPL_GridCompGetResource over manual
  ESMF_HConfigIsDefined+AsXxx... when the value is a simple scalar with a
  sensible default"). The 5 labels: `EDMF_NUMUP:` (default 10),
  `SCM_SL:` (default 0), plus 3 more not yet catalogued in Run bodies -
  **TODO: grep the other 3 `ESMF_ConfigGetAttribute` call sites (only 2
  are in SetServices per the parse below; the other 3 must be inside
  RUN1/RUN2/REFRESH/DIFFUSE/UPDATE - find and convert those too, they'll
  need their own `MAPL_GridCompGetResource(gc, ...)` using the `gc` dummy
  arg of whichever subroutine they're in, or host-associated `gc` if
  inside a contained procedure)**.
- `Iam`/`COMP_NAME` traceback pattern present in `SetServices` (dead code
  under `_RC`/`VERIFY_` macros per IRRAD's established finding - safe to
  strip, but low priority, do last).
- `dflt_false`/`dflt_q` module-level parameters (`.false.`/`'Q'`) used as
  `default=` args to `ESMF_InfoGet` inside `REFRESH`/`UPDATE`'s
  bundle-friendliness loops - plain ESMF, unaffected by MAPL3 porting,
  leave as-is.

## Counts (mechanical conversion sites)

- `MAPL_AddImportSpec`: 40, `MAPL_AddExportSpec`: 144 raw grep hits but
  only 136 are live calls (the other 8 are `!`-commented-out dead spec
  blocks: `SL2DIAG`/`QT2DIAG`/`SLQTDIAG`/`W3CANUTO` and re-commented
  `T`/`U`/`V`/`QV` re-exports around lines 968-1176 - confirmed by
  reading the source directly, not just counting; leave these commented
  out, don't resurrect them), `MAPL_AddInternalSpec`: 45. All 229 raw
  `MAPL_Add*Spec` occurrences live inside `SetServices` (213-2497), none
  in `RUN1`/`RUN2`/`REFRESH`/`DIFFUSE`/`UPDATE`.
- Of the 221 parsed static specs (40 Import + 136 Export + 45 Internal):
  - 6 are `DATATYPE = MAPL_BundleItem` (TR/TRG/DTG import,
    TRI/FSTAR/DFSTAR export) - genuinely need to stay as manual
    `MAPL_GridCompAddSpec(..., itemtype=MAPL_STATEITEM_FIELDBUNDLE, ...)`
    calls, NOT movable to the YAML (ACG's `ITEMTYPE` column only supports
    `F`/`V`, confirmed by reading `MAPL_GridCompSpecs_ACG.py` directly -
    same limitation IRRAD hit for its `AERO` nested-state import).
  - 11 are `FRIENDLYTO = trim(COMP_NAME)` / `'TURBULENCE'` /
    `FRIENDLIES_SHOC` INTERNAL specs, **none of which have a separate
    EXPORT row of the same short_name** (confirmed by cross-checking
    against every parsed EXPORT `SHORT_NAME`) - convert every one to
    `ADD2EXPORT: .TRUE.` in the YAML and drop the `FRIENDLYTO` column
    entirely, exactly matching the established IRRAD/Solar precedent (no
    `GEOS_Ocean_StateSpecs.rc`-style inert-column laziness).
    Names: `ZPBL_SC`, `cu_scm`, `ct_scm`, `ssurf_scm`, `qsurf_scm`,
    `ZPBL`, `PDF_A`, `TKESHOC`, `TKH`, `QT2`, `QT3`.
  - 10 are inside a runtime `if` block and need a `COND` column instead of
    being unconditional: `SHFX_SPRAY`/`LHFX_SPRAY` (import + the export
    pair at different lines) gated on
    `DO_WAVES/=0 .and. DO_SEA_SPRAY/=0`, and `cu_scm`/`ct_scm`/
    `ssurf_scm`/`qsurf_scm` (internal) gated on `SCM_SL /= 0`. Both
    `DO_WAVES`/`DO_SEA_SPRAY`/`SCM_SL` are computed via
    `MAPL_GridCompGetResource`/`ESMF_ConfigGetAttribute`-turned-
    `MAPL_GridCompGetResource` calls **before** the first `#include`
    point in `SetServices`, so a plain `COND: DO_WAVES /= 0 .and.
    DO_SEA_SPRAY /= 0` / `COND: SCM_SL /= 0` expression referencing those
    locals will resolve fine at the `#include` site (same pattern as
    IRRAD's `USE_RRTMG .or. USE_RRTMGP` COND expressions).
  - 3 are `UNGRIDDED_DIMS = (/NUMUP/)` (EDMF plume exports
    `EDMF_PLUMES_W`/`_THL`/`_QT`), `NUMUP` a plain SetServices-local
    integer (not `self%`-based, no private state exists here) set via
    `ESMF_ConfigGetAttribute(CF, NUMUP, Label="EDMF_NUMUP:", default=10,
    ...)` near the top of `SetServices`, well before the `#include`
    point - use `UNGRIDDED: NUMUP` column (see IRRAD's note that
    `UNGRIDDED` can reference any runtime variable in scope at the
    `#include` point, not just `self%`-prefixed ones).
  - 17 have a `DEFAULT = <value>` (16x `0.0`/`0.`, 1x `TKESHOC` at
    `1e-4`). **Verified (not just inspected) that IRRAD's notes' claim
    "FILL/fill_value column DOES work end-to-end" is WRONG for this ACG
    script as it exists today**: traced the actual Python
    (`MAPL_GridCompSpecs_ACG.py`) execution path for a minimal
    `FILL: 5.0` test spec - the value correctly survives dealiasing
    (`fill` -> `fill_value`) and `digest_spec`/`get_values`, but is
    silently dropped in `emit_args` because `is_printable`
    (`o if o else False` ... i.e. `not has_flags(...) if o else False`)
    treats the `FILL_VALUE` option's empty-dict spec-table entry (`{}`)
    as falsy and short-circuits to `False` - **any SPECIFICATIONS table
    entry defined as a bare `{}` with no `FLAGS`/`MAPPING` is silently
    unemittable**, not just `FILL_VALUE` (confirmed by testing
    `/tmp/mini.rc` with a real ACG invocation, and independently
    confirmed `Solar_StateSpecs.rc`'s own real `FILL` column - e.g.
    `TAULOPAR | ... | MAPL_UNDEF | ...` - also produces zero
    `fill_value=` output when run through the real ACG script). This is
    a genuine, currently-live ACG bug, not a usage mistake. **Decision**:
    drop the `FILL` column from `Turbulence_StateSpecs.rc` entirely for
    the 16 fields whose default is `0.0` (harmless - `0.0` is already
    the natural zero-initialized state of a freshly-allocated MAPL3
    field, so dropping this is a no-op in practice). For the one field
    with a *meaningful* nonzero default, `TKESHOC` (`1e-4`), pulled it
    OUT of the YAML entirely and added it as a 7th manual
    `MAPL_GridCompAddSpec(..., fill_value=1.0e-4, add_to_export=.true.,
    _RC)` call in `SetServices` (the `fill_value` *keyword argument* to
    `gridcomp_add_spec` itself is real and unaffected by the ACG bug -
    only the YAML/ACG *column* path to it is broken) - same treatment
    class as the 6 bundle specs (kept manual because ACG can't currently
    express what's needed). **This ACG bug is worth reporting/fixing
    upstream** (trivial one-line fix: `is_printable = lambda o: not
    has_flags(...) if o is not None else False`, i.e. check for `None`
    instead of falsy-emptiness) but that's out of scope for this
    component port - documented here so a future session doesn't have to
    rediscover it.
  - 2 (`SHFX_SPRAY`, `LHFX_SPRAY` exports at lines 1967/1974) use the old
    `__RC__` macro spelling instead of `RC=STATUS`/`_RC` - functionally
    identical (confirmed elsewhere in this repo), no special handling
    needed beyond normal parsing.
  - The remaining ~192 are plain static rows: `SHORT_NAME`, `LONG_NAME`
    (-> `LONG NAME`/`standard_name` column), `UNITS`, `DIMS`
    (`MAPL_DimsHorzOnly`/`MAPL_DimsHorzVert`/`MAPL_DimsVertOnly` - all
    accepted verbatim by the ACG `DIMS` mapping table, no need to
    convert to bare `xy`/`xyz`/`z` although that's the prettier style
    IRRAD's `.rc` uses), `VLOCATION`
    (`MAPL_VLocationCenter`/`None`/`Edge` -> `C`/`N`/`E` for the `VLOC`
    column), `RESTART` (`MAPL_RestartSkip`/`MAPL_RestartOptional` ->
    `SKIP`/`OPTIONAL`).
- **Second real ACG/MAPL3 gap found while validating** (via actually
  running `python3.12 MAPL_GridCompSpecs_ACG.py Turbulence_StateSpecs.rc
  -i ... -x ... -p ...` and inspecting the generated code, not just
  reading the script): the ACG script's `restart_mode` mapping table
  offers `OPTIONAL`/`BOOT`/`SKIP_INITIAL` values that emit
  `MAPL_RESTART_OPTIONAL`/`MAPL_RESTART_BOOT`/`MAPL_RESTART_SKIP_INITIAL`
  - but the real `mapl_RestartModes_mod`
  (`infrastructure/field/RestartModes.F90`) only defines two `RestartMode`
  parameters, `RESTART_REQUIRED`/`RESTART_SKIP`, and only those two are
  re-exported (as `MAPL_RESTART_REQUIRED`/`MAPL_RESTART_SKIP`) via
  `infrastructure/field/API.F90`. The 2 `RESTART = MAPL_RestartOptional`
  specs (`SHFX_SPRAY`/`LHFX_SPRAY`) would have generated a reference to
  the nonexistent `MAPL_RESTART_OPTIONAL` symbol - a real compile error,
  caught by actually running the generator and reading its output rather
  than trusting the mapping table's *documented* options. **Fix**: mapped
  `RESTART = OPTIONAL` -> `SKIP` in `Turbulence_StateSpecs.rc` (closest
  real semantic - "don't require this field in a restart file" - these
  two fields already have `FILL`->manual-fallback of `0.0` behavior
  covering the not-in-restart-file case). Worth reporting upstream
  alongside the `FILL`/`is_printable` bug above.
- `MAPL_GetPointer`: 237 call sites (scattered across `RUN1`/`REFRESH`/
  `DIFFUSE`/`RUN2`/`UPDATE`, not concentrated in one place like IRRAD's
  blanket `Run`-top fetch) - **TODO: check whether each site's local
  variable name matches its short_name** before deciding whether a
  blanket ACG `GET_POINTERS`/`DECLARE_POINTERS` conversion is safe per
  category (IRRAD's caution: only safe when names match, otherwise do a
  plain `MAPL_GetPointer` -> `MAPL_StateGetPointer` rename, deferring the
  ACG-generated-declaration optimization). Given the scattered-across-5-
  subroutines call pattern (vs. IRRAD's single blanket fetch in `Run`),
  the safe default here is the **rename-only** treatment
  (`MAPL_GetPointer(` -> `MAPL_StateGetPointer(`) at all 237 sites,
  matching IRRAD's treatment of its own dynamic/scattered fetches - do
  NOT attempt blanket DECLARE_POINTERS/GET_POINTERS in the same pass.
- `MAPL_GetResource`: 164 call sites, ALL should convert from
  `MAPL_GetResource(MAPL, var, Label="X:", default=..., RC=STATUS)` to
  `MAPL_GridCompGetResource(gc, "X", var, default=..., _RC)` - note the
  dummy arg is `GC`/`gc` in `SetServices`/`RUN1`/`RUN2` but **check
  `REFRESH`/`DIFFUSE`/`UPDATE`'s host-association** - they're nested
  inside `RUN1`/`RUN2` so `gc` should already be in scope via host
  association (RUN1/RUN2's own dummy arg name is `GC`), no signature
  changes needed for the contained procedures themselves.
- `MAPL_TimerOn`/`MAPL_TimerOff`: 21 each (42 total), `MAPL_TimerAdd`: 2
  (`"-RUN1"` at line 2452, `"-RUN2"` at line 2484, right before the
  `MAPL_GenericSetServices` call - all three of these lines should
  simply be **deleted**, not converted - MAPL3 doesn't need an explicit
  `Add` step, confirmed by reading `gridcomp_timer_start`/
  `gridcomp_timer_stop` in `MAPL_Generic.F90:1250-1276`, neither
  references any pre-registration/`Add` step). Convert every
  `MAPL_TimerOn(MAPL, "NAME", RC=STATUS)` /
  `MAPL_TimerOff(MAPL, "NAME", RC=STATUS)` to
  `MAPL_GridCompTimerStart(gc, "NAME", _RC)` /
  `MAPL_GridCompTimerStop(gc, "NAME", _RC)` - **note**: despite
  `docs/mapl3/api-changes.md` section 8 saying `MAPL_TimerOn`/`Off`/`Add`
  are "stubbed" and implying timers are now fully automatic with no
  user-facing API, `MAPL_GridCompTimerStart`/`Stop` ARE real, currently
  exported, working MAPL3 procedures (`MAPL_Public_API.md`'s
  `mapl_Generic_mod` export list, and IRRAD's own already-applied
  conversion) - the docs are stale/aspirational here, trust the IRRAD
  precedent + actual source over the markdown doc.
- `ESMF_UserCompSetInternalState`/`GetInternalState`: 0 (no private state
  to convert, see above).
- `ESMFL_BundleGetPointerToData`: used inside `REFRESH`'s tracer-bundle
  loop (and presumably `UPDATE`'s too - **TODO: count/verify all call
  sites**, only spotted via the `Friendly|Bundle` grep, not yet an exact
  count) to fetch a named field's data pointer out of an
  `ESMF_FieldBundle` (not an `ESMF_State`). This routine **no longer
  exists in MAPL3** - `base/ESMFL_Mod.F90` was deleted upstream (MAPL
  commit `20019e0e`, "Delete ESMFL_Mod.F90 and migrate all internal
  callers", with the note "no symbols used" for internal MAPL callers -
  i.e. nothing in MAPL3 itself replaced it, it's just gone). Its old body
  (recovered from `git show 20019e0e -- base/ESMFL_Mod.F90`) was exactly:
  ```fortran
  call ESMF_FieldBundleGet(BUNDLE, NAME, FIELD=FIELD, RC=STATUS)
  call ESMF_FieldGet(field, status=fieldStatus, rc=status)
  if (fieldStatus == ESMF_FIELDSTATUS_COMPLETE) then
     call ESMF_FieldGet(field, 0, Ptr, rc=status)
  else
     NULLIFY(ptr)
  end if
  ```
  i.e. plain ESMF, no MAPL-specific logic at all. **Fix**: inline this
  3-call sequence (or a small local helper subroutine, since there are
  multiple call sites with different pointer ranks (2D for TRG/DTG-style,
  3D for TR/TRI/FSTAR/DFSTAR-style) - a generic-interface local helper
  named e.g. `TURB_BundleGetPointerToData` with 2D/3D module procedures,
  declared `private` in this module, is the cleanest 1:1 replacement
  keeping every call site's syntax identical) - do NOT try to find a
  MAPL3-native bundle-pointer-fetch function, there isn't one; this is
  intentionally-dropped legacy surface area per the deletion commit.
- `ESMF_FieldBundleGet`/`ESMF_StateGet` (11x): plain ESMF, confirmed
  unaffected by MAPL3 (ESMF itself didn't remove these), no changes
  needed beyond whatever the `ESMFL_BundleGetPointerToData` replacement
  needs internally.

## Bundle import/export specs: MAPL3 replacement confirmed via a REAL working example

Unlike IRRAD/Solar (which never needed a bundle-typed spec and so never
had to solve this), **GOCART2G/DU2G/SS2G are already MAPL3-ported AND use
exactly this pattern** - `GOCART2G_GridCompMod.F90` is the concrete
template, not just inferred from `MAPL_Generic.F90` source-reading:

```fortran
use MAPL, only: MAPL_STATEITEM_FIELDBUNDLE
...
call MAPL_GridCompAddSpec( &
     gc, &
     state_intent=ESMF_STATEINTENT_IMPORT, &   ! or _EXPORT
     short_name="AERO_DP", &
     standard_name="aerosol_deposition_ng", &
     dims="xy", &                                ! or "xyz"
     vertical_stagger=MAPL_VERTICAL_STAGGER_NONE, &  ! or _CENTER
     units="kg m-2 s-1", &
     itemtype=MAPL_STATEITEM_FIELDBUNDLE, &
     _RC)
```

`MAPL_STATEITEM_FIELDBUNDLE = ESMF_STATEITEM_FIELDBUNDLE`
(`infrastructure/esmf/StateItem.F90`) is a real, fully-supported
`itemtype` value for `MAPL_GridCompAddSpec`/`gridcomp_add_spec`
(`superstructure/generic/MAPL_Generic.F90:622`, confirmed by reading the
full argument list and body - it passes straight through to
`make_VariableSpec(itemtype=itemtype, ...)`, no itemtype-specific
branching that would reject `FIELDBUNDLE`). Also confirmed generically
documented as a supported "Extended State Item Type" in
`docs/mapl3/diffs-from-mapl2.md` ("`MAPL_STATEITEM_FIELDBUNDLE` | Bundle
of fields"). The MAPL3 YAML/ACG schema's `ITEMTYPE` column does NOT
support this (`F`/`V` only) - these 6 specs stay as manual
`MAPL_GridCompAddSpec` calls in `SetServices`, same treatment as IRRAD's
`AERO` state.

Plan for Turbulence's 6 bundle specs (TR/TRG/DTG import,
TRI/FSTAR/DFSTAR export): convert each `MAPL_AddImportSpec`/
`AddExportSpec(..., DATATYPE=MAPL_BundleItem, ...)` 1:1 to
`MAPL_GridCompAddSpec(gc, state_intent=ESMF_STATEINTENT_IMPORT/EXPORT,
short_name=..., standard_name=..., units=..., dims=...,
vertical_stagger=..., itemtype=MAPL_STATEITEM_FIELDBUNDLE,
restart_mode=MAPL_RESTART_SKIP, _RC)` (drop the old `LONG_NAME`/`UNITS`
Fortran-string-literal style in favor of matching GOCART2G's
lowercase-keyword call style since these are staying as hand-written
calls anyway - use whichever style, functionally identical, but prefer
matching the modern GOCART2G idiom for readability since both are
"real" MAPL3 code side by side now).

## Friendly-bundle mechanism: NOT the MAPL3 "SERVICE" item type - already ESMF-native, needs ~zero porting

This was the scariest-looking unknown going in (IRRAD/Solar's own notes
don't cover it, and MAPL3's `diffs-from-mapl2.md` describes the "SERVICE
state item type" as *the* MAPL2-friendly-mechanism replacement, which
sounded like it might require a significant redesign here). **Turned out
not to apply.** Turbulence's own "friendly" reads (inside `REFRESH`
line ~5316 and `UPDATE`'s equivalent block) already do:
```fortran
call ESMF_InfoGet(infoh, key='FriendlyToTURBULENCE', value=FRIENDLY, default=dflt_false, RC=STATUS)
```
via `ESMF_InfoGetFromHost`/`ESMF_InfoGet` on each individual `ESMF_Field`
pulled out of the `TR` bundle at runtime - this is **plain ESMF Info-API
usage**, not a MAPL2-specific "friendly" framework call
(`MAPL_FriendlyGet`/`MAPL_GridCompGetFriendlies`, the ones actually
stubbed per `api-changes.md` section 1, are never called here). The
`FriendlyToTURBULENCE` attribute is presumably set by whichever upstream
component (Moist/Chem, not yet ported) constructs the `TR` bundle's
fields before passing them in as an import - entirely outside this
file's responsibility. **No change needed for the friendly-bundle read
logic itself.** The `FRIENDLYTO = trim(COMP_NAME)` occurrences found in
the *spec* calls (Import/Export/Internal declarations, not the bundle
read logic) are the unrelated MAPL2-`AddSpec`-level mechanism already
covered above (-> `ADD2EXPORT`).

## Two-phase Run: confirmed pattern from a real MAPL3 component (GOCART2G)

`GOCART2G_GridCompMod.F90`'s own `SetServices` (not just its use as a
parent calling children) registers:
```fortran
call MAPL_GridCompSetEntryPoint(gc, ESMF_Method_Initialize, Initialize, _RC)
call MAPL_GridCompSetEntryPoint(gc, ESMF_Method_Run, run1, phase_name="Run1", _RC)
call MAPL_GridCompSetEntryPoint(gc, ESMF_Method_Run, run2, phase_name="Run2", _RC)
```
i.e. **`phase_name` (a string), not positional-call-order phase
numbering** - confirmed the real signature
(`gridcomp_set_entry_point`, `MAPL_Generic.F90:587`) takes an optional
`character(len=*) :: phase_name`, and the framework
(`GenericGridComp.F90:124-128`) internally maps `phase_name`s to ESMF
integer phases via `outer_meta%get_phases(ESMF_METHOD_RUN)` - so the
*string* is the durable identifier now, not call order. Apply the same
to Turbulence:
```fortran
call MAPL_GridCompSetEntryPoint(GC, ESMF_METHOD_RUN, Run1, phase_name="Run1", _RC)
call MAPL_GridCompSetEntryPoint(GC, ESMF_METHOD_RUN, Run2, phase_name="Run2", _RC)
```
A parent component would later invoke a specific phase via
`MAPL_GridCompRunChild(gc, "TURBULENCE", phase_name="Run1", _RC)` /
`..."Run2"...` (real precedent: `GOCART2G_GridCompMod.F90` itself calling
its own DU2G/SS2G children this way, lines 526/612/648). **Not done in
this session**: `GEOS_PhysicsGridComp.F90` (the natural parent) is itself
still fully bypassed in `CMakeLists.txt` (commented out, "not yet ported
to MAPL3") because `GEOSmoist_GridComp`/`GEOSsurface_GridComp` - the
components that would need to run interleaved with Turbulence's two
phases - aren't ported either. So there is currently no live parent to
wire a `phase_name="Run1"`/`"Run2"` call into; that's blocked on Moist/
Surface, not on anything in this file. This session only makes
Turbulence itself buildable as a standalone library + correctly-shaped
child component, ready to be wired in later.

## Build verification

A real configured build exists at `GEOSgcm/build-debug` (ifx, Debug).
`GEOSirrad_GridComp` is already built there
(`build-debug/lib/libGEOSirrad_GridComp.a`,
`build-debug/include/GEOSirrad_GridComp/*.mod`) confirming the toolchain
works end-to-end for an already-ported MAPL3 component. Plan: after
wiring `CMakeLists.txt` + re-running cmake configure (needed since
`GEOSturbulence_GridComp` isn't currently in the parent's
`alldirs`/subdirectory list), do incremental `make GEOSturbulence_GridComp`
(or whatever the actual target name resolves to - check
`esma_set_this()`'s generated target name) builds after each major
conversion phase, not just at the very end - this project has the
infrastructure to actually catch real compile errors early, unlike the
IRRAD/Solar sessions which explicitly said "not build-verified" because
they lacked this.

## Session 1 result: component compiles clean under MAPL3 (real build-verified)

Unlike IRRAD/Solar's sessions, this one had a real configured MAPL3 build
(`GEOSgcm/build-debug`, ifx) available throughout, and was used
continuously - every API mapping below was confirmed against the real
compiler, not just by reading source. End state:
`make GEOSturbulence_GridComp -j4` succeeds with **zero errors**,
producing `build-debug/lib/libGEOSturbulence_GridComp.so` and the
expected `.mod` files. `GEOSirrad_GridComp` was rebuilt afterward to
confirm no regression from the shared `CMakeLists.txt`/build changes.

What was done, in order:
1. Parsed all 221 static `MAPL_Add*Spec` calls out of `SetServices`
   (which shrank from 2285 to ~180 lines) into `Turbulence_StateSpecs.rc`
   via an ad hoc Python script (see "Counts" section above), validated
   directly against the real `MAPL_GridCompSpecs_ACG.py` (not just by
   inspection) before wiring into the build.
2. Wired `CMakeLists.txt` (`mapl_acg()`, `TYPE SHARED`) and uncommented
   `GEOSturbulence_GridComp` in the parent `GEOSphysics_GridComp/
   CMakeLists.txt`'s `alldirs` (needed for `esma_add_subdirectories` to
   even see the component - the top-level `GEOS_PhysicsGridComp.F90`
   library itself stays bypassed, unrelated to this component, since
   Moist/Surface aren't ported yet).
3. Fixed `SetServices`/`RUN1`/`RUN2` entry-point signatures (no
   `intent`/`optional` on `gc`; `importState`/`exportState`/`clock` also
   plain-no-intent, matching `ESMF_Interfaces.F90`'s real `I_Run`
   abstract interface exactly - confirmed by compiler
   characteristics-mismatch errors on all 3 dummy args before this fix,
   not just IRRAD's precedent).
4. Added `phase_name="Run1"`/`"Run2"` to the two
   `MAPL_GridCompSetEntryPoint(..., ESMF_METHOD_RUN, ...)` calls, and a
   new `Initialize` subroutine (registered plain, no `phase_name`)
   creating a `"turbulence_run_alarm"` (heartbeat-cadence, since this
   component never configured a distinct run interval) - this component
   is naturally MAPL3's real, working precedent for a 2-Run-phase
   component (see `GOCART2G_GridCompMod.F90`'s own
   `run1`/`run2`/`phase_name="Run1"/"Run2"` registration, confirmed via
   source read, used as the template).
5. Added `MAPL_GridCompSet(gc, activate_all_exports=.true., _RC)` in
   `SetServices` - see the dedicated section below, this is the real
   fix for the ~150 now-unsupported MAPL2 `ALLOC=`/`ALLOC=PDFALLOC`
   fetch-time force-allocation keywords.
6. Converted, file-wide, via small parenthesis-balance-aware Python
   scripts (not `sed`, needed to handle multi-line `&`-continued calls
   and keyword-argument reordering correctly):
   - `MAPL_GetResource(MAPL, var, Label="X:"/"X:", default=..., RC=...)`
     (161 sites, several distinct call-shape variants) ->
     `MAPL_GridCompGetResource(gc, "X", var, default=..., RC=...)`
     (trailing `:` stripped from every label, including dynamic
     `trim(COMP_NAME)//"...:"` labels).
   - `MAPL_GetPointer(state, ptr, 'NAME', [ALLOC=...,] RC=...)` (237
     sites) -> `MAPL_StateGetPointer(state, ptr, 'NAME', RC=...)`,
     stripping every `ALLOC=`/`alloc=` keyword (150 raw occurrences;
     43 auto-stripped by the script's top-level-arg detection, 3 more
     multi-line `alloc=associated(...) .or. associated(...)` sites
     needed a manual follow-up fix - the script's comma-splitting missed
     these specifically for a reason not fully root-caused, worth
     re-checking if this script is reused elsewhere).
   - `MAPL_TimerOn`/`MAPL_TimerOff(MAPL, "NAME"[, RC=...])` (21+21
     sites, some using a `name=` keyword) ->
     `MAPL_GridCompTimerStart`/`Stop(gc, "NAME")`.
   - 2 `MAPL_TimerAdd` calls deleted outright (no MAPL3 `Add` step
     needed, confirmed by reading `gridcomp_timer_start`/`stop` - no
     pre-registration check).
   - `MAPL_GenericSetServices` call deleted outright (this component
     never had a custom `Initialize` before - MAPL3 needed one added
     purely for the run-alarm, see above; no other generic-lifecycle
     behavior needed porting).
7. Manually converted (not scriptable - each needed real API research,
   not just mechanical renaming):
   - `MAPL_MetaComp`/`MAPL_GetObjectFromGC` (2 sites, `RUN1`/`RUN2`
     headers) -> deleted; replaced `MAPL_Get(MAPL, IM=, JM=, LM=,
     RUNALARM=, INTERNAL_ESMF_STATE=, [LATS=])` with
     `MAPL_GridCompGet(gc, grid=esmfgrid, num_levels=LM, _RC)` +
     `MAPL_GridGet(esmfgrid, IM=, JM=, _RC)` +
     `MAPL_GridCompGetInternalState(gc, internal, _RC)` +ロ (`RUN2` only)
     `MAPL_GridGetCoordinates(esmfgrid, longitudes=lons, latitudes=lats,
     _RC)` (both mandatory outputs, `lats`/`lons` changed `pointer` ->
     `allocatable` to match) + `ESMF_ClockGetAlarm(clock,
     alarmname="turbulence_run_alarm", alarm=alarm, _RC)` (replacing
     `RUNALARM=`, fed by the new `Initialize`).
   - 5 `ESMF_ConfigGetAttribute(CF, var, Label="X:", ...)` calls (2 in
     `SetServices` for `NUMUP`/`SCM_SL`, 1 more in `SetServices` for
     `DO_SHOC`, 2 in `RUN1`/`RUN2` for `DT`/`RUN_DT`) ->
     `MAPL_GridCompGetResource(gc, "X", var, ...)`, same treatment as
     the `MAPL_GetResource` conversions above; the now-unused `type
     (ESMF_Config) :: CF` declarations and their
     `ESMF_GridCompGet(gc, CONFIG=CF, ...)` fetches were removed
     entirely (nothing else in this file used `CF`).
   - `shoc.F90`/`edmf.F90`/`scm_surface.F90` (the component's 3
     supporting source files, NOT `GEOS_TurbulenceGridComp.F90` itself)
     each had one dead MAPL2-era module-name `use` line that doesn't
     exist in MAPL3: `use MAPL_SatVaporMod, only: MAPL_EQsat` (real
     MAPL3 module is `mapl_SatVapor_mod`, lowercase-with-underscore,
     confirmed via `mp_utils/SatVapor.F90` - and NOT re-exported through
     any `mapl_*_api` aggregator, so the direct module-name `use` is
     required, not `use MAPL`) and (`edmf.F90` only)
     `use MAPL_Mod, only: mapl_undef` (no such module in MAPL3 at all;
     `MAPL_UNDEF` comes from `MAPL_Constants`, already `use`d
     elsewhere in the same file for other constants).
   - `ESMFL_BundleGetPointerToData` (10 call sites across `DIFFUSE`/
     `UPDATE`, fetching a named field's pointer out of an
     `ESMF_FieldBundle` rather than an `ESMF_State`) - this MAPL2-era
     routine was deleted upstream with no replacement (MAPL commit
     `20019e0e`). Added a small 2-procedure generic
     `TURB_BundleGetPointerToData` (rank-2/rank-3) at module scope,
     reproducing the old routine's exact (plain-ESMF, no MAPL-specific
     logic) body verbatim: `ESMF_FieldBundleGet` -> `ESMF_FieldGet`
     status check -> `ESMF_FieldGet(field, 0, ptr)` if
     `ESMF_FIELDSTATUS_COMPLETE`, else `nullify(ptr)`.
   - `MAPL_MaxMin(label, array)` (9 `DEBUG_TRB`-only diagnostic-print
     call sites) - MAPL3's real `MAPL_MaxMin` (`mp_utils/
     ArrayReductions.F90`) is a completely different API now: a pure
     parallel-reduction *function* returning a `[max,min]` array and
     requiring an explicit MPI communicator, with no label/print
     behavior at all - not a drop-in replacement for a debug
     print-with-label helper. Added a small local
     `TURB_MaxMin(label, array)` generic (rank-2/rank-3) doing a plain
     (non-globally-reduced, per-PET) `print *, label, minval(...),
     maxval(...)` instead - adequate for its actual use (interactive
     debugging only, gated behind `DEBUG_TRB=.true.`, never on by
     default), not worth wiring an MPI communicator through for.
8. Both new module-scope generics (`TURB_MaxMin`, `TURB_
   BundleGetPointerToData`) needed their `interface ... module
   procedure ... end interface` blocks moved to the module's
   specification part (right before the module's own `contains`,
   alongside `public SetServices`) - Fortran doesn't allow a generic
   interface block declaring dispatch to same-module procedures to sit
   *after* `contains`, interspersed among other procedure bodies (hit a
   real `#7111 interface-body incorrectly positioned` compiler error on
   the first attempt).

## The `ALLOC=`/`activate_all_exports` finding (real architecture question, resolved)

MAPL2's `MAPL_GetPointer(EXPORT, ptr, 'NAME', ALLOC=.TRUE.)` (150 raw
call sites in this file - far more than IRRAD/Solar, which had zero)
force-allocates an EXPORT field's backing memory even if no downstream
coupler/HISTORY requested it - needed here because much of this
component's own numerics writes into "diagnostic" export arrays
unconditionally, not just when something downstream asks for them.
`MAPL_StateGetPointer` (the real, confirmed-working MAPL3 replacement
for plain pointer fetches) has **no `alloc=` keyword at all** in its
real interface (`superstructure/state/StateGetPointer.F90`'s actual
`interface StateGetPointer` module-procedure list) - one stale
`python_fortran_bridge.F90` call site references `alloc=` but that
symbol doesn't exist on any real module procedure, so that call site is
itself presumably dead/untested code, not a counter-example.

Traced MAPL3's actual field-allocation lifecycle
(`superstructure/generic/specs/StateItemSpec.F90`'s `allocate`/
`is_active`, `OuterMetaComponent/advertise_variable.F90`'s
`set_default_activation`) to find the real mechanism: a field only gets
allocated if "active", and for EXPORT/IMPORT state intent, "active" only
becomes true if the component sets a *component-wide* misc flag
(`activate_all_exports`/`activate_all_imports`, default `.false.`,
comment: "used for testing in isolation") - INTERNAL fields are always
active unconditionally, importantly). This is exposed as a real,
public, working call: `MAPL_GridCompSet(gc, activate_all_exports=.true.,
rc=status)` (`gridcomp_set` in `MAPL_Generic.F90:367`, confirmed
via `MAPL_Public_API.md`'s `mapl_Generic_mod` export list which includes
`MAPL_GridCompSet`).

**Fix applied**: `call MAPL_GridCompSet(gc, activate_all_exports=.true.,
_RC)` once in `SetServices`, forcing every declared EXPORT to always
allocate - a safe, simple over-approximation of the ~40-ish individually
`ALLOC=`'d MAPL2 exports (extra memory for exports that were genuinely
meant to stay lazy/optional, in exchange for zero risk of a
null-pointer crash on the many exports this component's numerics write
into unconditionally). All `ALLOC=`/`ALLOC=PDFALLOC` keywords were
simply dropped during the `MAPL_GetPointer` -> `MAPL_StateGetPointer`
conversion; the `PDFALLOC` logical variable itself is still computed
(from `PDFSHAPE`) but may now be otherwise-unused - not cleaned up this
session, harmless (unused-variable warning at worst).

## KNOWN UNRESOLVED CORRECTNESS ISSUE (highest-priority follow-up): Edge (VLOC=E) 0-based indexing NOT yet audited/fixed

This is the single most important remaining risk before trusting this
component at runtime - flagging clearly and NOT attempting a rushed
fix, matching IRRAD's own real, confirmed bug of the same shape (see
IRRAD's `mapl3-porting-notes.md`, "Edge (VLOC=E) fields need an explicit
0-based bounds remap").

**The underlying MAPL3 fact** (confirmed once already for IRRAD, applies
identically here): field creation always gives the vertical dimension
Fortran bounds `1:num_levels` regardless of Center/Edge staggering, so
an Edge field (`num_levels = LM+1`) comes back from
`MAPL_StateGetPointer` bounded `1:LM+1`, **never** `0:LM` - there is no
automatic "Edge means 0-based" convention anywhere in MAPL3.

**Confirmed real usage evidence in THIS file**: a plain grep for
`NAME(:,:,0)`/`NAME(:,:,0:` 0-based-index patterns found genuine
0-based indexing on at least `PLE` (7 sites), `KM`/`KH` (2 each),
`EKM`/`EKH`/`KHRAD`/`KHSFC`/`TKE` (1 each), `QTFLXMF`/`SLFLXMF` (2
each), `UFLXTRB`/`VFLXTRB` (1 each), plus `ZL0` (a local var derived
from `ZLE`) - i.e. **this is real and affects fields this component
both imports (PLE) and internally computes/exports (KM/KH/EKM/EKH/
KHSFC/KHRAD/TKE/...)**.

**Why this was NOT fixed this session** (explicit scope decision, not an
oversight): there are 48 total `VLOC=E` fields in
`Turbulence_StateSpecs.rc` (`ALH, ALM, DU, EDMF_DRY_*(7), EDMF_MOIST_*
(7), EDMF_PLUMES_*(3), EDMF_TKE, EDMF_WQT, EDMF_WSL, EKH, EKM, ISOTROPY,
KH, KHLS, KHRAD, KHSFC, KM, KMLS, MFAW, PLE, PREF, QTFLXMF, QTFLXTRB, RI,
SHOCPRNUM, SLFLXMF, SLFLXTRB, TKE, TKH, UFLXTRB, VFLXTRB, WQT_DC, ZLE,
ZLES` - far more than IRRAD's 16+2). A reliable per-field
0-based-vs-1-based determination requires checking each field's *local
Fortran variable name* at its actual fetch site (which frequently
differs from the field's own `short_name` - e.g. `EDMF_MF` is fetched
into a local var literally named `edmf_mfx`, confirmed by grep) and
every downstream usage of that specific local variable, not just a
short_name-text grep (which undercounts real 0-based usage for any
field whose local alias differs from its short_name - most of the
"0 hits" results in this session's spot-check are almost certainly
false negatives for exactly this reason, not evidence those fields are
safely 1-based). Doing this correctly for all 48 fields, across a
~5300-line file with dozens of distinct fetch sites in `RUN1`/
`REFRESH`/`DIFFUSE`/`RUN2`/`UPDATE`, needs a full read of each
subroutine, not a mechanical script - the risk of a rushed/partial fix
(missing a real 0-based site, or wrongly remapping a field that's
actually used 1-based like IRRAD's `PREF`) is high enough that leaving
it clearly flagged is better than guessing.

**Recommended approach for the follow-up session** (once resumed):
1. For each of the 48 `VLOC=E` short names, find every
   `MAPL_StateGetPointer(state, LOCAL_VAR, 'SHORT_NAME', ...)` fetch
   site and note `LOCAL_VAR` (may differ from `SHORT_NAME`).
2. For each `LOCAL_VAR`, grep the enclosing subroutine (and anything it
   calls that receives `LOCAL_VAR` as an explicit argument, e.g.
   `LOUIS_KS`/`BELJAARS`/`ComputeMOScalingInterlevel` all take
   explicit-argument edge arrays) for `(...,0)`/`(...,0:`/
   `(...,LM+1)`-vs-`(...,LM)` patterns to determine 0-based vs 1-based
   usage - do NOT trust the `short_name` text alone.
3. For every field confirmed 0-based, apply the same fix IRRAD used:
   fetch through a local `real, pointer, contiguous, dimension(:,:,:)
   :: p3d` scratch variable (`p3d => LOCAL_VAR; LOCAL_VAR(1:IM,1:JM,
   0:LM) => p3d`) immediately after the fetch - plain self-remap
   (`LOCAL_VAR(bounds) => LOCAL_VAR`) is expected to fail to compile
   the same way it did for IRRAD (rank-remap requires a
   simply-contiguous, and apparently non-self, data-target) but this
   was NOT re-verified against `ifx` in this session (IRRAD's finding
   was against `gfortran-15`) - a quick isolated one-file `ifx` test of
   both the plain self-remap and the scratch-pointer forms should be
   done first before assuming IRRAD's gfortran-specific finding
   transfers exactly.
4. Also check the 5 explicit-argument module-level physics kernels
   (`LOUIS_KS`, `LOUIS_KS_OPTIMIZED`, `ComputeMOScalingInterlevel`,
   `BELJAARS`, plus the `VTRILU`/`VTRISOLVE*` solvers) for their OWN
   internal 0-based-vs-1-based assumptions on whichever edge arrays get
   passed into them - these were NOT audited at all this session (only
   confirmed compile-clean, never read for numerical/indexing
   correctness).

## Two other real ACG/MAPL3 API gaps found (beyond the FILL_VALUE/RESTART_OPTIONAL bugs above)

- **`ungridded_dim_array` YAML column emits the wrong type**: the ACG
  `UNGRIDDED`/`ungridded_dim_array` column (used here for
  `EDMF_PLUMES_W`/`_THL`/`_QT`, each needing `NUMUP` plume-count
  ungridded dims) emits a bare `[NUMUP]` integer-array literal, but the
  real `gridcomp_add_spec` dummy argument is `type(UngriddedDim),
  optional, intent(in) :: ungridded_dim_array(:)` - confirmed via a real
  compile error (`#6285 no matching specific subroutine`) before fixing.
  **Fix**: pulled these 3 specs out of the YAML entirely (same
  "genuinely can't be expressed via ACG" treatment as the 6 bundle
  specs and `TKESHOC`'s fill_value), added as manual
  `MAPL_GridCompAddSpec` calls using the real
  `MAPL_UngriddedDim(extent, name=, units=)` constructor (confirmed
  public via `infrastructure/esmf/API.F90`):
  `ungridded_dim_array=[MAPL_UngriddedDim(NUMUP, name='plume',
  units='1')]`.
- Combined with the `FILL`/`is_printable` bug and the
  `MAPL_RESTART_OPTIONAL`-doesn't-exist bug documented above, this
  session found **4 distinct real bugs/gaps** in
  `MAPL_GridCompSpecs_ACG.py`/its target API surface, all confirmed via
  actually running the generator and/or a real compile, not just code
  reading. Worth filing upstream as a batch:
  1. `is_printable()`'s `if o else False` treats a `{}`-valued
     SPECIFICATIONS option (e.g. `FILL_VALUE`) as unprintable regardless
     of whether a real value was supplied - breaks the `FILL` column
     entirely.
  2. `restart_mode` mapping table offers `OPTIONAL`/`BOOT`/
     `SKIP_INITIAL` values that don't exist as real exported
     `MAPL_RESTART_*` symbols (only `REQUIRED`/`SKIP` are real).
  3. `ungridded_dim_array` emission doesn't wrap the extent in
     `MAPL_UngriddedDim(...)`, so any *runtime-variable* (as opposed to
     literal-integer, untested either way) `UNGRIDDED` value produces a
     type mismatch.
  4. (Separately, in `emit_get_pointer`'s `ALLOC` handling) the
     `GET_POINTERS`-generated code path emits `alloc=` to
     `MAPL_StateGetPointer`, a keyword that doesn't exist on the real
     function - not hit directly this session (Turbulence never enabled
     `GET_POINTERS`/`DECLARE_POINTERS`) but would break any future
     component that both uses `GET_POINTERS` AND has an `ALLOC` column
     value in its `.rc`.

## Session 1 progress log

- Completed: full investigation above (this file).
- Completed: parsed all 221 static `MAPL_Add*Spec` calls out of
  `SetServices` into `/tmp/turb_specs.json` via a small ad hoc Python
  script (not committed anywhere permanent - re-derivable from the .F90
  source if needed; the script logic is described in this file's
  "Counts" section above in enough detail to reproduce it).
- Next: generate `Turbulence_StateSpecs.rc`, rewrite `SetServices`, wire
  `CMakeLists.txt`, then the mechanical `RUN1`/`RUN2`/`REFRESH`/
  `DIFFUSE`/`UPDATE` conversions, then attempt the first real build.
</content>
