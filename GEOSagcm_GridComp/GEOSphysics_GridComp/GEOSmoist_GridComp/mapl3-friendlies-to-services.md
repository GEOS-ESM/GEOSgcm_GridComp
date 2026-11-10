# Moist MAPL3 port: from `FRIENDLYTO` to services

This note records how Moist's MAPL2 "friendly" fields map onto MAPL3 services (2026-10-30; `CNV_TR` and `MTR` mechanics added 2026-11-04; updated 2026-11-09 for the BACM_1M port and the removal of `CONVECTION_TRACER == 1` support). It accompanies step 4 of `mapl3-porting-notes.md` and covers branch `feature/pchakrab/port-moist-to-mapl3`. MAPL2 references are to `origin/develop`.

File abbreviations:

| Short | File |
|---|---|
| `Moist` | `GEOS_MoistGridComp.F90` |
| `GFDL_1M` | `GEOS_GFDL_1M_InterfaceMod.F90` |
| `BACM_1M` | `GEOS_BACM_1M_InterfaceMod.F90` |
| `GF` | `GEOS_GF_InterfaceMod.F90` |
| `PL` | `Process_Library.F90` |
| `Turbulence` | `../GEOSturbulence_GridComp/GEOS_TurbulenceGridComp.F90` (MAPL2) |
| `Service` | `@MAPL/superstructure/generic/specs/ServiceClassAspect.F90` |

## Summary

| MAPL2 target | Scheme | Moist items | MAPL3 replacement | Status |
|---|---|---|---|---|
| `DYNAMICS` | GFDL_1M, BACM_1M | `Q`, `CLLS`, `CLCN`, `QLLS`, `QLCN`, `QILS`, `QICN`; GFDL_1M also `QRAIN`, `QSNOW`, `QGRAUPEL` | subscribe to DynCore's `TRADV` | done |
| `DYNAMICS` | GF | `CNV_TR` if `CONVECTION_TRACER == 1` | none: `CONVECTION_TRACER == 1` is not supported (`GF:54-59`) | dropped |
| `TURBULENCE` | GFDL_1M, BACM_1M | `Q`, `QLLS`, `QLCN`, `QILS`, `QICN`; GFDL_1M also `QRAIN`, `QSNOW`, `QGRAUPEL` | subscribe to `TRDIFF` (placeholder name; Turbulence not ported) | subscription done, no provider |
| `MOIST` | GF | `CNV_TR` (`CONVECTION_TRACER == 0`, the default and the only value MAPL3 accepts; see `CNV_TR` below) | Moist provides `MTR`; its own items go in a local bundle | done |
| `CHEMISTRY` | GFDL_1M, BACM_1M | `Q` | subscribe to `TRCHEM` (placeholder name; PChem not ported) | subscription done, no provider |
| `ANALYSIS` | GFDL_1M, BACM_1M | `Q` | subscribe to `TRANA` (MAPL2 bundle name; AGCM not ported) | subscription done, no provider |

Moist has one bundle in its states in both versions: `MTR`. In MAPL2 it was an import that Physics filled; in MAPL3 it is a service Moist provides. For `TRADV`, `TRDIFF`, `TRCHEM` and `TRANA`, Moist is a subscriber: it owns the fields it lists (its internals), and the provider acts on them.

## Moist in MAPL3

### Advection: `TRADV` (subscriber)

`GFDL_1M_Setup` and `BACM_1M_Setup` push their items onto `ADVECTED` (`GFDL_1M:96-105`, `BACM_1M:68-74`); `GF_Setup` no longer takes `ADVECTED`. SetServices adds one `TRADV` import with `service_items=ADVECTED` when the list is non-empty (`Moist:186-193`). Run never fetches `TRADV`; DynCore advects the fields in place.

> [!NOTE]
> **Possible bug in the original code.** `CNV_TR` with `CONVECTION_TRACER == 1`: MAPL2 builds `trim(COMP_NAME)//"DYNAMICS"`, i.e. `"MOISTDYNAMICS"` (MAPL2 `GF:68`), with no `:` separator, so the field was friendly to no component. MAPL3 does not support `CONVECTION_TRACER == 1` (not used, Bill 2026-11-06): `GF_Setup` asserts `CONVECTION_TRACER == 0` (`GF:54-59`), so `CNV_TR` is never advected. There is no behaviour change from MAPL2's default.

### Turbulent diffusion: `TRDIFF` (subscriber)

`GFDL_1M_Setup` and `BACM_1M_Setup` push their items onto `DIFFUSED` (`GFDL_1M:86-93`, `BACM_1M:61-65`). SetServices adds one `TRDIFF` import with `service_items=DIFFUSED` (`Moist:198-205`).

> [!NOTE]
> `DIFFUSED` is `ADVECTED` without the cloud fractions `CLLS` and `CLCN`. This is intentional and matches MAPL2: GFDL_1M and BACM_1M made them `FRIENDLYTO="DYNAMICS"` only (MAPL2 `GFDL_1M:89-90`, BACM_1M `:67-68`), and Physics filled Turbulence's `TR` with Moist's `TURBULENCE` friendlies only (MAPL2 `../GEOS_PhysicsGridComp.F90:2121`). So the cloud fractions are advected but not turbulently diffused. GF's `CNV_TR` is in neither list (in MAPL2 it was never friendly to `TURBULENCE`).

> [!NOTE]
> Turbulence is MAPL2 and provides no service yet, so `TRDIFF` is unconnected. The name is a placeholder; the parent YAML can map it to whatever Turbulence provides.

### Chemistry: `TRCHEM` (subscriber)

The only consumer of Moist's `FriendlyToCHEMISTRY` is PChem; no GOCART code checks it. PChem imports `Q` as a plain field (`../@GEOSchem_GridComp/GEOSpchem_GridComp/GEOS_PChemGridComp.F90:337-343`), which shares memory with Moist's internal. Only if `MAPL_VerifyFriendly(IMPORT,'Q','CHEMISTRY')` is true (`:1547`) does PChem apply its water production/loss to it. So the friendly flag is Moist's permission for PChem to modify `Q` in place, which is what a service does.

`GFDL_1M_Setup` and `BACM_1M_Setup` push `Q` onto `CHEM_UPDATED` (`GFDL_1M:80`, `BACM_1M:55`). SetServices adds one `TRCHEM` import with `service_items=CHEM_UPDATED` (`Moist:209-216`). PChem is MAPL2 (no `MAPL_GridCompAddSpec` on `release/MAPL-v3`), so `TRCHEM` is unconnected. The name is a placeholder. Not checked: whether a one-field service is the right tool, or MAPL3 has a lighter way to let an import be modified.

### Analysis: `TRANA` (subscriber)

There is no ANALYSIS component; the AGCM acts on the fields. In its Run it gets Physics' `TRANA` (`../../GEOS_AgcmGridComp.F90:2022`), asserts it has exactly two members (`:2027`), applies the analysis increments to them in place (`:2171`ff; `OX` gets its own coefficients, `:2183`), and fills negative `Q` in place (`FILL_Friendly`, `:2774-2790`). The two members are Moist's `Q` and `OX`, from PChem by default (`ANALYSIS_OX_PROVIDER`, `GEOS_PChemGridComp.F90:372-375`) or e.g. `GOCART::OX` (AGCM `:2183`).

`GFDL_1M_Setup` and `BACM_1M_Setup` push `Q` onto `ANALYZED` (`GFDL_1M:83`, `BACM_1M:58`). SetServices adds one `TRANA` import with `service_items=ANALYZED` (`Moist:220-227`). The AGCM is MAPL2, so `TRANA` is unconnected; the name is MAPL2's bundle name. Open for the AGCM port: the two-member assert and the lookup of fields by name (`'Q'`, and `'OX'` after `MAPL_RmQualifier`); not checked whether service-bundle field names are qualified.

### Convective transport: `MTR` (provider)

MAPL2's `MTR` import (MAPL2 `Moist:488-495`) is now an export with `itemtype=MAPL_STATEITEM_SERVICE` (`Moist:248-252`), so Moist is the provider. Components that were `FRIENDLYTO=MOIST` (e.g. the GOCART tracers) must subscribe to `MTR`. None do yet, so the bundle is empty in practice.

Run reads `MTR` from both states and asserts they share data (`MAPL_FieldBundleSameData`). Because Moist cannot subscribe to its own `MTR` (see Self-subscription below), each Run builds a local bundle `TR_CNV` from MTR's fields plus `CNV_TR` (added whenever `CONVPAR_OPTION` is `GF`, `Moist:704-707`), passes it to `CNV_Tracers_Init`, and destroys it (`Moist:708-710`; the fields are not destroyed). `CNV_Tracers_Init` no longer checks `FriendlyToMOIST`: every field in the service bundle comes from a subscriber.

How the bundle gets filled and used:

1. Provider side: `create` makes an empty service bundle (`Service:111`); `allocate` is a no-op, since fields are added at connection time (`Service:160`).
2. Each subscriber declares an `MTR` import of class `service` listing its own internals, and the parent YAML connects `MOIST`'s `MTR` to it (e.g. `../regression/moist-gfdl1m-gf-uw-sa/data-tracers.yaml:18-20`, `../regression/moist-gfdl1m-gf-uw-sa/root.yaml:286-289`).
3. On connection, `connect_to_import` adds each listed subscriber internal to Moist's bundle (`Service:268-304`, add at `:294`, via `FieldClassAspect.F90:525`). The bundle holds the subscriber's own fields, not copies.
4. `add_to_state` puts the bundle in both Moist's import and export states (`Service:182`).
5. On the first Run, `CNV_Tracers_Init` allocates `CNV_Tracers`, reads each field's `ESMF_Info` metadata, and points `CNV_Tracers(F)%Q` at the field data (`PL:507`). Later Runs only assert the count is unchanged (`PL:382-383`).
6. GF (`ConvPar_GF2020.F90`, `ConvPar_GF_Shared.F90`) and UW (`uwshcu.F90:279-352`) update `CNV_Tracers(:)%Q` in place, so the results land directly in the subscribers' internals. Nothing is copied back and Run does not touch `MTR` again.

Why `CNV_TR` has to go through `CNV_Tracers` at all: neither GF nor UW has a separate convective transport routine. GF2020 copies `CNV_Tracers(:)%Q` into `se_chem` inside the column loop (`ConvPar_GF2020.F90:967`), computes in-cloud updraft/downdraft values (`ConvPar_GF_Shared.F90:101`, `:350`), accumulates tendencies (`ConvPar_GF2020.F90:1223`, `:3164`), and applies them (`:523-530`); UW reads `CNV_Tracers(m)%Q` directly. Convecting `CNV_TR` separately would need a second tracer path through both schemes.

> [!TIP]
> Possible simplification: the bundle only carries a field list into `CNV_Tracers_Init`. Passing `TR_FIELDS` (or `CNV_TR` as an extra optional field) directly would remove the per-step create/destroy of `TR_CNV`, with the same behaviour. Small signature change in `PL`.

### Self-subscription

In MAPL2 a component could mark its own internals friendly to itself. In MAPL3 that means subscribing to its own service. We assume a component cannot (connections are made between components in the parent YAML; not traced in MAPL), so each case needs a workaround:

- Moist: GF's `CNV_TR` is friendly to `MOIST` (MAPL2 with `CONVECTION_TRACER == 0`, the only value MAPL3 accepts), and Moist provides `MTR`. The local bundle `TR_CNV` (see `MTR` above) stands in for the self-subscription.
- Moist, MGB2_2M with `MGVERSION <= 1`: its MOIST-only items need the same treatment (see Out-of-scope schemes).
- Turbulence: its own internals are friendly to itself (`Turbulence:1827`, `:2170-2424`); left for the Turbulence port.

### `CNV_TR`

`CNV_TR` is GF's convection tracer (cold-pool memory). `GF_Setup` registers it as a Moist internal that starts at 0 on every run (`MAPL_RestartSkip`, default 0.0) and, in MAPL2, sets its `FRIENDLYTO` from the `CONVECTION_TRACER` resource (default 0): `"MOIST"` by default, `"MOIST"//"DYNAMICS"` = `"MOISTDYNAMICS"` when it is 1 (MAPL2 `GF:66-81`).

MAPL2 mechanics (references to the v11.11.0 checkout, MAPL v2.71.0):

- `FRIENDLYTO` is split on `:`, and each piece becomes a boolean attribute `FriendlyTo<piece>` (`MAPL_Generic.F90:6888-6915`). `"MOISTDYNAMICS"` gives only `FriendlyToMOISTDYNAMICS`.
- A `FRIENDLYTO` internal is also added to the component's exports (`MAPL_Generic.F90:3919-3946`), so `CNV_TR` is a Moist export in both cases.
- Physics fills `TRADV` with its children's `DYNAMICS` friendlies (`GEOS_PhysicsGridComp.F90:2425-2431`) and Moist's `MTR` with their `MOIST` friendlies (`:2477-2484`), Moist included. Nothing asks for `MOISTDYNAMICS`.

What happens in MAPL2:

| `CONVECTION_TRACER` | Bundles | What GF2020 does | Net effect |
|---|---|---|---|
| 0 (default) | `MTR` only | Nothing: all its `CNV_TR` code is under `CONVECTION_TRACER == 1` (`ConvPar_GF2020.F90:462`, `:614`, `:1254`, `:1680`) | Stays 0 for the whole run; convective transport of a zero field. Harmless. |
| 1 | None (missing `:`) | Uses it as buoyancy excess (`:462`ff, `:1680`ff) and updates it per column with a downdraft source and ~1 h decay (`:614-630`); the comment reads `ADV + SRC - SINK` (`:626`) | Never advected, never convectively transported. Defeats the cold-pool memory, with no error or warning. |

Only GF2020 reads `CNV_TR`; `GF_GEOS5_Interface` is not passed it (MAPL2 `GF:862-879`). `USE_GF2020` defaults to 0 at 72 levels and 1 otherwise (MAPL2 `GF:158`, `:167`), and nothing in GEOSgcm_App sets it or `CONVECTION_TRACER`, so the `CONVECTION_TRACER: 1` path is not exercised by any standard setup. The line dates from `8bf8341d` (2022-09-23) and is still on GEOSgcm_GridComp `main` and `develop`.

MAPL3 supports only `CONVECTION_TRACER == 0` (2026-11-06; Bill: `CONVECTION_TRACER: 1` is not used). `GF_Setup` asserts it (`GF:54-59`), and `CNV_TR` is only in the local `TR_CNV` bundle (convective transport, see `MTR` above), never in `ADVECTED`. This matches MAPL2's default. (Until 2026-11-06, `CONVECTION_TRACER == 1` also added `CNV_TR` to `ADVECTED`, as `MOIST:DYNAMICS` intended.)

Bug report: drafted (2026-11-04) for GEOSgcm_GridComp, not yet filed. It concerns MAPL2 only; per Bill the fix there is just the friendly statement (`":DYNAMICS"`, giving `MOIST:DYNAMICS`). If MAPL3 ever needs `CONVECTION_TRACER == 1`, restoring it means passing `ADVECTED` back to `GF_Setup` and pushing `CNV_TR` onto it.

### Not replaced

- Physics' `MTRI`/`MCHEMTRI` increment bundles read `MTR` (`../MBundle_IncrementMod.F90:102`); not ported.

## Out-of-scope schemes

What the unported schemes will need (MAPL2 `FRIENDLIES%...` settings). BACM_1M, ported, is listed for reference: its MAPL3 lists are above.

| Scheme | DYNAMICS and TURBULENCE | DYNAMICS only | MOIST only |
|---|---|---|---|
| BACM_1M (`:66-72`; ported) | `Q`, `QLLS`, `QLCN`, `QILS`, `QICN` | `CLLS`, `CLCN` | |
| THOM_1M (`:78-90`), NSSL_2M (`:78-90`) | as BACM_1M plus `QRAIN`, `QSNOW`, `QGRAUPEL`, `NCPL`, `NCPI`, `NRAIN` | `CLLS`, `CLCN` | |
| MGB2_2M (`:117-143`) | as BACM_1M plus `NCPI`, `NCPL`; with `MGVERSION > 1` also `QRAIN`, `QSNOW`, `QGRAUPEL`, `NRAIN`, `NSNOW`, `NGRAUPEL` | `CLLS`, `CLCN` | with `MGVERSION <= 1`: `QRAIN`, `QSNOW`, `QGRAUPEL`, `NRAIN`, `NSNOW`, `NGRAUPEL` |

`Q` is also `CHEMISTRY:ANALYSIS` in every scheme (`TRCHEM`, `TRANA`). MGB2's MOIST-only items need the same local-bundle treatment as `CNV_TR`. UW and RAS declare no friendlies.

## Not checked

- Unconnected service imports: accepted for `TRADV` (SU2G notes step 4); assumed the same for `TRDIFF`, `TRCHEM` and `TRANA`.
- Several service subscriptions from one component (Moist has four): no precedent among the ported components. Nothing in `Service` keys on the component, but the registry path is not traced. The three Moist regression cases run with all four unconnected; several connected subscriptions are untested.
- Whether a component can subscribe to its own service (assumed not; the `TR_CNV` workaround rests on it).
- What MAPL does with a provider service that has no subscribers (`MTR` in `moist-gfdl1m-gf-uw-sa`).
- Cost of creating and destroying `TR_CNV` every step.
- That subscribers' internal arrays never move after the first Run: `CNV_Tracers(:)%Q` pointers are set only once (`PL:507`).
- One subscriber connected to several providers: allowed (adv-dyn connects `DataMoist`'s `TRADV` to both DYN and ADV, `root.yaml:41-49`), but `add_to_bundle` uses `multiflag=.true.` (`FieldClassAspect.F90:525`), so the subscriber's own bundle probably holds its items twice. Moist never reads its subscriber bundles.
