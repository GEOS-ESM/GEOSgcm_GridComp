# GEOSphysics_GridComp (Physics): architecture

**Audience:** coding agents and developers working on the Physics parent gridded component (the driver that runs gravity wave drag, moist processes, turbulence, surface, chemistry and radiation).
**Describes:** the MAPL3 line, branch `release/MAPL-v3`.
**Verified against:** `GEOSgcm_GridComp` @ `84f0aae3` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `f50622d`; repository `GEOSgcm_GridComp`, mounted at `src/Components/@GEOSgcm_GridComp`): `src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp`.
Other fixtures (for example GEOSadas, GEOSldas) may place this repository elsewhere; paths inside this repository do not change.

Paths below are relative to `GEOSagcm_GridComp/GEOSphysics_GridComp/` unless they start with `GEOSagcm_GridComp/`.
Children: [GWD](GEOSgwd_GridComp/ARCHITECTURE.md), [Moist](GEOSmoist_GridComp/ARCHITECTURE.md), [Turbulence](GEOSturbulence_GridComp/ARCHITECTURE.md), [Surface](GEOSsurface_GridComp/ARCHITECTURE.md). Chemistry (architecture (`@GEOSchem_GridComp/docs/ARCHITECTURE.md` in the fixture)) and Radiation (architecture (`@GEOSradiation_GridComp/docs/ARCHITECTURE.md` in the fixture)) are separate repositories (section 2). Repository overview: [../../docs/ARCHITECTURE.md](../../docs/ARCHITECTURE.md). Conventions for agents: [../../AGENTS.md](../../AGENTS.md). Decision log: [../../docs/decisions.md](../../docs/decisions.md).

> **Build status on this line:** the top-level library is not compiled, not
> ported. On `release/MAPL-v3`, `GEOSphysics_GridComp/CMakeLists.txt` calls
> `esma_add_subdirectories` on `GEOSchem_GridComp` ("restricted to GOCART only"),
> `GEOSgwd_GridComp` and `GEOSradiation_GridComp` ("restricted to IRRAD only"; the nested `@GEOSradiation_GridComp/CMakeLists.txt` builds a library with `GEOS_RadiationShared`, Solar, IRRAD and SatSim);
> `GEOSmoist_GridComp`, `GEOSsurface_GridComp` and `GEOSturbulence_GridComp` appear
> only in commented lines ("not yet ported to MAPL3 - restore when ported"), and
> the block that builds `GEOS_PhysicsGridComp.F90` and `MBundle_IncrementMod.F90` is
> commented out ("USE-s non-ported components"). `GEOS_PhysicsGridComp.F90` is identical to
> the MAPL2 line; `MBundle_IncrementMod.F90` differs only in its attribute calls
> (`ESMF_InfoSet` / `ESMF_InfoGet` instead of `ESMF_AttributeSet` / `ESMF_AttributeGet`
> for `DIMS`, `VLOCATION`, `UNITS`, `LONG_NAME`). This document
> describes the MAPL2-style composite, kept for a later port; there is no MAPL3 Physics
> composite on this line. Built children: GWD (ported; standalone regression
> `GEOSgwd_GridComp/regression`, see the GWD document), Chemistry and Radiation (nested
> repositories, restricted as stated). Statements below about the other children
> (Moist, Turbulence, Surface) describe code that is not built.

## 1. Purpose
Physics is the parent of every column-physics process in the atmospheric model. It owns no physical scheme itself. It does five things (`GEOS_PhysicsGridComp.F90`, header comment, SetServices, Run):
1. Creates the six children and wires their imports and exports together (`MAPL_AddConnectivity`).
2. Runs the children in a fixed order, once per Physics step.
3. Hands updated temperature, winds and moisture from one child to the next ("synchronization", sections 3.2 and 7).
4. Collects the children's tendencies into the five quantities Dynamics needs: `DUDT`, `DVDT`, `DWDT`, `DTDT`, `DPEDT` (section 3.3).
5. Builds the lists of "friendly" tracers (fields that declare which components may transform them) that Turbulence, Moist and Dynamics use (section 4.3).

## 2. Place in the hierarchy
- Parent is the Agcm (atmospheric general circulation model) component, which adds Physics as `PHYS` (`GEOSagcm_GridComp/GEOS_AgcmGridComp.F90`, SetServices); that Agcm file is itself not compiled on this line ([Agcm document](../ARCHITECTURE.md)). In a checkout that includes the Held-Suarez files (the fixture's sparse file excludes them), the `PHYS` child of the Held-Suarez composite is the Held-Suarez physics, not this composite (`GEOSagcm_GridComp/GEOS_AgcmSimpleGridComp.F90`, `SetServices`). Agcm calls Physics with `PHASE=1` and then applies the five tendencies to Dynamics (`GEOSagcm_GridComp/GEOS_AgcmGridComp.F90`, Run, `DO_UPDATE_PHY`).
- Children, added in this order with `MAPL_AddChild` (`GEOS_PhysicsGridComp.F90`, SetServices):

| Module variable | Name | `NAME=` (`MAPL_AddChild`) | Set-services | Where |
|---|---|---|---|---|
| `GWD` | gravity wave drag | `'GWD'` | `GwdSetServices` | `GEOSgwd_GridComp/` (this repository) |
| `MOIST` | moist processes | `'MOIST'` | `MoistSetServices` | `GEOSmoist_GridComp/` (this repository) |
| `TURBL` | turbulence | `'TURBULENCE'` | `TurblSetServices` | `GEOSturbulence_GridComp/` (this repository) |
| `CHEM` | chemistry, aerosols | `'CHEMISTRY'` | `AChemSetServices` | `@GEOSchem_GridComp/`, nested repository `GEOSchem_GridComp` (fixture HEAD `d6ad829`; architecture (`@GEOSchem_GridComp/docs/ARCHITECTURE.md` in the fixture)) |
| `SURF` | surface | `'SURFACE'` | `SurfSetServices` | `GEOSsurface_GridComp/` (this repository) |
| `RAD` | radiation | `'RADIATION'` | `RadiationSetServices` | `@GEOSradiation_GridComp/`, nested repository `GEOSradiation_GridComp` (fixture HEAD `0b402f6`; architecture (`@GEOSradiation_GridComp/docs/ARCHITECTURE.md` in the fixture)) |

The first column is the module variable holding the child's index (`GEOS_PhysicsGridComp.F90`, SetServices); for `TURBL`, `CHEM`, `SURF` and `RAD` this differs from the string passed as `NAME=` (third column), unlike `GWD` and `MOIST` where the two match (`GEOS_PhysicsGridComp.F90`, SetServices, `MAPL_AddChild` calls).

- Chemistry is added before Surface on purpose, so that Chemistry can fill the `AERO_DP` bundle before Surface reads its import restart (`GEOS_PhysicsGridComp.F90`, SetServices, comment at `MAPL_AddChild(... 'CHEM' ...)`).
- The child order in `AddChild` does **not** set the run order. `Run` calls the children in an order written in the code (section 3.2).
- Chemistry and Radiation are nested repositories that this repository's `.gitignore` excludes (`.gitignore`). This document only describes them as children. Their own structure is documented in those repositories: Chemistry (`@GEOSchem_GridComp/docs/ARCHITECTURE.md` in the fixture) (`@GEOSchem_GridComp/`; aerosols: GOCART (`@GEOSchem_GridComp/@GOCART/docs/ARCHITECTURE.md` in the fixture)) and Radiation (`@GEOSradiation_GridComp/docs/ARCHITECTURE.md` in the fixture) (`@GEOSradiation_GridComp/`), in the fixture.
- On the MAPL2 line `CMakeLists.txt` builds `GEOS_PhysicsGridComp.F90` and `MBundle_IncrementMod.F90`, with `alldirs` listing the child directories, linking `MAPL`, `GMAO_mpeu`, `GMAO_stoch` and ESMF. On this line that `esma_add_library` block is commented out and only `esma_add_subdirectories (${alldirs})` runs, with `alldirs` as in the build status above (`CMakeLists.txt`).

## 3. Entry points
Registered with `MAPL_GridCompSetEntryPoint` (`GEOS_PhysicsGridComp.F90`, SetServices): `INITIALIZE` = `Initialize`, and one `RUN` = `Run`. There is **one run phase** and **no Finalize** entry point (MAPL generic default).

### 3.1 SetServices
Reads `USE_OCEANOBIOGEOCHEM`, `USE_WAVES`, `USE_SEA_SPRAY`, `SYNCUV`, `SYNCTQ` (both default 1.0), `CO2_PROVIDER` and `CO2` (from `AGCM.rc`), `ATM_CO2` (from the Surface resource file `SURFRC`) and `USE_CNNEE` (from `CO2_GridComp.rc`) (`GEOS_PhysicsGridComp.F90`, SetServices). These choose which connections and `MAPL_TerminateImport` calls are made (section 5).
Then: adds the children, declares its own specs, adds exports promoted from children (`MAPL_AddExportSpec` with `CHILD_ID`), adds the connectivities, terminates the imports that Physics fills by hand, and ends with `MAPL_GenericSetServices`.

### 3.2 Run (`GEOS_PhysicsGridComp.F90`, Run)
Everything below happens when `DOPHYSICS==1` (resource `DOPHYSICS`). Between most steps Physics calls `MAPL_GenericRunCouplers` for the child just run. This only copies the exports of connections for which `MAPL_AddConnectivity` had to create a coupler component in `SetServices`, because the two specs differed; a connection whose specs match exactly is wired as a direct, shared field at `SetServices` time (`MAPL_Reconnect`) and needs no run-time copy (detail and citations: [Agcm architecture](../ARCHITECTURE.md), section 3.1). So "couplers run" below means only that any coupler created for that child's connections is invoked, not that every connected export is copied at that point.

1. `FILLQ2ZERO` "PHYS0" on the water species; `QW` is set to their sum (section 3.4).
2. **GWD** runs first. A comment says it must see the `Q` from Dynamics.
3. Reset the increment bundles for Moist (`Initialize_IncBundle_run`, `Initialize_IncMBundle_run`).
4. **MOIST** runs. Then `FILLQ2ZERO` "AFMST" and `Compute_IncBundle` / `Compute_IncMBundle`.
5. **Synchronization, stage 1**: Moist's after-Moist state `UAFMOIST`, `VAFMOIST`, `TAFMOIST`, `THAFMOIST`, `QAFMOIST`, `SAFMOIST` is copied into Surface, Turbulence and (when `SYNCTQ==1`) Chemistry imports. `DTDT_BL` and `DQDT_BL` (Moist imports) are set from the state after Moist.
6. **SURF phase 1** (no coupler run), **CHEM phase 1** (coupler run), **TURBL phase 1** (no coupler run) (`GEOS_PhysicsGridComp.F90`, Run, "Surface Stage 1", "Aerosol/Chemistry Stage 1", "Turbulence Stage 1"). Then `FILLQ2ZERO` "AFSTG1".
7. **Synchronization, stage 2**: Turbulence's `UAFDIFFUSE`, `VAFDIFFUSE`, `SAFDIFFUSE` replace the winds and the dry static energy `S` in the Turbulence and Surface imports. `T` is derived from `S` as `(S - g*z_mid)/cp`, so all change in `S` is assumed to be a change in `T`.
8. **SURF phase 2**, **TURBL phase 2** (couplers run). Then `FILLQ2ZERO` "AFSTG2".
9. **Synchronization, stage 3**: Radiation `T` is derived from Turbulence's `SAFUPDATE`; the Turbulence `S` is reset to `SAFUPDATE`; Chemistry `T`, `TH` are reset when `SYNCTQ==1`. `DTDT_BL`, `DQDT_BL` become the change over the step divided by `DT`.
10. **CHEM phase 2**, **RAD** (couplers run).
11. Stochastic pattern generation (`SPPT`, `SKEB`), then the tendency block (3.3).

Child timers `TIME_IN_GWD`, `TIME_IN_MST`, `TIME_IN_TRB`, `TIME_IN_SFC`, `TIME_IN_CHM`, `TIME_IN_RAD` use `MPI_Wtime`; Surface, Turbulence and Chemistry add their two phases together (`GEOS_PhysicsGridComp.F90`, Run).
Phases: Surface and Chemistry register `Run1`/`Run2`, Turbulence `Run1`/`Run2`; Moist, GWD and Radiation have one run phase each.

### 3.3 Tendency block (end of Run)
- Winds: `DUDT = UIM + UIT + UIG`, `DVDT = VIM + VIT + VIG`, `DWDT = WIM`. `UIM`, `WIM` come from Moist exports, `UIT` from the Physics export `UIT` (Turbulence increment), `UIG` from GWD `DUDT`. `SPPT` multiplies the wind tendencies by a random pattern; `SKEB` adds backscatter terms (`GEOS_PhysicsGridComp.F90`, Run).
- Temperature (mass-weighted, `Pa K s-1`): `TOT = TIR + STN + TTN + FRI + TIG + TICU`. These are Radiation `DTDT`, the Turbulence increment of dry static energy `SIT/CP`, Moist `DTDT`, Turbulence `INTDIS + TOPDIS` (frictional heating), GWD `DTDT`, and Moist `DTDTFRIC`. `NEED_TOT` is hard-coded `.TRUE.`; a missing pointer aborts with `status=99`. `SCM_NO_RAD` drops the radiation term (`GEOS_PhysicsGridComp.F90`, Run).
- Dry-mass step (resource `DPEDT_PHYS`, default true): the water change `DQ` of the whole column is accumulated downward into the pressure tendency `DPEDT`. Physics then builds a new pressure `ple_new`, **scales every `TRADV` tracer by the ratio of old to new layer mass** (cloud fractions `CLCN`, `CLLS` excluded), and re-weights `DTDT` for the new layer mass. Water rescale increments go to the `H2ORTRI` bundle. With `DPEDT_PHYS` false the pressure tendency is zero (`GEOS_PhysicsGridComp.F90`, Run).
- Under `SPPT`, `QV` is perturbed and all water is scaled by an area-mean ratio to conserve dry mass (`GEOS_PhysicsGridComp.F90`, Run).
- Diagnostic exports are computed here: `PERAD`, `PETRB`, `PEMST`, `PEFRI`, `PEGWD`, `PECUF`, `PEPHY` (potential-energy tendencies), `DTDT`, `DTDTTOT`, `DTDTRAD`, `DQ*DT{MST,PHY}INT`, `DQVDTTRBINT`, and others (`GEOS_PhysicsGridComp.F90`, SetServices exports, Run).

### 3.4 `FILLQ2ZERO` (`GEOS_PhysicsGridComp.F90`, FILLQ2ZERO)
Per-column fill of negative water. It clips negatives to zero and then spreads the lost column water over the positive layers. Called at five points: `PHYS0` (before Physics), `AFMST`, `AFSTG1`, `AFSTG2`, `PHYS1` (end). Applies to `QV`, `QLLS`, `QLCN`, `QILS`, `QICN`, `QRAIN`, `QSNOW`, `QGRAUPEL`, all pointers into the `TRADV` bundle. Optional exports `DQ*DT_FILL_<point>` hold the change. A warning with the global count of negatives is printed from rank 0. Species absent from `TRADV` are replaced by a zero array.
`QW` (4.1) holds the total water at the start of the step, so the total-water change `DQ` can be formed after the last fill (`PHYS1`) (`GEOS_PhysicsGridComp.F90`, Run).

### 3.5 Initialize (`GEOS_PhysicsGridComp.F90`, Initialize)
Calls `MAPL_GenericInitialize` first (children are initialized), then builds the tracer bundles in section 4.3.

## 4. State
### 4.1 Internal state and restarts
- One internal spec, `QW` (passive total water), declared `RESTART=MAPL_RestartSkip`, `FRIENDLYTO='TURBULENCE:MOIST'` (`GEOS_PhysicsGridComp.F90`, SetServices). It is reset at the start of every `Run`.
- Physics has no other internal state. Child restart files are named in `AGCM.rc` (`MODELATMGWD_INTERNAL_RESTART_FILE`, `MODELATMTURBULENCE_..._RESTART_FILE`, and similar keys); see each child's document (`src/Applications/@GEOSgcm_App/AGCM.rc.tmpl` in the fixture).

### 4.2 Imports and exports
- Physics declares few imports of its own: `U`, `V`, `W`, `DZ`, `T`, `S`, `ZLE`, `PLE`, `AREA`, `DTDTDYN`, `DQVDTDYN`, and the `*_DYN_IN` fields (`GEOS_PhysicsGridComp.F90`, SetServices). Dynamics supplies them through Agcm connections.
- `PHIS`, `SGH`, `VARFLT` are not declared in Physics, but Agcm reads them from Physics's import state and fills them by hand, and terminates their connectivity (`GEOSagcm_GridComp/GEOS_AgcmGridComp.F90`, `MAPL_TerminateImport` with `CHILD = PHYS`; Initialize). *Inferred:* MAPL generic promotes unsatisfied child imports (GWD, Turbulence) to the parent, which is why they exist in Physics's import state at all.
- Exports: tendencies and energy diagnostics (3.3), fill diagnostics, stochastic-physics fields, timers, and the bundles in 4.3. Exports promoted from children with `CHILD_ID`: Chemistry `O3PPMV`, `OX`; Turbulence `PPBL`; Moist `Q`, `QCTOT`; and fifteen Surface fields (`U10M`, `V10M`, `U10N`, `V10N`, `SNOMAS`, `WET1`, `TSOIL1`, `LWI`, `TS`, `FRLAND`, `FRLANDICE`, `FRLAKE`, `FROCEAN`, `FRACI`, `Z0`) (`GEOS_PhysicsGridComp.F90`, SetServices).
- Increment exports `<name>IT` are created from the table `TRI_increments::` in `AGCM.rc` (`ESMF_ConfigGetDim`), one per diffused tracer (`GEOS_PhysicsGridComp.F90`, SetServices). `AOADAYS` gets the unit `days s-1` there.

### 4.3 Tracer bundles (built in Initialize)
Fields carry attributes `FriendlyTo<NAME>`, set from the `FRIENDLYTO` string (colon-separated) on the spec. `MAPL_GridCompGetFriendlies` collects the fields that carry a given name (`GEOS_PhysicsGridComp.F90`, Initialize; MAPL2 (MAPL `release/v2`, `generic/MAPL_Generic.F90`; not present in this fixture)).

| Bundle | Where | Content |
|---|---|---|
| `TR` | Turbulence import | Dynamics `S` (`DiffuseLike='S'`, `WeightedTendency`), `U`, `V`; `TKESHOC`, `QT2` from Turbulence's exports; friendlies `TURBULENCE` from Physics (`QW`), Moist and Chemistry (Chemistry's share comes through its internal bundle `CHEM_TRACERS`: Chemistry document (`@GEOSchem_GridComp/docs/ARCHITECTURE.md` in the fixture)) |
| `TRI` | Turbulence export | increments `<NAME>IT` from Physics exports, `QVIT` aliased as `QIT`; missing fields become empty |
| `TRADV` | Physics export | friendlies `DYNAMICS` from Physics and all children; linked to Dynamics by Agcm |
| `TRANA` | Physics export | friendlies `ANALYSIS` |
| `MTR` | Moist import | friendlies `MOIST` |
| `H2ORTRI`, `MTRI`, `MCHEMTRI` | Physics exports | water-rescale, Moist-tracer and 2D mass-weighted Chemistry increments |

Increment bundles come from two helpers. `Bundle_IncrementMod` lives in the Chemistry repository (`@GEOSchem_GridComp/Shared/Chem_Shared/Bundle_IncrementMod.F90`: `Initialize_IncBundle_init`, `_run`, `Compute_IncBundle`). `MBundle_IncrementMod.F90` (this repository) is its 2D mass-weighted variant. Each reads a table `<bundle>_increments::` from `AGCM.rc`.

## 5. Configuration that changes structure
| Key | File | Effect |
|---|---|---|
| `SYNCUV`, `SYNCTQ` | `AGCM.rc` (default 1.0) | whether winds (`SYNCUV`), temperature and moisture (`SYNCTQ`) are handed forward between children. `SYNCTQ==1` also feeds Chemistry; `SYNCTQ>=1` terminates the Surface, Turbulence and Radiation `T`/`TH` imports, so Physics fills them (section 7). |
| `DOPHYSICS` | `AGCM.rc` | if not 1, `Run` skips all children |
| `DPEDT_PHYS` | `AGCM.rc` (default true) | dry-mass rescale (3.3) |
| `HGT_SURFACE` | `AGCM.rc` (default 50; 0 if `LM==72`) | height of the near-surface level reduction; `PK` and `HGT` are allocated only when this is positive and `SYNC*` is on |
| `SPPT`, `SKEB` | `AGCM.rc` | stochastic physics; need `stoch_module` (`GMAO_stoch`) |
| `SCM_NO_RAD` | `AGCM.rc` | single-column mode, drops the radiation heating term |
| `USE_WAVES`, `USE_SEA_SPRAY`, `USE_OCEANOBIOGEOCHEM` | `AGCM.rc` | extra imports are terminated or connected |
| `CO2_PROVIDER`, `CO2` | `AGCM.rc` | `CO2` is connected to Radiation when the provider is `GOCART` or `RRG` and `CO2: -2` |
| `ATM_CO2` | Surface resource file | `CO2SC001` to Surface `CO2SC` when `ATM_CO2==4` |

Child scheme choices (Moist, Turbulence SHOC and EDMF, GWD) are in the child documents.

## 6. Dependencies
- Libraries: MAPL, `GMAO_mpeu`, `GMAO_stoch`, ESMF (`CMakeLists.txt`).
- Modules: `GEOS_UtilsMod` (`GEOS_Qsat`), `stoch_module` (`sppt_pattern`, `skeb_pattern`), `Bundle_IncrementMod` (Chemistry repository), `MBundle_IncrementMod`.
- Needs the six children's `SetServices` by name; on this line three of the child directories are not built (build status above), which is why the library is bypassed. Chemistry and Radiation come from nested repositories, so a checkout without them does not build Physics.
- Agcm expects from Physics: the five tendencies, `TRADV`, and the exports promoted in 4.2.

## 7. Known traps
1. **Order is in `Run`, not in `AddChild`.** Changing the order of children there changes nothing; changing the call order in `Run` changes the model. Each child's imports are valid only if the child that feeds them has already run (section 3.2).
2. **Synchronization writes into other components' imports.** Because `SYNCUV` and `SYNCTQ` terminate the imports at set-services, a changed value of `SYNCTQ` or `SYNCUV` also changes which connections exist. Test both on and off.
3. **`T` is reconstructed from `S`.** After Turbulence, all `S` change is treated as a change of `T` (stages 2 and 3). Any new source of `S` change that is not a temperature change breaks this.
4. **Stage 3 reuses the variable `BUNDLE`.** The `S` pointer into the Turbulence `TR` bundle at stage 3 is taken from `BUNDLE`, which was set to `TR` at stage 2 (`GEOS_PhysicsGridComp.F90`, Run, `ESMF_StateGet(GIM(TURBL),'TR',BUNDLE)`). Nothing re-assigns `BUNDLE` between stage 2 and stage 3 today. A new `ESMF_StateGet(...,BUNDLE)` inserted in between would silently change what stage 3 writes.
5. **Fixed set of water species in `FILLQ2ZERO` and the dry-mass rescale.** A new water species must be added in both and in the `TRADV` friendly string.
6. **Tracer order in `TR` and `TRI` must match.** Turbulence asserts the counts are equal and pairs fields by name (`GEOSturbulence_GridComp/ARCHITECTURE.md`, section 7).
7. **Pressure-weighted versus plain tendencies.** `DTDT` handed to Dynamics is mass-weighted and re-weighted by the dry-mass step; Moist and GWD `DTDT` are also mass-weighted but `TTMGW` and `DTDT_*` diagnostics are not. Check the unit string before adding terms.
8. **Lagged import.** Moist reads `DTDT_BL` and `DQDT_BL` computed in the previous Physics step (5 and 9 of 3.2).
9. **Comments say "five children"; the code has six** (`GEOS_PhysicsGridComp.F90`, header).

## 8. Decisions
### One run phase with children called in a fixed sequence
- Pattern: a single `Run` calls each child (and its phases) in an order written in code, running couplers by hand between them (`GEOS_PhysicsGridComp.F90`, Run). [Verified]
- Likely reason: Surface and Turbulence exchange fluxes and sensitivities in two phases, and the order of Moist, Surface, Turbulence and Radiation matters for energy consistency; a framework-driven order cannot interleave phases. [Inferred]

### Explicit synchronization of `T`, `S`, winds and moisture (`SYNCUV`, `SYNCTQ`)
- Pattern: Physics copies the state after each child into the next child's imports and terminates the framework connection (`GEOS_PhysicsGridComp.F90`, SetServices `MAPL_TerminateImport`, Run). [Verified]
- Likely reason: the parent treats the physics as a split-step sequence in which each child sees the result of the previous one; the default connectivity would give every child the start-of-step state. [Inferred]
- Question for maintainers: is `SYNCTQ==0` still a supported configuration? It keeps a separate branch for `DQDT_BL` and `DTDT_BL`.

### Dry-mass rescale at the end of Physics
- Pattern: Physics integrates the column water change into a pressure tendency, then rescales all `TRADV` tracers so mixing ratios stay consistent with the new layer mass (`GEOS_PhysicsGridComp.F90`, Run, `DPEDT_PHYS`). [Verified]
- Likely reason: keep Dynamics' dry-air mass constant while Physics adds or removes water. [Inferred]

### Friendly-attribute tracer lists built in Initialize
- Pattern: components mark fields `FRIENDLYTO` names; Physics assembles the bundles for Turbulence, Moist, Dynamics and analysis (4.3). [Verified]
- Likely reason: Chemistry and Moist can add tracers without changes in Physics or Dynamics. [Inferred]

## 9. Observations (noted, not fixed) and open questions
1. **`DQIMST` fetched with the liquid label.** `MAPL_GetPointer(GIM(MOIST), DQIMST, 'DQLDTTRB', ...)` should be `'DQIDTTRB'`; both are then terminated imports (`GEOS_PhysicsGridComp.F90`, Run, "Pass Turbulence QL, QI tendencies to Moist"). Moist never reads either import, so there is no effect today (see the Moist document, Observation 2). The `STATUS` is not checked either.
2. **`FILLQ2ZERO` redistribution is dimensionally inconsistent.** It adds `d_tpw * (MASS(l) / positive_mass)` to the mixing ratio `Q(l)`, where `d_tpw` is a column water mass (mixing ratio times mass). The change in column water is then `d_tpw * sum(MASS^2) / positive_mass`, not `d_tpw` (`GEOS_PhysicsGridComp.F90`, FILLQ2ZERO). A conservative form would add `d_tpw / positive_mass`. The arithmetic is verified; the intent is *Inferred:*. The safety clip afterwards can also lose water.
3. **`DWDT` is declared in Moist but not computed there.** Physics reads it through `alloc=.true.`, so `DWDT` from Physics is zero (`GEOS_PhysicsGridComp.F90`, Run, `WIM`). *Inferred:* intended.
4. **`NEED_TOT` is hard-coded `.TRUE.`**, so the "not needed" branches are dead code.
5. **The stage-3 radiation temperature check only prints** a message for `T > 333` or out-of-range `QV` (`GEOS_PhysicsGridComp.F90`, Run, "SYNTQ for Radiation error").
6. **`OX` is connected to Radiation, but Radiation has no `OX` import.** Physics connects `OX, O3, CH4, N2O, CFC11, CFC12, HCFC22` from `CHEM` to `RAD` (`GEOS_PhysicsGridComp.F90`, SetServices, `MAPL_AddConnectivity` with `DST_ID=RAD, SRC_ID=CHEM`); no Radiation child imports `OX` (see the Radiation document (`@GEOSradiation_GridComp/docs/ARCHITECTURE.md` in the fixture), Observations, which searched `RADIATION`, `SOLAR`, `IRRAD`, `SATSIM`). *Inferred:* the connection therefore has no effect; Radiation imports `O3`. Question for maintainers: is the `OX` connection to Radiation used?
7. **The Turbulence `TKESHOC` and `QT2` fields are added to `TR` directly** from Turbulence's exports, rather than through friendly attributes (4.3). *Unverified:* the reason.

### Open questions for maintainers
- **Conservation in `FILLQ2ZERO`.** The fill adds a mass-weighted amount to each mixing ratio, so column water changes slightly at each of the five call sites (section 9, item 2). Question for maintainers: is this non-conservation intended as a deliberate smoothing?
- **Run order as a contract.** The child run order lives only in code and in comments, not in a named sequence checked anywhere (section 3.2). Question for maintainers: is an implicit, code-only order intended to remain the contract between `SetServices` connectivity and `Run`?
