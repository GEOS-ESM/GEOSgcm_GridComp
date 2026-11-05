# GEOSturbulence_GridComp (Turbulence): architecture

**Audience:** coding agents and developers working on the Turbulence gridded component (vertical diffusion, boundary-layer mixing, shallow-convection mass flux, orographic form drag).
**Describes:** the MAPL2 line, branch `develop`.
**Verified against:** `GEOSgcm_GridComp` @ `45e1f86f` (2026-11-04).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `0f2e0fb`; repository `GEOSgcm_GridComp`, mounted at `src/Components/@GEOSgcm_GridComp`): `src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSturbulence_GridComp`.
Other fixtures (for example GEOSadas, GEOSldas) may place this repository elsewhere; paths inside this repository do not change.

Paths below are relative to `GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSturbulence_GridComp/` unless they start with `GEOSagcm_GridComp/`.
Parent: [Physics](../ARCHITECTURE.md). Siblings: [Moist](../GEOSmoist_GridComp/ARCHITECTURE.md), [GWD](../GEOSgwd_GridComp/ARCHITECTURE.md), [Surface](../GEOSsurface_GridComp/ARCHITECTURE.md). Repository overview: [../../../docs/ARCHITECTURE.md](../../../docs/ARCHITECTURE.md).

## 1. Purpose
Turbulence mixes heat, moisture, momentum and every friendly tracer in the vertical. It computes diffusivities (a Louis-type scheme, optionally the Lock entrainment scheme, or the Simplified Higher-Order Closure, SHOC), optionally adds a shallow-convection mass flux (EDMF, eddy-diffusivity mass-flux), adds orographic form drag (the Beljaars scheme), and solves a tridiagonal implicit system for each field (`GEOS_TurbulenceGridComp.F90`, header comment, REFRESH, DIFFUSE).
It is the vertical-diffusion half of the surface-to-atmosphere exchange. Surface supplies the lower boundary condition; Turbulence supplies fluxes and their sensitivities back to Surface, which in turn feeds back an increment (section 3.2).
The implementation is one large file, `GEOS_TurbulenceGridComp.F90` (about 7000 lines), plus `edmf.F90`, `shoc.F90`, `LockEntrain.F90` and `scm_surface.F90`.

## 2. Place in the hierarchy
- Parent is Physics, which adds Turbulence as a leaf child `TURBL` (`GEOSagcm_GridComp/GEOSphysics_GridComp/GEOS_PhysicsGridComp.F90`, SetServices).
- No children. SHOC, EDMF and Lock are Fortran modules called directly (`GEOS_TurbulenceGridComp.F90`, `use shoc`, `use edmf_mod`, `use LockEntrain`).
- `CMakeLists.txt` builds one library with all four extra sources.

## 3. Entry points
`SetServices` registers two run phases and nothing else: `MAPL_GridCompSetEntryPoint(..., ESMF_METHOD_RUN, Run1, ...)` and `Run2` (`GEOS_TurbulenceGridComp.F90`, SetServices). Initialize and Finalize are the MAPL generic defaults; the header says so.

### 3.1 SetServices
- Reads `DEBUG_TRB`, `USE_WAVES`, `USE_SEA_SPRAY`, `SCM_SL` (idealized single-column surface layer) and `EDMF_NUMUP` (default 10, the number of plumes; sizes the ungridded dimension of the EDMF plume exports).
- Reads `TURBULENCE_DO_SHOC:` with `ESMF_ConfigGetAttribute`. If non-zero, the friendly string of `PDF_A`, `TKESHOC`, `QT2` and `QT3` becomes `'DYNAMICS:TURBULENCE'` (so Dynamics advects them); otherwise it is the component name only (`GEOS_TurbulenceGridComp.F90`, SetServices, `FRIENDLIES_SHOC`).
- Declares imports, internal specs, exports and timers, then `MAPL_GenericSetServices`.

### 3.2 Run1 (`GEOS_TurbulenceGridComp.F90`, RUN1) and Run2 (RUN2 → UPDATE)
Physics calls the phases around Surface (`GEOSagcm_GridComp/GEOSphysics_GridComp/GEOS_PhysicsGridComp.F90`, Run): Surface phase 1, Chemistry phase 1, **Turbulence phase 1**, then Surface phase 2, **Turbulence phase 2**.

**Run1**
1. Copy sea-spray fluxes to exports when `USE_WAVES` and `USE_SEA_SPRAY` are on.
2. Fetch internal-state pointers (matrix diagonals `AK*`, `BK*`, `CK*`, `DK*`, the SHOC/EDMF terms, `SINC`).
3. If the component's `RUNALARM` rings, turn the ringer off and call `REFRESH`, which recomputes the diffusivities and the matrices.
4. Always call `DIFFUSE`.

`REFRESH` order, each step optional by resource:
1. Boundary-layer height options.
2. EDMF mass flux (`RUN_EDMF`) if `EDMF_DOMF` is non-zero.
3. SHOC (`RUN_SHOC`) if `TURBULENCE_DO_SHOC` is non-zero; this also sets `LOCK_ON=0`.
4. Louis diffusivities: `LOUIS_KS` if `JASON_LOUIS`, else `LOUIS_KS_OPTIMIZED` (only when SHOC is off).
5. Lock entrain scheme (`ENTRAIN`) if `LOCK_ON==1`, and the TKE diagnostics.
6. Beljaars orographic form drag if `C_B` is non-zero; it writes `BKV`, `BKUU` and `FKV`.
7. LU factorization of the tridiagonal systems (`VTRILU`) for `S`, `Q`, `V` and the EDMF variants; `VTRISOLVESURF` gives `DKS`, `DKQ`, `DKV` (the response to a unit change of the surface value).

`DIFFUSE` loops over the fields of the import bundle `TR` and, for each, solves the implicit system with the surface value from Surface held fixed:
- Per field it reads the attributes `DiffuseLike` (`'U'`, `'S'` or `'Q'`, default `'Q'`), `FriendlyToTURBULENCE` (default false) and `WeightedTendency` (default false). The coefficient set is chosen from `DiffuseLike`; the exchange coefficient is `CM` (momentum), `CT` (heat) or `CQ` (moisture).
- Only fields named `S`, `U`, `V`, `Q` have a surface value (`<N>HAT` import), a flux (`<N>FLX` export) and a sensitivity (`<N>DFL` export). All other tracers get zero surface flux.
- The names `S`, `Q`, `QLLS`, `QILS`, `U`, `V` use the EDMF/SHOC matrices (`AKSS`/`YS`, `AKQQ`/`YQV`/`YQL`/`YQI`, `AKUU`/`YU`/`YV`). Other tracers get optional EDMF tracer transport (`MFPARAMS%DOTRACERS`) and then `VTRISOLVE`. `QLLS` and `QILS` skip the mass-flux part when SHOC is on.
- The flux is `SF = CX*(SG - SX(LM))`. The tendency is `(SX - S)/DT`, multiplied by layer thickness `DP` when `WeightedTendency` is true, and is written into the field `<NAME>IT` of the export bundle `TRI`.
- If the field is friendly to Turbulence, `S = SX` is written in place, so later stages see the diffused value.
- Exports `UAFDIFFUSE`, `VAFDIFFUSE`, `SAFDIFFUSE`, `QAFDIFFUSE` hold the state after diffusion. Physics uses them for stage-2 synchronization (Physics document, 3.2).

**Run2** (`UPDATE`)
1. Takes the surface increments from Surface's second phase (imports `SDEL`, `UDEL`, `VDEL`, `QDEL`) and updates `SX = S + DKX*DSG`, where `DKX` is the sensitivity stored in Run1.
2. Adds `DKX*DSG/DT` to the `TRI` tendency (times `DP` if weighted), and updates the flux `SF += DSG*SDF`.
3. Computes dissipation exports `INTDIS` (shear, from `EKV`), `SRFDIS`, `TOPDIS` (from `FKV`). Each is capped by `CAP_INTDIS` and `CAP_TOPDIS` (default 1 K per step, 0 if `LM==72`). With SHOC, `INTDIS = -TKEDISS*DP/CP`. Kinetic-energy diagnostics `KETRB`, `KESRF`, `KETOP`, `KEINT`.
4. Optional SHVC parameterization when `SHVC_EFFECT>0` (default 0; the source does not expand the acronym), which uses `VARFLT`; `S` is processed first (`GEOS_TurbulenceGridComp.F90`, UPDATE, `SHVC_EFFECT:`).
5. Friendly fields get `S = SX`; `SAFUPDATE` is the dry static energy after the update. Physics uses it for stage-3 synchronization.
6. Flux diagnostics `QTFLXTRB`, `SLFLXTRB`, `UFLXTRB`, `VFLXTRB`, `WQT`, `WSL`, the exports `QT` and `SL`, and `DPDTTRB = g*SF` for `Q` at the surface (`GEOS_TurbulenceGridComp.F90`, UPDATE).
7. Dissipation heating is **not** added to the temperature here. It is exported, and Physics adds it (`INTDIS + TOPDIS`, Physics document, 3.3).

## 4. State
### 4.1 Imports
Declared in `GEOS_TurbulenceGridComp.F90`, SetServices. Restart column: *skip* means `RESTART=MAPL_RestartSkip`.

| Group | Names | Source (via Physics) |
|---|---|---|
| Tracer bundle (skip) | `TR` | filled by Physics Initialize (Physics document, 4.3) |
| State (skip) | `PLE`, `ZLE`, `T`, `TH`, `U`, `V`, `AREA`, `PREF`, `VARFLT` | Dynamics; `T`, `TH`, `U`, `V` overwritten by Physics when `SYNCTQ`/`SYNCUV` are on |
| Moist (skip) | `QV`, `QLTOT`, `QITOT`, `QRTOT`, `QSTOT`, `QGTOT`, `FCLD`, `EIS` | Moist |
| Moist (restarted) | `WTHV2`, `WQT_DC` | Moist |
| Surface (skip) | `CT`, `CM`, `CQ`, `BSTAR`, `USTAR`, `Z0`, `Z0H`, `FRLAND`, `EVAP`, `SHAT`, `UHAT`, `VHAT`, `QHAT`, `SDEL`, `UDEL`, `VDEL`, `QDEL` | Surface |
| Surface (restarted) | `SH`, `SHFX_SPRAY`, `LHFX_SPRAY` (optional restart) | Surface |
| Radiation (restarted) | `RADLW`, `RADLWC` | Radiation |
| Other (restarted) | `PHIS`, `OMEGA`, `SHOBS`, `LHOBS` | Agcm, observation inputs for single-column runs |

### 4.2 Exports
- `TRI` bundle of tendencies `<NAME>IT`, one per `TR` field (Physics fills it, section 3.2).
- To Surface: `SFLX`, `UFLX`, `VFLX`, `QFLX` (fluxes) and `SDFL`, `UDFL`, `VDFL`, `QDFL` (sensitivities) (`GEOS_PhysicsGridComp.F90`, SetServices, "Surface Imports").
- To Moist: `KH`, `KM`, `RI`, `KPBL`, `KPBL_SC`, `TKE`, `TKESHOC`, `PDF_A`, `SL2`, `SL3`, `W2`, `W3`, `SLQT`, `WQT`, `WSL`, `QT2`, `QT3`, `ALH`. To Chemistry: `ZPBL`, `PPBL`.
- To Physics: `INTDIS`, `TOPDIS`, `SRFDIS`, `DPDTTRB`, `UAFDIFFUSE`, `VAFDIFFUSE`, `SAFDIFFUSE`, `QAFDIFFUSE`, `SAFUPDATE`.
- EDMF and SHOC diagnostics (`EDMF_MF`, `EDMF_ENTR`, many others). Plume exports have `EDMF_NUMUP` levels in the ungridded dimension.

### 4.3 Internal state and restarts
- **Restarted:** `ZPBL`, `ZPBL_SC`, `PDF_A`, `TKESHOC`, `TKH`, `QT2`, `QT3`, `MFAW`, `SLFLXMF`, `QTFLXMF`, and `cu_scm`, `ct_scm`, `ssurf_scm`, `qsurf_scm` (the last four are declared only when `SCM_SL` is non-zero).
- **Not restarted** (skip): the matrix coefficients `AKS`..`DKS`, `AKQ`..`DKQ`, `AKV`..`FKV`, the SHOC/EDMF matrices `AKSS`.., `AKQQ`.., `AKUU`.., right-hand sides `YS`, `YQV`, `YQL`, `YQI`, `YU`, `YV`, `SINC` and `TKEDISS`.
- The unrestarted coefficients are recomputed only when `RUNALARM` rings (3.2). *Inferred:* after a restart in the middle of a diffusivity interval, `DIFFUSE` uses coefficients that are zero or stale until the next ring. *Unverified:* whether the alarm rings on the first step after a restart.
- Restart file names are set in `AGCM.rc`: `MODELATMTURBULENCE_INTERNAL_RESTART_FILE: turb_internal_rst` and `MODELATMTURBULENCE_IMPORT_RESTART_FILE: turb_import_rst` (`src/Applications/@GEOSgcm_App/AGCM.rc.tmpl` in the fixture).

## 5. Configuration that changes structure
Resources are read with the component prefix `TURBULENCE_` as well as bare, in that order of preference (`src/Shared/@MAPL/base/MAPL_Resource.F90` in the fixture, label lookup with prefix). Many defaults depend on the number of levels: `LM==72` selects an older tuning (`JASON_*` switches, `PBLHT_OPTION`, `SMTH_HGT`, `C_B`, `CAP_*`) (`GEOS_TurbulenceGridComp.F90`, REFRESH).

| Key | Effect |
|---|---|
| `TURBULENCE_DO_SHOC` | selects SHOC; turns off Lock and the optimized Louis path; changes which fields Dynamics advects (3.1) |
| `EDMF_DOMF`, `EDMF_NUMUP`, `EDMF_*` | shallow-convection mass flux on or off; plume count and tuning |
| `LOCK_ON`, `JASON_LOUIS`, `PBLHT_OPTION`, `C_B`, `SHVC_EFFECT` | choice of diffusivity scheme, boundary-layer height definition, orographic form drag, SHVC parameterization |
| `DT:` / `RUN_DT:` | interval at which `REFRESH` runs (via `RUNALARM`); `RUN_DT` is the time step used in `DIFFUSE` |
| `USE_WAVES`, `USE_SEA_SPRAY` | add sea-spray flux imports and terms |
| `SCM_SL`, `SCM_SL_FLUX`, `SCM_SH`, `SCM_LH` | idealized surface layer for single-column runs (`scm_surface.F90`) |
| `CAP_INTDIS`, `CAP_TOPDIS` | limits on frictional heating |

The `TURBULENCE_*` and `EDMF_*`, `SHC_*` values are set in `AGCM.rc.tmpl` in the fixture.

## 6. Dependencies
- Libraries: `GEOS_Shared`, MAPL, ESMF, NetCDF Fortran (`CMakeLists.txt`, `DEPENDENCIES`).
- Modules: `LockEntrain`, `shoc`, `edmf_mod` (`run_edmf`, `mfparams`), `scm_surface`, `GEOS_Mod`.
- `_CUDA` preprocessor sections exist in the file. *Unverified:* whether any build still uses them.
- `SOURCE.txt` records a stale origin path (`GEOS-LM_v12-rc3`).
- Physics supplies everything through connectivity; see the Physics document, 3.2 and 4.3.

## 7. Known traps
1. **The `TR` and `TRI` bundles must have the same fields in the same order.** `DIFFUSE` and `UPDATE` index them in parallel and assert equal counts. Physics builds both (Physics document, 4.3). A tracer added to `TR` needs an `<NAME>IT` export in Physics; if the export is missing, the field is empty and the tendency is skipped.
2. **`Q` is called `Q` in `TR` but `QV` in the Moist state.** The tendency name is `QIT` in the bundle and `QVIT` in the Physics export; the special case is in Physics Initialize (`GEOS_PhysicsGridComp.F90`, Initialize, `case ('Q')`). `UPDATE` compares the name with `'QV'` as well, which seems redundant (`GEOS_TurbulenceGridComp.F90`, UPDATE). *Inferred:* a leftover from an earlier naming.
3. **Surface values exist for four names only.** A new tracer that needs a surface flux must be added in `DIFFUSE` and `UPDATE` (the `select case ('S','U','V','Q')`) and in Physics.
4. **Run2 must follow Surface phase 2.** The `*DEL` imports are valid only after Surface has run its second phase; Physics couples them with `MAPL_AddConnectivity` (`GEOS_PhysicsGridComp.F90`, SetServices).
5. **Imports `T`, `TH`, `U`, `V` are overwritten by Physics.** Physics terminates their connectivity and fills them with the state after Moist (stage 1) and after diffusion (stage 2) (Physics document, 3.2). Turbulence reads them unchanged during a step.
6. **Frictional heating is not applied here.** `INTDIS`, `TOPDIS` are exports; if Turbulence runs outside Physics, heating is lost.
7. **Diffusivity coefficients are cached across steps** and recomputed only when the alarm rings. Changing a resource at run time has no effect until then. The surface sensitivities (`DKS`, `DKQ`, `DKV`) are recomputed at the same time.
8. **Run2 uses `RUN_DT:`, and the alarm uses `DT:`** (default `RUN_DT`). If `DT:` is set to a multiple of `RUN_DT`, only `REFRESH` is supercycled (`GEOS_TurbulenceGridComp.F90`, RUN2, `ESMF_ConfigGetAttribute ... RUN_DT`).
9. **`UPDATE` returns failure for an unknown `DiffuseLike`** (`RETURN_(ESMF_FAILURE)`).

## 8. Decisions
### Two run phases around Surface's second phase
- Pattern: Turbulence solves with a fixed surface value (Run1), exports the sensitivity, and corrects after Surface supplies an increment (Run2) (`GEOS_TurbulenceGridComp.F90`, header comment, RUN1, RUN2). [Verified]
- Likely reason: Surface (land and sea ice) depends on the atmospheric fluxes while Turbulence depends on the surface value; the split gives an implicit coupling without iteration. [Inferred]

### Fields to diffuse are discovered from attributes in the `TR` bundle
- Pattern: any field in `TR` is diffused; `DiffuseLike`, `FriendlyToTURBULENCE`, `WeightedTendency` set the treatment (`GEOS_TurbulenceGridComp.F90`, DIFFUSE). [Verified]
- Likely reason: Chemistry and Moist can add tracers with no change in Turbulence. [Inferred]

### Matrices stored in the internal state and refreshed on an alarm
- Pattern: tridiagonal coefficients and the surface sensitivities live in unrestarted internal specs and are rebuilt by `REFRESH` on `RUNALARM` (`GEOS_TurbulenceGridComp.F90`, RUN1, internal specs). [Verified]
- Likely reason: factoring is the expensive step and the diffusivities can be held over several steps, while the solve is cheap and uses current fields. [Inferred]
- Question for maintainers: is the first step after a restart guaranteed to call `REFRESH` (section 4.3)?

### Scheme choice by flags inside `REFRESH`
- Pattern: Louis, Lock, SHOC and EDMF are selected by run-time flags tested in sequence, not by separate components (`GEOS_TurbulenceGridComp.F90`, REFRESH). [Verified]
- Likely reason: the schemes share work arrays and the `DIFFUSE` matrix form. [Inferred]

## 9. Observations (noted, not fixed) and open questions
1. **Stale header comment.** The header describes `MY_STEP` and `ACCUMINT`, but no code uses them (`GEOS_TurbulenceGridComp.F90`, header comment). *Inferred:* the refresh control is now the generic `RUNALARM`.
2. **Descriptions say Initialize and Finalize are inherited from the generic component**, which is still true. The header's Run1 and Run2 descriptions predate the SHOC and EDMF additions.
3. **`QT3` is commented out of the `TR` additions** while `TKESHOC` and `QT2` are included (Physics document, 4.3). *Unverified:* intent.
4. **Redundant `'QV'` test in `UPDATE`** (trap 2).
5. **`SOURCE.txt` refers to an old path.**
6. **Different defaults for `LM==72`.** The tuning for 72 levels is chosen from `LM`, not from a named configuration, so a 72-level run with a different intended tuning needs explicit keys.

### Open questions for maintainers
- **Frictional heating.** Turbulence exports `INTDIS` and `TOPDIS`; Physics adds them to the temperature tendency (trap 6). Question for maintainers: is the split between computing dissipation (Turbulence) and applying it as heating (Physics) intentional?
