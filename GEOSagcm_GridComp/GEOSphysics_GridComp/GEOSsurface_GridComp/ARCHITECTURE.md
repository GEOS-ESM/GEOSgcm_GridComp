# Surface gridded component (parent only)

**Audience:** coding agents and developers working on the Surface gridded component and its coupling to Turbulence.
**Describes:** the MAPL3 line, branch `release/MAPL-v3`.
**Verified against:** `GEOSgcm_GridComp` @ `84f0aae3` (checked 2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `f50622d`): this directory is
`src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/`.
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

Paths below are relative to the `GEOSgcm_GridComp` root. `SURF/` abbreviates
`GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/`.
Related documents: [Physics](../ARCHITECTURE.md),
[Turbulence](../GEOSturbulence_GridComp/ARCHITECTURE.md),
[repository architecture](../../../docs/ARCHITECTURE.md), [repository AGENTS.md](../../../AGENTS.md).

> **Build status on this line:** not built, not ported. On `release/MAPL-v3`, `GEOSagcm_GridComp/GEOSphysics_GridComp/CMakeLists.txt` builds `GEOSchem_GridComp`, `GEOSgwd_GridComp` and `GEOSradiation_GridComp`; `GEOSsurface_GridComp` appears only in a commented-out list ("not yet ported to MAPL3 - restore when ported"). The Physics top-level library (`GEOS_PhysicsGridComp.F90`, which adds this component) is present but not compiled ("Temporarily bypass GEOSphysics top-level library"). `GEOSsurface_GridComp/CMakeLists.txt` itself is unchanged in structure: it builds `GEOS_SurfaceGridComp.F90` with the five child directories as subcomponents and `Shared`, and adds `Utils`; nothing in it is reached from the built tree. Of the files this document cites, `GEOSland_GridComp/GEOScatch_GridComp/GEOS_CatchGridComp.F90` (gFTL iterator calls `ftn_begin`, `ftn_end`, `first`) and the `Utils/mk_restarts` sources and `CMakeLists.txt` (gFTL v2: `gFTL2_StringVector`, `gftl-shared-v2`) differ from the MAPL2 line. All other claims are the MAPL2 claims; the code is MAPL2 code, kept for a later port.

## 1. Purpose
Surface (`GEOS_SurfaceGridCompMod`, instance name `SURFACE`) is a light-weight parent that
connects the atmosphere grid to tiled surface components. Its own imports and exports are on the
atmosphere grid; it moves data to and from each child's tile list and runs the children
(`SURF/GEOS_SurfaceGridComp.F90`, module header comment). This document covers the parent only; the children get one line each
(section 2). A *tile* is one piece of the exchange grid: the intersection of an atmosphere grid
cell with one surface type (and, for ocean, an ocean cell). *Inferred:* from the exchange-grid use in `Initialize`; the tile file format is defined by the boundary-condition tools in section 8.

## 2. Place in the hierarchy
Two parents add Surface as child `SURFACE`:
- `PHYSICS`, in the atmosphere (`GEOSagcm_GridComp/GEOSphysics_GridComp/GEOS_PhysicsGridComp.F90`, `SetServices`, `MAPL_AddChild`).
- `DATAATM`, the data atmosphere used when `USE_DATA_ATM4OCN` is true (`GEOSdataatm_GridComp/GEOS_DataAtmGridComp.F90`, `SetServices`, `MAPL_AddChild(NAME='SURFACE')`). Its only child is `SURFACE`; there is no Turbulence. In `Initialize` it replaces Surface's exchange grid by an ocean-only stream (`MAPL_LocStreamCreate(..., MASK=(/MAPL_OCEAN/))`, `MAPL_ExchangeGridSet`), and in `Run` it calls Surface phase 1 and phase 2 itself (`ESMF_GridCompRun(..., PHASE=1)`, `PHASE=2`). *Inferred:* LAKE, LANDICE and LAND then get empty tile lists, because their streams are subsets of that ocean-only stream.

Surface adds its children in this order (`SURF/GEOS_SurfaceGridComp.F90`, `SetServices`, `MAPL_AddChild`); the MAPL2 `MAPL_AddChild` numbers a child `get_num_children() + 1`, so the indexes are OCEAN=1, LAKE=2, LANDICE=3, LAND=4, ROUTE=5 (MAPL2 behavior, `MAPL_Generic.F90`, `AddChildFromMeta`; not re-verified on this line):

| Child name | Directory (under `SURF/`) | One line |
|---|---|---|
| `SALTWATER` | `GEOSsaltwater_GridComp/` | Ocean and salt-water tiles. Own children, in this order: `SEAICETHERMO` (`USE_CICE_Thermo` 0 = SimpleSeaice, 1 = CICE4ColumnPhys, other = SeaiceInterface), `OPENWATER`, optional `OBIOIMPORTS` (`USE_OCEANOBIOGEOCHEM`) (`GEOSsaltwater_GridComp/GEOS_SaltWaterGridComp.F90`, `SetServices`; its comment says sea ice must precede open water because sub-tile indexes depend on it). |
| `LAKE` | `GEOSlake_GridComp/` | Slab lake tiles; no children (`GEOSlake_GridComp/GEOS_LakeGridComp.F90`, module comment). |
| `LANDICE` | `GEOSlandice_GridComp/` | Permanent land ice tiles; optional child `ISSM` (ice sheet model; directory `GEOSlandice_GridComp/GEOSissm_GridComp/`), real only if CMake finds ISSM, otherwise a stub; added only when `DO_ISSM`=1 and built with `HAVE_ISSM` (`GEOSlandice_GridComp/CMakeLists.txt`; `GEOS_LandIceGridComp.F90`, `SetServices`). |
| `LAND` | `GEOSland_GridComp/` | All other land. Children: `VEGDYN` (`GEOSland_GridComp/GEOSvegdyn_GridComp/`), then `CATCH` (`GEOSland_GridComp/GEOScatch_GridComp/`) for `LSM_CHOICE`=1 or `CATCHCN` (`GEOSland_GridComp/GEOScatchCN_GridComp/`) for 2, optional `IGNI` (`GEOSland_GridComp/GEOSigni_GridComp/`, fire danger) (`GEOSland_GridComp/GEOS_LandGridComp.F90`, `SetServices`). With `NUM_CATCH_ENSEMBLES`>1 there are several Catchment children named `ensNNN:CATCH`. |
| `ROUTE` | `GEOSroute_GridComp/` | River routing in Pfafstetter catchment space; no children; always added, inactive unless `RUN_ROUTE`>0 (`GEOSroute_GridComp/GEOS_RouteGridComp.F90`, `Initialize`). |

Also in this directory: `Shared/` (library `GEOS_SurfaceShared`: `StieglitzSnow.F90`, `SurfParams.F90`;
and the installed resource file `GEOS_SurfaceGridComp.rc`) (`Shared/CMakeLists.txt`);
`Utils/` (section 8). `CMakeLists.txt` here builds `GEOS_SurfaceGridComp.F90` with the five child
directories as subcomponents.

`CATCHCN` itself has one child, `CATCHCNCLM40`, loaded by shared-object name (`libGEOScatchCNCLM40_GridComp.so`); it asserts `LSM_CHOICE`==2 (`GEOSland_GridComp/GEOScatchCN_GridComp/GEOS_CatchCNGridComp.F90`, `SetServices`). Land accepts `LSM_CHOICE`=3 but Surface `SetServices` rejects it ("unknown LSM_CHOICE"), so 3 is unreachable from Surface.

## 3. Tiles, tile space and the location stream
(Routine names without a path are in `SURF/GEOS_SurfaceGridComp.F90`.)
- The exchange grid is read once by the top-level GCM from the tile file (`TILING_FILE`, default `tile.data`; `src/Applications/@GEOSgcm_App/AGCM.rc.tmpl`, `TILING_FILE`, sets `tile.bin`)
  with `MAPL_LocStreamCreate` and attached to the atmosphere component by `MAPL_ExchangeGridSet`
  (`GEOS_GcmGridComp.F90`, `Initialize`); under `DATAATM` the stream is cut down to ocean tiles (section 2). Surface itself does not read the tile file; its
  `Initialize` comment says it does (Observation 5).
- A *location stream* (`MAPL_LocStream`) is MAPL's ordered list of tiles, distributed over processors, each tile carrying its type, area weight and atmosphere-grid cell
  (fixture path `src/Shared/@MAPL/base/MAPL_LocStreamMod.F90`, `MAPL_LocStreamType`). "Tile space" means arrays of length `NT_LOCAL`
  indexed by position in a stream.
- Surface takes the exchange grid (`MAPL_Get(..., ExchangeGrid=)`), makes it its own stream, then makes one
  child stream per child as a subset by tile type (`MAPL_LocStreamCreate(..., MASK=)`) and stores
  it in the child's MAPL state with `MAPL_Set(..., LOCSTREAM=)` (`Initialize`).
  Type masks: SALTWATER=`MAPL_OCEAN`, LAKE=`MAPL_LAKE`, LANDICE=`MAPL_LANDICE`, LAND=`MAPL_LAND` (`SetServices`, `CHILD_MASK`). The four streams
  do not overlap and together cover all tiles (module header comment).
- ROUTE gets no stream of its own from Surface: it receives LAND's stream, and `XFORM_IN(ROUTE)`/`XFORM_OUT(ROUTE)` are built from it
  (the out-transform masks all tiles off: ROUTE has no grid exports) (`Initialize`, `I == LAND` block). ROUTE then replaces its own stream by a
  Pfafstetter-catchment stream (`GEOSroute_GridComp/GEOS_RouteGridComp.F90`, `Initialize`, `MAPL_Set(... locstream=)`).
- Two kinds of transform, both `MAPL_LocStreamTransform`:
  - grid <-> tile on Surface's stream: tile to grid is an area-weighted average of the tile values (tiles with `MAPL_UNDEF` skipped; cells with no tile get `MAPL_UNDEF`; optional `SAMPLE` takes the largest-weight tile instead); grid to tile copies the cell value, with optional bilinear `INTERP` (weights `D` documented as "Bilinear weights" in `MAPL_LocStreamType`) (`src/Shared/@MAPL/base/MAPL_LocStreamMod.F90`, `MAPL_LocStreamTransformT2G`; `MAPL_LocStreamTransformG2T`).
  - Surface-stream tile <-> child-stream tile, via `MAPL_LocStreamXFORM` objects `XFORM_IN(child)` and `XFORM_OUT(child)` created in `Initialize` and kept in the private internal state `SURF_state` (`T_SURFACE_STATE`). Helpers `FILLIN_TILE` / `FILLOUT_TILE` apply them by field name and silently skip fields a child does not have.
- Static exports `FRLAND`, `FRLAKE`, `FRLANDICE`, `FROCEAN` are filled in `Initialize` with `MAPL_LocStreamFracArea`.

## 4. Entry points
All in `SURF/GEOS_SurfaceGridComp.F90`.
`SetServices`: registers `Initialize`, `Run1`, `Run2` (two `ESMF_METHOD_RUN` entry points, so phases 1 and 2), adds imports, exports and internal specs,
adds the five children, and calls `MAPL_GenericSetServices` (`SetServices`).
Surface registers no Finalize of its own (no `ESMF_METHOD_FINALIZE` entry point in the file); the MAPL generic default applies.

`Initialize` (`Initialize`): create streams (section 3); `MAPL_GenericInitialize` (initializes children, so child streams must exist first);
fill `FR*` exports; create transforms; mark those four exports `MAPL_InitialRestart` so connected consumers do not overwrite them by default;
`SurfParams_init(LAND_PARAMS, LSM_CHOICE)`; if `ROUTING_FILE` is set, read the instantaneous runoff-to-outlet table (`InitializeRiverRoutingTeleport`) and optionally `DISCHARGE_ADJUST_FILE`.

`Run1` (called as phase 1): per child except ROUTE, in order OCEAN, LAKE, LANDICE, LAND, internal routine `DOCDS`:
fill the child's import tiles from grid fields (`FILLIN_TILE`), `ESMF_GridCompRun(PHASE=1)`, copy the child's tile exports back (`FILLOUT_TILE`).
Then aggregate to the grid (section 6) the exchange coefficients and the effective surface values (section 5) (`Run1`, `DOCDS`).

`Run2` (phase 2): see section 5. Per child, internal routine `DOTYPE`: fill child imports, `ESMF_GridCompRun(PHASE=2)` (`PHASE=1` for ROUTE, which has one phase), copy exports back.
Continental children LAKE, LANDICE, LAND, ROUTE run first, in that order; river discharge is then built; OCEAN runs last (`Run2`, loops `if (I == OCEAN) cycle`, `call DOTYPE(TYPE=OCEAN)`).

Who calls the phases: Physics runs, per step, GWD, MOIST, `SURFACE` phase 1, CHEMISTRY phase 1, `TURBULENCE` phase 1, `SURFACE` phase 2, `TURBULENCE` phase 2, CHEMISTRY phase 2, RADIATION
(`GEOSagcm_GridComp/GEOSphysics_GridComp/GEOS_PhysicsGridComp.F90`, `Run`). Under `DATAATM` the same two phases are called without Turbulence in between (section 2).

**Routines in the file** (about 11,000 lines): `SetServices` (contains `OBIO_setServices`); `Initialize`; `InitializeRiverRoutingTeleport` (contains `FindTile`, `Tile2Index`); `Run1` (contains `DOCDS`); `Run2` (contains `DOTYPE`); helpers `MKTILE` (`MKTILE_1D`, `MKTILE_UNGRIDDED`), `FILLIN_TILE` (`FILLIN_TILE1D`, `FILLIN_TILE2D`), `FILLOUT_TILE` (`FILLOUT_TILE1D`, `FILLOUT_UNGRIDDED`); `RouteRunoffTeleport`; `OBIO_fillExports`.

**To add a Surface export fed by a child** (follow an existing field such as `Z0H` or `HLATN`): (1) `MAPL_AddExportSpec` in `SetServices`; (2) declare the grid pointer and tile array in `Run1` or `Run2`, and `MAPL_GetPointer(EXPORT, ...)`; (3) `MKTILE(grid, tile, NT)`; (4) in `DOCDS` or `DOTYPE`, `MAPL_GetPointer(GEX(type), dum, name, ALLOC=associated(tile), notFoundOK=.true.)`; (5) after `ESMF_GridCompRun`, `FILLOUT_TILE(GEX(type), name, tile, XFORM)`; (6) after the child loop, `MAPL_LocStreamTransform(LOCSTREAM, grid, tile)`; (7) deallocate the tile array in the clean-up block; (8) add the export to each child that should produce it. For an import to a child use `FILLIN_TILE` in the same routines.

## 5. Coupling with Turbulence: the split implicit flux solve
(`Run1`, `Run2`: `SURF/GEOS_SurfaceGridComp.F90`; `DIFFUSE`, `UPDATE`: `GEOSturbulence_GridComp/GEOS_TurbulenceGridComp.F90`; connectivities: `GEOSagcm_GridComp/GEOSphysics_GridComp/GEOS_PhysicsGridComp.F90`, `SetServices`.)
Surface and Turbulence together solve a semi-implicit surface flux. The solve is split over the two phases of each:

1. **Surface phase 1** (`Run1`): children compute exchange coefficients `CT`, `CQ`, `CM`, `CN` and effective surface values (temperature-like `TH`, `QH`, `UH`, `VH`) per tile
   (`GEOS_SurfaceGridComp.F90`, `Run1`, `DOCDS` getting `CHT CQT CMT CNT TH QH UH VH`). Surface multiplies tile values by the coefficient, averages to the grid,
   and divides by the averaged coefficient: so each grid "HAT" value is a coefficient-weighted average of tile values
   (`Run1`, `THTILE = THTILE*CTTILE` ... `TH = TH/CT`). Results live in Surface's internal state (`CT CQ CM CN THAT QHAT UHAT VHAT TS QS RHOS D0`),
   and the export `SHAT = cp*THAT + PHIS` (`Run1`).
2. **Turbulence phase 1** (`DIFFUSE` in Turbulence `RUN1`): with surface values `SHAT UHAT VHAT QHAT` and `CT CM CQ` as imports, it solves the column and exports
   a preliminary surface flux with the surface value held fixed (`SFLX UFLX VFLX QFLX`) and its sensitivity to the surface value (`SDFL UDFL VDFL QDFL`):
   `SF = CX*(SG - SX(LM))`, `SDF = CX*(1-DX(LM))` (`GEOSturbulence_GridComp/GEOS_TurbulenceGridComp.F90`, `DIFFUSE`).
3. **Physics connectivity** moves them (`GEOS_PhysicsGridComp.F90`, `SetServices`, `MAPL_AddConnectivity`):
   - Surface to Turbulence: `CT CM CQ BSTAR USTAR Z0 Z0H FRLAND EVAP SH`; `SHAT UHAT VHAT QHAT`; `DELSS->SDEL`, `DELUS->UDEL`, `DELVS->VDEL`, `DELQS->QDEL`.
   - Turbulence to Surface: `SFLX->SH`, `UFLX->TAUX`, `VFLX->TAUY`, `QFLX->EVAP`, `SDFL->DSH`, `UDFL->DFU`, `VDFL->DFV`, `QDFL->DEVAP`.
4. **Surface phase 2** (`Run2`): puts the fluxes and sensitivities on tiles (`SH EVAP TAUX TAUY DSH DEVAP`) and runs the children, which update their state (`DELTS`, `DELQS`, ...)
   and return final fluxes. Surface averages `DELTS`, `DELQS` to the grid as `DTS`, `DQS`, updates `TH += DTS`, `QH += DQS`, and exports
   `DELTS=DTS`, `DELQS=DQS`, `DELSS=cp*DTS`, `DELUS=0`, `DELVS=0`; final fluxes are `EVAP + DEVAP*DQS` and `SH + DSH*DTS*cp` (exports `EVAP`, `SH`, `TAUX`, `TAUY`, with sign conventions in `Run2`).
5. **Turbulence phase 2** (`UPDATE` in `RUN2`): adds the surface change into the diffused column, `SX = S + DKX*DSG`, and corrects the surface flux, `SF = SF + DSG*SDF`.

Consequences for an agent:
- Surface owns the coefficients, Turbulence the flux solve; neither phase pair may be reordered or skipped. Physics sets the order.
- Momentum has no implicit feedback: `DELUS` and `DELVS` are always zero (`Run2`).
- Turbulence phase 1 reads Surface's exports `SH` and `EVAP` (for SHOC and EDMF) before Surface phase 2 of the same step has refreshed them. *Inferred:* they are therefore the previous step's final values.
- `SYNCTQ` / `SYNCUV` (Physics resource, default 1.0): Physics writes the freshest `TA`, `QA`, `UA`, `VA`, `SPEED` straight into Surface's import state before each Surface phase
  (from MOIST before phase 1, from TURBULENCE phase 1 before phase 2, interpolated to `HGT_SURFACE` m when that is > 0), and terminates the corresponding connectivities (`GEOS_PhysicsGridComp.F90`, `Run`; `SetServices`, `MAPL_TerminateImport`).
  Unwiring the import in Surface's `SetServices` is therefore not enough to change which values it sees.
- Tile-resolved (not grid-resolved) stress inputs are made in `Run2` from `UUA = -TAUX/CM + UH` (`INTERPOLATE_ATMTAU`); `INTERPOLATE_ATMLM` applies `INTERP` in `Run1` to `PS TA QA DZ UU` and the winds (not winds only).
- `Run1` also forms the Monin-Obukhov diagnostics when `CHOOSEMOSFC`=1 (children return `MOT2M` ... `MOV50M`), otherwise `Run2` forms 2, 10, 50 m values from the grid fluxes.

## 6. Aggregation of child exports to the atmosphere grid
(Citations in `SURF/GEOS_SurfaceGridComp.F90`.)
- Pattern: for each exported field, allocate tile array `*TILE(NT)` only if the grid export is requested (`MKTILE`, filled with `MAPL_UNDEF`),
  ask the child to allocate the export (`MAPL_GetPointer(GEX(type), ..., ALLOC=associated(tile), notFoundOK=.true.)`), run the child, `FILLOUT_TILE` into the Surface-stream array
  (only tiles of that child's type are filled by `XFORM_OUT(MASK_OUT = TILETYPE == CHILD_MASK)`), repeat for all children, then
  one `MAPL_LocStreamTransform(LOCSTREAM, gridfield, TILE)` area-weights over all types (`Run1`; `Run2`, "Transform other exports from exchange grid to agcm grid").
- `FRACI` / `OFRACI` (ice fraction) and `LWI` (land/water/ice flag, derived from `FROCEAN`, `FRLAKE`, `FRI`, `TS`) are built in `Run2`.
- Land fields exist only for land-like children; unused ones stay `MAPL_UNDEF` on other tiles, and the area average skips them. For example `LAI`, `GRN`, `ROOTL` are only produced by LAND (`Run1`, `notFoundOK`).
- Precipitation seen by children (`PCU PLS SNO ICE FRZR`) are Surface exports. They equal the MOIST imports unless `PRECIP_FILE` or `PRECIP_FILE_CLIMSCALE` is set; the two options are mutually exclusive (`Run2`).
- Radiation: normalized beam and diffuse fluxes (`DRPARN` ...) from RADIATION are multiplied by insolation and solar constant in `Run2` before being given to children.
- Export count is large: about 296 `MAPL_AddExportSpec` and 43 `MAPL_AddImportSpec` calls in `SetServices`; do not copy them here, read `SetServices`.

### River runoff and discharge
- LAND export `RUNOFF` (tile space) goes to ROUTE import `RUNOFF` by `MAPL_AddConnectivity(SRC_ID=LAND, DST_ID=ROUTE)` inside Surface. This is the only child-to-child connection. Each child's imports are otherwise cut from automatic resolution by `MAPL_TerminateImport(CHILD=...)`; the one exception is SALTWATER when `USE_CICE_Thermo`=2, where only its `SURFSTATE` import is terminated (`SetServices`). So Surface alone moves data between children. *Inferred:* the `RUNOFF` specs of the two components are identical, so MAPL connects by pointer sharing; Surface never calls `MAPL_GenericRunCouplers`.
- Two routing schemes coexist (code comment in `InitializeRiverRoutingTeleport`):
  - "teleport" (per the code comment): runoff of all land tiles is sent instantaneously to ocean outlet tiles (`ROUTING_FILE`, `RouteRunoffTeleport`), producing `DISCHARGE`, then given to SALTWATER as import `DISCHARGE`. Used for coupled ocean (the template sets `ROUTING_FILE` only for coupled runs: `AGCM.rc.tmpl`, `@MODELATM@COUPLED ROUTING_FILE`). `DISCHARGE_ADJUST_FILE` (with a precipitation replacement) scales it by internal `DISCHARGE_ADJUST`.
  - `ROUTE`: hydraulic routing, switched by `RUN_ROUTE` in `GEOS_SurfaceGridComp.rc`; per the same comment it does not feed the ocean yet.
- Requesting the `DISCHARGE` export without `ROUTING_FILE` fails an assertion (`Run2`). With `USE_DATA_ATM4OCN`, `DISCHARGE` is an import transformed from grid to tile instead.

## 7. State and configuration
**State.** Imports, exports and internal specs are declared in `SetServices` of `SURF/GEOS_SurfaceGridComp.F90` (no ACG file here).
Internal spec items `QS TS CT CQ CM CN THAT QHAT UHAT VHAT RHOS D0` (grid) and `DISCHARGE_ADJUST` (tile) pass data from `Run1` to `Run2` and, being `FriendlyTO` the component, are also exported (`SetServices`).
`SURF_state` holds `XFORM_IN/OUT` and the teleport routing table (`T_SURFACE_STATE`).
**Restarts.** Surface's internal state is not checkpointed: the module header says it need not be saved, and `src/Applications/@GEOSgcm_App/AGCM.rc.tmpl` has no `SURFACE_INTERNAL_RESTART_FILE` or `SURFACE_INTERNAL_CHECKPOINT_FILE` key, which MAPL needs before it reads or writes an internal state (MAPL2 (MAPL `release/v2`, `generic/MAPL_Generic.F90`), `INTERNAL_RESTART_FILE:`). The parent's real restart is its import file, `SURFACE_IMPORT_RESTART_FILE` (`surf_import_rst`, `AGCM.rc.tmpl`); imports not marked `MAPL_RestartSkip` are written to it (`SetServices`).
Restart state of the children (keys `*_INTERNAL_RESTART_FILE` in `AGCM.rc.tmpl`): `lake_internal_rst`, `landice_internal_rst`, `catch_internal_rst`, `catchcn_internal_rst`, `route_internal_rst`, `vegdyn.data` (`VEGDYN_INTERNAL_RESTART_FILE`), `openwater_internal_rst`, `seaicethermo_internal_rst`.

**Configuration that changes structure** (read in `SetServices` / `Initialize`):
| Key | File read from | Effect |
|---|---|---|
| `SURFRC` (default `GEOS_SurfaceGridComp.rc`) | AGCM.rc | Names the surface resource file read by Surface, Physics, Land, Route (`SURF/Shared/GEOS_SurfaceGridComp.rc`). Its header says it is used jointly by the GCM and GEOSldas, with defaults edited by the setup scripts. |
| `LSM_CHOICE` (1 Catchment, 2 CatchmentCN-CLM40) | AGCM.rc | Selects `CATCH` or `CATCHCN` child; Surface asserts "unknown LSM_CHOICE" otherwise (`SetServices`). |
| `NUM_CATCH_ENSEMBLES`, `DO_ISSM`, `FIRE_DANGER` | AGCM.rc; AGCM.rc; `GEOS_SurfaceGridComp.rc` | Number of Catchment children (`ensNNN:CATCH`); ISSM child of LANDICE; `IGNI` child of LAND (`GEOS_LandGridComp.F90`, `SetServices`; `GEOS_LandIceGridComp.F90`, `SetServices`). |
| `CHOOSEMOSFC`, `LAND_PARAMS`, `ATM_CO2`, `N_CONST_LAND4SNWALB`, `N_CONST_LANDICE4SNWALB` | `GEOS_SurfaceGridComp.rc` | Surface layer scheme; land parameter set; extra imports/exports (OBIO CO2 bands when `ATM_CO2`=4, aerosol deposition tiles when SWIM constants or OBIO are non-zero, fire index exports) (`SetServices`). |
| `RUN_ROUTE`, `RRM_DT`, `TILING_NC_FILE` | `GEOS_SurfaceGridComp.rc` | ROUTE active or not (`GEOS_RouteGridComp.F90`, `Initialize`). Route also needs `ROUTE_INTERNAL_RESTART_FILE` (AGCM.rc). |
| `USE_CICE_Thermo`, `USE_OCEANOBIOGEOCHEM`, `USE_WAVES`, `USE_DATA_ATM4OCN` | AGCM.rc | Sea-ice scheme; OBIO imports; `CHARNOCK` import; with `USE_DATA_ATM4OCN` the parent becomes `DATAATM` (section 2) and `DISCHARGE` becomes an import (`SetServices`, `Run2`). With `USE_CICE_Thermo`=2 only `SURFSTATE` is terminated for SALTWATER; for other values all SALTWATER imports are (`SetServices`, `MAPL_TerminateImport`). *Inferred:* SALTWATER's other imports then stay open as Surface imports. |
| `ROUTING_FILE`, `DISCHARGE_ADJUST_FILE`, `PRECIP_FILE`, `PRECIP_FILE_CLIMSCALE` | AGCM.rc | See section 6. |
| `INTERPOLATE_ATMLM`, `INTERPOLATE_ATMTAU` | AGCM.rc | Bilinear rather than cell-constant grid-to-tile transfer: `ATMLM` for `PS TA QA DZ UU` and the winds in `Run1`; `ATMTAU` for the effective winds in `Run2`. The template sets 0 and 1 only for coupled runs (`AGCM.rc.tmpl`, `@COUPLED INTERPOLATE_ATMLM`). |

Build switches: `AQUA_PLANET` (preprocessor) reduces the children to SALTWATER only (`NUM_CHILDREN`=1); `HAVE_ISSM` (LandIce).

## 8. Utils (one line each)
- `Utils/Raster/makebcs/` : `make_bcs.py` and Fortran programs that make tile files (`.til`), rasters and land parameter files (boundary conditions) for a combination of atmosphere, ocean and land resolutions (`Utils/Raster/makebcs/make_bcs.py`, usage text).
- `Utils/Raster/preproc/` : programs that make *inputs* for `make_bcs` (soil, outlets, routing model, topography, albedo, lake, ISSM) (`Utils/Raster/preproc/README`).
- `Utils/Raster/misc/` : IDL (`.pro`) and shell plotting scripts and one Fortran parameter-check program (`chk_clsm_params.F90`).
- `Utils/mk_restarts/` : programs that create restarts for all tile components (`Utils/mk_restarts/README`); called by `regrid.pl` in another repository (*Unverified:* from README).

## 9. Nested repositories under this directory
One: `Topo` (NCAR topography tool), at `Utils/Raster/preproc/topography/@ncar_topo`, tag `geos/v2.2.1` (`components.yaml`, entry `Topo`); *Unverified:* used by the topography preprocessing.

## 10. Dependencies
`GEOS_SurfaceGridComp.F90` uses `ESMF`, `MAPL` (including `MAPL_LocStream*`), `GEOS_UtilsMod`, `StieglitzSnow` (aerosol-in-snow bin counts), `SurfParams`, `CATCH_CONSTANTS`, `m_mpif90`
(module `use` list). CMake: `GEOS_SurfaceShared MAPL GMAO_mpeu ESMF::ESMF` (`SURF/CMakeLists.txt`). Child libraries are `SUBCOMPONENTS`; a missing one becomes a stub (`@cmake/esma_support/esma_add_library.cmake` at the GEOSgcm fixture root).

## 11. Known traps
(Citations in `SURF/GEOS_SurfaceGridComp.F90` unless a path is given.)
- Child index order is fixed by `MAPL_AddChild` order, and loops in `Run1`/`Run2` rely on it (ocean first in `Run1`, last in `Run2`). Adding a child means updating `NUM_CHILDREN`, `CHILD_MASK`, `TerminateImport` and both run loops (timers loop over `GCNames`). SALTWATER's own children are order-sensitive too (section 2).
- Surface calls children with `ESMF_GridCompRun`, not the MAPL generic run, and connects them by hand. Ordinary MAPL import/export resolution among the children is turned off by `MAPL_TerminateImport` for each child (only `SURFSTATE` for SALTWATER when `USE_CICE_Thermo`=2); new child imports must be filled in `DOCDS`/`DOTYPE` by name.
- `FILLIN_TILE`/`FILLOUT_TILE` skip a missing field without error; a misspelled name makes a silent zero or `MAPL_UNDEF`.
- Run1 internal values are "weighted by coefficient", not area; area-weighted averages are only used for ordinary exports.
- The Run1 to Run2 hand-off is through Surface's internal state; between the two phases Chemistry and Turbulence run, and Physics may rewrite `TA`, `QA`, `UA`, `VA`, `SPEED` (section 5).
- Surface has an import and an export with the same short name: `SH`, `EVAP`, `TAUX`, `TAUY`. The import is Turbulence's preliminary flux (`SFLX->SH`, `UFLX->TAUX`, ...); the export is the final flux, and `TAUX`/`TAUY` exports have the sign flipped (`TAUXO = -TAUX`, `Run2`). Physics connects Surface's exports `SH`, `EVAP` to Turbulence and Moist (`GEOS_PhysicsGridComp.F90`, `SetServices`).
- Turbulence looks up the surface fields by built name, `trim(NAME)//'HAT'`, `'FLX'`, `'DFL'`, `'DEL'` for `NAME` in S, U, V, Q (`GEOS_TurbulenceGridComp.F90`, `DIFFUSE`, `UPDATE`); grepping for `SDEL` or `SFLX` finds only the spec.
- `AERO_DP` bin counts are fixed in Surface. `Run2` counts the `AERO_DP` fields it finds by name and stops (`STOP 'NUM_DUDP mismatch with AERO_DP'` and the like) if a count differs from the `NUM_*` parameters in `Shared/StieglitzSnow.F90` (DU 5, BC 2, OC 2, SU 1, SS 5 for each of the four suffixes `DP`, `SV`, `WT`, `SD`). A change of the number of aerosol bins in Chemistry therefore needs a change here; see GOCART document (`../@GEOSchem_GridComp/@GOCART/docs/ARCHITECTURE.md` in the fixture), section 5.6.
- ROUTE is always created and runs every step, returning immediately if `RUN_ROUTE`=0.
- `DELUS`/`DELVS` are zero by construction (no implicit momentum feedback).

## 12. Decisions
### Two-phase surface with implicit coupling to Turbulence
- Pattern: coefficients and effective surface values first, fluxes and state updates after Turbulence's solve; the surface change is fed back to Turbulence (`SURF/GEOS_SurfaceGridComp.F90`, header comment, `Run1`, `Run2`).
- Likely reason: stable, large-step coupling without iterating. [Inferred]

### Children run on tile subsets with Surface as the only mediator
- Pattern: separate non-overlapping streams per child, all data movement by Surface, `MAPL_TerminateImport` on every child (`SetServices`, `Initialize`).
- Likely reason: the children share no state or grid; the code comment says MAPL's automatic resolution would not be wanted. [Inferred]
- Question for maintainers: the module header says the stored transforms "should be done more elegantly once ESMF has exchange grid support"; is that still the plan?

### Surface does not own the tile file
- Pattern: exchange grid created by the GCM parent and inherited (`GEOS_GcmGridComp.F90`, `Initialize`; Surface `Initialize`, `MAPL_Get(ExchangeGrid=)`).
- Likely reason: atmosphere and ocean must share one tile file. [Inferred]

### River routing added beside teleport routing
- Pattern: `ROUTE` always a child, active by resource; legacy teleport kept in Surface (`InitializeRiverRoutingTeleport`, comment block).
- Likely reason: keep coupled-ocean behavior unchanged while ROUTE is introduced. [Inferred]
- Question for maintainers: the comment says the two will be unified so ROUTE feeds the ocean; is that planned for this component?

### Physics overwrites some Surface imports
- Pattern: with `SYNCTQ`/`SYNCUV` >= 1, Physics writes `TA QA UA VA SPEED` into Surface's import state before each phase and terminates the connectivities (`GEOS_PhysicsGridComp.F90`, `Run`, `SetServices`). [Verified]
- Likely reason: Surface should see the state left by the latest Moist/Turbulence stage. [Inferred]

## 13. Observations (noted, not fixed)
1. Under `AQUA_PLANET`, `LAKE`, `LANDICE`, `LAND`, `ROUTE` are never assigned (assignments are inside `#ifndef AQUA_PLANET`), yet `Initialize`, `Run1` and `Run2` compare `I == ROUTE` and `I == LAND`. *Unverified:* whether this matters; no CMake definition of `AQUA_PLANET` was found in the fixture.
2. Imports `DFU` and `DFV` are declared and connected from Turbulence (`UDFL`, `VDFL`) but never read in `Run2`.
3. `GEOSland_GridComp/CMakeLists.txt` lists `GEOSlana_GridComp`, which is not in this repository; `esma_add_library` builds a stub for it.
4. Surface's `RUN2` has an empty description block; section 5 here is read from code.
5. The `Initialize` description comment says it "reads the tiling file"; the code takes the exchange grid from MAPL (`MAPL_Get(..., ExchangeGrid=)`).
6. Several `_ASSERT` calls use the message `'needs informative message'`.
