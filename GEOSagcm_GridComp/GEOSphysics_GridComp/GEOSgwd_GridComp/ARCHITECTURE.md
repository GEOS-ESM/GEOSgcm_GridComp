# GEOSgwd_GridComp (GWD): architecture

**Audience:** coding agents and developers working on the gravity wave drag (GWD) gridded component (orographic, convective and background gravity-wave drag, Rayleigh friction).
**Describes:** the MAPL2 line, branch `develop`.
**Verified against:** `GEOSgcm_GridComp` @ `45e1f86f` (2026-11-04).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `0f2e0fb`; repository `GEOSgcm_GridComp`, mounted at `src/Components/@GEOSgcm_GridComp`): `src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSgwd_GridComp`.
Other fixtures (for example GEOSadas, GEOSldas) may place this repository elsewhere; paths inside this repository do not change.

Paths below are relative to `GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSgwd_GridComp/` unless they start with `GEOSagcm_GridComp/`.
Parent: [Physics](../ARCHITECTURE.md). Siblings: [Moist](../GEOSmoist_GridComp/ARCHITECTURE.md), [Turbulence](../GEOSturbulence_GridComp/ARCHITECTURE.md), [Surface](../GEOSsurface_GridComp/ARCHITECTURE.md). Repository overview: [../../../docs/ARCHITECTURE.md](../../../docs/ARCHITECTURE.md).

## 1. Purpose
GWD computes the momentum and heat tendencies caused by sub-grid gravity waves that the model does not resolve. It has two independent implementations whose tendencies are added (`GEOS_GwdGridComp.F90`, Gwd_Driver):
- The **GEOS** scheme (`gw_drag.F90`, module `gw_drag`, entry `gw_intr`): an older scheme with a background (non-orographic) spectrum and orographic drag.
- The **NCAR** scheme (`ncar_gwd/`, module `gw_drag_ncar`, entry `gw_intr_ncar`): convective (Beres) and frontal/background sources, plus either an isotropic orographic scheme (`gw_oro`) or a ridge-based scheme (`gw_rdg`).

A Rayleigh friction layer (`RAYLEIGH_TAU1`) can be added in `postintr`. Which sources are active is chosen at run time (section 5). The GEOS and NCAR schemes are meant for different sources: the GEOS scheme for the extratropical background, the NCAR scheme for convective and orographic sources (`GEOS_GwdGridComp.F90`, Initialize, comments at the resource defaults).

## 2. Place in the hierarchy
- Parent is Physics, which adds GWD as a leaf child `GWD` (`GEOSagcm_GridComp/GEOSphysics_GridComp/GEOS_PhysicsGridComp.F90`, SetServices). It is the **first** child to run in each Physics step.
- No children. The schemes are Fortran modules (`gw_drag`, `gw_drag_ncar`) called from the driver.
- `CMakeLists.txt` builds one shared library from `GEOS_GwdGridComp.F90`, `gw_drag.F90`, `machine.F90`, `gwdps.f`, `gwdc.f` and the `ncar_gwd/` sources. It runs the MAPL automatic code generator (`mapl_acg`) on `GWD_StateSpecs.rc` to generate the import, export and internal specs and the pointer declarations. It installs `GWD_GridComp.rc` into `etc`.

## 3. Entry points
Registered in `SetServices` (`GEOS_GwdGridComp.F90`): `INITIALIZE` = `Initialize`, `RUN` = `Run`. There is one run phase. Finalize is the MAPL generic default.

### 3.1 SetServices
- Allocates a private state of derived type `GEOS_GwdGridComp` and stores it with `ESMF_UserCompSetInternalState` under the name `'GEOS_GwdGridComp'` (`GEOS_GwdGridComp.F90`, SetServices, `wrap_`).
- Reads `GWD_GridComp.rc` (key `use_threads:`, default `.FALSE.`), tells MAPL whether the component may use OpenMP threads, and allocates per-thread workspaces.
- Reads `NCAR_NRDG` first (default 0) because the generated internal-spec header needs it to decide whether the ridge fields exist (`GEOS_GwdGridComp.F90`, SetServices, comment "We need to get NCAR_NRDG").
- Includes the generated headers `GWD_Import___.h`, `GWD_Export___.h`, `GWD_Internal___.h`, then `MAPL_GenericSetServices`.

### 3.2 Initialize
Calls `MAPL_GenericInitialize` first, then (`GEOS_GwdGridComp.F90`, Initialize):
1. Reads `AGCM.GRIDNAME` and `AGCM.STRETCH_FACTOR` to get a horizontal size `imsize`, and from it a tuning factor `sigma = 1 - 0.9839*exp(-0.09835*4e7*0.9/imsize/3000)`. For a cubed-sphere name, `imsize` is the digits between the third character and the `x`, times 4 when the last two characters are `CF`.
2. Reads the scheme resources and defaults (section 5).
3. Sets up the Beres convective source (`gw_common_init`, `gw_beres_init`) from the table in `BERES_FILE_NAME` (default `ExtData/g5gcm/gwd/newmfspectra40_dc25.nc`) and the `NCAR_*` resources.
4. Sets up the orographic scheme: `gw_rdg_init` if `NCAR_NRDG>0`, else `gw_oro_init`.
5. Sets up the Newtonian cooling profile with `gw_newtonian_set(LM, PREF, alpha)`.

### 3.3 Run
The body runs only when the component's `RUNALARM` rings, and then calls the internal routine `Gwd_Driver`. The time step used is the alarm's ring interval (`ESMF_AlarmGet ... ringInterval`) (`GEOS_GwdGridComp.F90`, Run, Gwd_Driver). `ESMF_AlarmRingerOff` is commented out; the generic code handles it.
`Gwd_Driver` does:
1. Get the vertical grid: pressure at layer midpoints (`PMID`), thickness (`PDEL`), height (`ZM`) from `pregeo` and `geopotential`.
2. **NCAR path** (`gw_intr_ncar`), if `NCAR_EFFGWORO` or `NCAR_EFFGWBKG` is non-zero. The convective source uses `HT_dc` and `(1 - CNV_FRC)^4 * HT_mi` and `WSPD_STABLE300M` when `NCAR_EFFGWBKG>0`. The orographic source needs `NCAR_EFFGWORO>0`: ridge scheme if `NCAR_NRDG>0`, isotropic otherwise.
3. **GEOS path** (`gw_intr`), if `GEOS_EFFGWORO` or `GEOS_EFFGWBKG` is non-zero, using `GEOS_PGWV` and `GEOS_BGSTRESS`.
4. Add the two sets of tendencies.
5. `postintr`: Rayleigh friction when `RAYLEIGH_TAU1>0` (with `RAYLEIGH_Z1`, `H0`, `HH`), and the potential- and kinetic-energy diagnostics.
6. Export the results. `DTDT` is multiplied by `PDEL` (mass-weighted, `Pa K s-1`, for Physics). `TTMGW` is the unweighted heating. The `*_EXP` aliases hold the winds and temperature after GWD (section 4).

## 4. State
All specs come from `GWD_StateSpecs.rc`; that file is the source of truth for the list. The generated headers are not checked in.

### 4.1 Imports (`GWD_StateSpecs.rc`, category IMPORT)
| Names | Restart | Source (via Physics) |
|---|---|---|
| `PLE`, `T`, `Q`, `U`, `V`, `AREA`, `PREF` | skip | `Q` from Moist; the others from Dynamics through the Agcm (`GEOS_PhysicsGridComp.F90`, SetServices, "Imports for GWD") |
| `PHIS`, `SGH`, `VARFLT` | skip | topography, set by Agcm in Initialize (`GEOSagcm_GridComp/GEOS_AgcmGridComp.F90`, Initialize) |
| `HT_dc`, `HT_mi` | restarted | Moist `DTDT_DC`, `DTDT_micro` |
| `CNV_FRC` | restarted | Moist |
| `WSPD_STABLE300M` | skip | Dynamics export of the same name (`GEOSagcm_GridComp/GEOSsuperdyn_GridComp/@FVdycoreCubed_GridComp/DynCore_GridCompMod.F90`, export spec). *Unverified:* the Agcm-level connection that carries it to GWD. |

`Q` is the Moist internal humidity, exported by Moist and connected to GWD (`MAPL_AddConnectivity`, `SRC_ID=MOIST`). Because GWD runs first, *Inferred:* `Q` is the value left by Dynamics and the previous Moist step. `HT_dc`, `HT_mi` and `CNV_FRC` are exports of Moist's previous run, copied by the coupler after Moist ran. *Inferred:* this is why these three imports are restarted and the others are not.

## 4.2 Exports
- Aliased copies of the state after GWD: `PLE_EXP`, `T_EXP`, `Q_EXP`, `U_EXP`, `V_EXP`, `SGH_EXP`, `PREF_EXP` (exported under the names `PLE`, `T`, `Q`, `U`, `V`, `SGH`, `PREF`). Physics connects `U`, `V`, `T` to Moist (`GEOS_PhysicsGridComp.F90`, SetServices, "Moist Imports").
- Tendencies: `DUDT`, `DVDT`, `DTDT` (mass-weighted), `TTMGW`; split by source `*_ORO`, `*_BKG`, `*_RAY`; stresses `TAUOROX/Y`, `TAUBKGX/Y`, `TAUGWX/Y`, `BKG_TAU_*`, `TAUGW*_TOT/EAST/WEST`; energy diagnostics such as `PEGWD`, `FEGW`; ridge diagnostics `RDG1_*`.
- Physics reads `DUDT`, `DVDT`, `DTDT` (as `UIG`, `VIG`, `TIG`) to build the total tendencies, and `TAUOROX`, `TAUOROY` are connected to Moist (`GEOS_PhysicsGridComp.F90`, Run, SetServices).

### 4.3 Internal state and restarts
- The ridge data `KWVRDG`, `EFFRDG`, `GBXAR`, `HWDTH`, `CLNGT`, `MXDIS`, `ANGLL`, `ANIXY` exist only if `NCAR_NRDG /= 0` (condition in `GWD_StateSpecs.rc`, category INTERNAL). All have an ungridded dimension of 16 ridges except `GBXAR`.
- They are read from the internal restart file. Restart names in `AGCM.rc`: `MODELATMGWD_INTERNAL_RESTART_FILE: gwd_internal_rst`, `MODELATMGWD_IMPORT_RESTART_FILE: gwd_import_rst` (`src/Applications/@GEOSgcm_App/AGCM.rc.tmpl`, fixture). The run script handles per-face ridge files (`gwd_internal_c${AGCM_IM}_face_${n}`) (`src/Applications/@GEOSgcm_App/gcm_run.j`, fixture).
- The private derived-type state (section 3.1) holds the scheme parameters and per-thread workspaces; it is not restarted.

## 5. Configuration that changes structure
Resources are read in Initialize with `MAPL_GetResource` (`GEOS_GwdGridComp.F90`, Initialize). Defaults depend on `LM`.

| Key | Default | Effect |
|---|---|---|
| `JASON_BKG` | `LM==72` | selects the legacy background-drag defaults |
| `GEOS_PGWV`, `GEOS_BGSTRESS`, `GEOS_EFFGWBKG`, `NCAR_EFFGWBKG`, `RAYLEIGH_TAU1` | with `JASON_BKG`: `4`, `0.9`, `0.125`, `0`, `0`. Without: `NINT(32*LM/181)`, `0`, `0`, `0.375`, `0` | which background scheme is on and with what strength |
| `JASON_ORO` | `LM==72` | selects the legacy orographic defaults |
| `GEOS_EFFGWORO`, `NCAR_EFFGWORO`, `NCAR_NRDG` | with `JASON_ORO`: `0.25`, `0`, `0`. Without: `0`, `0.25` (`1.0` if `NCAR_NRDG==16`), `0` | orographic scheme: GEOS, NCAR isotropic, or NCAR ridge |
| `NCAR_NRDG` | 0 (read in SetServices and Initialize) | 0 disables the ridge scheme; 1 to 16 enables it and adds the ridge internal state |
| `BERES_FILE_NAME`, `NCAR_DC_BERES` (default true), `NCAR_ET_FAC_DTDTM` (0.75), `NCAR_ET_FAC_WS300` (1.25), `NCAR_BKG_TAU` | | Beres convective source |
| `RAYLEIGH_TAU1`, `RAYLEIGH_Z1`, `H0`, `HH` | | Rayleigh friction |
| `use_threads:` in `GWD_GridComp.rc` | `.FALSE.` | OpenMP over horizontal blocks |
| `DEBUG_GWD`, `DEBUG_TQ_ERRORS` | | print diagnostics |

`NCAR_NRDG` is substituted into `AGCM.rc` by the setup script (`AGCM.rc.tmpl`, "TRB/GWD Parameterization" block). Changing it also changes the restart contents.
The names `JASON_BKG` and `JASON_ORO` are legacy tuning switches and are kept as is.

## 6. Dependencies
- Libraries: `GEOS_Shared`, MAPL, ESMF, NetCDF Fortran, and OpenMP (not linked for NAG) (`CMakeLists.txt`).
- Uses `MAPL_get_current_thread` and `MAPL_find_bounds` from MAPL for threading.
- Inputs: the Beres spectrum file (3.2); topography from Agcm; Moist heating and convective fraction; `WSPD_STABLE300M` from Dynamics.
- Consumed by: Physics (tendencies), Moist (`U`, `V`, `T`, `TAUOROX/Y`), and Physics's stochastic pattern code, which reads `PREF` from the GWD import (`GEOS_PhysicsGridComp.F90`, Run).

## 7. Known traps
1. **Two schemes are summed.** If both `GEOS_EFFGWBKG` and `NCAR_EFFGWBKG` (or both orographic ones) are non-zero, both act. Check all four keys when tuning.
2. **`NCAR_NRDG` changes the restart.** Turning the ridge scheme on without ridge data in `gwd_internal_rst` leaves the fields at their default values. *Unverified:* what the default is; the run script must supply the ridge files.
3. **`Initialize` parses `AGCM.GRIDNAME`.** It assumes the third character starts the digits and an `x` follows, as in `PEx...`/`CF` cubed-sphere names. Other formats make the `read` fail (`GEOS_GwdGridComp.F90`, Initialize, `imchar`). *Inferred:* the code is not tested with other grid-name forms.
4. **Run happens only on the alarm**, and the time step is the alarm interval. If `DT:` is longer than `RUN_DT`, GWD is supercycled and the tendencies are applied only when it runs. *Unverified:* whether the held exports are re-applied on the in-between steps.
5. **`DTDT` is mass-weighted; `TTMGW` is not.** Use the right one.
6. **Generated headers.** `GWD_StateSpecs.rc` is the only source of the spec lists. Editing the Fortran for specs has no effect.
7. **Run order.** GWD runs before Moist. Its `CNV_FRC`, `HT_dc`, `HT_mi` are therefore from the previous step (4.1).
8. **Per-thread workspaces.** When `use_threads` is on, any new module-level scratch array in the schemes is a race.

## 8. Decisions
### Two scheme families under one component
- Pattern: both the GEOS and the NCAR scheme are called from `Gwd_Driver` and added together; resources decide which is active (`GEOS_GwdGridComp.F90`, Gwd_Driver). [Verified]
- Likely reason: the NCAR convective and ridge sources were added alongside the older scheme without removing it, and tuning moves between them by resolution. [Inferred]
- Question for maintainers: is the long-term plan to keep both?

### State specs generated from `GWD_StateSpecs.rc`
- Pattern: imports, exports and the optional ridge internal specs are declared in a table, and CMake generates the Fortran includes (`CMakeLists.txt`, `mapl_acg`). [Verified]
- Likely reason: removes boilerplate and keeps the list of fields in one reviewable file. [Inferred]

### Ridge internal state conditional on a resource read before specs
- Pattern: `NCAR_NRDG` is read in `SetServices` so that the table's `COND` can create or omit the ridge fields (`GWD_StateSpecs.rc`, `COND`). [Verified]
- Likely reason: avoids a ~16-level ungridded set of unused fields and restart variables when the ridge scheme is off. [Inferred]

### Run first, from Dynamics' state
- Pattern: GWD is the first child in Physics `Run`, with a comment that it must use `Q` from Dynamics (`GEOS_PhysicsGridComp.F90`, Run). [Verified]
- Likely reason: gravity-wave drag is tied to the dynamical state and not to the moist or turbulent adjustments; Moist then sees GWD-updated winds and temperature. [Inferred]

## 9. Observations (noted, not fixed) and open questions
1. **`gwdps.f` and `gwdc.f` appear unused.** They are in the CMake source list, but a search of the `.F90` files finds no call to `GWDPS` or `GWDC` (only their own definitions). *Inferred:* dead legacy code.
2. **`NCAR_NRDG` is read twice** with the same default (SetServices and Initialize, twice in Initialize's branches) (`GEOS_GwdGridComp.F90`).
3. **Copy-pasted description.** The `Run` description text refers to a "composite Moist Gridded Component" (`GEOS_GwdGridComp.F90`, Run, `DESCRIPTION`).
4. **A commented-out `ESMF_AlarmRingerOff`** in `Run` (3.3).
5. **A commented-out alternative default** for `RAYLEIGH_TAU1` (`172800.`) remains in Initialize.
6. **`NCAR_BKG_TAU`** is set to 3 times `sigma` if either NCAR source is non-zero, else 6.4 (`GEOS_GwdGridComp.F90`, Initialize). *Unverified:* the physical meaning of the constant.

### Open questions for maintainers
- **Grid-name parsing.** `Initialize` derives the cubed-sphere size `imsize` by parsing the text of `AGCM.GRIDNAME` and fails on names it cannot parse (trap 3). Question for maintainers: is parsing the grid name, rather than reading the size from the grid itself, intended to remain the source of `imsize`?
