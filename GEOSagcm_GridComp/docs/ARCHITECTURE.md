# GEOSagcm_GridComp architecture

A map of the atmosphere subtree of `GEOSgcm_GridComp`: the Agcm composite, Superdyn (dynamics) and Physics with its children Gravity Wave Drag, Moist and Turbulence. FVdycoreCubed, Chemistry and Radiation are separate repositories with their own documents; Surface is not covered here. Paths are relative to `GEOSagcm_GridComp`. It names modules and routines rather than linking to lines; use them to search the code. It describes only things that change rarely; status, rules for changes, known traps and to-dos are in [`AGENTS.md`](../AGENTS.md).

## Bird's-eye view

The Agcm is the atmospheric model of GEOS. Each step it advances the resolved flow with the dynamical core (Superdyn), computes the unresolved processes with Physics, and applies the physics tendencies to the dynamics state. Physics runs its children (gravity wave drag, moist processes, turbulence and surface exchange, chemistry and aerosols, radiation) in a fixed, staged order within the step, keeps temperature, humidity and winds consistent between stages, and returns the summed tendencies.

## Codemap

- `GEOS_AgcmGridComp.F90`: the Agcm composite. Adds `SUPERDYNAMICS` (or `SCMDYNAMICS` for single-column runs), `PHYSICS` and `ORBIT`, applies analysis increments or bias tendencies, and connects the children.
- `GEOSsuperdyn_GridComp/`: Superdyn (`GEOS_SuperdynGridComp.F90`), the wrapper around the dynamical core. `@FVdycoreCubed_GridComp` is FV3, a separate nested repository; `GEOSdatmodyn_GridComp/` is a single-column (data) dynamics driver.
- `GEOSphysics_GridComp/GEOS_PhysicsGridComp.F90`: the Physics composite. Adds six children, connects them, runs them, and sums their tendencies; `MBundle_IncrementMod.F90` computes tracer increments into a bundle.
- `GEOSphysics_GridComp/GEOSgwd_GridComp/`: Gravity Wave Drag (`GEOS_GwdGridComp.F90`, `GWD_StateSpecs.rc`), with the GEOS scheme (`gw_drag.F90`) and the NCAR scheme (`ncar_gwd/`).
- `GEOSphysics_GridComp/GEOSmoist_GridComp/`: Moist (`GEOS_MoistGridComp.F90`, `Moist_StateSpecs.rc`). Each scheme sits behind an interface module `GEOS_<scheme>_InterfaceMod.F90` (`<scheme>_Setup`, `_Initialize`, `_Run`, `_Finalize`): microphysics GFDL_1M, BACM_1M, MGB2_2M, THOM_1M, NSSL_2M; deep convection GF, RAS; shallow convection UW. Shared code is in `Process_Library.F90`; aerosol activation in `aer_actv_single_moment.F90`. `pyMoist/` is a Python port with its own tests; `mapl3-porting-notes.md` and its two companions record the MAPL3 port.
- `GEOSphysics_GridComp/GEOSturbulence_GridComp/`: Turbulence (`GEOS_TurbulenceGridComp.F90`), with Lock entrainment (`LockEntrain.F90`), EDMF (`edmf.F90`) and SHOC (`shoc.F90`).
- `GEOSphysics_GridComp/@GEOSchem_GridComp`, `GEOSphysics_GridComp/@GEOSradiation_GridComp`: Chemistry and Radiation, separate nested repositories. `GEOSphysics_GridComp/GEOSsurface_GridComp/`: Surface (not covered here).
- `GEOSphysics_GridComp/regression/`, `GEOSphysics_GridComp/GEOSgwd_GridComp/regression/`: standalone regression cases for Moist and GWD.

## How the Agcm is put together

- **Agcm run order within one step:** `ORBIT`; `SUPERDYNAMICS` phase `run` (dynamics, with analysis increments or bias tendencies); `PHYSICS`; `SUPERDYNAMICS` phase `run_add_incs` (adds the physics tendencies to the dynamics state).
- **Superdyn:** adds the dynamical core as child `DYN`, chosen by `DYCORE` (default `FV3`), with configuration `dyn.yaml`. It runs `DYN`'s two phases under the same names and re-exports most of `DYN`'s exports (the state, `PLE`, `ZLE`, `PREF`, `DELP`, `OMEGA`, the analysis and physics increment diagnostics, `WSPD_STABLE300M`, ...) as its own, including the advection service `TRADV`.
- **Physics children**, added in this order: `GWD`, `MOIST`, `TURBULENCE`, `CHEMISTRY`, `SURFACE`, `RADIATION`. Chemistry is added before Surface so that it fills `AERO_DP` before Surface reads its restart.
- **Physics run order within one step:** `GWD`; `MOIST`; stage 1: `SURFACE` phase 1, `CHEMISTRY` phase 1, `TURBULENCE` phase 1; stage 2: `SURFACE` phase 2, `TURBULENCE` phase 2; then `CHEMISTRY` phase 2 and `RADIATION`. Between stages the Physics composite synchronizes temperature, humidity and winds (`SYNCTQ`, `SYNCUV`) and clips negative water species (`FILLQ2ZERO`).
- **Tracers that others change** (advected by the dynamics, convected by Moist, diffused by Turbulence, updated by chemistry or analysis) were `FRIENDLYTO` bundles in MAPL2 and are services in MAPL3: `TRADV`, `MTR`, `TRDIFF`, `TRCHEM`, `TRANA`. A component subscribes the fields it owns; the provider updates them in place.
- **Gravity Wave Drag:** one phase; the GEOS and NCAR schemes both act, scaled by `GEOS_EFFGWBKG`, `GEOS_EFFGWORO`, `NCAR_EFFGWBKG`, `NCAR_EFFGWORO`, and their tendencies are summed.
- **Moist:** one scheme of each kind, chosen by `CLDMICR_OPTION` (default `BACM_1M`), `CONVPAR_OPTION` (default `GF`) and `SHALLOW_OPTION` (default `UW`) in `if` chains in `SetServices`, `Initialize`, `Run` and `Finalize`. `Run` sets up the convective tracers, computes aerosol activation, runs deep and shallow convection (UW first when `SH_MD_DP` is set), combines their outputs, then runs the microphysics.
- **Turbulence:** two phases, in stages 1 and 2; it diffuses heat, moisture, momentum and the tracers subscribed to it.

## Contracts

Sources and consumers are the connections in the composites' `SetServices`; other imports come from ExtData.

- **Agcm and Superdyn.** Superdyn exports the dynamics state that Physics and its children import (`PLE`, `ZLE`, `T`, `U`, `V`, `OMEGA`, `PREF`, `PHIS`, `AREA`, the pre-dynamics state and dynamics tendencies, `WSPD_STABLE300M` for GWD) and provides `TRADV`. Physics returns its summed tendencies, which Superdyn's `run_add_incs` applies. The Agcm connects Superdyn's `TRADV` to Physics so that the Physics children's subscriptions reach FV3.
- **Gravity Wave Drag.**
  - Imports: the state and subgrid orography from the parent; from Moist, `CNV_FRC` and its convective and microphysical heating (`DTDT_DC` and `DTDT_micro`, imported as `HT_dc` and `HT_mi`).
  - Exports: drag tendencies and diagnostics; Moist takes `TAUOROX`, `TAUOROY` and the updated winds and temperature.
  - Services: none.
- **Moist.**
  - Imports: the state and dynamics tendencies from the parent; PBL quantities and turbulence moments from Turbulence; surface fields from Surface; orographic stress from GWD; radiative heating (`RADLW`, `RADSW`) from Radiation; the `AERO` state from Chemistry; subgrid vertical-velocity and droplet-number climatologies from ExtData.
  - Exports: water species, cloud fraction and effective radii to Radiation; precipitation to Surface; total water and buoyancy fluxes to Turbulence; condensates, precipitation fluxes and convective mass fluxes to Chemistry; `CNV_FRC` and heating to GWD.
  - Services: provides `MTR` (convective transport and scavenging); subscribes its own fields to `TRADV`, `TRDIFF`, `TRCHEM` and `TRANA`.
- **Turbulence.**
  - Imports: the state from the parent; total water and buoyancy fluxes from Moist; exchange coefficients and surface fluxes from Surface; `RADLW`, `RADLWC` from Radiation.
  - Exports: diffusivities, PBL depth and turbulence moments to Moist and Chemistry; implicit-diffusion fields that the Physics composite uses between stages.
  - Services: provides `TRDIFF` (turbulent diffusion; MAPL2's `TR` import bundle); subscribes `PDF_A`, `TKESHOC`, `QT2`, `QT3` to `TRADV` when SHOC is on. Its own `S`, `U`, `V`, `TKESHOC` and `QT2` are diffused directly, not through the bundle.

## Invariants

- The Agcm and Physics run orders are part of the model; changing either changes results.
- Because of the run order, GWD uses the previous step's convection, and Moist and Turbulence use the previous step's radiative heating.
- Physics tendencies reach the dynamics state only through Superdyn's `run_add_incs`.

## Cross-cutting concerns

- **Generated specs:** GWD, Moist and Turbulence declare their specs in `*_StateSpecs.rc` (ACG); Moist's scheme fields are added by the `_Setup` routines. Superdyn declares none of its own: it re-exports `DYN`'s.
- **Threading:** GWD keeps per-thread workspaces and the Moist schemes use OpenMP loops; module-level scratch data and loop-carried dependencies are races.
- **Exports as work arrays:** Moist uses some exports as scratch space and reads others back by name, so which exports are requested matters (`GEOSphysics_GridComp/GEOSmoist_GridComp/mapl3-alloc-exports.md`).
- **Testing:** each regression case runs one component (GWD, Moist, Turbulence, or in FVdycoreCubed the dynamical core) for one step under a configurable Root, with imports from a data component restarted from a MAPL2 capture, and compares checkpoints exactly against a stored MAPL3 baseline validated against that capture.

## Where to go next

- [`AGENTS.md`](../AGENTS.md): rules, status on this line, to-dos, known traps, build and test for the Agcm subtree.
- `GEOSphysics_GridComp/GEOSmoist_GridComp/mapl3-porting-notes.md`: the Moist port, its decisions and their reasons.
- FVdycoreCubed, Chemistry and Radiation: their `AGENTS.md` and `docs/ARCHITECTURE.md`.
