# Moist MAPL3 port: the 197 `ALLOC=.TRUE.` exports

This is the survey for step 6 of `mapl3-porting-notes.md` (2026-10-29). It covers the driver and the in-scope schemes (GFDL_1M, GF, UW) on branch `feature/pchakrab/port-moist-to-mapl3` at `8bde1bc9`. The out-of-scope schemes (BACM_1M, THOM_1M, MGB2_2M, NSSL_2M, RAS) and pyMoist are not covered.

File abbreviations used throughout:

| Short | File |
|---|---|
| `Moist` | `GEOS_MoistGridComp.F90` |
| `GFDL_1M` | `GEOS_GFDL_1M_InterfaceMod.F90` |
| `GF` | `GEOS_GF_InterfaceMod.F90` |
| `UW` | `GEOS_UW_InterfaceMod.F90` |
| `PL` | `Process_Library.F90` |
| `Physics` | `../GEOS_PhysicsGridComp.F90` (MAPL2) |

## Summary

Moist registers 559 exports. **197 of them are fetched with `ALLOC=.TRUE.`**, at 226 call sites (driver 45, GFDL_1M 70, GF 60, UW 51). In MAPL2 that call makes the export a buffer that is zero-filled once, persists across steps, and is shared by name, and the code depends on all three properties. MAPL3 has no equivalent: an export nobody connects is never allocated, and a field is zero-filled only when its spec says so.

Sorted by what each export needs:

| Group | Count | Needs | Exports |
|---|---|---|---|
| **A. Shared within a step** | 49 | Storage visible to several routines within one step | Driver inputs to the schemes, GFDL_1M tendencies, precipitation/flux outputs read by the driver epilogue, and 4 fields the driver overwrites |
| **B. Persisting across steps** | 25 | Storage that also keeps its value until a later step | `CNPCPRATE` (every default run), plus 24 GF/UW results when the convection timestep is longer than Moist's (production c1080 and above). 21 feed model state; 4 only reach diagnostics |
| **C. Local to one routine** | 123 | No shared storage | 91 plain work arrays, 22 that rely on the zero fill, 2 re-fetched by name, 8 only written under `associated()` |

Main consequences:

- **The schemes read A1 and B without guards** (`CNV_FRC`, `SRF_TYPE`, `EIS`, `BYNCY`, `CNPCPRATE`, ...). Dropping `ALLOC` leaves them unallocated, and the first read is a null dereference.
- **The driver epilogue guards its terms**, so with unallocated A2/A3/B2 fields it still runs but silently leaves terms out of `TPREC`, `RAIN`, `PCU`, `PLS` and `DTDT`, ... `TPREC` feeds Surface.
- **22 exports (C1) are never written in some configurations** and need `FILL: 0.0`, connected or not. That includes `WQT_DC` and `WTHV2`, which Turbulence imports.
- **The 74 in A and B need storage that outlives one routine call.** Internals re-exported with `ADD2EXPORT` reproduce MAPL2 directly. A private-state derived type is the alternative. See "Options for step 6". **Decided 2026-10-29/30: A and B both become internals** (porting notes, decisions 1b and 1c).

## Background

### What `ALLOC=.TRUE.` did in MAPL2

- `MAPL_GetPointer(..., ALLOC=.TRUE.)` allocates the export through `MAPL_AllocateCoupling` and `MAPL_FieldAllocCommit` (MAPL2 `base/Base/Base_Base_implementation.F90:40-107`). On the first fetch only, the field is initialized to the spec's `DEFAULT`, or to 0.0 when there is none (`:160-164`).
- After that the field lives in the export state. Any later `MAPL_GetPointer` on the same name, from any routine and on any step, returns the same buffer with its last values.

### What MAPL3 does instead

- `MAPL_StateGetPointer` has no `alloc` argument. An export is allocated only when something connects to it or History requests it. Internals are always allocated (porting notes, issue 2).
- A field is filled at allocation only when its spec has a fill value (`@MAPL/superstructure/generic/specs/FieldClassAspect.F90:254-255`). Otherwise its contents are undefined.
- An internal can be re-exported under the same name with `add_to_export` (`@MAPL/superstructure/generic/MAPL_Generic.F90:740-750`), which ACG exposes as the `ADD2EXPORT` column (`@MAPL/apps/MAPL_GridCompSpecs_ACG.py:249`). Re-exporting from `<self>` is exercised in `@MAPL/gridcomps/configurable/tests/scenarios/scenario_1/child_A.yaml:22-28`. GEOS precedents:
  - **DynCore, through ACG.** `DynCore_StateSpecs.rc:239-252` (`@FVdycoreCubed_GridComp`, on `release/MAPL-v3`) re-exports internals through the `ADD2EXPORT` column: `AK`, `BK` and the 3D edge field `PE`, which goes out as `PLE` via `EXPORT_NAME`. All three use `RESTART: REQUIRED`. The same commit (`68b4705`) replaced DynCore's `alloc=.true.` export fetches (`PREF`, `DXC`, `DYC`, `AREA`) with `MAPL_StateGetPointer` plus `if (associated(...))` copies.
  - **AgcmSimple, hand-written.** `GEOS_AgcmSimpleGridComp.F90:96-111` (#1455, `071ca608`) passes `add_to_export=.true.` for the 2D internals `PHIS` and `VARFLT`. The held-suarez regression loads it (`GEOShs_GridComp/regression/held-suarez/agcm-simple.yaml:3-5`). The file is excluded from this checkout by sparse-checkout, so use `git show HEAD:GEOSagcm_GridComp/GEOS_AgcmSimpleGridComp.F90`.

### One Moist step

1. **Driver prologue** (`Moist:419-690`) computes:
   - `SRF_TYPE` (`:458-467`);
   - `LTS`/`EIS` (`FIND_EIS`, `:535`);
   - `BYNCY`, `CAPE`, ... (`:618-622`);
   - `CNV_FRC` (`:624-645`);
   - aerosol activation (`:655-684`).
2. **Convection** (`Moist:694-725`). With the defaults (`CONVPAR_OPTION = "GF"`, `SHALLOW_OPTION = "UW"`, `Moist:113, 123`; `SH_MD_DP = .TRUE.`, `GFDL_1M:315`), `UW_Run` runs before `GF_Run`.
   - Each scheme computes only when its own alarm rings (`GF:705-949`, `UW:409-636`), but applies its stored tendencies on every Moist step (`GF:955-975`, `UW:641-680`).
   - `GF_DT`/`UW_DT` default to Moist's timestep (`GF:116`, `UW:75`). The GCM sets both to `CONV_DT` (`@GEOSgcm_App/AGCM.rc.tmpl:64-65`), which is longer than `DT` at the highest resolutions (group B2).
3. **Microphysics** (`Moist:792`): `GFDL_1M_Run`.
4. **Driver epilogue** (`Moist:742-1330`): mass-flux, tendency and precipitation totals, the precipitation-type diagnostic, and the `*AFMOIST` copies.
   - When Moist's own alarm does not ring, the else branch (`Moist:1279-1330`) only recomputes the `*AFMOIST` copies from the current state.

### Method

- A script found every `MAPL_GetPointer`/`MAPL_StateGetPointer` on `EXPORT` in the four files, joining continuation lines, and grouped the fetches by short name and by routine.
- For each fetch it followed the pointer variable until it was re-fetched or re-pointed, and marked each use as guarded when an enclosing `if` tested `associated()` on that variable.
- Every unguarded use was read by hand. Kernel outputs were traced into the kernel to check whether every element is assigned on every call.
- This is static analysis only: nothing was built or run. The script missed statements after `;` on a fetch line; the only ones are at `GFDL_1M:717-722` and were read by hand.

## A. Shared within a step (49)

These are fetched in more than one routine and read in the same step they are written.

### A1. Driver-computed inputs to the schemes (5)

The driver computes them in its prologue. The schemes fetch them again, with `ALLOC`, and read them without guards.

| Export | Written | Read |
|---|---|---|
| `SRF_TYPE` | `Moist:458-467` | `GFDL_1M:839, 914, 950`; `GF:827, 878, 919`; `UW:467` |
| `EIS` | `Moist:535` (`FIND_EIS`, `PL:3842-3873`) | `GFDL_1M:839, 1158, 1239`; `UW:467` |
| `LTS` | `Moist:535` | Nowhere. `GFDL_1M:724` and `UW:281` fetch it but never use it |
| `CNV_FRC` | `Moist:624-645` (from `CAPE`) | `GF:791, 797, 827, 878, 919`; `GFDL_1M:914, 950, 960, 979, 1158, 1239, 1387-1388, 1395`. `UW:423` fetches it but never uses it |
| `BYNCY` | `Moist:618-622` (`BUOYANCY2`, `BUOYANCY`) | `GF:814` (`compute_sgs_vvel`), `GF:821` (GF2020) |

If the schemes got these from local arrays, deep convection, macrophysics and microphysics would run on null or uninitialized data. Outside Moist, `CNV_FRC` goes to GWD and Chem, `EIS` to Turbulence, and `BYNCY` to Chem (`Physics:1542-1546, 1693, 1716-1728`).

### A2. GFDL_1M macro- and microphysics tendencies (20)

These are `DUDT`, `DVDT`, `DTDT`, `DQVDT`, `DQLDT`, `DQIDT`, `DQADT`, `DQRDT`, `DQSDT` and `DQGDT`, each with a `_macro` and a `_micro` version.

- Every step, each is rewritten as an after-minus-before difference (`GFDL_1M:777-786` then `:1038-1047`; `:1106-1115` then `:1438-1447`).
- GFDL_1M uses them itself in `dissipative_ke_heating` (`:1490`) and in the `DQRDT` sum (`:1485`).
- The driver sums them under `associated()` (`Moist:882-999`), so a missing field silently drops out of the total. `DTDT_micro` also goes to GWD as `HT_mi` (`Physics:1698-1703`).

### A3. Precipitation, flux and mass-flux outputs read by the driver epilogue (20)

| Subgroup | Exports | Written | Read by the driver |
|---|---|---|---|
| Large-scale precipitation | `LS_PRCP`, `LS_SNR`, `ICE`, `FRZR` | `GFDL_1M:1355-1358` | `LS_PRCP` unguarded: `Moist:854-866` adds `ER_PRCP`. `ICE`/`FRZR` unguarded: `Moist:1161-1169`, where `DIAGNOSE_PRECIP_TYPE` takes them `inout`. The rest guarded (`:1061-1148`) |
| Precipitation by type | `CN_PRCP`, `AN_PRCP`, `SC_PRCP`, `CN_SNR`, `AN_SNR`, `SC_SNR` | Only zeroed, by GFDL_1M (`GFDL_1M:717-722`) | The driver accumulates into them (`Moist:1101-1134`) and sums them (`:1048-1148`), guarded. GFDL_1M's fetch exists only to allocate them and reset them each step |
| Precipitation fluxes (edge) | `PFL_AN`, `PFL_LS`, `PFI_AN`, `PFI_LS`, `PFL_CN`, `PFI_CN` | `GFDL_1M:1133-1134, 1168-1183, 1378-1383`; `GF:894-895, 934-935` | Guarded. They build `DPDTMST`, `PFL_LSAN` and `PFI_LSAN` (`Moist:1005-1039`) |
| Mass fluxes | `MFD_DC`, `MFD_SC`, `UMF_DC` (edge), `UMF_SC` (edge) | GF and UW kernels | Guarded (`Moist:742-746, 782-786`). Also used inside the schemes: `GF:887, 928`, `UW:561-570` |

Outside Moist, `ICE`, `FRZR` and `CN_PRCP` go to Surface, and `PFL_CN` and `PFI_CN` go to Chem (`Physics:1639-1643, 1716-1728`). The totals built from them go outside Moist as well: `TPREC` to Surface and Chem, `PCU`/`PLS` to Surface, and `PFL_LSAN`/`PFI_LSAN` to Chem.

### A4. Written in a scheme and overwritten by the driver (4)

- **`QLTOT`, `QITOT`.** UW sets them (`UW:507-508`) and uses them as kernel inputs (`UW:510-519`). The driver overwrites them at the end of the step (`Moist:1195-1199`, guarded). They go to Turbulence and Chem.
- **`RL`, `RI`** (pointer `CLDREFFI`). The driver fills them with `MAPL_UNDEF` (`Moist:498-501`, guarded), then GFDL_1M writes every cell through `RADCOUPLE_SCALE_AWARE` (`GFDL_1M:1365-1400`, `PL:1454-1533`). They go to Radiation and Chem.

Strictly, these four need no shared storage, since no routine reads a value another routine wrote. Only the final value has to reach the export.

## B. Persisting across steps (25)

These must keep their value until a later step.

### B1. Every step, default configuration (1)

**`CNPCPRATE`.**
- GF writes it, in a kernel that zeroes it on entry (`ConvPar_GF2020.F90:303`, `ConvPar_GF_GEOS5.F90:235`).
- UW passes it to its kernel as `cnvtr` (`UW:519`, `uwshcu.F90:86, 128`) and uses it over ocean (`uwshcu.F90:323-325`).
- Because UW runs before GF by default, UW always reads the value GF wrote on the previous step, or the zero fill on the first step.
- The driver also reads it, guarded (`Moist:1052, 1080, 1103`).

### B2. When the convection timestep is longer than Moist's (24)

Both alarms are created to ring on the first step (`RingTime = currentTime`) and then every `GF_DT`/`UW_DT` (`GF:123-128`, `UW:82-87`). Nothing outside the ringing blocks rewrites these fields (`GF:587-704`, `UW:280-408`), except `DQLDT_SC`/`DQIDT_SC` (below). So on a non-ringing step, every field below still holds the value from the last ring.

| Exports | Computed when the alarm rings | Read on every Moist step | Effect |
|---|---|---|---|
| GF tendencies (7): `DUDT_DC`, `DVDT_DC`, `DTDT_DC`, `DQVDT_DC`, `DQLDT_DC`, `DQIDT_DC`, `DQADT_DC` | `GF:705-949` | Added to `U`, `V`, `Q`, `T`, `QLCN`, `QICN`, `CLCN` (`GF:955-975`) | Model state |
| UW tendencies (5): `DUDT_SC`, `DVDT_SC`, `DTDT_SC`, `DQVDT_SC`, `DQADT_SC` | `UW:409-636` | Added to `Q`, `T`, `U`, `V`, `CLCN` (`UW:661-666`) | Model state |
| UW condensate sources (6): `QLDET_SC`, `QIDET_SC`, `QLSUB_SC`, `QISUB_SC`, `QLENT_SC`, `QIENT_SC` | `UW:519`, `:573-588` | Fetched again (`UW:641-646`) and added to `QLCN`, `QICN`, `QLLS`, `QILS` (`UW:668-676`) | Model state |
| UW precipitation for GFDL_1M (2): `SHLW_PRC3`, `SHLW_SNO3` | `UW:596-599` (copies of `DQRDT_SC`/`DQSDT_SC`) | Added to `QRAIN`/`QSNOW` by GFDL_1M (`GFDL_1M:797-819`, guarded) | Model state |
| UW precipitation tendencies (2): `DQRDT_SC`, `DQSDT_SC` | `UW:519` | Only in the driver's tendency sums (`Moist:961-977`) | Diagnostics only |
| UW condensate tendencies (2): `DQLDT_SC`, `DQIDT_SC` | `UW:510-511` (before state), `UW:519` | Rewritten from their own previous value (`UW:679-680`) | Diagnostics, wrong on non-ringing steps |

So **20 of the 24 feed model state**. With B1, that makes 21 of the 25 in group B. The other 4 persist but only reach diagnostics. The driver sums all 16 GF/UW tendencies under `associated()` (`Moist:878-999`), and `DTDT_DC` also goes to GWD as `HT_dc` (`Physics:1705-1710`).

**`DQLDT_SC`/`DQIDT_SC` are computed wrongly on non-ringing steps.**
- On a ringing step, UW stores the condensate before convection in them (`UW:510-511`). The apply loop then turns that into a tendency: `DQLDT_SC = (QLLS + QLCN - DQLDT_SC) / MOIST_DT` (`UW:679-680`).
- The apply loop runs on every Moist step. On a non-ringing step the field already holds last step's tendency, so the formula subtracts a rate from a mixing ratio.
- Model state is unaffected. The result does flow into the `DQLDT`/`DQIDT` totals (`Moist:935-955`), which Chem imports (`Physics:1716-1728`).
- Reproducing MAPL2 bit-for-bit still requires this value to persist.

When `GF_DT`/`UW_DT` equals Moist's timestep, the alarms ring every step and nothing carries over. They are set from `CONV_DT`, which the setup script chooses per resolution (`@GEOSgcm_App/set_resolution_params.csh`):

| Resolution | `DT` (s) | `CONV_DT` (s) | Moist steps per convection call | Lines |
|---|---|---|---|---|
| c12, c24, c48 | 1200 | 1200 | 1 | `:29-30, 52-53, 70-71` |
| c90 | 900 | 900 | 1 | `:88-89` |
| c180, c270 | 600 | 600 | 1 | `:121-122, 283-284` |
| c360 | 450 | 450 | 1 | `:157-158` |
| c540, c720, c1120 | 300 | 300 | 1 | `:305-306, 176-177, 196-197` |
| c1080, c1440 | 150 | 300 | 2 | `:327-328, 216-217` |
| c1536, c2160, c2880, c4320, c5760 | 75 | 300 | 4 | `:349-350, 371-372, 236-237, 393-394, 258-259` |

Two overrides keep the two equal: coupled c180 runs (`:138-141`) and BACM_1M runs (`:435-437`) both use 450 s. A regression at c1120 or below never exercises B2. At those resolutions the 24 behave like group A: shared within a step, with only `CNPCPRATE` carried over to the next step.

Moist itself runs every heartbeat, since `AGCM.rc.tmpl` sets no Moist-specific run interval. This relies on MAPL2's `RUN_DT` default, which was not re-checked.

## C. Local to one routine (123)

Each of these is fetched in only one routine and needs no storage beyond that routine. The four subgroups differ in what a local array must still do.

### C1. Not written, so they rely on the zero fill (22)

| Exports | When they are not written | In the default config? | Evidence |
|---|---|---|---|
| `LCL_AGL` | Never: `BUOYANCY2` takes `LCL` but never assigns it | Yes | `PL:1268-1409`, `Moist:615-619` |
| `WQT_DC` | On the `GF_GEOS5_Interface` path, which does not take it | Yes when `LM = 72` (`USE_GF2020` defaults to 0, `GF:143-144`) | `GF:846-862` |
| `LFR_GF`, `CNV_TOPP_DP`, `CNV_TOPP_MD`, `CNV_TOPP_SH`, `MDNDP`, `ENTR_DP`, `ENTR_MD`, `ENTR_SH`, `SGS_VVEL_DP`, `SGS_VVEL_MD`, `SGS_VVEL_SH` | Same path. `compute_sgs_vvel` runs only on the GF2020 path (`GF:814`) | Yes when `LM = 72` | `GF:846-862` |
| `WTHV2`, `WQL` | `hystpdf` sets them only when `PDFSHAPE == 6` | Yes (`PDFSHAPE` defaults to 1, `GFDL_1M:339`) | `PL:2923-3013` |
| `VFALL_ICE`, `VFALL_SNOW`, `VFALL_GRAUPEL`, `VFALL_RAIN` | `gfdl_mp_driver` path: never passed. Old driver path: only levels `ktop:kbot` | Yes (`GFDL_MP3` defaults to `.TRUE.`, `GFDL_1M:290`) | `GFDL_1M:1152-1161`, `gfdl_cloud_microphys.F90:829` |
| `NWFA`, `NIFA` | Unless `USE_AEROSOL_NN` is on and `USE_NCLOUD_CLIM` is off. When written, zeroed on entry (`aer_actv_single_moment.F90:89-90`) | No (defaults `.TRUE.`/`.false.`, `PL:85`, `Moist:142`) | `Moist:655-684` |
| `PRCP_WATER` | Old GFDL driver path: not passed, though `GFDL_1M:1349` reads and rescales it every step | No | `GFDL_1M:1230-1242` |

A local array must be zeroed, and the spec row needs `FILL: 0.0`. MAPL3 leaves even a connected export undefined without a fill value. `WQT_DC` and `WTHV2` are Turbulence imports (`Physics:1542-1546`), and in a default 72-level MAPL2 run Turbulence receives zeros for both.

### C2. Re-fetched by name later in the same routine (2)

`DTDT_ER` and `DQVDT_ER` are written through the ALLOC'd pointer (`Moist:848-871`), then fetched again without ALLOC in the epilogue (`:912`, `:928`) behind `associated()`. With a local array, those re-fetches must use the local variable instead.

### C3. Only written under `associated()` (8)

- **`UAFMOIST`, `VAFMOIST`, `TAFMOIST`, `QAFMOIST`, `THAFMOIST`, `SAFMOIST`** (`Moist:1233-1251` and `:1311-1329`, one site in each branch of Moist's alarm). Moist itself doesn't need ALLOC for these. It exists because the MAPL2 Physics parent reads them from Moist's export state without a guard (`Physics:3118-3176`, e.g. `UFORTURB = UAFMOIST`). That is a step 12 wiring question.
- **`PREC_STRAT`, `PREC_CONV`** (`Moist:1153-1158`). Dropping ALLOC is enough.

### C4. Fully written every step (91)

Every element is assigned on every call before it is read. Most are kernel outputs that the kernel zeroes on entry.

| Where | Exports | Evidence |
|---|---|---|
| Driver: buoyancy | `CAPE`, `INHB`, `SBCAPE`, `SBCIN`, `MLCAPE`, `MLCIN`, `MUCAPE`, `MUCIN`, `ZLFC`, `ZLNB` | `BUOYANCY` assigns all its outputs (`PL:1411-1452`). `RETURN_CAPE_CIN` initializes `LFC`/`LNB` (`PL:1179-1183`). `BUOYANCY2` checks `associated()` on the `ML*`/`MU*` pointer dummies (`PL:1286-1287`) |
| Driver: other | `UMST0`, `VMST0`, `CFICE`, `ER_PRCP`, `RAIN`, `SNO`, `TPREC`, `PTYPE` | Whole-array assignments (`Moist:492-495, 807-813, 855-866, 1069-1148`). `PTYPE = 0` on entry (`PL:3629`) |
| GFDL_1M | `QV`, `QL`, `QI`, `QR`, `QS`, `QG`, `FCLD` (pointers `RAD_*`), `EVAPC`, `SUBLC`, `PDFITERS`, `RHX`, `RHCRIT`, `PRCP_RAIN`, `PRCP_SNOW`, `PRCP_ICE`, `PRCP_GRAUPEL`, `REV_LS`, `RSU_LS`, `PFR_LS`, `PFS_LS`, `PFG_LS` | `RAD_*` assigned before the call (`GFDL_1M:1138-1145`). `EVAPC`, `SUBLC`, `PDFITERS`, `RHX` zeroed each step (`:788-791`). `RHCRIT` set per cell (`:910`). `gfdl_mp_driver` zeroes all its outputs (`gfdl_mp.F90:758-776`) |
| GF | `CNV_MF0`, `CNV_PRC3`, `CNV_DQCDT`, `CNV_UPDF`, `CNV_CVW`, `CNV_QC`, `ENTLAM`, `SIGMA_DEEP`, `SIGMA_MID`, `MUPDP`, `MUPSH`, `MUPMD`, `MFDP`, `MFSH`, `MFMD`, `ERRDP`, `ERRSH`, `ERRMD`, `AA0`, `AA1`, `AA2`, `AA3`, `AA1_BL`, `AA1_CIN`, `TAU_BL`, `TAU_DP`, `TAU_MD`, `REV_CN`, `RSU_CN`, `TPWI`, `TPWI_star`, `DNDCNV`, `DNICNV` | Both kernels zero their outputs (`ConvPar_GF2020.F90:281-292`, `ConvPar_GF_GEOS5.F90:226-263`). `REV_CN`/`RSU_CN` assigned as whole arrays (`GF:888-893, 926-932`). `TPWI*` computed in `GF_Run` (`GF:732-733`). `DNDCNV`/`DNICNV` zeroed (`GF:805-806`) |
| UW | `PLCL_SC`, `PLFC_SC`, `PINV_SC`, `PREL_SC`, `PBUP_SC`, `CBMF_SC`, `CLDTOP_SC`, `CUFRC_SC`, `DCM_SC`, `ENTR_SC`, `DETR_SC`, `SC_NDROP`, `SC_NICE`, `TPERT_SC`, `QPERT_SC`, `QTFLX_SC`, `SLFLX_SC`, `UFLX_SC`, `VFLX_SC` | `compute_uwshcu_inv` copies every column and level out of the kernel, with no `cycle` (`uwshcu.F90:316-434`) |

Some of these feed prognostic state:

- `RAD_*` are written back into `Q`, `QRAIN`, `QSNOW` and `QGRAUPEL` (`GFDL_1M:1300-1303`).
- `RHCRIT` feeds `EVAP3`/`SUBL3` (`GFDL_1M:960, 979`).
- `CAPE` feeds `CNV_FRC` (`Moist:628-640`).
- `PRCP_*` become `LS_PRCP`, `LS_SNR` and `ICE` (`GFDL_1M:1355-1357`).
- `DCM_SC`/`DETR_SC` become `MFD_SC` (`UW:561-567`).
- `CNV_DQCDT` feeds the GF ice/liquid split (`GF:881, 920`).

A local array is still correct for all of them, because each is assigned before use.

Outside Moist:
- **Radiation:** `QV`, `QL`, `QI`, `QR`, `QS`, `QG` and `FCLD`.
- **Turbulence:** `QV` and `FCLD`.
- **Chem:** `QL`, `QI`, `FCLD`, `CNV_QC`, `CNV_CVW`, `REV_LS`, `REV_CN`, `TPREC`, `CAPE`, `INHB` and `ZLFC`.
- **Surface:** `SNO` and `TPREC`.

These connections are at `Physics:1542-1728`.

## Other findings

- **Fill values.** Anything that keeps a value without writing it (C1), or is read before its first write (B, on the first step), needs `FILL: 0.0` to match MAPL2 (`FieldClassAspect.F90:254-255`).
- **Edge fields** need the `0:LM` remap (porting notes, step 5 follow-ups): `PFL_AN`, `PFL_LS`, `PFI_AN`, `PFI_LS`, `PFR_LS`, `PFS_LS`, `PFG_LS`, `UMF_DC`, `PFI_CN`, `PFL_CN`, `WQT_DC`, `UMF_SC`, `QTFLX_SC`, `SLFLX_SC`, `UFLX_SC`, `VFLX_SC`.
- **Totals allocated by the parent.** In MAPL2 the Physics parent allocates Moist's epilogue totals with `alloc=.true.` from outside Moist (`Physics:2914-2986`): `DUDT`, `DVDT`, `DWDT`, `DTDT`, `DTDTFRIC`, `DPDTMST`, `DQDT`, `DQLDT`, `DQIDT`, `DQRDT`, `DQSDT`, `DQGDT`. In the full model the epilogue's `if (associated(...))` blocks therefore always run. These are not among the 197, but a standalone MAPL3 Moist will differ unless they are connected or activated.

## Suspected upstream bugs

Record these, but keep them unchanged in the port so results stay bit-for-bit with MAPL2:

- **`LCL_AGL` is never computed.** `BUOYANCY2` ignores its `LCL` argument (`PL:1268-1409`), so the export is always 0.
- **`WQT_DC` is always 0 at `LM = 72`.** The `GF_GEOS5_Interface` path never computes it, and Turbulence imports it.
- **`WTHV2` and `WQL` are 0 unless `PDFSHAPE == 6`.** Turbulence imports `WTHV2`. This may be intended.
- **UW reads the previous step's GF precipitation** (B1). This may be intended.
- **Potential bug, unverified: `CNPCPRATE` may break restart reproducibility.**
  - MAPL2 does not checkpoint exports, so after a restart UW's first step reads `CNPCPRATE = 0` (the ALLOC zero fill). In a continuous run, the same step reads GF's value from the step before.
  - UW uses the value over ocean (`uwshcu.F90:323-325`, capped at `1e-5`). A run split at a restart could therefore differ from a continuous run.
  - Not checked: whether `cnvtr` changes UW's output in practice, and whether GEOS start/stop regression tests compare bit-for-bit and would have caught it.
  - The port keeps MAPL2's behaviour (`RESTART: SKIP`). Making it a restarted internal would fix it, but would change results relative to MAPL2.
- **`DQLDT_SC`/`DQIDT_SC` are wrong on non-ringing UW steps.** When `UW_DT` is longer than Moist's timestep, `UW:679-680` subtracts last step's tendency from the current condensate (B2). This reaches Chem through `DQLDT`/`DQIDT`.
- **`CLDTOP_SC` is passed twice under `#ifdef UWDIAG`** (`UW:531`, `:536`), aliasing two dummies of `compute_uwshcu_inv`. The default build is unaffected.

## Options for step 6

What each group needs under the two approaches discussed:

| Group | Count | Internals with `ADD2EXPORT` | Local arrays copied to the export when associated |
|---|---|---|---|
| A1 | 5 | Works | Breaks: the schemes need the driver's values |
| A2, A3 | 40 | Works | The driver loses terms unless the arrays reach it. `ICE`/`FRZR`/`LS_PRCP` are read unguarded (null dereference) |
| A4 | 4 | Works | Works, if the scheme copies to the export |
| B | 25 | Works | Breaks: needs storage that lasts across steps |
| C1 | 22 | Needs `FILL: 0.0` | Needs zero-initialization, plus `FILL: 0.0` on the connected exports |
| C2 | 2 | Works | Works, if the epilogue re-fetch uses the local |
| C3 | 8 | Not needed; drop `ALLOC` | Same. `*AFMOIST` waits for step 12 |
| C4 | 91 | Works | Works (zero-initialize anyway; it costs little) |

So the C groups can use local arrays, and the 74 in A and B need storage that outlives one routine call. There are two ways to provide it:

- **Internals with `ADD2EXPORT`, `RESTART: SKIP` and `FILL: 0.0`.**
  - This reproduces MAPL2 at the same memory cost and with no copy. The export is the internal field itself, re-exported from `<self>` under the same name (`@MAPL/superstructure/generic/specs/ComponentSpec.F90:136-152`).
  - The edits are small: the spec rows move from `EXPORT` to `INTERNAL`, and each fetch switches from `EXPORT` to `INTERNAL`.
  - Cost: MAPL3 honours `RESTART_SKIP` only when reading a restart. On write it drops only incomplete fields (`@MAPL/superstructure/generic/RestartHandler.F90:36-37, 93, 204-220`). Checkpoints would therefore carry about 59 extra 3D fields, against roughly 16 internals in Moist today. Restarts stay correct. Avoiding this needs a MAPL change to `RestartHandler%write` (porting notes, decision 1c).
  - Precedent: DynCore re-exports internals, including a 3D edge field, through the same ACG column, and AgcmSimple does it by hand (see Background). What would be new for Moist is `RESTART: SKIP` combined with `ADD2EXPORT`: every existing case uses `REQUIRED` or leaves `RESTART` empty.
  - Internal `RI` needs `ALIAS: CLDREFFI`, because ACG declares one pointer per short name and import `RI` exists.
- **A private derived-type state** (`_SET_NAMED_PRIVATE_STATE`/`_GET_NAMED_PRIVATE_STATE`, `@MAPL/include/MAPL_private_state.h:39-52`), the pattern Solar, IRRAD, SU2G and GWD use.
  - The 74 become allocatable arrays, allocated and zeroed in Initialize. Each writer copies its array to the export when the export is associated.
  - Nothing is added to the checkpoints.
  - Cost: a copy and a second buffer for each connected export, and more edits at every fetch site and in the driver epilogue.

**Decision:** C uses local arrays (2026-10-29). B (2026-10-29) and A (2026-10-30) use internals with `ADD2EXPORT`. See decisions 1, 1b and 1c in `mapl3-porting-notes.md`.

## Caveats

- This is static analysis: nothing was built or run, and the regex scan can miss unusual syntax. One such case was found and handled by hand.
- The out-of-scope schemes also fetch several of these names with ALLOC, for example `CN_PRCP` in BACM_1M, THOM_1M, MGB2_2M and NSSL_2M. Revisit them when those schemes are ported (step 11).
- That MAPL3 does not allocate unconnected exports and always allocates internals is taken from the porting notes (issue 2), not re-read for this report. The fill-value, `add_to_export` and restart behaviour were checked in the MAPL source.
- "Moist's timestep" is assumed to equal the model `DT`. Moist's own `RUN_DT` setting was not checked.
- Line numbers refer to `8bde1bc9` in `@GEOSgcm_GridComp`. `@GEOSgcm_App` line numbers refer to its current checkout.

## Appendix: all 197 exports

Fetch sites use the file abbreviations above, and groups are as defined above.
- **Other fetches** are fetches without `ALLOC`. An `IMPORT` fetch of the same name is not included.
- **MAPL2 consumers** are the components Physics connects to the export (`Physics:1542-1728`), plus the Physics parent's direct reads (`Physics:3118-3176`).

| Export | Group | `ALLOC=.TRUE.` fetches | Other fetches | MAPL2 consumers |
|---|---|---|---|---|
| `BYNCY` | A1 | Moist:604, GF:724 | - | Chem |
| `CNV_FRC` | A1 | Moist:603, GFDL_1M:684, GF:721, UW:423 | - | GWD, Chem |
| `EIS` | A1 | Moist:524, GFDL_1M:725, UW:282 | - | Turbulence |
| `LTS` | A1 | Moist:523, GFDL_1M:724, UW:281 | - | - |
| `SRF_TYPE` | A1 | Moist:458, GFDL_1M:685, GF:722, UW:424 | - | - |
| `DQADT_MACRO` | A2 | GFDL_1M:760 | Moist:996 | - |
| `DQADT_MICRO` | A2 | GFDL_1M:1086 | Moist:998 | - |
| `DQGDT_MACRO` | A2 | GFDL_1M:763 | Moist:983 | - |
| `DQGDT_MICRO` | A2 | GFDL_1M:1089 | Moist:985 | - |
| `DQIDT_MACRO` | A2 | GFDL_1M:758 | Moist:952 | - |
| `DQIDT_MICRO` | A2 | GFDL_1M:1084 | Moist:954 | - |
| `DQLDT_MACRO` | A2 | GFDL_1M:759 | Moist:939 | - |
| `DQLDT_MICRO` | A2 | GFDL_1M:1085 | Moist:941 | - |
| `DQRDT_MACRO` | A2 | GFDL_1M:761 | Moist:963 | - |
| `DQRDT_MICRO` | A2 | GFDL_1M:1087 | Moist:965 | - |
| `DQSDT_MACRO` | A2 | GFDL_1M:762 | Moist:974 | - |
| `DQSDT_MICRO` | A2 | GFDL_1M:1088 | Moist:976 | - |
| `DQVDT_MACRO` | A2 | GFDL_1M:757 | Moist:924 | - |
| `DQVDT_MICRO` | A2 | GFDL_1M:1083 | Moist:926 | - |
| `DTDT_MACRO` | A2 | GFDL_1M:766 | Moist:908 | - |
| `DTDT_MICRO` | A2 | GFDL_1M:1092 | Moist:910 | GWD |
| `DUDT_MACRO` | A2 | GFDL_1M:764 | Moist:882 | - |
| `DUDT_MICRO` | A2 | GFDL_1M:1090 | Moist:884 | - |
| `DVDT_MACRO` | A2 | GFDL_1M:765 | Moist:895 | - |
| `DVDT_MICRO` | A2 | GFDL_1M:1091 | Moist:897 | - |
| `AN_PRCP` | A3 | GFDL_1M:718 | Moist:1063, Moist:1074, Moist:1123, Moist:1143 | - |
| `AN_SNR` | A3 | GFDL_1M:721 | Moist:1091, Moist:1125 | - |
| `CN_PRCP` | A3 | GFDL_1M:717 | Moist:1048, Moist:1076, Moist:1101, Moist:1145 | Surface |
| `CN_SNR` | A3 | GFDL_1M:720 | Moist:1093, Moist:1105 | - |
| `FRZR` | A3 | Moist:1162, GFDL_1M:698 | Moist:1117 | Surface |
| `ICE` | A3 | Moist:1161, GFDL_1M:697 | Moist:1115 | Surface |
| `LS_PRCP` | A3 | Moist:854, GFDL_1M:695 | Moist:1061, Moist:1072, Moist:1111, Moist:1141 | - |
| `LS_SNR` | A3 | GFDL_1M:696 | Moist:1089, Moist:1113 | - |
| `MFD_DC` | A3 | GF:714 | Moist:782 | - |
| `MFD_SC` | A3 | UW:418 | Moist:783 | - |
| `PFI_AN` | A3 | GFDL_1M:704 | Moist:1009, Moist:1036 | - |
| `PFI_CN` | A3 | GF:747 | Moist:1005 | Chem |
| `PFI_LS` | A3 | GFDL_1M:705 | Moist:1011, Moist:1038 | - |
| `PFL_AN` | A3 | GFDL_1M:702 | Moist:1017, Moist:1027 | - |
| `PFL_CN` | A3 | GF:748 | Moist:1013 | Chem |
| `PFL_LS` | A3 | GFDL_1M:703 | Moist:1019, Moist:1029 | - |
| `SC_PRCP` | A3 | GFDL_1M:719 | Moist:1050, Moist:1078, Moist:1130, Moist:1147 | - |
| `SC_SNR` | A3 | GFDL_1M:722 | Moist:1095, Moist:1132 | - |
| `UMF_DC` | A3 | GF:715 | Moist:742 | - |
| `UMF_SC` | A3 | UW:426 | Moist:743 | - |
| `QITOT` | A4 | UW:445 | Moist:1198 | Turbulence, Chem |
| `QLTOT` | A4 | UW:444 | Moist:1195 | Turbulence, Chem |
| `RI` | A4 | GFDL_1M:682 | Moist:500 | Radiation, Chem |
| `RL` | A4 | GFDL_1M:681 | Moist:498 | Radiation, Chem |
| `CNPCPRATE` | B1 | GF:726, UW:422 | Moist:1052, Moist:1080, Moist:1103 | - |
| `DQADT_DC` | B2 | GF:586 | Moist:992 | - |
| `DQADT_SC` | B2 | UW:279 | Moist:994 | - |
| `DQIDT_DC` | B2 | GF:585 | Moist:948 | - |
| `DQIDT_SC` | B2 | UW:275 | Moist:950 | - |
| `DQLDT_DC` | B2 | GF:584 | Moist:935 | - |
| `DQLDT_SC` | B2 | UW:276 | Moist:937 | - |
| `DQRDT_SC` | B2 | UW:277 | Moist:961 | - |
| `DQSDT_SC` | B2 | UW:278 | Moist:972 | - |
| `DQVDT_DC` | B2 | GF:583 | Moist:920 | - |
| `DQVDT_SC` | B2 | UW:274 | Moist:922 | - |
| `DTDT_DC` | B2 | GF:582 | Moist:904 | GWD |
| `DTDT_SC` | B2 | UW:273 | Moist:906 | - |
| `DUDT_DC` | B2 | GF:580 | Moist:878 | - |
| `DUDT_SC` | B2 | UW:271 | Moist:880 | - |
| `DVDT_DC` | B2 | GF:581 | Moist:891 | - |
| `DVDT_SC` | B2 | UW:272 | Moist:893 | - |
| `QIDET_SC` | B2 | UW:420, UW:642 | - | - |
| `QIENT_SC` | B2 | UW:545, UW:646 | - | - |
| `QISUB_SC` | B2 | UW:431, UW:645 | - | - |
| `QLDET_SC` | B2 | UW:419, UW:641 | - | - |
| `QLENT_SC` | B2 | UW:544, UW:644 | - | - |
| `QLSUB_SC` | B2 | UW:430, UW:643 | - | - |
| `SHLW_PRC3` | B2 | GFDL_1M:797, UW:596 | - | - |
| `SHLW_SNO3` | B2 | GFDL_1M:811, UW:598 | - | - |
| `CNV_TOPP_DP` | C1 | GF:767 | - | - |
| `CNV_TOPP_MD` | C1 | GF:768 | - | - |
| `CNV_TOPP_SH` | C1 | GF:769 | - | - |
| `ENTR_DP` | C1 | GF:759 | - | - |
| `ENTR_MD` | C1 | GF:760 | - | - |
| `ENTR_SH` | C1 | GF:761 | - | - |
| `LCL_AGL` | C1 | Moist:615 | - | - |
| `LFR_GF` | C1 | GF:734 | - | - |
| `MDNDP` | C1 | GF:736 | - | - |
| `NIFA` | C1 | Moist:602 | - | - |
| `NWFA` | C1 | Moist:601 | - | - |
| `PRCP_WATER` | C1 | GFDL_1M:689 | - | - |
| `SGS_VVEL_DP` | C1 | GF:763 | - | - |
| `SGS_VVEL_MD` | C1 | GF:764 | - | - |
| `SGS_VVEL_SH` | C1 | GF:765 | - | - |
| `VFALL_GRAUPEL` | C1 | GFDL_1M:714 | - | - |
| `VFALL_ICE` | C1 | GFDL_1M:712 | - | - |
| `VFALL_RAIN` | C1 | GFDL_1M:715 | - | - |
| `VFALL_SNOW` | C1 | GFDL_1M:713 | - | - |
| `WQL` | C1 | GFDL_1M:710 | - | - |
| `WQT_DC` | C1 | GF:758 | - | Turbulence |
| `WTHV2` | C1 | GFDL_1M:709 | - | Turbulence |
| `DQVDT_ER` | C2 | Moist:849 | Moist:928 | - |
| `DTDT_ER` | C2 | Moist:848 | Moist:912 | - |
| `PREC_CONV` | C3 | Moist:1157 | - | - |
| `PREC_STRAT` | C3 | Moist:1153 | - | - |
| `QAFMOIST` | C3 | Moist:1242, Moist:1320 | - | Physics (parent) |
| `SAFMOIST` | C3 | Moist:1248, Moist:1326 | - | Physics (parent) |
| `TAFMOIST` | C3 | Moist:1239, Moist:1317 | - | Physics (parent) |
| `THAFMOIST` | C3 | Moist:1245, Moist:1323 | - | Physics (parent) |
| `UAFMOIST` | C3 | Moist:1233, Moist:1311 | - | Physics (parent) |
| `VAFMOIST` | C3 | Moist:1236, Moist:1314 | - | Physics (parent) |
| `AA0` | C4 | GF:749 | - | - |
| `AA1` | C4 | GF:750 | - | - |
| `AA1_BL` | C4 | GF:753 | - | - |
| `AA1_CIN` | C4 | GF:754 | - | - |
| `AA2` | C4 | GF:751 | - | - |
| `AA3` | C4 | GF:752 | - | - |
| `CAPE` | C4 | Moist:605 | - | Chem |
| `CBMF_SC` | C4 | UW:289 | - | - |
| `CFICE` | C4 | Moist:807 | - | - |
| `CLDTOP_SC` | C4 | UW:290 | - | - |
| `CNV_CVW` | C4 | GF:716 | - | Chem |
| `CNV_DQCDT` | C4 | GF:718 | - | - |
| `CNV_MF0` | C4 | GF:725 | - | - |
| `CNV_PRC3` | C4 | GF:719 | - | - |
| `CNV_QC` | C4 | GF:717 | - | Chem |
| `CNV_UPDF` | C4 | GF:720 | - | - |
| `CUFRC_SC` | C4 | UW:421 | - | - |
| `DCM_SC` | C4 | UW:427 | - | - |
| `DETR_SC` | C4 | UW:429 | - | - |
| `DNDCNV` | C4 | GF:771 | - | - |
| `DNICNV` | C4 | GF:772 | - | - |
| `ENTLAM` | C4 | GF:727 | - | - |
| `ENTR_SC` | C4 | UW:428 | - | - |
| `ERRDP` | C4 | GF:742 | - | - |
| `ERRMD` | C4 | GF:744 | - | - |
| `ERRSH` | C4 | GF:743 | - | - |
| `ER_PRCP` | C4 | Moist:855 | - | - |
| `EVAPC` | C4 | GFDL_1M:687 | - | - |
| `FCLD` | C4 | GFDL_1M:674 | - | Turbulence, Radiation, Chem |
| `INHB` | C4 | Moist:606 | - | Chem |
| `MFDP` | C4 | GF:739 | - | - |
| `MFMD` | C4 | GF:741 | - | - |
| `MFSH` | C4 | GF:740 | - | - |
| `MLCAPE` | C4 | Moist:609 | - | - |
| `MLCIN` | C4 | Moist:610 | - | - |
| `MUCAPE` | C4 | Moist:611 | - | - |
| `MUCIN` | C4 | Moist:612 | - | - |
| `MUPDP` | C4 | GF:735 | - | - |
| `MUPMD` | C4 | GF:738 | - | - |
| `MUPSH` | C4 | GF:737 | - | - |
| `PBUP_SC` | C4 | UW:288 | - | - |
| `PDFITERS` | C4 | GFDL_1M:711 | - | - |
| `PFG_LS` | C4 | GFDL_1M:708 | - | - |
| `PFR_LS` | C4 | GFDL_1M:706 | - | - |
| `PFS_LS` | C4 | GFDL_1M:707 | - | - |
| `PINV_SC` | C4 | UW:286 | - | - |
| `PLCL_SC` | C4 | UW:284 | - | - |
| `PLFC_SC` | C4 | UW:285 | - | - |
| `PRCP_GRAUPEL` | C4 | GFDL_1M:693 | - | - |
| `PRCP_ICE` | C4 | GFDL_1M:692 | - | - |
| `PRCP_RAIN` | C4 | GFDL_1M:690 | - | - |
| `PRCP_SNOW` | C4 | GFDL_1M:691 | - | - |
| `PREL_SC` | C4 | UW:287 | - | - |
| `PTYPE` | C4 | Moist:1164 | - | - |
| `QG` | C4 | GFDL_1M:680 | - | Radiation |
| `QI` | C4 | GFDL_1M:677 | - | Radiation, Chem |
| `QL` | C4 | GFDL_1M:676 | - | Radiation, Chem |
| `QPERT_SC` | C4 | UW:435 | - | - |
| `QR` | C4 | GFDL_1M:678 | - | Radiation |
| `QS` | C4 | GFDL_1M:679 | - | Radiation |
| `QTFLX_SC` | C4 | UW:436 | - | - |
| `QV` | C4 | GFDL_1M:675 | - | Turbulence, Radiation |
| `RAIN` | C4 | Moist:1069 | - | - |
| `REV_CN` | C4 | GF:746 | - | Chem |
| `REV_LS` | C4 | GFDL_1M:700 | - | Chem |
| `RHCRIT` | C4 | GFDL_1M:855 | - | - |
| `RHX` | C4 | GFDL_1M:699 | - | - |
| `RSU_CN` | C4 | GF:745 | - | - |
| `RSU_LS` | C4 | GFDL_1M:701 | - | - |
| `SBCAPE` | C4 | Moist:607 | - | - |
| `SBCIN` | C4 | Moist:608 | - | - |
| `SC_NDROP` | C4 | UW:432 | - | - |
| `SC_NICE` | C4 | UW:433 | - | - |
| `SIGMA_DEEP` | C4 | GF:728 | - | - |
| `SIGMA_MID` | C4 | GF:729 | - | - |
| `SLFLX_SC` | C4 | UW:437 | - | - |
| `SNO` | C4 | Moist:1086 | - | Surface |
| `SUBLC` | C4 | GFDL_1M:688 | - | - |
| `TAU_BL` | C4 | GF:755 | - | - |
| `TAU_DP` | C4 | GF:756 | - | - |
| `TAU_MD` | C4 | GF:757 | - | - |
| `TPERT_SC` | C4 | UW:434 | - | - |
| `TPREC` | C4 | Moist:1138 | - | Surface, Chem |
| `TPWI` | C4 | GF:730 | - | - |
| `TPWI_STAR` | C4 | GF:731 | - | - |
| `UFLX_SC` | C4 | UW:438 | - | - |
| `UMST0` | C4 | Moist:492 | - | - |
| `VFLX_SC` | C4 | UW:439 | - | - |
| `VMST0` | C4 | Moist:494 | - | - |
| `ZLFC` | C4 | Moist:613 | - | Chem |
| `ZLNB` | C4 | Moist:614 | - | - |
