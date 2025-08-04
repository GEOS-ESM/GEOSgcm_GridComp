# Catchment-CN5.1

Jana Kolassa, Rolf H. Reichle, and Randal D. Koster

August 2025

## 1 Introduction

Catchment-CN is a hybrid land surface model that incorporates routines from the Community Terrestrial System Model (CTSM) into the NASA Catchment land surface model. While water and energy balance calculations are handled by Catchment, CTSM routines are used for carbon and nitrogen dynamics, which include processes such as photosynthesis, phenology, decomposition, and wildfire.

## 2 Catchment-CN5.1

Catchment-CN5.1 was constructed by merging the Catchment model with a pre-release of CTSM version 5.1 (tag: *branch_tags/PPE.n08_ctsm5.1.dev023*), which is available through the ESCOMP GitHub repository:

- [https://github.com/ESCOMP/CTSM/tags?after=ctsm5.1.dev030](https://github.com/ESCOMP/CTSM/tags?after=ctsm5.1.dev030)
- [https://github.com/ESCOMP/CTSM/releases/tag/branch\_tags/PPE.n08\_ctsm5.1.dev023](https://github.com/ESCOMP/CTSM/releases/tag/branch_tags/PPE.n08_ctsm5.1.dev023)

Information on the CTSM subroutines that are mentioned in this documentation is available through the [CLM50\_Tech\_Note](https://escomp.github.io/CTSM/release-clm5.0/tech_note/index.html).

The initial version of Catchment-CN5.1 is consistent with the Catchment model of release v2.7.5 of the GEOSgcm\_GridComp repository:

- [https://github.com/GEOS-ESM/GEOSgcm\_GridComp/releases/tag/v2.7.5](https://github.com/GEOS-ESM/GEOSgcm\_GridComp/releases/tag/v2.7.5)

Catchment-CN5.1 will be updated along with the original version of Catchment with the continued development of the GEOSgcm\_GridComp repository and associated repositories.

## 3 Structural implementation

In Catchment-CN5.1, the 'split PFTs' used by previous versions of Catchment-CN have been removed. This means there are now a total of 15 PFTs (corresponding to the 15 CTSM non-crop-model PFTs) instead of the 19 PFTs used in previous Catchment-CN versions. This also means that each model tile only has up to 2 PFTs instead of up to 4 as before to accommodate PFT splitting.

Catchment and CTSM are linked within the *GEOS_CatchCNCLM51GridComp.F90* module through three key connection points.

1. *CN_init* initializes the CTSM variables used in Catchment-CN5.1.

2. *catchcn_calc_rc* connects to the CTSM photosynthesis modules, executing at the GEOS model "heartbeat" time step (currently 5 or 7.5 minutes). The heartbeat time step allows for the provision of important heartbeat-scale variations in the canopy conductance to the Catchment energy balance calculations.

3. *CN_Driver* links to all other CTSM modules (phenology, decomposition, wildfire, etc.) on a 90-minute cycle.

Each connection operates through specialized "bridging modules" (section 3.1) that translate between the Catchment and CTSM spaces, handling variable mapping, unit conversions, and spatial translations between the Catchment tile space and the CTSM column/patch space (section 3.2). These bridging modules include *CNCLM_init_mod.F90*, *CNCLM_Photosynthesis.F90*, and *CNCLM_Driver.F90*, each responsible for different aspects of the Catchment-CN integration. Each of the bridging modules is discussed in detail below. 

### 3.1 Bridging Routines

**CNCLM_init_mod.F90:** This module handles the initialization of the CTSM variables, the reading of the CTSM configuration from the *CN_CLM51.nml* namelist file, and the reading of the CTSM parameter file. It contains the subroutine *CN_init*, which calls the CTSM initialization routines. These routines have in some cases been adapted to pass in the Catchment-CN restart variables CNCOL and CNPFT in order to initialize with values from the restart file. Modules that have been modified in this manner are listed in Table B1 in Appendix B. This module and the *CN_init* subroutine are called once during the initialization phase.  

**CNCLM_Photosynthesis.F90:** This module handles the call to the CTSM Photosynthesis module. It contains the subroutine *catchcn\_calc\_rc*, which is organized as follows: 

1. Perform the mapping of Catchment variables needed for the photosynthesis calculations to their CTSM equivalents. 
2. Call CTSM subroutines that calculate the sunlit and shaded fraction of the canopy (*TwoStream* and *CanopySunShadeFracs*), information that is needed in the photosynthesis calculations.
3. Call the CTSM *PhotosynthesisHydraulicStress* routine, which performs the main photosynthesis calculations.
4. Map the CTSM photosynthesis outputs back into the Catchment space.

The above steps are performed at the model "heartbeat" time step (currently 5 or 7.5 minutes).

**CNCLM_DriverMod.F90:** This module handles the call to the non-photosynthesis CTSM routines and is engaged every 90 minutes. It contains the subroutine *CN_Driver*, which:

1. maps Catchment variables to their CTSM equivalents,
2. calls the CTSM routines for the non-photosynthesis calculations (phenology, decomposition, wildfire, etc.) including the carbon and nitrogen balance checks, and
3. maps the CTSM outputs back to Catchment space.

The module also contains the subroutine *CN_exit*, which maps CTSM variables into the arrays CNCOL and CNPFT that are written into the CNCLM restart file.

**CN2CLMType.F90:** Defines a special type used to pass forcing data from Catchment into the nested CTSM fire data types.

**update_model_para4cn.F90:** Utility module for updating current time. Introduced in previous versions of Catchment-CN.

### 3.2 Catchment to CTSM mapping

In addition to assigning the values of Catchment variables to their CTSM equivalents, three additional mapping steps are needed to use Catchment inputs within the CTSM modules.

1. The first step only applies to the variables passed to *CN_Driver* and consists of computing 90-minute average states and fluxes from the 5- or 7.5-minute Catchment states and fluxes. (At present, these calculations are done before the call to subroutine *catchmentCN*, and as a result, extra model diagnostics calculations and restart variables are needed.  To avoid this overhead and duplication of code segments, it is planned to calculate the averages after the call to catchment in future versions of Catchment-CN5.1.)  
 
2. The second step involves mapping hydrologic and temperature variables from the three dynamic hydrologic zones used in Catchment to the three static “carbon zones” used in Catchment-CN, which correspond to valley bottoms (10%), hillslopes (45%) and hilltops (45%). This mapping happens prior to the call to *CN_Driver* and uses an areal weighting approach (see Figure 2 of Koster et al. 2014). 

3. The third mapping step relates to the different organization of array variables in Catchment(-CN) and CTSM. PFT-level variables in Catchment-CN typically have the dimensions (ntiles, nveg, nzone), where ntiles is the number of tiles, nveg=2 is the maximum number of PFTs per tile, and nzone=3 is the number of carbon zones. In CTSM, the equivalent variable would have the dimension (1:ntiles\*maxpft\*nzones), where maxpft=15 is the total number of PFTs in the model. In a CTSM array, only the entries corresponding to PFTs that are present in a location have data values. Computations are only performed for the array entries that have data values, which is implemented through indices (called 'filters' in CTSM). This mapping from the Catchment-CN order to the CTSM order is handled by the bridging routines, which also handle the mapping of the CTSM outputs back into the Catchment-CN order. The layout of a CTSM variable array is schematically illustrated in Table B2 in Appendix B.

## 4 Experiment setup and model configuration

Configuring and setting up a Catchment-CN5.1 experiment is handled through ldas\_setup or gcm\_setup (to be implemented), the basic approach also used when running Catchment or previous versions of Catchment-CN. To choose Catchment-CN5.1, set LSM\_CHOICE = 4 and LAND\_PARAMS = CN\_CLM51. As with previous versions of Catchment-CN, it is also necessary to choose an option for the source of atmospheric CO2 by modifying the ATM\_CO2 setting.

CTSM also offers a choice of multiple schemes for different parts of the model code, such as the wildfire calculations, the stomatal conductance scheme, or the decomposition model. These elements of the model configuration are controlled through the *CN_CLM51.nml* file, which can be modified in the ./run directory after setting up an experiment.  Note, however, that configuration changes may require additional spinup of the hydrological or carbon states.  For more permanent changes to the default configuration, the original file located in the ./GEOScatchCNCLM51\_GridComp directory needs to be modified. Changes made here will only take effect after re-compiling the model as the compilation moves the files to the ./install/etc/ directory, from where *ldas_setup* moves it to the ./run directory.  

Finally, the model configuration can be changed by modifying the CTSM parameters. These can be found in the *ctsm51_params.c210923_forCNCLM.nc* file in the ./input directory. For Catchment-CN5.1, a few parameters were changed from their original CTSM default values (section 5).  

A set of spun-up restart files for Catchment-CN5.1 will be available, and these can be used to initialize an experiment. If no restart file is available at the required resolution, use the GEOS *remap_restarts* functionality to map an existing restart file to the desired resolution.  Note that remapping a Catchment[CN] restart to a different resolution generally requires additional model spin-up.   

## 5 Science changes from the original CTSM version

A few science changes have been made to the original CTSM implementation to adapt the CTSM routines for use with Catchment. These include:

**'fcur' parameter:** 'fcur' is the fraction of carbon allocation that goes to currently displayed growth. In previous versions of Catchment-CN, fcur was set to 0.5 for the PFTs that, in CLM4 or CLM4.5, have this parameter set to 0. This change was made to address problems with LAI overshoot in years following years with no freeze-related offset as well as to bring temporal phenological variations more in line with satellite-based estimates. This same modification has been adopted here. This change has been made directly in the parameter file (*ctsm51_params.c210923_forCNCLM.nc*).

**'leaf_long' parameter:** The leaf longevity was set to be at least 1 year in previous versions of Catchment-CN. This modification has been adopted here for consistency.

**Water stress threshold:** In the original CTSM, the water stress threshold is set through a global constant. In Catchment-CN5.1 the water stress threshold at a given location is set equal to the wilting point at that location.

**Solar Induced Fluorescence:** The original CTSM does not compute solar induced fluorescence (SIF). In Catchment-CN, the capability to compute SIF was added by including the fluorescence routine developed by Jung-Eun Lee (Brown) using the approach of van der Tol and Berry (2012). This fluorescence routine was initially implemented for Catchment-CN by Greg Walker and has been adapted by Jana Kolassa for Catchment-CN5.1. The adaption to Catchment-CN5.1 accounts for the fact that photosynthesis calculations in CTSM5.1 are separated for the sunlit and shaded parts of the canopy and the number of canopy layers above snow.

**Heterotrophic respiration calculations:** In the original CTSM tag that was used to develop Catchment-CN5.1 (section 2), the heterotrophic respiration (HR) is calculated as the sum of the litter HR and the soil organic matter HR only. That is, HR does not include the coarse woody debris HR, because the original model parameters were such that the coarse woody debris HR was always 0. The parameter file used for Catchment-CN5.1 (*ctsm51_params.c210923_forCNCLM.nc*) corresponds to a slightly newer version of CTSM, with the parameter values set such that the coarse woody debris HR is not always 0. Consequently, HR from coarse woody debris was added as an additional input to the total HR calculation, consistent with newer versions of CTSM.

**Number of soil layers:** All versions of Catchment-CN assume only one soil layer for CLM/CTSM calculations. This contrasts with the standard CLM/CTSM framework, which supports multiple soil layer configurations and indeed requires a minimum of 5 layers.

**Root-zone wetness in fire code:** Instead of calculating the root-zone wetness inside the fire code as is done in the original CTSM, we instead use the root-zone wetness calculated by Catchment.

## Appendix A: Restarting Catchment-CN5.1 from a Catchment-CN4.5 restart file

**IMPORTANT: The instructions in Appendix A are for developers only! They do NOT apply for setting up a science experiment.

If for some reason it becomes necessary to initialize an experiment from a Catchment-CN4.5 restart file, a few extra steps must be taken.

First, the Catchment-CN4.5 restart file needs to be modified, such that it reflects the PFT-distribution of Catchment-CN5.1 and the fact that in Catchment-CN5.1 each tile only has up to 2 PFTs as opposed to 4 in previous Catchment-CN versions. The variables affected by this are the PFT (ITY), the PFT fraction (FVG), the array of PFT-level variables (CNPFT), and the sunlit and shaded photosynthesis (PSNSUNM and PSNSHAM). Additionally, the dimension on the CNPFT array needs to be padded with zeroes to reflect the larger number of variables saved to it in Catchment-CN5.1.

Next, the toggles *init_accum* and *init_accum_365* need to be set to *.true.* in *GEOS_CatchCNCLM51GridComp.F90*. These toggles control how the multi-day average variables are computed before their nominal accumulation period is reached. *init_accum* needs to be set to *.false.* after one month and *init_accum_365* needs to be set to *.false.* after 365 days. Additionally, the toggle *no_cn51_rst* needs to be set to *.true.* when starting the experiment. This toggle allows initialization of the nitrate and ammonia variables from the total nitrogen content when using a Catchment-CN4.5 file. This toggle needs to be switched back to *.false.* after the first time that a restart file is written.

## Appendix B: CN\_CLM51 files

Files that are specific to Catchment-CN5.1 are in the ./GEOScatchCNCLM51\_GridComp/ directory. Within this directory, the CTSM files needed to build Catchment-CN5.1 are in the ./CLM51/ directory, which also contains the CLM-to-CN bridging modules (section 3.1). Most of the CTSM files used in Catchment-CN5.1 have been modified from their original CTSM version. For reference, the original, unaltered CTSM files are available in the ./CLM51\_orig\_files/ directory. The modifications can be grouped into the categories outlined below. Table B1 lists all the files in ./CLM51/ and the associated categories.

### Modification Categories

**C1:** The initialization of variables was modified to use the CNCOL or CNPFT for saving CTSM variables to the CNCLM restart files. And/or the initialization was simplified to combine the variable allocation and initialization in a single Init subroutine. And/or the subroutines to initialize from CTSM restart files were removed. These changes have no science impact.

**C2:** The module was modified to load MAPL and/or ESMF types or variables (e.g., "MAPL\_R8" instead of CTSM's "r8"). These changes have no science impact.

**C3:** Code changes were made that have scientific impact.

**C4:** Unused parts of the CTSM code were commented out or removed.

**C5:** Changes were made to use MAPL routines instead of CTSM routines (often for reading files), but there is no science impact from this change.

**C6:** A new custom CLM-to-CN bridging module (section 3.1).

**C7:** Code changes were made to match Catchment-CN requirements or to be consistent with previous versions of Catchment-CN.

**C8:** External files.

### Table B1: List of modules in the ./CLM51 directory

| Module | Category |
|--------|----------|
| abortutils.F90 | C4 |
| ActiveLayerMod.F90 | C1, C4 |
| AnnualFluxDribbler.F90 | C4 |
| atm2lndType.F90 | C1, C2, C4 |
| CanopyStateType.F90 | C1, C2 |
| ch4Mod.F90 | C1, C2, C4 |
| clm\_time\_manager.F90 | C2, C4, C7 |
| clm\_varcon.F90 | C2, C3, C4, C7 |
| clm\_varctl.F90 | C2, C4, C7 |
| clm\_varpar.F90 | C2, C4, C7 |
| cmake/genf90\_utils.cmake | C8 (from GitHub CESM-Development repository; tag: CMake\_Fortran\_utils\_150308) |
| CN2CLMType.F90 | C6 |
| CNAnnualUpdateMod.F90 | No changes |
| CNBalanceCheckMod.F90 | C1, C5 |
| CNCLM\_DriverMod.F90 | C6 |
| CNCLM\_init\_mod.F90 | C6 |
| CNCLM\_Photosynthesis.F90 | C6 |
| CNCStateUpdate1Mod.F90 | C4 |
| CNCStateUpdate2Mod.F90 | C4 |
| CNCStateUpdate3Mod.F90 | C4 |
| CNDriverMod.F90 | C4 |
| CNDVType.F90 | C1, C4 |
| CNFireBaseMod.F90 | C2, C4 |
| CNFireEmissionsMod.F90 | C1 |
| CNFireFactoryMod.F90 | No changes |
| CNFireLi2014Mod.F90 | C3, C4, C7 |
| CNFireLi2016Mod.F90 | C3, C4, C7 |
| CNFireLi2021Mod.F90 | C3, C4, C7 |
| CNFireNoFireMod.F90 | C4 |
| CNFUNMod.F90 | C2, C4 |
| CNGapMortalityMod.F90 | C4 |
| CNGRespMod.F90 | No changes |
| CNMRespMod.F90 | C4 |
| CNNDynamicsMod.F90 | C4 |
| CNNStateUpdate1Mod.F90 | C4 |
| CNNStateUpdate2Mod.F90 | C4 |
| CNNStateUpdate3Mod.F90 | C4 |
| CNPhenologyMod.F90 | C3, C4 |
| CNPrecisionControlMod.F90 | No changes |
| CNProductsMod.F90 | C1, C2, C4 |
| CNRootDynMod.F90 | No changes |
| CNSharedParamsMod.F90 | C2, C5, C7 |
| CNVegCarbonFluxType.F90 | C1, C2, C4 |
| CNVegCarbonStateType.F90 | C1, C2, C4 |
| CNVegetationFacade.F90 | C1, C2, C4 |
| CNVegNitrogenFluxType.F90 | C1, C2, C4 |
| CNVegNitrogenStateType.F90 | C1, C2, C4 |
| CNVegStateType.F90 | C1, C2, C4 |
| CNVegStructUpdateMod.F90 | No changes |
| ColumnType.F90 | C1, C2, C4 |
| column\_varcon.F90 | C4 |
| CropType.F90 | C1, C2, C4 |
| decompMod.F90 | C1, C2, C4 |
| dynSubgridControlMod.F90 | C1, C4 |
| EnergyFluxType.F90 | C1, C4 |
| fileutils.F90 | No changes |
| filterColMod.F90 | No changes |
| filterMod.F90 | C1, C2, C4 |
| FireDataBaseType.F90 | C2, C4 |
| FireMethodType.F90 | C4 |
| FrictionVelocityMod.F90 | C1, C2, C4 |
| genf90.pl | C8 (from GitHub PARALLELIO repository; tag: unknown) |
| GridcellType.F90 | C1, C2, C4 |
| initSubgridMod.F90 | No changes |
| initVerticalMod.F90 | C4 |
| LandunitType.F90 | C1, C2, C4 |
| landunit\_varcon.F90 | No changes |
| ncdio\_pio.F90 | C2, C5 (CTSM original file name: ncdio\_pio.F90.in) |
| NutrientCompetitionCLM45defaultMod.F90 | C4 |
| NutrientCompetitionFactoryMod.F90 | No changes |
| NutrientCompetitionFlexibleCNMod.F90 | C1, C2, C4 |
| NutrientCompetitionMethodMod.F90 | No changes |
| OzoneBaseMod.F90 | C1, C2, C4 |
| paramUtilMod.F90 | C4, C5 |
| PatchType.F90 | C1, C2 |
| perf\_mod.F90 | C4 |
| pftconMod.F90 | C1, C2, C7 |
| PhotosynthesisMod.F90 | C1, C2, C7 |
| QSatMod.F90 | No changes |
| quadraticMod.F90 | No changes |
| RootBiophysMod.F90 | C4 |
| SaturatedExcessRunoffMod.F90 | C2, C4 |
| shr\_abort\_mod.F90 | C2, C4, C5 |
| shr\_assert.h | C7 |
| shr\_assert\_mod.F90.in | C4 |
| shr\_const\_mod.F90 | C2, C5 |
| shr\_file\_mod.F90 | No changes |
| shr\_fire\_emis\_mod.F90 | C4 |
| shr\_infnan\_mod.F90.in | C4 |
| shr\_kind\_mod.F90 | C2, C5 |
| shr\_log\_mod.F90 | C7 |
| shr\_mpi\_mod.F90 | C4 |
| shr\_nl\_mod.F90 | C4 |
| shr\_sys\_mod.F90 | C4 |
| SoilBiogeochemCarbonFluxType.F90 | C2, C3, C4 |
| SoilBiogeochemCarbonStateType.F90 | C1, C2, C4 |
| SoilBiogeochemCompetitionMod.F90 | C5 |
| SoilBiogeochemDecompCascadeBGCMod.F90 | C4 |
| SoilBiogeochemDecompCascadeCNMod.F90 | C4 |
| SoilBiogeochemDecompCascadeConType.F90 | C2, C4 |
| SoilBiogeochemDecompMod.F90 | C4 |
| SoilBiogeochemLittVertTranspMod.F90 | C4 |
| SoilBiogeochemNitrifDenitrifMod.F90 | No changes |
| SoilBiogeochemNitrogenFluxType.F90 | C1, C2, C4 |
| SoilBiogeochemNitrogenStateType.F90 | C1, C2, C4 |
| SoilBiogeochemNLeachingMod.F90 | No Changes |
| SoilBiogeochemNStateUpdate1Mod.F90 | No changes |
| SoilBiogeochemPotentialMod.F90 | No changes |
| SoilBiogeochemPrecisionControlMod.F90 | No changes |
| SoilBiogeochemStateType.F90 | C1, C2 |
| SoilBiogeochemVerticalProfileMod.F90 | No changes |
| SoilStateInitTimeConstMod.F90 | C4 |
| SoilStateType.F90 | C2, C4 |
| SoilWaterRetentionCurveMod.F90 | No changes |
| SolarAbsorbedType.F90 | C2, C4 |
| spmdMod.F90 | C7 |
| subgridAveMod.F90 | C4 |
| SurfaceAlbedoMod.F90 | C4 |
| SurfaceAlbedoType.F90 | C1, C2, C4 |
| SurfaceRadiationMod.F90 | C4 |
| TemperatureType.F90 | C1, C2, C4 |
| TridiagonalMod.F90 | No changes |
| update\_model\_para4cn.F90 | C6 |
| Wateratm2lndBulkType.F90 | C4 |
| Wateratm2lndType.F90 | C1, C4 |
| WaterDiagnosticBulkType.F90 | C1, C2, C4 |
| WaterDiagnosticType.F90 | C4 |
| WaterFluxBulkType.F90 | C1, C2, C4 |
| WaterFluxType.F90 | C1, C2, C4 |
| WaterStateBulkType.F90 | C1, C5 |
| WaterStateType.F90 | C1, C4 |
| WaterType.F90 | C1, C4 |

### Table B2: Illustration of the layout of a CTSM variable array

In this example, tile 1 contains PFTs 7 and 10 and has data values for the corresponding entries. Only the gray shaded data is present in the CTSM variable array, the other columns of the table serve to illustrate the layout.

| Array Index | Tile | Zone | PFT | Data |
|-------------|------|------|-----|------|
| 1 | 1 | 1 | 1 | NaN |
| 2 | 1 | 1 | 2 | NaN |
| 3 | 1 | 1 | 3 | NaN |
| 4 | 1 | 1 | 4 | NaN |
| 5 | 1 | 1 | 5 | NaN |
| 6 | 1 | 1 | 6 | NaN |
| **7** | **1** | **1** | **7** | **Some data** |
| 8 | 1 | 1 | 8 | NaN |
| 9 | 1 | 1 | 9 | NaN |
| **10** | **1** | **1** | **10** | **Some data** |
| 11 | 1 | 1 | 11 | NaN |
| 12 | 1 | 1 | 12 | NaN |
| 13 | 1 | 1 | 13 | NaN |
| 14 | 1 | 1 | 14 | NaN |
| 15 | 1 | 1 | 15 | NaN |
| 16 | 1 | 2 | 1 | NaN |
| 17 | 1 | 2 | 2 | NaN |
| 18 | 1 | 2 | 3 | NaN |
| 19 | 1 | 2 | 4 | NaN |
| 20 | 1 | 2 | 5 | NaN |
| 21 | 1 | 2 | 6 | NaN |
| **22** | **1** | **2** | **7** | **Some data** |
| 23 | 1 | 2 | 8 | NaN |
| 24 | 1 | 2 | 9 | NaN |
| **25** | **1** | **2** | **10** | **Some data** |
| 26 | 1 | 2 | 11 | NaN |
| 27 | 1 | 2 | 12 | NaN |
| 28 | 1 | 2 | 13 | NaN |
| 29 | 1 | 2 | 14 | NaN |
| 30 | 1 | 2 | 15 | NaN |
| ... | ... | ... | ... | ... |
