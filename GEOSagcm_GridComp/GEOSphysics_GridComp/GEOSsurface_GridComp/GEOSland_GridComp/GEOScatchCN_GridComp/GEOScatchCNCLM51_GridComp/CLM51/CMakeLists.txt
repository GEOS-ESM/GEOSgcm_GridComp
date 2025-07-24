esma_set_this ()
string (REPLACE GEOScatchCN_GridComp_ "" is_openmp ${this})

list (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
include(genf90_utils)

find_program(GENF90 genf90.pl PATHS ${CMAKE_CURRENT_LIST_DIR})

set(genf90_files shr_infnan_mod.F90.in shr_assert_mod.F90.in)

process_genf90_source_list("${genf90_files}" ${CMAKE_CURRENT_BINARY_DIR}
  share_genf90_sources)

set (srcs
  abortutils.F90
  AnnualFluxDribbler.F90
  clm_time_manager.F90
  clm_varcon.F90
  clm_varctl.F90
  clm_varpar.F90
  CN2CLMType.F90
  CNAnnualUpdateMod.F90
  CNBalanceCheckMod.F90
  CNCLM_Photosynthesis.F90
  ActiveLayerMod.F90
  atm2lndType.F90
  CanopyStateType.F90
  ch4Mod.F90
  CNDVType.F90
  CNFireBaseMod.F90
  CNProductsMod.F90
  CNVegCarbonFluxType.F90
  CNVegCarbonStateType.F90
  CNVegNitrogenFluxType.F90
  CNVegNitrogenStateType.F90
  CNVegStateType.F90
  ColumnType.F90
  CropType.F90
  decompMod.F90
  CNCLM_DriverMod.F90
  dynSubgridControlMod.F90
  EnergyFluxType.F90
  filterMod.F90
  FireDataBaseType.F90
  FrictionVelocityMod.F90
  GridcellType.F90
  initVerticalMod.F90
  LandunitType.F90
  OzoneBaseMod.F90
  PatchType.F90
  pftconMod.F90
  SaturatedExcessRunoffMod.F90
  SoilBiogeochemCarbonFluxType.F90
  SoilBiogeochemCarbonStateType.F90
  SoilBiogeochemCompetitionMod.F90
  SoilBiogeochemDecompCascadeConType.F90
  SoilBiogeochemNitrogenFluxType.F90
  SoilBiogeochemNitrogenStateType.F90
  SoilBiogeochemStateType.F90
  SoilStateType.F90
  SolarAbsorbedType.F90
  SurfaceAlbedoType.F90
  TemperatureType.F90
  Wateratm2lndBulkType.F90
  Wateratm2lndType.F90
  WaterDiagnosticType.F90
  WaterDiagnosticBulkType.F90
  WaterFluxBulkType.F90
  WaterFluxType.F90
  WaterStateBulkType.F90
  WaterStateType.F90
  WaterType.F90
  CNCStateUpdate1Mod.F90
  CNCStateUpdate2Mod.F90
  CNCStateUpdate3Mod.F90
  CNDriverMod.F90
  CNFireEmissionsMod.F90
  CNFireFactoryMod.F90
  CNFireLi2014Mod.F90
  CNFireLi2016Mod.F90
  CNFireLi2021Mod.F90
  CNFireNoFireMod.F90
  CNFUNMod.F90
  CNGapMortalityMod.F90
  CNGRespMod.F90
  CNCLM_init_mod.F90
  CNMRespMod.F90
  CNNDynamicsMod.F90
  CNNStateUpdate1Mod.F90
  CNNStateUpdate2Mod.F90
  CNNStateUpdate3Mod.F90
  CNPhenologyMod.F90
  CNPrecisionControlMod.F90
  CNRootDynMod.F90
  CNSharedParamsMod.F90
  CNVegetationFacade.F90
  CNVegStructUpdateMod.F90
  column_varcon.F90
  fileutils.F90
  filterColMod.F90
  FireMethodType.F90
  initSubgridMod.F90
  landunit_varcon.F90
  ncdio_pio.F90
  NutrientCompetitionCLM45defaultMod.F90
  NutrientCompetitionFactoryMod.F90
  NutrientCompetitionFlexibleCNMod.F90
  NutrientCompetitionMethodMod.F90
  paramUtilMod.F90
  perf_mod.F90
  PhotosynthesisMod.F90
  QSatMod.F90
  quadraticMod.F90
  RootBiophysMod.F90
  shr_abort_mod.F90
  shr_assert.h
  shr_const_mod.F90
  shr_file_mod.F90
  shr_fire_emis_mod.F90
  shr_kind_mod.F90
  shr_log_mod.F90
  shr_mpi_mod.F90
  shr_nl_mod.F90
  shr_sys_mod.F90
  SoilBiogeochemDecompCascadeBGCMod.F90
  SoilBiogeochemDecompCascadeCNMod.F90
  SoilBiogeochemDecompMod.F90
  SoilBiogeochemLittVertTranspMod.F90
  SoilBiogeochemNitrifDenitrifMod.F90
  SoilBiogeochemNLeachingMod.F90
  SoilBiogeochemNStateUpdate1Mod.F90
  SoilBiogeochemPotentialMod.F90
  SoilBiogeochemPrecisionControlMod.F90
  SoilBiogeochemVerticalProfileMod.F90
  SoilStateInitTimeConstMod.F90
  SoilWaterRetentionCurveMod.F90
  spmdMod.F90
  subgridAveMod.F90
  SurfaceAlbedoMod.F90
  SurfaceRadiationMod.F90
  TridiagonalMod.F90
  )

list(APPEND srcs
    ${share_genf90_sources})

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES MAPL GEOS_LandShared GEOS_CatchCNShared esmf NetCDF::NetCDF_Fortran
  TYPE SHARED)
target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})

if (is_openmp)
  target_compile_options(${this} PRIVATE ${OpenMP_Fortran_FLAGS})
endif ()

