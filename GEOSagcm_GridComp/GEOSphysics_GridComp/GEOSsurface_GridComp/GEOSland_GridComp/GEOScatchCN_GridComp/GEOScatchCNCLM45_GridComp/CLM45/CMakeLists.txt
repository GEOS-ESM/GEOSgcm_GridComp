esma_set_this ()
string (REPLACE GEOScatchCN_GridComp_ "" is_openmp ${this})

set (srcs
  update_model_para4cn.F90
  pso_params.F90
  clmtype.F90
  clm_time_manager.F90
  shr_const_mod.F90
  clm_varpar.F90
  clm_varcon.F90
  clm_varctl.F90
  subgridAveMod.F90
  CNAllocationMod.F90
  CNAnnualUpdateMod.F90
  CNCStateUpdate1Mod.F90
  CNCStateUpdate2Mod.F90
  CNCStateUpdate3Mod.F90
  CNNStateUpdate1Mod.F90
  CNNStateUpdate2Mod.F90
  CNNStateUpdate3Mod.F90
  CNDecompCascadeMod_BGC.F90
  CNDecompMod.F90
  CNFireMod.F90
  CNGapMortalityMod.F90
  pftvarcon.F90
  CNMRespMod.F90
  CNGRespMod.F90
  CNNDynamicsMod.F90
  CNPhenologyMod.F90
  CNPrecisionControlMod.F90
  CNSetValueMod.F90
  CNVegStructUpdateMod.F90
  CNVerticalProfileMod.F90    
  CNSoilLittVertTranspMod.F90 
  CNWoodProductsMod.F90
  CNSummaryMod.F90
  CNEcosystemDynMod.F90
  clmtypeInitMod.F90
  CNiniTimeVar.F90
  CNBalanceCheckMod.F90
  CN_DriverMod.F90
  compute_rc.F90
  TridiagonalMod.F90
  getco2.F90
  pft_clm_to_pso.F90
  read_env.F90
  )

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES MAPL GEOS_LandShared GEOS_CatchCNShared
  TYPE SHARED)
#target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})

if (is_openmp)
  target_link_libraries(${this} PRIVATE OpenMP::OpenMP_Fortran)
endif ()

#if (is_openmp)
#  target_compile_options(${this} PRIVATE ${OpenMP_Fortran_FLAGS})
#endif ()

