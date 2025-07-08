esma_set_this ()

string (REPLACE GEOScatchCNCLM51_GridComp_ "" is_openmp ${this})

esma_add_subdirectories (CLM51)

set (srcs
  GEOS_CatchCNCLM51GridComp.F90
  clm_varcon.F90
  clm_varctl.F90
  clm_varpar.F90
  CN2CLMType.F90
  CNCLM_DriverMod.F90
  CNCLM_init_mod.F90
  CNCLM_Photosynthesis.F90
  CNFireLi2014Mod.F90
  CNFireLi2016Mod.F90
  CNFireLi2021Mod.F90
  CNPhenologyMod.F90
  SoilBiogeochemNitrogenStateType.F90
  SoilStateType.F90
  )

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES MAPL GEOS_Shared GEOS_LandShared CLM51 GEOS_CatchCNShared ESMF::ESMF NetCDF::NetCDF_Fortran
  TYPE SHARED)

if (is_openmp)
  target_link_libraries(${this} PRIVATE OpenMP::OpenMP_Fortran)
endif ()

file(GLOB nml_files CN_CLM51.nml      )

install(
   FILES ${nml_files}
   DESTINATION etc
   )

