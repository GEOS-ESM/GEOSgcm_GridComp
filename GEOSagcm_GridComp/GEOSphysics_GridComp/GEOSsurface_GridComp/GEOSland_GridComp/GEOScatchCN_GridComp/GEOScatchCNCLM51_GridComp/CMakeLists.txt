esma_set_this ()

string (REPLACE GEOScatchCNCLM51_GridComp_ "" is_openmp ${this})

esma_add_subdirectories (CLM51)

set (srcs
  GEOS_CatchCNCLM51GridComp.F90
  update_model_para4cn.F90
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

