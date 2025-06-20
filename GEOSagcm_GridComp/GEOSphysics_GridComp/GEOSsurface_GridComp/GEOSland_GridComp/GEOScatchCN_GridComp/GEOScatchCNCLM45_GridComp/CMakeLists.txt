esma_set_this ()

string (REPLACE GEOScatchCNCLM45_GridComp_ "" is_openmp ${this})

esma_add_subdirectories (CLM45)

set (srcs
  GEOS_CatchCNCLM45GridComp.F90
  )

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES MAPL GEOS_Shared GEOS_LandShared CLM45 GEOS_CatchCNShared ESMF::ESMF NetCDF::NetCDF_Fortran
  TYPE SHARED)

if (is_openmp)
  target_link_libraries(${this} PRIVATE OpenMP::OpenMP_Fortran)
endif ()

