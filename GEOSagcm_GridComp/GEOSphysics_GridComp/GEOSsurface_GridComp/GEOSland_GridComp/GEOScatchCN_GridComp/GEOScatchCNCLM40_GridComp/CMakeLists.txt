esma_set_this ()

string (REPLACE GEOScatchCNCLM40_GridComp_ "" is_openmp ${this})

esma_add_subdirectories (CLM40)

set (srcs
  GEOS_CatchCNCLM40GridComp.F90
  )

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES MAPL GEOS_Shared GEOS_LandShared CLM40 GEOS_CatchCNShared
  TYPE SHARED)
target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})

if (is_openmp)
  target_compile_options(${this} PRIVATE ${OpenMP_Fortran_FLAGS})
endif ()

