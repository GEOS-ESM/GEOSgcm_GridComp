esma_set_this ()
string (REPLACE GEOScatchCN_GridComp_ "" is_openmp ${this})

set (srcs
  CNCLM_Photosynthesis.F90
  PhotosynthesisMod.F90
  )

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES MAPL GEOS_LandShared GEOS_CatchCNShared 
  TYPE SHARED)
target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})

if (is_openmp)
  target_compile_options(${this} PRIVATE ${OpenMP_Fortran_FLAGS})
endif ()

