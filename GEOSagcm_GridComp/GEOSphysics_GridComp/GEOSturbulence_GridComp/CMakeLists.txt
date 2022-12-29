esma_set_this()

set (srcs
  GEOS_TurbulenceGridComp.F90
  LockEntrain.F90
  shoc.F90
  edmf.F90
  scm_surface.F90
  shocparams.F90
  edmfparams.F90
  )

set (resource_files
  GEOS_TurbulenceGridComp.rc
  )

install( FILES ${resource_files}
  DESTINATION etc
  )

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES GEOS_Shared MAPL esmf NetCDF::NetCDF_Fortran)

if (NOT CMAKE_Fortran_COMPILER_ID MATCHES "NAG")
  target_link_libraries(${this} PRIVATE OpenMP::OpenMP_Fortran)
endif ()
