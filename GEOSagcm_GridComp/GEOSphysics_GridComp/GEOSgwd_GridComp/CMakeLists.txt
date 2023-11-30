esma_set_this ()

set (subdirs
  ncar_gwd
  )

set (srcs
  gw_drag.F90
  GEOS_GwdGridComp.F90
  machine.F90
  gwdps.f
  gwdc.f
  )

set (resource_files
   GWD_GridComp.rc
   )

install( FILES ${resource_files}
   DESTINATION etc
   )

esma_add_library (${this}
  SRCS ${srcs}
  SUBCOMPONENTS ${subdirs}
  DEPENDENCIES GEOS_Shared MAPL esmf NetCDF::NetCDF_Fortran
  )

# CMake has an OpenMP issue with NAG Fortran: https://gitlab.kitware.com/cmake/cmake/-/issues/21280
if (NOT CMAKE_Fortran_COMPILER_ID MATCHES "NAG")
   target_link_libraries(${this} PRIVATE OpenMP::OpenMP_Fortran)
endif ()
