esma_set_this ()

set (srcs
  gw_drag.F90
  GEOS_GwdGridComp.F90
  machine.F90
  gwdps.f
  gwdc.f
  ncar_gwd/linear_1d_operators.F90
  ncar_gwd/interpolate_data.F90
  ncar_gwd/vdiff_lu_solver.F90
  ncar_gwd/coords_1d.F90
  ncar_gwd/gw_utils.F90
  ncar_gwd/gw_common.F90
  ncar_gwd/gw_convect.F90
  ncar_gwd/gw_rdg.F90
  ncar_gwd/gw_oro.F90
  ncar_gwd/gw_drag.F90
  )

set (resource_files
   GWD_GridComp.rc
   )



install( FILES ${resource_files}
   DESTINATION etc
   )

esma_add_library (${this} SRCS ${srcs} DEPENDENCIES GEOS_Shared MAPL esmf NetCDF::NetCDF_Fortran)

# CMake has an OpenMP issue with NAG Fortran: https://gitlab.kitware.com/cmake/cmake/-/issues/21280
if (NOT CMAKE_Fortran_COMPILER_ID MATCHES "NAG")
   target_link_libraries(${this} PRIVATE OpenMP::OpenMP_Fortran)
endif ()
