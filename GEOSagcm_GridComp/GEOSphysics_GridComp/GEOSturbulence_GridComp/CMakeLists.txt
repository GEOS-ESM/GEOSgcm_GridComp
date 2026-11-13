esma_set_this()

set (srcs
  GEOS_TurbulenceGridComp.F90
  LockEntrain.F90
  shoc.F90
  edmf.F90
  scm_surface.F90
  )

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES GEOS_Shared MAPL ESMF::ESMF NetCDF::NetCDF_Fortran TYPE SHARED)

mapl_acg (${this} Turbulence_StateSpecs.rc
  IMPORT_SPECS EXPORT_SPECS INTERNAL_SPECS GET_POINTERS DECLARE_POINTERS)

# CMake has an OpenMP issue with NAG Fortran: https://gitlab.kitware.com/cmake/cmake/-/issues/21280
if (NOT CMAKE_Fortran_COMPILER_ID MATCHES "NAG")
   target_link_libraries(${this} PRIVATE OpenMP::OpenMP_Fortran)
endif ()

add_subdirectory(regression)
