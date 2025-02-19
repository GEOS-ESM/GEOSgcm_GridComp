esma_set_this()

option(BUILD_PYMLINC_INTERFACE "Build pyMLINC interface" OFF)

set (srcs
  IAU_GridCompMod.F90
  GEOS_mkiauGridComp.F90
  DFI_GridComp.F90
  mkiau_specmod.F90
  DynVec_GridComp.F90
  )

set(dependencies
  MAPL_cfio_r4
  NCEP_sp_r4i4
  GEOS_Shared
  GMAO_mpeu MAPL
  FVdycoreCubed_GridComp
  ESMF::ESMF
  NetCDF::NetCDF_Fortran)  

if (BUILD_PYMLINC_INTERFACE)
  include (pyMLINC.cmake)
  set(dependencies pyMLINC_interface_py ${dependencies})
endif ()

esma_add_library (${this} SRCS ${srcs} DEPENDENCIES ${dependencies})
