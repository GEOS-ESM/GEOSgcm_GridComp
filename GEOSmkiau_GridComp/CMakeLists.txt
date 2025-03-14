esma_set_this()

set(srcs
  IAU_GridCompMod.F90
  GEOS_mkiauGridComp.F90
  DFI_GridComp.F90
  mkiau_specmod.F90
  DynVec_GridComp.F90)

if (BUILD_WITH_PYMLINC)
  add_compile_definitions(-DHAS_PYMLINC)
  set(pymlinc_srcs
    pyMLINC/interface/interface.f90
    pyMLINC/interface/interface.c)
  list(APPEND srcs ${pymlinc_srcs})
  include(pyMLINC.cmake)
  set(dependencies pyMLINC_interface_py)
  ecbuild_add_executable(
    TARGET pymlinc-driver
    SOURCES pyMLINC/driver/driver.F90 ${pymlinc_srcs}
    LIBS ${PYMLINC_INTERFACE_LIBRARY})
endif ()

set(dependencies
  ${dependencies}
  MAPL_cfio_r4
  NCEP_sp_r4i4
  GEOS_Shared
  GMAO_mpeu
  MAPL
  FVdycoreCubed_GridComp
  ESMF::ESMF
  NetCDF::NetCDF_Fortran)  

esma_add_library(${this} SRCS ${srcs} DEPENDENCIES ${dependencies})

# install(
#   DIRECTORY pyMLINC/pyMLINC pyMLINC/geos_state_bias
#   DESTINATION lib/Python
#   USE_SOURCE_PERMISSIONS)
