esma_set_this()

set (srcs
  GEOS_SaltWaterGridComp.F90
  GEOS_OpenWaterGridComp.F90
  GEOS_SimpleSeaiceGridComp.F90
  GEOS_ObioGridComp.F90)

if (GEOS_SEAICE STREQUAL CICE4)
  list( APPEND srcs GEOS_CICE4ColumnPhysGridComp.F90)
  set (CICE_LIBRARY CICE4)
elseif (GEOS_SEAICE STREQUAL CICE6)
  list( APPEND srcs GEOS_SeaiceInterfaceGridComp.F90)
elseif (GEOS_SEAICE STREQUAL dataseaice)
  # no real sea ice model, just datasea
else()
  message (FATAL_ERROR "Invalid sea ice model selection: ${GEOS_SEAICE}")
endif()

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES GEOS_Shared MAPL ${CICE_LIBRARY} esmf NetCDF::NetCDF_Fortran)
