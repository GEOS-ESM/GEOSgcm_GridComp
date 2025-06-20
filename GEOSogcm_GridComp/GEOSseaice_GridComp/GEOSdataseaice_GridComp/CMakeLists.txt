esma_set_this ()

# This is for selecting the sea ice data handle.
## We default to VERSION_1 and look if something is passed in.
### VERSION_1 is what was used in MERRAA, MERRA-2 reanalyses,
### VERSION_2 is a "test" with ice thickness from ExtData.

set (DEFAULT_DATA_SEAICE "VERSION_1")
if (NOT DATA_SEAICE)
  set (DATA_SEAICE ${DEFAULT_DATA_SEAICE})
endif ()

## ... and make sure it is allowed
set (ALLOWED_DATA_SEAICE "VERSION_1" "VERSION_2")
if (NOT DATA_SEAICE IN_LIST ALLOWED_DATA_SEAICE)
  message (FATAL_ERROR "The only allowed values for DATA_SEAICE are ${ALLOWED_DATA_SEAICE}")
endif ()

message(STATUS "DATA_SEAICE set to: ${DATA_SEAICE}")

if (DATA_SEAICE STREQUAL "VERSION_1")
  esma_add_library (${this} SRCS GEOS_DataSeaIceGridComp.F90         DEPENDENCIES MAPL CICE4 ESMF::ESMF NetCDF::NetCDF_Fortran)
elseif (DATA_SEAICE STREQUAL "VERSION_2")
  esma_add_library (${this} SRCS GEOS_DataSeaIceGridComp_ExtData.F90 DEPENDENCIES MAPL ESMF::ESMF NetCDF::NetCDF_Fortran)
else ()
  message (FATAL_ERROR "The only allowed values for DATA_SEAICE are ${ALLOWED_DATA_SEAICE}")
endif ()

