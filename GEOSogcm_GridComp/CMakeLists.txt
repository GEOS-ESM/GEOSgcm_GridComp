esma_set_this ()

set (GEOS_OCEAN MOM6 CACHE STRING "Options for ocean model: MOM5, MOM6, MIT")
message (STATUS "Selected ocean model: " ${GEOS_OCEAN})

if (GEOS_OCEAN STREQUAL MIT)
  if (GEOS_SEAICE NOT STREQUAL CICE4)
    message(FATAL_ERROR "CANNOT USE ${GEOS_SEAICE} SEA ICE MODEL with MIT OCEAN MODEL.")
  endif()
else()
  if ( BUILD_MIT_OCEAN )
    message(FATAL_ERROR "INCONSISTENT CHOICE OF OCEAN MODEL.")
  endif()
endif()

set (alldirs
  GEOS_OceanBioGeoChemGridComp
  GEOS_OradGridComp
  GEOS_OradBioGridComp
  GEOS_OceanGridComp
  )

set (dependencies GEOS_Seaice_GridComp MAPL esmf)

if ( BUILD_MIT_OCEAN )
  add_compile_definitions(BUILD_MIT_OCEAN)
endif()

esma_add_library (${this}
  SRCS GEOS_OgcmGridComp.F90
  SUBCOMPONENTS ${alldirs}
  SUBDIRS GEOSseaice_GridComp
  DEPENDENCIES ${dependencies})
