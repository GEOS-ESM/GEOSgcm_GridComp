esma_set_this ()

set (GEOS_OCEAN MOM6 CACHE STRING "Options for ocean model: MOM5, MOM6, MIT")
message (STATUS "Selected ocean model: " ${GEOS_OCEAN})

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
