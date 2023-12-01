esma_set_this (OVERRIDE GEOS_Seaice_GridComp)

option (BUILD_MIT_OCEAN "BUILD_MIT_OCEAN is turned off by default" OFF) #OFF by default

set (subcomponents_seaice GEOSdataseaice_GridComp)

if (GEOS_SEAICE STREQUAL CICE4)
  list( APPEND subcomponents_seaice GEOSCICEDyna_GridComp)
  set (CICE_LIBRARY CICE4)
elseif (GEOS_SEAICE STREQUAL CICE6)
  list( APPEND subcomponents_seaice CICE_GEOSPlug)
 elseif (GEOS_SEAICE STREQUAL dataseaice)
   # no real sea ice model, just dataseaice
 else()
   message (FATAL_ERROR "Invalid sea ice model selection: ${GEOS_SEAICE}")
 endif()

if ( NOT BUILD_MIT_OCEAN )
 esma_add_library (${this}
   SRCS GEOS_SeaIceGridComp.F90
   SUBCOMPONENTS ${subcomponents_seaice}
   DEPENDENCIES MAPL ${CICE_LIBRARY} esmf)
else()

  add_compile_definitions(BUILD_MIT_OCEAN)

  esma_add_library (${this}
    SRCS GEOS_SeaIceGridComp.F90
    SUBCOMPONENTS GEOSdataseaice_GridComp GEOSMITDyna_GridComp
    DEPENDENCIES MAPL
    INCLUDES ${INC_ESMF})
endif()
