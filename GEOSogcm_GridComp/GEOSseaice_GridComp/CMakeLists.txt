esma_set_this (OVERRIDE GEOS_Seaice_GridComp)

option (BUILD_MIT_OCEAN "BUILD_MIT_OCEAN is turned off by default" OFF) #OFF by default

if ( NOT BUILD_MIT_OCEAN )
 esma_add_library (${this}
   SRCS GEOS_SeaIceGridComp.F90
   SUBCOMPONENTS GEOSdataseaice_GridComp GEOSCICEDyna_GridComp CICE_GEOSPlug 
   DEPENDENCIES MAPL CICE4 esmf)
else()

  add_compile_definitions(BUILD_MIT_OCEAN)

  esma_add_library (${this}
    SRCS GEOS_SeaIceGridComp.F90
    SUBCOMPONENTS GEOSdataseaice_GridComp GEOSMITDyna_GridComp
    DEPENDENCIES MAPL
    INCLUDES ${INC_ESMF})
endif()
