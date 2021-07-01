esma_set_this (OVERRIDE GEOS_Seaice_GridComp)

esma_add_library (${this}
  SRCS GEOS_SeaIceGridComp.F90
if ( BUILD_MIT_OCEAN )
   add_definitions(-DBUILD_MIT_OCEAN=1)
   SUBCOMPONENTS GEOSdataseaice_GridComp GEOSMITDyna_GridComp
else()
  SUBCOMPONENTS GEOSdataseaice_GridComp GEOSCICEDyna_GridComp
endif()
  DEPENDENCIES MAPL CICE4
  INCLUDES ${INC_ESMF})
