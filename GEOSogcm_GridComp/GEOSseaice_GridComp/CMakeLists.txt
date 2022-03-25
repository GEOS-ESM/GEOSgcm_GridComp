esma_set_this (OVERRIDE GEOS_Seaice_GridComp)

esma_add_library (${this}
  SRCS GEOS_SeaIceGridComp.F90
  SUBCOMPONENTS GEOSdataseaice_GridComp GEOSCICEDyna_GridComp CICE_GEOSPlug 
  DEPENDENCIES MAPL CICE4
  INCLUDES ${INC_ESMF})
