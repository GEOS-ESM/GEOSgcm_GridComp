esma_set_this ()

add_subdirectory(GEOSdataseaice_GridComp)
add_subdirectory(GEOSCICEDyna_GridComp)

esma_add_library (${this}
  SRCS GEOS_SeaIceGridComp.F90
  DEPENDENCIES MAPL LANL_cice
  INCLUDES ${INC_ESMF})
