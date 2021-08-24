esma_set_this ()

if ( NOT BUILD_MIT_OCEAN )
  esma_add_subdirectories(
    GEOSdataseaice_GridComp
    GEOSCICEDyna_GridComp
    GEOSMITDyna_GridComp
  )
else()
  esma_add_subdirectories(
    GEOSdataseaice_GridComp
    GEOSMITDyna_GridComp
  )
endif()

