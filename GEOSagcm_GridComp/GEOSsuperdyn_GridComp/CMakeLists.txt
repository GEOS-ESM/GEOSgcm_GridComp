esma_set_this ()

set (alldirs
  FVdycoreCubed_GridComp
  FVdycore_GridComp
  GEOSdatmodyn_GridComp
  ARIESg3_GridComp
  )

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/GEOS_SuperdynGridComp.F90)

  esma_add_library (${this}
    SRCS GEOS_SuperdynGridComp.F90
    SUBCOMPONENTS ${alldirs}
    DEPENDENCIES MAPL GEOS_Shared
    INCLUDES ${INC_ESMF})

else ()

  esma_add_subdirectories (${alldirs})

endif ()
  



