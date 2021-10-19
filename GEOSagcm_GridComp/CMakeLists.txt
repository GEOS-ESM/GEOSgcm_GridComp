esma_set_this()

set (alldirs
  GEOSsuperdyn_GridComp
  GEOShs_GridComp
  )
if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/GEOS_AgcmSimpleGridComp.F90)

  esma_add_library (${this}
    SRCS GEOS_AgcmSimpleGridComp.F90
    SUBCOMPONENTS ${alldirs}
#    DEPENDENCIES MAPL GEOS_Shared Chem_Shared
    DEPENDENCIES MAPL
    INCLUDES ${INC_ESMF})

else ()

  esma_add_subdirectories (${alldirs})

endif()
