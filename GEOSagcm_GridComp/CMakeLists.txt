esma_set_this()

set (alldirs
  GEOSsuperdyn_GridComp
  GEOSphysics_GridComp
  )

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/GEOS_AgcmGridComp.F90)

  esma_add_library (${this}
    SRCS GEOS_AgcmGridComp.F90
    SUBCOMPONENTS ${alldirs}
    DEPENDENCIES MAPL GEOS_Shared Chem_Shared
    INCLUDES ${INC_ESMF})

else ()

  esma_add_subdirectories (${alldirs})

endif()
