esma_set_this()

set (alldirs
  GEOSsuperdyn_GridComp
  GEOSphysics_GridComp
  )

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/GEOS_AgcmSimpleGridComp.F90)

  message(STATUS "Building Held-Suarez GridComp")
  esma_add_library (${this}
    SRCS GEOS_AgcmSimpleGridComp.F90
    SUBCOMPONENTS GEOSsuperdyn_GridComp GEOShs_GridComp
    DEPENDENCIES MAPL ESMF::ESMF)

elseif (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/GEOS_AgcmGridComp.F90)

  esma_add_library (${this}
    SRCS GEOS_AgcmGridComp.F90
    SUBCOMPONENTS ${alldirs}
    DEPENDENCIES MAPL GEOS_Shared Chem_Shared ESMF::ESMF)

  target_compile_definitions (${this} PRIVATE $<$<BOOL:${BUILD_WITH_GIGATRAJ}>:HAS_GIGATRAJ>)

else ()

  esma_add_subdirectories (${alldirs})

endif()
