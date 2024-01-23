esma_set_this ()

set (alldirs
  GEOSchem_GridComp
  GEOSmoist_GridComp
  GEOSsurface_GridComp
  GEOSturbulence_GridComp
  GEOSgwd_GridComp
  )


if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/GEOS_PhysicsGridComp.F90)

  esma_add_subdirectory(GEOSradiation_GridComp)

  esma_add_library (${this}
    SRCS GEOS_PhysicsGridComp.F90 MBundle_IncrementMod.F90
    SUBCOMPONENTS ${alldirs}
    DEPENDENCIES MAPL GMAO_mpeu GMAO_stoch esmf GEOS_RadiationGridComp)

else ()

  esma_add_subdirectories (${alldirs} GEOSradiation_GridComp)

endif()

