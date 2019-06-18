esma_set_this ()

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/GEOS_PhysicsGridComp.F90)

  set (alldirs
    @GEOSchem_GridComp
    GEOSmoist_GridComp
    GEOSsurface_GridComp
    GEOSturbulence_GridComp
    GEOSgwd_GridComp
    GEOSradiation_GridComp
    )

  esma_add_library (${this}
    SRCS GEOS_PhysicsGridComp.F90
    SUBCOMPONENTS ${alldirs}
    DEPENDENCIES MAPL_Base GMAO_mpeu GMAO_stoch
    INCLUDES ${INC_ESMF})

else ()

  esma_add_subdirectories (${alldirs})

endif()

