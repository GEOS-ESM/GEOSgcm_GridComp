esma_set_this ()

set (alldirs
  GEOSchem_GridComp
  GEOSmoist_GridComp
  GEOSsurface_GridComp
  GEOSturbulence_GridComp
  GEOSgwd_GridComp
  GEOSradiation_GridComp
)

# LIS_GridComp integrates LISF's `LIS` library from source
# opt-in via -DBUILD_LIS=ON
if (BUILD_LIS)
  list (APPEND alldirs LIS_GridComp)
endif ()

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/GEOS_PhysicsGridComp.F90)

  esma_add_library (${this}
    SRCS GEOS_PhysicsGridComp.F90 MBundle_IncrementMod.F90
    SUBCOMPONENTS ${alldirs}
    DEPENDENCIES MAPL GMAO_mpeu GMAO_stoch ESMF::ESMF)

else ()

  esma_add_subdirectories (${alldirs})

endif()

