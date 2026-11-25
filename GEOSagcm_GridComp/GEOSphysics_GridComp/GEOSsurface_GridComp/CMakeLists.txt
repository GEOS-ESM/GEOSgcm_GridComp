esma_set_this ()

set (alldirs
   GEOSlandice_GridComp
   GEOSlake_GridComp
   GEOSland_GridComp
   GEOSroute_GridComp
   GEOSsaltwater_GridComp
   )

# LIS_GridComp integrates LISF's `LIS` library from source
# opt-in via -DBUILD_LIS=ON
if (BUILD_LIS)
  list (APPEND alldirs LIS_GridComp)
endif ()

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/GEOS_SurfaceGridComp.F90)

  esma_add_library (${this}
    SRCS GEOS_SurfaceGridComp.F90
    SUBCOMPONENTS ${alldirs}
    SUBDIRS Shared
    DEPENDENCIES GEOS_SurfaceShared MAPL GMAO_mpeu ESMF::ESMF)

else ()

   esma_add_subdirectories (${alldirs} Shared)

endif()

esma_add_subdirectories(Utils)


