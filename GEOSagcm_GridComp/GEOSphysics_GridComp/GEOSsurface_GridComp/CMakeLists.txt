esma_set_this ()

set (alldirs
   GEOSlandice_GridComp
   GEOSlake_GridComp
   GEOSland_GridComp
   GEOSsaltwater_GridComp
   )

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/GEOS_SurfaceGridComp.F90)

  option(AQUAPLANET "Build for Aquaplanet simulation" OFF)
  if (AQUAPLANET)
    add_compile_definitions(AQUA_PLANET)
  endif()

  esma_add_library (${this}
    SRCS GEOS_SurfaceGridComp.F90
    SUBCOMPONENTS ${alldirs}
    SUBDIRS Shared
    DEPENDENCIES GEOS_SurfaceShared MAPL GMAO_mpeu esmf)

else ()

   esma_add_subdirectories (${alldirs} Shared)

endif()

esma_add_subdirectories(Utils)


