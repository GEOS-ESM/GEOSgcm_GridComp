esma_set_this ()

set (alldirs
   GEOSlandice_GridComp
   GEOSlake_GridComp
   GEOSland_GridComp
   GEOSsaltwater_GridComp
   GEOShetprecip_GridComp
   )

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/GEOS_SurfaceGridComp.F90)

  esma_add_library (${this}
    SRCS GEOS_SurfaceGridComp.F90
    SUBCOMPONENTS ${alldirs}
    SUBDIRS Shared
    DEPENDENCIES GEOS_SurfaceShared MAPL GMAO_mpeu esmf)

else ()

   esma_add_subdirectories (${alldirs} Shared)

endif()

esma_add_subdirectories(Utils)


