esma_set_this ()

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/GEOS_surfaceGridComp.F90)

  set (alldirs
    GEOSlandice_GridComp
    GEOSlake_GridComp
    GEOSland_GridComp
    GEOSsaltwater_GridComp
    )
  
  esma_add_library (${this}
    SRCS GEOS_SurfaceGridComp.F90
    SUBCOMPONENTS ${alldirs}
    SUBDIRS Shared
    DEPENDENCIES GEOS_SurfaceShared MAPL_Base GMAO_mpeu
    INCLUDES ${INC_ESMF})

else ()

  esma_add_subdirectories (${alldirs})

endif()



