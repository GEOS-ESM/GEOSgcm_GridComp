esma_set_this()

set (alldirs
  GEOSumwm_GridComp
  GEOSwavewatch_GridComp
  )

set (resource_files
  WGCM.rc
  )

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/GEOS_WgcmGridComp.F90)

  set_property(
    SOURCE bl_seaspray.F90
    PROPERTY COMPILE_DEFINITIONS GEOS
    )

  esma_add_library (${this}
    SRCS bl_seaspray.F90 GEOS_WgcmGridComp.F90
    SUBCOMPONENTS ${alldirs}
    DEPENDENCIES MAPL esmf
    )

else ()

  esma_add_subdirectories (${alldirs})

endif()


install (FILES ${resource_files}
  DESTINATION etc
  )
