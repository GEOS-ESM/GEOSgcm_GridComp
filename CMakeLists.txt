esma_set_this()

set (alldirs
  GEOSagcm_GridComp
  GEOSdataatm_GridComp
  GEOSmkiau_GridComp
  GEOSogcm_GridComp
  )


if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/GEOS_GcmGridComp.F90)

  esma_add_library(${this}
    SRCS GEOS_GcmGridComp.F90
    SUBCOMPONENTS ${alldirs}
    DEPENDENCIES MAPL
    INCLUDES ${INC_ESMF})

else ()

  esma_add_subdirectories (${alldirs})

endif()

  


