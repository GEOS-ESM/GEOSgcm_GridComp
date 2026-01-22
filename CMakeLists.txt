esma_set_this()

set (alldirs
  GEOSagcm_GridComp
  GEOSdataatm_GridComp
  GEOSmkiau_GridComp
  GEOSogcm_GridComp
  GEOSwgcm_GridComp
  )

option(BUILD_WITH_GIGATRAJ "Build GEOSgcm with Gigatraj" OFF)
if (BUILD_WITH_GIGATRAJ)
  list(APPEND alldirs GEOSgigatraj_GridComp)
endif()

option(BUILD_WITH_PYMLINC "Build pyMLINC interface" OFF)

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/GEOS_GcmGridComp.F90)

 ecbuild_declare_project()
  esma_add_library(${this}
    SRCS GEOS_GcmGridComp.F90
    SUBCOMPONENTS ${alldirs}
    DEPENDENCIES MAPL ESMF::ESMF)

  target_compile_definitions (${this} PRIVATE $<$<BOOL:${BUILD_WITH_GIGATRAJ}>:HAS_GIGATRAJ>)
  target_compile_definitions (${this} PRIVATE $<$<BOOL:${BUILD_WITH_PYMLINC}>:HAS_PYMLINC>)

  ecbuild_install_project( NAME GEOSgcm_GridComp)

else ()

  esma_add_subdirectories (${alldirs})

endif()




