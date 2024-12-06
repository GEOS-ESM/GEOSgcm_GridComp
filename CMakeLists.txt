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
  add_definitions (-DHAS_GIGATRAJ)
  set (alldirs ${alldirs} GEOSgigatraj_GridComp) 
endif()


if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/GEOS_GcmGridComp.F90)

 ecbuild_declare_project()
  esma_add_library(${this}
    SRCS GEOS_GcmGridComp.F90
    SUBCOMPONENTS ${alldirs}
    DEPENDENCIES MAPL ESMF::ESMF)

  ecbuild_install_project( NAME GEOSgcm_GridComp)

else ()

  esma_add_subdirectories (${alldirs})

endif()




