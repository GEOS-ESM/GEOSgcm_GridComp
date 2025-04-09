esma_set_this()

# This set(rename_) is for the automatic stubber
# If directory is not available, this renames the
# module in the stubbed routine to be
# name which is expected by the parent.
set (rename_GEOSwgcm_GridComp GEOS_WgcmGridComp)

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

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/GEOS_GcmGridComp.F90)

 ecbuild_declare_project()
  esma_add_library(${this}
    SRCS GEOS_GcmGridComp.F90
    SUBCOMPONENTS ${alldirs}
    DEPENDENCIES MAPL ESMF::ESMF)

  target_compile_definitions (${this} PRIVATE $<$<BOOL:${BUILD_WITH_GIGATRAJ}>:HAS_GIGATRAJ>)

  ecbuild_install_project( NAME GEOSgcm_GridComp)

else ()

  esma_add_subdirectories (${alldirs})

endif()




