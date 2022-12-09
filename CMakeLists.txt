esma_set_this()

set (alldirs
  GEOSmkiau_GridComp
  )


if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/GEOS_GcmGridComp.F90)

 ecbuild_declare_project()
  esma_add_library(${this}
    #SRCS GEOS_GcmGridComp.F90
    SRCS gcm_stub.F90
    SUBCOMPONENTS ${alldirs}
    DEPENDENCIES MAPL esmf)

  ecbuild_install_project( NAME GEOSgcm_GridComp)

else ()

  esma_add_subdirectories (${alldirs})

endif()

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/include/modelE)
include(ExternalProject)
message(INFO "MODELE TEMPLATE: $ENV{ModelE_Template}")
message(INFO "MODELE RUN: $ENV{ModelE_Run}")
message(INFO "MODELERC:  $ENV{MODELERC}")
ExternalProject_Add(modelE
   PREFIX "modelE"
   SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/@Model5Eagcm_GridComp/decks
   DOWNLOAD_COMMAND ""
   CONFIGURE_COMMAND ""
   UPDATE_COMMAND ""
   BUILD_ALWAYS 1
   BUILD_COMMAND make clean OVERWRITE=YES && make rundeck RUN=$ENV{ModelE_Run} RUNSRC=$ENV{ModelE_Template} OVERWRITE=YES && $(MAKE) gcm RUN=$ENV{ModelE_Run} MAPL=YES OVERWRITE=YES GEOS_BINARY_DIR=${CMAKE_BINARY_DIR}  VERBOSE_OUTPUT=YES MODELERC=$ENV{MODELERC} MPI=YES
   BUILD_IN_SOURCE 1
   INSTALL_COMMAND find ../model/ -type f -name "*.a" -exec /bin/cp {} ${CMAKE_BINARY_DIR}/lib $<SEMICOLON> &&  find ../model/ -type f -name "*.mod" -exec /bin/cp {} ${CMAKE_BINARY_DIR}/include/modelE $<SEMICOLON>
   DEPENDS MAPL
)

