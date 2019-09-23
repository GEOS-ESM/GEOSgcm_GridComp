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
    DEPENDENCIES MAPL_Base
    INCLUDES ${INC_ESMF} ${CMAKE_BINARY_DIR}/include/modelE)

else ()

  esma_add_subdirectories (${alldirs})

endif()

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/include/modelE)
include(ExternalProject)

ExternalProject_Add(modelE
	PREFIX "modelE"
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/@Model5Eagcm_GridComp/decks
	DOWNLOAD_COMMAND ""
	CONFIGURE_COMMAND ""
	UPDATE_COMMAND ""
	BUILD_COMMAND make setup RUN=$ENV{ModelE_Run} MAPL=YES OVERWRITE=YES GEOS_BINARY_DIR=${CMAKE_BINARY_DIR}
	BUILD_IN_SOURCE 1
	INSTALL_COMMAND find ../model/ -type f -name "*.a" -exec /bin/cp {} ${CMAKE_BINARY_DIR}/lib $<SEMICOLON> &&  find ../model/ -type f -name "*.mod" -exec /bin/cp {} ${CMAKE_BINARY_DIR}/include/modelE $<SEMICOLON>
	DEPENDS MAPL_Base
)

  


