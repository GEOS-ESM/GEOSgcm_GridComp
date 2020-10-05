esma_set_this ()

set (alldirs
  umwm
  )


if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/GEOS_UMWMGridComp.F90)

  esma_add_library (${this}
    SRCS GEOS_UMWMGridComp.F90
    SUBCOMPONENTS ${alldirs}
    DEPENDENCIES MAPL
    INCLUDES ${INC_ESMF})

else ()

  esma_add_subdirectories (${alldirs})

endif()


file (GLOB_RECURSE rc_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.rc)
foreach ( file ${rc_files} )
   get_filename_component( dir ${file} DIRECTORY )
   install( FILES ${file} DESTINATION etc/${dir} )
endforeach()
