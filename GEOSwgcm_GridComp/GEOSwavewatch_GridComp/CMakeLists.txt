esma_set_this ()

set (alldirs
    ww3_multi_esmf
  )

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/GEOS_WaveWatchGridComp.F90)

  esma_add_library (${this}
    SRCS GEOS_WaveWatchGridComp.F90
    SUBCOMPONENTS ${alldirs}
    DEPENDENCIES MAPL esmf
    )

else ()

  esma_add_subdirectories (${alldirs})

endif()


file (GLOB_RECURSE rc_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.rc)
foreach ( file ${rc_files} )
   get_filename_component( dir ${file} DIRECTORY )
   install( FILES ${file} DESTINATION etc/${dir} )
endforeach()
