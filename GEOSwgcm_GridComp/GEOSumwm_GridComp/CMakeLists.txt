esma_set_this ()


esma_add_subdirectories(
  umwm_cmake
  )

esma_add_library (${this}
  SRCS GEOS_UMWMGridComp.F90
  DEPENDENCIES MAPL umwm
  INCLUDES ${INC_ESMF}
  )


file (GLOB_RECURSE rc_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.rc)
foreach ( file ${rc_files} )
   get_filename_component( dir ${file} DIRECTORY )
   install( FILES ${file} DESTINATION etc/${dir} )
endforeach()


install (
    FILES UMWM_ExtData.yaml
    DESTINATION etc
)
