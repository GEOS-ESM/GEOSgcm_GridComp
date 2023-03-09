esma_set_this ()

set (SRCS
  CICE_GEOSPlug.F90
  )

esma_add_subdirectories(
  cice6_cmake
  )

esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES MAPL cice6
  TYPE SHARED
  )


install(
   DIRECTORY cice6_app
   DESTINATION etc/CICE6
   USE_SOURCE_PERMISSIONS
   )
