esma_set_this ()

set (SRCS
  MOM6_GEOSPlug.F90
  )

esma_add_subdirectories(
  mom6_cmake
  )

esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES MAPL fms_r8 mom6
  TYPE SHARED
  )

# Specs for r8 version
string(REPLACE " " ";" tmp ${FREAL8})
foreach (flag ${tmp})
  target_compile_options (${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${flag}>)
endforeach ()

install(
   DIRECTORY mom6_app
   DESTINATION etc/MOM6
   USE_SOURCE_PERMISSIONS
   )
