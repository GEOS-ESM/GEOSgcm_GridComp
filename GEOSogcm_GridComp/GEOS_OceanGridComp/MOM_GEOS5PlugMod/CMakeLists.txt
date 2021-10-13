esma_set_this ()

set (SRCS
  MOM_GEOS5PlugMod.F90
  )

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/seawater.F90)
  list (APPEND SRCS seawater.F90)
endif ()

esma_add_library (${this}
  SRCS ${SRCS}
  SUBCOMPONENTS @mom
  DEPENDENCIES MAPL
  TYPE SHARED
  )

# Specs for r8 version
string(REPLACE " " ";" tmp ${FREAL8})
foreach (flag ${tmp})
  target_compile_options (${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${flag}>)
endforeach ()
