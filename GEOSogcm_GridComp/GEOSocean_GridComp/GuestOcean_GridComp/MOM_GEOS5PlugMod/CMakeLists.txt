esma_set_this ()

set (SRCS
  MOM_GEOS5PlugMod.F90
  )

esma_add_library (${this}
  SRCS ${SRCS}
  SUBCOMPONENTS @mom @mom6
  DEPENDENCIES MAPL_Base
  )

# Specs for r8 version
string(REPLACE " " ";" tmp ${FREAL8})
foreach (flag ${tmp})
  target_compile_options (${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${flag}>)
endforeach ()
