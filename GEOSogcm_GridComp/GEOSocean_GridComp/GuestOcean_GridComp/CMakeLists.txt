esma_set_this ()

add_subdirectory(MOM6_GEOSPlug)
add_subdirectory(MOM_GEOS5PlugMod)

esma_add_library (${this}
  SRCS GuestOcean_GridComp.F90
  DEPENDENCIES GEOSdatasea_GridComp MAPL
  INCLUDES ${INC_ESMF})
