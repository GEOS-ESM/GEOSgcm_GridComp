esma_set_this ()

add_subdirectory(MOM6_GEOSPlug)
add_subdirectory(MOM_GEOS5PlugMod)
add_subdirectory(GEOSdatasea_GridComp)

esma_add_library (${this}
  SRCS GEOS_OceanGridComp.F90
  DEPENDENCIES GEOSdatasea_GridComp MAPL
  INCLUDES ${INC_ESMF})
