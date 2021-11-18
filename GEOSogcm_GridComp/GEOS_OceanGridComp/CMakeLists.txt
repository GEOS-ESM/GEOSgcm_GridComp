esma_set_this ()

esma_add_library (${this}
  SRCS GEOS_OceanGridComp.F90
  SUBCOMPONENTS MOM6_GEOSPlug MOM_GEOS5PlugMod GEOSdatasea_GridComp
  DEPENDENCIES GEOSdatasea_GridComp MAPL
  INCLUDES ${INC_ESMF})
