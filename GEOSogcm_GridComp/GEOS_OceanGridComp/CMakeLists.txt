esma_set_this ()

#target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})

option (BUILD_MIT_OCEAN "BUILD_MIT_OCEAN is turned off by default" OFF) #OFF by default
esma_add_library (${this}
  SRCS GEOS_OceanGridComp.F90
  if ( NOT BUILD_MIT_OCEAN )
   SUBCOMPONENTS MOM6_GEOSPlug MOM_GEOS5PlugMod GEOSdatasea_GridComp
  else()
   SUBCOMPONENTS MIT_GEOS5PlugMod
  endif()
  DEPENDENCIES GEOSdatasea_GridComp MAPL
  INCLUDES ${INC_ESMF})
