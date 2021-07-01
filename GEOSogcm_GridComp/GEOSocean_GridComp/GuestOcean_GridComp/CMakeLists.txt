esma_set_this ()

option (BUILD_MIT_OCEAN "BUILD_MIT_OCEAN is turned off by default" OFF) #OFF by default
if ( NOT BUILD_MIT_OCEAN )
  add_subdirectory(MOM6_GEOSPlug)
  add_subdirectory(MOM_GEOS5PlugMod)
else()
  add_subdirectory(MIT_GEOS5PlugMod)
endif()

esma_add_library (${this}
  SRCS GuestOcean_GridComp.F90
  DEPENDENCIES GEOSdatasea_GridComp MAPL
  INCLUDES ${INC_ESMF})

