esma_set_this ()

set (rename_MOM_GEOS5PlugMod MOM_GEOS5Plug)
set (rename_MIT_GEOS5PlugMod MIT_GEOS5Plug)

option (BUILD_MIT_OCEAN "BUILD_MIT_OCEAN is turned off by default" OFF) #OFF by default
if ( NOT BUILD_MIT_OCEAN )
  esma_add_library (${this}
    SRCS GuestOcean_GridComp.F90
    SUBCOMPONENTS MOM_GEOS5PlugMod MOM6_GEOSPlug
    DEPENDENCIES GEOSdatasea_GridComp MAPL
    INCLUDES ${INC_ESMF} ${INC_NETCDF})
else()
  esma_add_library (${this}
    SRCS GuestOcean_GridComp.F90
    SUBCOMPONENTS MIT_GEOS5PlugMod
    DEPENDENCIES GEOSdatasea_GridComp MAPL
    INCLUDES ${INC_ESMF} ${INC_NETCDF})
endif()




