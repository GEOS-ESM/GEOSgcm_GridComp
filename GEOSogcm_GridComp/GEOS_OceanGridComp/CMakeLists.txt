esma_set_this ()

#target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})

option (BUILD_MIT_OCEAN "BUILD_MIT_OCEAN is turned off by default" OFF) #OFF by default
if ( NOT BUILD_MIT_OCEAN )
  esma_add_library (${this}
    SRCS GEOS_OceanGridComp.F90
    SUBCOMPONENTS MOM6_GEOSPlug MOM_GEOS5PlugMod GEOSdatasea_GridComp
    DEPENDENCIES GEOSdatasea_GridComp MAPL
    INCLUDES ${INC_ESMF})
else()
  esma_add_library (${this}
    SRCS GEOS_OceanGridComp.F90
    SUBCOMPONENTS MIT_GEOS5PlugMod
    DEPENDENCIES GEOSdatasea_GridComp MAPL
    INCLUDES ${INC_ESMF})
endif()
