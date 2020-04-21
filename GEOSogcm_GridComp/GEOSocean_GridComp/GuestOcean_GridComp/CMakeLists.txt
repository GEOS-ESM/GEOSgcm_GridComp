esma_set_this ()

set (rename_MOM_GEOS5PlugMod MOM_GEOS5Plug)

esma_add_library (${this}
  SRCS GuestOcean_GridComp.F90
  SUBCOMPONENTS MOM_GEOS5PlugMod MOM6_GEOSPlug
  DEPENDENCIES GEOSdatasea_GridComp MAPL
  INCLUDES ${INC_ESMF} ${INC_NETCDF})




