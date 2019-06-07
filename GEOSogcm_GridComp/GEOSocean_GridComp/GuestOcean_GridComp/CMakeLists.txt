esma_set_this ()


esma_add_library (${this}
  SRCS GuestOcean_GridComp.F90
  SUBCOMPONENTS MOM_GEOS5PlugMod
  DEPENDENCIES GEOSdatasea_GridComp MAPL_Base
  INCLUDES ${INC_ESMF} ${INC_NETCDF})




