esma_set_this ()

set (srcs
  gw_drag.F90
  GEOS_GwdGridComp.F90
  )

esma_add_library (${this} SRCS ${srcs} DEPENDENCIES GEOS_Shared MAPL)
target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})
