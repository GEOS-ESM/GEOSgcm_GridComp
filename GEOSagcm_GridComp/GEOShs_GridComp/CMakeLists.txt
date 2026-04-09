esma_set_this()

set (srcs
  GEOS_HsGridComp.F90
  )

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES GEOS_Shared MAPL)
target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})
