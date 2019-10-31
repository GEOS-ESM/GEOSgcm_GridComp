set (this GEOSturbulence_GridComp)

set (srcs
  GEOS_TurbulenceGridComp.F90
  LockEntrain.F90
  shoc.F90
  sl3.F90
  )

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES GEOS_Shared MAPL_Base)
target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})
