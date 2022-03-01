esma_set_this()

set (srcs
  GEOS_TurbulenceGridComp.F90
  LockEntrain.F90
  shoc.F90
  edmf.F90
  scm_surface.F90
  shocparams.F90
  edmfparams.F90
  )

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES GEOS_Shared MAPL)
target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})
