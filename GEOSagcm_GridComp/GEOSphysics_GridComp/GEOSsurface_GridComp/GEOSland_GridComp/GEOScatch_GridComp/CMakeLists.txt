esma_set_this ()

set (srcs
  GEOS_CatchGridComp.F90
  catchment.F90
  catch_incr.F90         
  )

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES MAPL GEOS_LandShared)
target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})

