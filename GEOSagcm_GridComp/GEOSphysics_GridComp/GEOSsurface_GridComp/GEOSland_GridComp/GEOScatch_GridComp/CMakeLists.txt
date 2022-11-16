esma_set_this ()

set (srcs
  GEOS_CatchGridComp.F90
  catchment.F90
  catch_incr.F90
  )

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES MAPL GEOS_LandShared esmf NetCDF::NetCDF_Fortran)

