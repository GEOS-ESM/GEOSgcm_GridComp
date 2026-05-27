esma_set_this ()

esma_add_library (${this}
  SRCS GEOS_LandIceGridComp.F90
  DEPENDENCIES MAPL GEOS_Shared GEOS_SurfaceShared ESMF::ESMF NetCDF::NetCDF_Fortran
  )
