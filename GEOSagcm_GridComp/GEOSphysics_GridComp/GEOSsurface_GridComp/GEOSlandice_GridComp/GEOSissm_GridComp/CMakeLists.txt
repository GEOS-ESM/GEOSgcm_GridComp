find_package (ISSM REQUIRED COMPONENTS Core)

esma_set_this ()

esma_add_library (${this}
  SRCS GEOS_ISSMGridComp.F90
  DEPENDENCIES MAPL GEOS_Shared GEOS_SurfaceShared ESMF::ESMF NetCDF::NetCDF_Fortran ISSM::Core
  )

