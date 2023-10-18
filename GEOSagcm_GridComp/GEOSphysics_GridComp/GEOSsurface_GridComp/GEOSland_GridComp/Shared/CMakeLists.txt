esma_set_this(OVERRIDE GEOS_LandShared)

set (srcs
  catch_constants.f90 lsm_routines.F90 sibalb_coeff.f90 catch_wrap_state.F90
  )

set_source_files_properties (catch_constants.f90 PROPERTIES COMPILE_FLAGS ${PP})

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES GEOS_SurfaceShared GEOS_Shared MAPL esmf
  TYPE SHARED)
