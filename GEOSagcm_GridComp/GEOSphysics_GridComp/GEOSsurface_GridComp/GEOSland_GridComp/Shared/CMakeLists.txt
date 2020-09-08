esma_set_this(OVERRIDE GEOS_LandShared)

set (resource_files
    GEOS_LandGridComp_ExtData.rc
   ) 

install( FILES ${resource_files} 
   DESTINATION etc
   )

set (srcs
  catch_constants.f90 lsm_routines.F90 sibalb_coeff.f90 GEOSland_modules.F90
  )

set_source_files_properties (catch_constants.f90 PROPERTIES COMPILE_FLAGS ${PP})

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES GEOS_SurfaceShared GEOS_Shared MAPL)
target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF} ${INC_HDF5} LIBS hdf5hl_fortran hdf5_fortran MAPL)
