set (this GEOS_LandShared)

set (resource_files
    GEOS_LandGridComp.rc
   ) 

install( FILES ${resource_files} 
   DESTINATION etc
   )
set (srcs
  catch_constants.f90 lsm_routines.F90 sibalb_coeff.f90
  )

set_source_files_properties (catch_constants.f90 PROPERTIES COMPILE_FLAGS ${PP})

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES GEOS_SurfaceShared GEOS_Shared MAPL_Base)
target_include_directories (${this} PUBLIC ${INC_ESMF})
