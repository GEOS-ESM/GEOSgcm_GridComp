esma_set_this ()

set (srcs
   GEOS_RouteGridComp.F90
   routing_model.F90
   reservoir.F90
  )

if (CMAKE_Fortran_COMPILER_ID MATCHES GNU AND CMAKE_BUILD_TYPE MATCHES Release)
  set_source_files_properties(routing_model.F90 PROPERTIES COMPILE_OPTIONS ${FOPT2})
endif ()

esma_add_library (${this} SRCS ${srcs} DEPENDENCIES MAPL GEOS_LandShared ESMF::ESMF NetCDF::NetCDF_Fortran)
