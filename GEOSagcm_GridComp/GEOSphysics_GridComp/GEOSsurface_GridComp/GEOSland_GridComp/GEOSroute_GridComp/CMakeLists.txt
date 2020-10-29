esma_set_this ()

set (srcs
   GEOS_RouteGridComp.F90
   routing_model.F90
  )

include_directories (${INC_ESMF} ${INC_NETCDF})

esma_add_library (${this} SRCS ${srcs} DEPENDENCIES MAPL)
target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})

ecbuild_add_executable (TARGET GEOSroute.x SOURCES GEOSroute.F90 INCLUDES ${INC_ESMF} ${INC_NETCDF} ${INC_HDF5} LIBS hdf5hl_fortran hdf5_fortran GEOS_LandShared MAPL ${this})

install(PROGRAMS build_rivernetwork.py run_rrm.x DESTINATION bin)

