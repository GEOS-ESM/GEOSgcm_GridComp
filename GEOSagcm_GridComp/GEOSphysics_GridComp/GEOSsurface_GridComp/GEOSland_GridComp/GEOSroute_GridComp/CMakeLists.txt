esma_set_this ()

set (srcs
   #GEOS_RouteGridComp.F90
   routing_model.F90
  )

esma_add_library (${this} SRCS ${srcs} DEPENDENCIES MAPL esmf NetCDF::NetCDF_Fortran)

install(PROGRAMS build_rivernetwork.py DESTINATION bin)
