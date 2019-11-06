esma_set_this ()

set (srcs
   #GEOS_RouteGridComp.F90
   routing_model.F90
  )

include_directories (${INC_ESMF} ${INC_NETCDF})

esma_add_library (${this} SRCS ${srcs} DEPENDENCIES MAPL_Base)
target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})

install(PROGRAMS build_rivernetwork.py DESTINATION bin)