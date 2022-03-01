esma_set_this ()

esma_add_library (${this} SRCS GEOS_DataSeaGridComp.F90 DEPENDENCIES MAPL)
target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})


