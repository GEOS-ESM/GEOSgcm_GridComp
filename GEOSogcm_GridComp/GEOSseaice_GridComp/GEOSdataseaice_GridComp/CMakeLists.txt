esma_set_this ()

esma_add_library (${this} SRCS GEOS_DataSeaIceGridComp.F90 DEPENDENCIES MAPL_Base LANL_cice)
target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})

