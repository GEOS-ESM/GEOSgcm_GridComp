esma_set_this ()

esma_add_library (${this} 
   SRCS GEOS_MITDynaGridComp.F90 
   DEPENDENCIES MAPL
   INCLUDES ${INC_ESMF} ${INC_NETCDF})

target_compile_definitions(${this} PRIVATE DIAGOUT MODIFY_TOPOGRAPHY USE_R8)



