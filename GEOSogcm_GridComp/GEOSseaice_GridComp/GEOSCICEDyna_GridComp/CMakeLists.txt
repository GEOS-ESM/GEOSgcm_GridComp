esma_set_this ()

esma_add_library (${this}
   SRCS GEOS_CICEDynaGridComp.F90
   DEPENDENCIES MAPL CICE4 ESMF::ESMF)

target_compile_definitions(${this} PRIVATE DIAGOUT MODIFY_TOPOGRAPHY USE_R8)



