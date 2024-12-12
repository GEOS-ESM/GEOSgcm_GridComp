esma_set_this ()

set (srcs
   GEOS_DataAtmGridComp.F90
  )

esma_add_library (${this}
   SRCS ${srcs}
   DEPENDENCIES MAPL GEOSsurface_GridComp GEOS_Shared ESMF::ESMF)

target_compile_definitions (${this} PRIVATE USE_CICE USE_R8)

install(
   FILES JRA55-DO_DataAtm_Forcings_ExtData.yaml CORE_NYF_Data_AtmForcings_ExtData.yaml
   DESTINATION etc
   )
