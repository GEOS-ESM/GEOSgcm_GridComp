esma_set_this ()

set (srcs
   GEOS_DataAtmGridComp.F90
  )

esma_add_library (${this}
   SRCS ${srcs}
   DEPENDENCIES MAPL GEOSsurface_GridComp GEOS_Shared esmf)

target_compile_definitions (${this} PRIVATE USE_CICE USE_R8)

