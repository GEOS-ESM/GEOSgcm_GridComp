esma_set_this ()

set (SRCS
   GEOS_SeaiceInterfaceGridComp.F90
  )

esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES GEOS_Shared MAPL
  TYPE SHARED
  )

