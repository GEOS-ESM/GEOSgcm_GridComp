esma_set_this ()

set (SRCS
  GEOS_CICE4ColumnPhysGridComp.F90
  )

esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES GEOS_Shared MAPL CICE4
  TYPE SHARED
  )

