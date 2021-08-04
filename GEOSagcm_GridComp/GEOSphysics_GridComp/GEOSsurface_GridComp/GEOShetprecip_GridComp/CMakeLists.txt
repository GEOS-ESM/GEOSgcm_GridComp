esma_set_this ()

set (SRCS
  random.F90
  GEOS_HetPrecipGridComp.F90
  )

esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES MAPL
  TYPE SHARED
  )

