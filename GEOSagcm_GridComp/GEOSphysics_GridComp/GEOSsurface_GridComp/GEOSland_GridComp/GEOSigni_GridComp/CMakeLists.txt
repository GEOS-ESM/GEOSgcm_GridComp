esma_set_this ()

set (srcs
  GEOS_IgniGridComp.F90
  cffwi.F90
  )


esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES MAPL GEOS_Shared esmf
  )

