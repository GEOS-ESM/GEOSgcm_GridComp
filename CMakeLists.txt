esma_set_this()

set (alldirs
  GEOSagcm_GridComp
  GEOSdataatm_GridComp
  GEOSmkiau_GridComp
  GEOSogcm_GridComp
  )

esma_add_library(${this}
  SRCS GEOS_GcmGridComp.F90
  SUBCOMPONENTS ${alldirs}
  DEPENDENCIES MAPL_Base
  INCLUDES ${INC_ESMF})


