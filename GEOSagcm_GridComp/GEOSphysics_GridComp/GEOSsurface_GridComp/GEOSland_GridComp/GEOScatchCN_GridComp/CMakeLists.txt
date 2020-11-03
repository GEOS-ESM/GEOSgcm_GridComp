esma_set_this ()

esma_add_subdirectories (Shared)
esma_add_subdirectories (CLM40 CLM45)

set (srcs
  GEOS_CatchCNGridComp.F90) 

esma_add_library (${this}
    SRCS ${srcs}
    DEPENDENCIES MAPL GEOS_Shared GEOS_LandShared CLM40 CLM45 GEOS_CatchCNShared)

