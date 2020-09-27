# Note that directory and library name do not correspond!
esma_set_this (OVERRIDE GEOS_CatchCNRoot)
esma_add_subdirectories (CLM45)

set (srcs
  catchmentCN.F90
  catchcn_iau.F90)

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES MAPL GEOS_Shared GEOS_LandShared)
