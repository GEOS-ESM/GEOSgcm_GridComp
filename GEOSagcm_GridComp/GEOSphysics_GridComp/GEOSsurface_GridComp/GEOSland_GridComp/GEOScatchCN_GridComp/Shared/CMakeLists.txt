esma_set_this(OVERRIDE GEOS_CatchCNShared)

set (srcs
  nanMod.F90
  catchmentCN.F90
  catchcn_iau.F90
  clm_varpar_shared.F90
  )

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES MAPL GEOS_Shared GEOS_LandShared
  TYPE SHARED)
