esma_set_this()

esma_add_library (${this}
  SRCS GEOS_LakeGridComp.F90
  DEPENDENCIES GEOS_Shared MAPL
  INCLUDES ${INC_ESMF})

