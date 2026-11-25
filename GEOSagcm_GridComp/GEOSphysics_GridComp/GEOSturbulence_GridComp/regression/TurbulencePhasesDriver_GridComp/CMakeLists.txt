esma_set_this (OVERRIDE TurbulencePhasesDriver_GridComp)

esma_add_library(${this}
  SRCS TurbulencePhasesDriver_GridCompMod.F90
  DEPENDENCIES MAPL.generic MAPL ESMF::ESMF
  TYPE SHARED)
