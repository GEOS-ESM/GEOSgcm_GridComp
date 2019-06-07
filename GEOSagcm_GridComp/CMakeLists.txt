esma_set_this()

set (alldirs
  GEOSsuperdyn_GridComp
  GEOSphysics_GridComp
  )

esma_add_library (${this}
  SRCS GEOS_AgcmGridComp.F90
  SUBCOMPONENTS ${alldirs}
  DEPENDENCIES MAPL_Base GEOS_Shared Chem_Shared
  INCLUDES ${INC_ESMF})


