esma_set_this ()

set (alldirs
  FVdycore_GridComp
  FVdycoreCubed_GridComp
  GEOSdatmodyn_GridComp
  ARIESg3_GridComp
  )

esma_add_library (${this}
  SRCS GEOS_SuperdynGridComp.F90
  SUBCOMPONENTS ${alldirs}
  DEPENDENCIES MAPL_Base GEOS_Shared
  INCLUDES ${INC_ESMF})



