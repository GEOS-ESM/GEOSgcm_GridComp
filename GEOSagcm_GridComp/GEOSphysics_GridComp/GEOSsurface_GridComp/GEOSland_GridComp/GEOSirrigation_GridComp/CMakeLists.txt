esma_set_this (OVERRIDE GEOS_Irrigation)

set (srcs
   GEOS_IrrigationGridComp.F90
   irrigation_model.F90
  )

include_directories (${INC_ESMF})

esma_add_library (${this} SRCS ${srcs} DEPENDENCIES MAPL)
target_include_directories (${this} PUBLIC ${INC_ESMF})


