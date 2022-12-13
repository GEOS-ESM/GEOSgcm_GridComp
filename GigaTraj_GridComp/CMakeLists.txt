esma_set_this()

set (dependencies MAPL esmf geos_giga metsources filters gigatraj)

esma_add_library (${this}
  SRCS GEOS_Giga_InterOp.F90 GigaTraj_GridComp.F90
  DEPENDENCIES ${dependencies})

esma_add_subdirectories( @GigaTraj)
