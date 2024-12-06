esma_set_this()

set (dependencies MAPL ESMF::ESMF geos_giga metsources filters gigatraj)

esma_add_library (${this}
  SRCS GEOS_Giga_InterOp.F90 Gigatraj_Utils.F90 GEOS_GigatrajGridComp.F90 
  DEPENDENCIES ${dependencies})

esma_add_subdirectories( @GigaTraj)
