esma_set_this()

set (srcs
  GEOS_IrradGridComp.F90
  irradconstants.F90
  irrad.F90
  )

esma_add_library (${this}
  SRCS ${srcs}
  SUBDIRS RRTMG
  DEPENDENCIES GEOS_Shared MAPL GEOS_RadiationShared RRTMG RRTMG_mods RRTMGP)

if(MKL_FOUND)
   target_compile_definitions(${this} PRIVATE HAVE_MKL)
endif ()
