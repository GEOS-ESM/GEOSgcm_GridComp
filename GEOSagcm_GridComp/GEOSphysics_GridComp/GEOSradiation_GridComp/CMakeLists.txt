esma_set_this()

set (alldirs
  GEOS_RadiationShared
  GEOSsolar_GridComp
  GEOSirrad_GridComp
  GEOSsatsim_GridComp
  )

esma_add_library (${this}
  SRCS GEOS_RadiationGridComp.F90
  SUBCOMPONENTS ${alldirs}
  DEPENDENCIES MAPL_Base GEOS_Shared
  INCLUDES ${INC_ESMF})

