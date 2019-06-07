esma_set_this ()

set (alldirs
#  Shared
  GEOSlandice_GridComp
  GEOSlake_GridComp
  GEOSland_GridComp
  GEOSsaltwater_GridComp
  )

esma_add_library (${this}
  SRCS GEOS_SurfaceGridComp.F90
  SUBCOMPONENTS ${alldirs}
  SUBDIRS Shared
  DEPENDENCIES GEOS_SurfaceShared MAPL_Base GMAO_mpeu
  INCLUDES ${INC_ESMF})




