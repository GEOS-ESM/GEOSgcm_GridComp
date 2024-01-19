esma_set_this ()

set (alldirs
  GEOSirrigation_GridComp
  GEOSvegdyn_GridComp
  GEOScatch_GridComp
  GEOScatchCN_GridComp
  GEOSlana_GridComp
  GEOSroute_GridComp
  GEOSigni_GridComp
  )

esma_add_library (${this}
  SRCS GEOS_LandGridComp.F90
  SUBCOMPONENTS ${alldirs}
  SUBDIRS Shared
  DEPENDENCIES GEOS_SurfaceShared GEOS_LandShared MAPL esmf)

