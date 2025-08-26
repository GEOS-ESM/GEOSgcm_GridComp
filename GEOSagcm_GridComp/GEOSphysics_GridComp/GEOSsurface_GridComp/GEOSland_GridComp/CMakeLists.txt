esma_set_this ()

set (alldirs
  GEOSvegdyn_GridComp
  GEOScatch_GridComp
  GEOScatchCN_GridComp
  GEOSlana_GridComp
  GEOSroute_GridComp
  GEOSigni_GridComp
  GEOSirrigation_GridComp
  )

esma_add_library (${this}
  SRCS GEOS_LandGridComp.F90
  SUBCOMPONENTS ${alldirs}
  SUBDIRS Shared
  DEPENDENCIES GEOS_SurfaceShared GEOS_LandShared MAPL ESMF::ESMF)

