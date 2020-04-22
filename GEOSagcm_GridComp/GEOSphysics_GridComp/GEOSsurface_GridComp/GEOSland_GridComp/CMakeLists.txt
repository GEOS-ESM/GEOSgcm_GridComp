esma_set_this ()

set (alldirs
  GEOSvegdyn_GridComp
  GEOScatch_GridComp
  GEOScatchCN_GridComp
  GEOSlana_GridComp
  GEOSroute_GridComp
  )

esma_add_library (${this}
  SRCS GEOS_LandGridComp.F90
  SUBCOMPONENTS ${alldirs}
  SUBDIRS Shared
  DEPENDENCIES GEOS_SurfaceShared GEOS_LandShared MAPL
  INCLUDES ${INC_ESMF})

# Special case for GEOScatch and GEOScatchCN components - build with / without openmp
if (USE_OPENMP)
  add_subdirectory (GEOScatch_GridComp GEOScatch_GridComp_openmp)
  add_subdirectory (GEOScatchCN_GridComp GEOScatchCN_GridComp_openmp)
endif ()

