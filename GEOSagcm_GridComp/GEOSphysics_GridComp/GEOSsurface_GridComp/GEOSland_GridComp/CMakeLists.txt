esma_set_this ()

set (alldirs
  GEOSvegdyn_GridComp
  GEOScatch_GridComp
  GEOSlana_GridComp
  GEOSroute_GridComp
  )

# For CatchCN, the actual directory has no library, thus it cannot be a "SUBCOMPONENT".
# Instead, we add it as a SUBDIR and then add the library inside it (called CLM40 because of
# esma_set_this() makes a library named for the directory. Add that as a DEPENDENCIES

esma_add_library (${this}
  SRCS GEOS_LandGridComp.F90
  SUBCOMPONENTS ${alldirs}
  SUBDIRS Shared GEOScatchCN_GridComp
  DEPENDENCIES GEOS_SurfaceShared CLM40 GEOS_LandShared MAPL_Base
  INCLUDES ${INC_ESMF})

# Special case for GEOScatch and GEOScatchCN components - build with / without openmp
if (USE_OPENMP)
  add_subdirectory (GEOScatch_GridComp GEOScatch_GridComp_openmp)
  add_subdirectory (GEOScatchCN_GridComp GEOScatchCN_GridComp_openmp)
endif ()

