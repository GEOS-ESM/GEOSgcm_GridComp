esma_set_this ()

set (alldirs
  GEOScatchCNCLM40_GridComp
  GEOScatchCNCLM45_GridComp
  )

esma_add_library (${this}
  SRCS GEOS_CatchCNGridComp.F90 
  SUBCOMPONENTS ${alldirs}
  SUBDIRS Shared 
  DEPENDENCIES GEOS_SurfaceShared GEOS_LandShared GEOS_CatchCNShared MAPL
  INCLUDES ${INC_ESMF})

# Special case for GEOScatch and GEOScatchCN components - build with / without openmp
if (USE_OPENMP)
  add_subdirectory (GEOScatchCN_GridComp GEOScatchCN_GridComp_openmp)
endif ()
