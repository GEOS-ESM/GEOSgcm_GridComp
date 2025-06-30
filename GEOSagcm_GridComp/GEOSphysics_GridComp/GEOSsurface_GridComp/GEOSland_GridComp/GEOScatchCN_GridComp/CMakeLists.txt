esma_set_this ()

add_subdirectory (GEOScatchCNCLM40_GridComp)
add_subdirectory (GEOScatchCNCLM45_GridComp)
add_subdirectory (GEOScatchCNCLM51_GridComp)

esma_add_library (${this}
  SRCS GEOS_CatchCNGridComp.F90
  SUBDIRS Shared
  DEPENDENCIES GEOS_SurfaceShared GEOS_LandShared GEOS_CatchCNShared MAPL ESMF::ESMF)

# Special case for GEOScatch and GEOScatchCN components - build with / without openmp
if (USE_OPENMP)
  add_subdirectory (GEOScatchCN_GridComp GEOScatchCN_GridComp_openmp)
endif ()
