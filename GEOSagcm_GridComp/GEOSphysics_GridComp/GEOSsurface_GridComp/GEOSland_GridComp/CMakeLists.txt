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

set(CATCHCN_VERSION "4.0" CACHE STRING "CatchCN Version")

if ("${CATCHCN_VERSION}" STREQUAL "4.0")
   set (CATCHCN_COMPONENT CLM40)
elseif("${CATCHCN_VERSION}" STREQUAL "4.5")
   set (CATCHCN_COMPONENT CLM45)
else ()
   ecbuild_error ("Unrecognized CatchCN Version ${CATCHCN_VERSION}. Only 4.0 and 4.5 allowed")
endif ()
ecbuild_info ("Building CatchCN Version ${CATCHCN_VERSION}")
ecbuild_info ("Setting  CATCHCN_COMPONENT to ${CATCHCN_COMPONENT}")

esma_add_library (${this}
  SRCS GEOS_LandGridComp.F90
  SUBCOMPONENTS ${alldirs}
  SUBDIRS Shared GEOScatchCN_GridComp
  DEPENDENCIES GEOS_SurfaceShared ${CATCHCN_COMPONENT} GEOS_LandShared MAPL
  INCLUDES ${INC_ESMF})

# Special case for GEOScatch and GEOScatchCN components - build with / without openmp
if (USE_OPENMP)
  add_subdirectory (GEOScatch_GridComp GEOScatch_GridComp_openmp)
  add_subdirectory (GEOScatchCN_GridComp GEOScatchCN_GridComp_openmp)
endif ()

