esma_set_this()

set (alldirs
  GEOSsuperdyn_GridComp
  GEOSphysics_GridComp
  )

# Temporarily bypass GEOSagcm top-level library: GEOS_AgcmGridComp.F90
# USE-s non-ported components (e.g. GEOS_superdynGridCompMod). Restore
# the block below when all SUBCOMPONENTS are ported to MAPL3.
esma_add_subdirectories (${alldirs})

# if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/GEOS_AgcmSimpleGridComp.F90)
#
#   message(STATUS "Building Held-Suarez GridComp")
#   esma_add_library (${this}
#     SRCS GEOS_AgcmSimpleGridComp.F90
#     SUBCOMPONENTS GEOSsuperdyn_GridComp GEOShs_GridComp
#     DEPENDENCIES MAPL ESMF::ESMF)
#
# elseif (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/GEOS_AgcmGridComp.F90)
#
#   esma_add_library (${this}
#     SRCS GEOS_AgcmGridComp.F90
#     SUBCOMPONENTS ${alldirs}
#     DEPENDENCIES MAPL GEOS_Shared Chem_Shared ESMF::ESMF)
#
#   target_compile_definitions (${this} PRIVATE $<$<BOOL:${BUILD_WITH_GIGATRAJ}>:HAS_GIGATRAJ>)
#
# else ()
#
#   esma_add_subdirectories (${alldirs})
#
# endif()
