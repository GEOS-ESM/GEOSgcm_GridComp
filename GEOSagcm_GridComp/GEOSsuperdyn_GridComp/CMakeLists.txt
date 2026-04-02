esma_set_this ()

set (alldirs
  FVdycoreCubed_GridComp
  # FVdycore_GridComp       # not yet ported to MAPL3 - restore when ported
  # GEOSdatmodyn_GridComp   # not yet ported to MAPL3 - restore when ported
  # ARIESg3_GridComp        # not yet ported to MAPL3 - restore when ported
  )

# if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/GEOSdatmodyn_GridComp)
#
#   set (rename_GEOSdatmodyn_GridComp GEOS_DatmodynGridComp)
#
# endif ()

# Temporarily bypass GEOSsuperdyn top-level library: GEOS_SuperdynGridComp.F90
# USE-s non-ported components. Restore the block below when all SUBCOMPONENTS
# are ported to MAPL3.
esma_add_subdirectories (${alldirs})

# if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/GEOS_SuperdynGridComp.F90)
#
#   esma_add_library (${this}
#     SRCS GEOS_SuperdynGridComp.F90
#     SUBCOMPONENTS ${alldirs}
#     DEPENDENCIES MAPL GEOS_Shared ESMF::ESMF)
#
#   target_compile_definitions (${this} PRIVATE $<$<BOOL:${BUILD_WITH_GIGATRAJ}>:HAS_GIGATRAJ>)
#
# else ()
#
#   esma_add_subdirectories (${alldirs})
#
# endif ()
