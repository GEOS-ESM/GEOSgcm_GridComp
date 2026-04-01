esma_set_this ()

set (alldirs
  # GEOSchem_GridComp      # not yet ported to MAPL3 - restore when ported
  # GEOSmoist_GridComp     # not yet ported to MAPL3 - restore when ported
  # GEOSsurface_GridComp   # not yet ported to MAPL3 - restore when ported
  # GEOSturbulence_GridComp # not yet ported to MAPL3 - restore when ported
  GEOSgwd_GridComp
  # GEOSradiation_GridComp # not yet ported to MAPL3 - restore when ported
  )

# Temporarily bypass GEOSphysics top-level library: GEOS_PhysicsGridComp.F90
# USE-s non-ported components. Restore the block below when all SUBCOMPONENTS
# are ported to MAPL3.
esma_add_subdirectories (${alldirs})

# if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/GEOS_PhysicsGridComp.F90)
#
#   esma_add_library (${this}
#     SRCS GEOS_PhysicsGridComp.F90 MBundle_IncrementMod.F90
#     SUBCOMPONENTS ${alldirs}
#     DEPENDENCIES MAPL GMAO_mpeu GMAO_stoch ESMF::ESMF)
#
# else ()
#
#   esma_add_subdirectories (${alldirs})
#
# endif()

