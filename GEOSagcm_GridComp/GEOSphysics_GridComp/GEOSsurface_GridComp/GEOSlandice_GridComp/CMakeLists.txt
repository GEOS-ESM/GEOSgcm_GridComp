esma_set_this ()

# First set alldirs to empty
set(alldirs)
# ... and srcs to the LandIceGridComp
set (srcs GEOS_LandIceGridComp.F90)

# Then try to find ISSM
find_package (ISSM QUIET COMPONENTS Core)

if(ISSM_FOUND)
  # If ISSM is found, append the directory to alldirs
  message(STATUS "ISSM found, building gridded component")
  list (APPEND alldirs GEOSissm_GridComp)
else()
  # If ISSM is not found, append the stub component to srcs
  message(STATUS "ISSM not found, using stub component")
  # For esma_create_stub_component, we need the name of the *module*
  # that the landice gc is expecting, minus the mod. Also we
  # have to use the variable srcs but not the ${srcs} because
  # esma_create_stub_component is expecting the list itself, not the
  # contents of the list.
  esma_create_stub_component(srcs GEOS_IssmGridComp)
endif()

# We have to do the above way because what esma_create_stub_component
# does is to append the stub component it creates in the build tree to a
# list of sources srcs.  So if ISSM is found, we go into the subdir and
# build the real component, and if ISSM is not found, we create the stub
# component and append it to srcs.

esma_add_library (${this}
  SRCS ${srcs}
  SUBCOMPONENTS ${alldirs}
  DEPENDENCIES MAPL GEOS_Shared GEOS_SurfaceShared ESMF::ESMF NetCDF::NetCDF_Fortran)
