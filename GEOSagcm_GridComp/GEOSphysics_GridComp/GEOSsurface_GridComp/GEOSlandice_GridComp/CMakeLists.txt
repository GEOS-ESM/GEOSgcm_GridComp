esma_set_this ()

set (alldirs
   GEOSissm_GridComp
   )

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/GEOS_LandIceGridComp.F90)
  esma_add_library (${this}
    SRCS GEOS_LandIceGridComp.F90
    SUBCOMPONENTS ${alldirs}
    DEPENDENCIES MAPL GEOS_Shared GEOS_SurfaceShared ESMF::ESMF NetCDF::NetCDF_Fortran)
else ()
  esma_add_subdirectories (${alldirs})
 
endif()