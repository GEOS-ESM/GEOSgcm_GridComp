esma_set_this ()

set(srcs
  outlets_constants.f90
)

set (exe_srcs
   get_finalID_msk.f90
   get_landocean_Greenland_real.f90
   get_outlets_land_allcat.f90
   get_sinkxy_land.f90
   get_outlets_catchindex.f90
   get_outlets_land.f90
   Pfaf_to_2d_30s_land.f90
)

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES GEOS_LandShared NetCDF::NetCDF_Fortran
)
foreach (src ${exe_srcs})
  string (REGEX REPLACE ".f90" ".x" exe ${src})
  ecbuild_add_executable (
    TARGET ${exe}
    SOURCES ${src}
    NOINSTALL
    LIBS    ${this} NetCDF::NetCDF_Fortran)
endforeach ()

# copy to the build directory
configure_file(run_routing_raster.py run_routing_raster.py)
