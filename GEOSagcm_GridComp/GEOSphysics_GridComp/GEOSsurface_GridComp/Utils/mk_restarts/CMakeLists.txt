esma_set_this ()

set(srcs
  CatchmentRst.F90
  CatchmentCNRst.F90
  getids.F90
)

set (exe_srcs
  cv_SaltRestart.F90
  SaltIntSplitter.F90
  SaltImpConverter.F90
  mk_CICERestart.F90
  mk_LakeLandiceSaltRestarts.F90
  mk_RouteRestarts.F90
  mk_catchANDcnRestarts.F90
)

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES MAPL GEOS_SurfaceShared GEOS_LandShared ESMF::ESMF NetCDF::NetCDF_Fortran GEOS_CatchCNShared)

foreach (src ${exe_srcs})
  string (REGEX REPLACE ".F90" ".x" exe ${src})
  ecbuild_add_executable (
    TARGET ${exe}
    SOURCES ${src}
    LIBS MAPL GFTL_SHARED::gftl-shared GEOS_SurfaceShared GEOSroute_GridComp GEOS_LandShared GEOS_CatchCNShared ${this})
endforeach ()

foreach (src ${exe_srcs})
  string (REGEX REPLACE ".F90" ".x" exe ${src})
  string (REGEX REPLACE ".F90" ""  lname ${src})
  install(CODE "execute_process( \
    COMMAND ${CMAKE_COMMAND} -E create_symlink \
    ${CMAKE_INSTALL_PREFIX}/bin/${exe} \
    ${CMAKE_INSTALL_PREFIX}/bin/${lname}  \
    )")
endforeach ()
