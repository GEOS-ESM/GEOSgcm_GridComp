esma_set_this ()

set(srcs
  CatchmentRst.F90
  CatchmentCNRst.F90
  getids.F90
)

set (exe_srcs
  Scale_Catch.F90
  Scale_CatchCN.F90
  cv_SaltRestart.F90
  SaltIntSplitter.F90
  SaltImpConverter.F90
  mk_CICERestart.F90
  mk_CatchCNRestarts.F90
  mk_CatchRestarts.F90
  mk_LakeLandiceSaltRestarts.F90
  mk_GEOSldasRestarts.F90
  mk_catchANDcnRestarts.F90
)

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES MAPL GEOS_SurfaceShared GEOS_LandShared ESMF::ESMF NetCDF::NetCDF_Fortran GEOS_CatchCNShared)

# NAG: MPI_ISend is intentionally called with varying Fortran types
if (CMAKE_Fortran_COMPILER_ID MATCHES "NAG")
  target_compile_options(${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-wmismatch=MPI_ISend>)
endif()

foreach (src ${exe_srcs})
  string (REGEX REPLACE ".F90" ".x" exe ${src})
  ecbuild_add_executable (
    TARGET ${exe}
    SOURCES ${src}
    LIBS MAPL GFTL_SHARED::gftl-shared GEOS_SurfaceShared GEOSroute_GridComp GEOS_LandShared GEOS_CatchCNShared ${this})
  # NAG: executable sources also call MPI_ISend with varying Fortran types
  if (CMAKE_Fortran_COMPILER_ID MATCHES "NAG")
    target_compile_options(${exe} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-wmismatch=MPI_ISend,NF_DEF_VAR>)
  endif()
endforeach ()

install(PROGRAMS mk_Restarts DESTINATION bin)
foreach (src ${exe_srcs})
  string (REGEX REPLACE ".F90" ".x" exe ${src})
  string (REGEX REPLACE ".F90" ""  lname ${src})
  install(CODE "execute_process( \
    COMMAND ${CMAKE_COMMAND} -E create_symlink \
    ${CMAKE_INSTALL_PREFIX}/bin/${exe} \
    ${CMAKE_INSTALL_PREFIX}/bin/${lname}  \
    )")
endforeach ()
