esma_set_this ()

set(srcs
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
  mk_RouteRestarts.F90
  mk_GEOSldasRestarts.F90
)  

include_directories(${INC_ESMF})
include_directories(${INC_NETCDF})

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES MAPL)

foreach (src ${exe_srcs})
  string (REGEX REPLACE ".F90" "" exe ${src})
  if (ESMA_USE_GFE_NAMESPACE)
    ecbuild_add_executable (
      TARGET ${exe}
      SOURCES ${src}
      LIBS MAPL GFTL_SHARED::gftl-shared GEOS_SurfaceShared GEOSroute_GridComp GEOS_LandShared ${this})
  else ()
    ecbuild_add_executable (
      TARGET ${exe}
      SOURCES ${src}
      LIBS MAPL gftl-shared GEOS_SurfaceShared GEOSroute_GridComp GEOS_LandShared ${this})
  endif ()
endforeach ()

install(PROGRAMS mk_Restarts DESTINATION bin)
