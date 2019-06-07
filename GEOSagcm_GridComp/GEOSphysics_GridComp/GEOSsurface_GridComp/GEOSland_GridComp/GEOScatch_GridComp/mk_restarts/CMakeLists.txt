# No library
set (srcs
  Scale_Catch.F90			    
  Scale_CatchCN.F90		    
  cv_SaltRestart.F90		    
  SaltIntSplitter.F90
  SaltImpConverter.F90		    
  mk_CICERestart.F90		    
  mk_CatchCNRestarts.F90		    
  mk_CatchRestarts.F90		    
  mk_LDASsaRestarts.F90		    
  mk_LakeLandiceSaltRestarts.F90	    
  mk_RouteRestarts.F90		    
#  mk_catch_restart.F90		    
#  mk_vegdyn_restart.F90		    
  newcatch.F90			    
#  newvegdyn.f90			    
#  replace_params.F90		    
  strip_vegdyn.F90                    
)  


include_directories(${INC_ESMF})
include_directories(${INC_NETCDF})

foreach (src ${srcs})
  string (REGEX REPLACE ".F90" "" exe ${src})
  add_executable (${exe} ${src})
  target_link_libraries (${exe} PUBLIC MAPL_Base GEOS_SurfaceShared GEOSroute_GridComp GEOS_LandShared)
  install (
    TARGETS ${exe}
    DESTINATION ${CMAKE_BINARY_DIR}/bin)
endforeach ()

