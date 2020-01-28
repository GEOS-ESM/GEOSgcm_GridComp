esma_set_this(OVERRIDE raster)
set (srcs
date_time_util.F90 
leap_year.F90 
easeV1_conv.F90   
mod_process_hres_data.F90
easeV2_conv.F90  
rasterize.F90
read_riveroutlet.F90
CubedSphere_GridMod.F90 
rmTinyCatchParaMod.F90
comp_CATCHCN_AlbScale_parameters.F90
sibalb.f
zenith.f
zip.c
util.c
)  

if(NOT FORTRAN_COMPILER_SUPPORTS_FINDLOC)
   list(APPEND srcs findloc.F90)
endif ()

include_directories(${INC_ESMF})
include_directories(${INC_NETCDF})

esma_add_library(${this} SRCS ${srcs} DEPENDENCIES MAPL_Base)
target_link_libraries(${this} PRIVATE OpenMP::OpenMP_Fortran OpenMP::OpenMP_C)

if(NOT FORTRAN_COMPILER_SUPPORTS_FINDLOC)
   target_compile_definitions(${this} PRIVATE USE_EXTERNAL_FINDLOC)
endif ()

# MAT NOTE This should use find_package(ZLIB) but Baselibs currently
#     confuses find_package(). This is a hack until Baselibs is
#     reorganized.
set (INC_ZLIB ${BASEDIR}/include/zlib)
target_include_directories(${this} PRIVATE ${INC_ZLIB})

ecbuild_add_executable (TARGET chk_clsm_params.x SOURCES chk_clsm_params.F90  LIBS MAPL_Base ${this})
ecbuild_add_executable (TARGET CombineRasters.x SOURCES CombineRasters.F90 LIBS MAPL_Base ${this})
ecbuild_add_executable (TARGET mkCatchParam.x SOURCES mkCatchParam.F90 LIBS MAPL_Base ${this} ${OpenMP_Fortran_LIBRARIES})
set_target_properties(mkCatchParam.x PROPERTIES COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}")
set_target_properties(mkCatchParam.x PROPERTIES LINK_FLAGS "${OpenMP_Fortran_FLAGS}")

ecbuild_add_executable (TARGET mkCubeFVRaster.x SOURCES mkCubeFVRaster.F90 LIBS MAPL_Base ${this})
ecbuild_add_executable (TARGET mkLandRaster.x SOURCES mkLandRaster.F90 LIBS MAPL_Base ${this})
ecbuild_add_executable (TARGET mkLatLonRaster.x SOURCES mkLatLonRaster.F90 LIBS MAPL_Base ${this})
ecbuild_add_executable (TARGET mkLISTilesPara.x SOURCES mkLISTilesPara.F90 LIBS MAPL_Base ${this})
ecbuild_add_executable (TARGET mkMITAquaRaster.x SOURCES mkMITAquaRaster.F90 LIBS MAPL_Base ${this})
ecbuild_add_executable (TARGET mkMOMAquaRaster.x SOURCES mkMOMAquaRaster.F90 LIBS MAPL_Base ${this})
ecbuild_add_executable (TARGET mk_runofftbl.x SOURCES mk_runofftbl.F90 LIBS MAPL_Base ${this})
ecbuild_add_executable (TARGET mkSMAPTilesPara.x SOURCES mkSMAPTilesPara.F90 LIBS MAPL_Base ${this})
ecbuild_add_executable (TARGET mkSMAPTilesPara_v2.x SOURCES mkSMAPTilesPara_v2.F90 LIBS MAPL_Base ${this})

install(PROGRAMS make_bcs clsm_plots.pro plot_curves.pro create_README.csh plot_curves.csh DESTINATION bin)
