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
sibalb.f
zenith.f
zip.c
util.c
)  


include_directories(${INC_ESMF})
include_directories(${INC_NETCDF})

esma_add_library(${this} SRCS ${srcs} DEPENDENCIES MAPL_Base)

ecbuild_add_executable (TARGET chk_clsm_params SOURCES chk_clsm_params.F90  LIBS MAPL_Base ${this})
ecbuild_add_executable (TARGET CombineRasters SOURCES CombineRasters.F90 LIBS MAPL_Base ${this})
ecbuild_add_executable (TARGET mkCatchParam SOURCES mkCatchParam.F90 LIBS MAPL_Base ${this})
ecbuild_add_executable (TARGET mkCatchParam_openmp SOURCES mkCatchParam.F90 LIBS MAPL_Base ${this})
ecbuild_add_executable (TARGET mkCubeFVRaster SOURCES mkCubeFVRaster.F90 LIBS MAPL_Base ${this})
ecbuild_add_executable (TARGET mkLandRaster SOURCES mkLandRaster.F90 LIBS MAPL_Base ${this})
ecbuild_add_executable (TARGET mkLatLonRaster SOURCES mkLatLonRaster.F90 LIBS MAPL_Base ${this})
ecbuild_add_executable (TARGET mkLISTilesPara SOURCES mkLISTilesPara.F90 LIBS MAPL_Base ${this})
ecbuild_add_executable (TARGET mkMITAquaRaster SOURCES mkMITAquaRaster.F90 LIBS MAPL_Base ${this})
ecbuild_add_executable (TARGET mkMOMAquaRaster SOURCES mkMOMAquaRaster.F90 LIBS MAPL_Base ${this})
ecbuild_add_executable (TARGET mk_runofftbl SOURCES mk_runofftbl.F90 LIBS MAPL_Base ${this})
ecbuild_add_executable (TARGET mkSMAPTilesPara SOURCES mkSMAPTilesPara.F90 LIBS MAPL_Base ${this})
ecbuild_add_executable (TARGET mkSMAPTilesPara_v2 SOURCES mkSMAPTilesPara_v2.F90 LIBS MAPL_Base ${this})
ecbuild_add_executable (TARGET rmTinyTiles SOURCES rmTinyTiles.F90 LIBS MAPL_Base ${this})
