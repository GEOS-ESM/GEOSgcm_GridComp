esma_set_this()

set (srcs
LDAS_DateTimeMod.F90
EASE_conv.F90
mod_process_hres_data.F90
rasterize.F90
read_riveroutlet.F90
CubedSphere_GridMod.F90
rmTinyCatchParaMod.F90
zip.c
util.c
)

if(NOT FORTRAN_COMPILER_SUPPORTS_FINDLOC)
   list(APPEND srcs findloc.F90)
endif ()

set_source_files_properties(mkMITAquaRaster.F90 PROPERTIES COMPILE_FLAGS "${BYTERECLEN}")

esma_add_library(${this} SRCS ${srcs} DEPENDENCIES MAPL GEOS_SurfaceShared GEOS_LandShared ESMF::ESMF NetCDF::NetCDF_Fortran OpenMP::OpenMP_Fortran)

if(NOT FORTRAN_COMPILER_SUPPORTS_FINDLOC)
   target_compile_definitions(${this} PRIVATE USE_EXTERNAL_FINDLOC)
endif ()

# MAT NOTE This should use find_package(ZLIB) but Baselibs currently
#     confuses find_package(). This is a hack until Baselibs is
#     reorganized.
if (Baselibs_FOUND)
  set (INC_ZLIB ${BASEDIR}/include/zlib)
  target_include_directories(${this} PRIVATE ${INC_ZLIB})
else ()
  find_package(ZLIB)
  target_link_libraries(${this} PRIVATE ZLIB::zlib)
  # Spack testing shows we need stdlib.h, so this ifdef will
  # enable it in zip.c
  target_compile_definitions(${this} PRIVATE STDC)
endif ()

ecbuild_add_executable (TARGET CombineRasters.x SOURCES CombineRasters.F90 LIBS MAPL ${this})
ecbuild_add_executable (TARGET mkCatchParam.x SOURCES mkCatchParam.F90 LIBS MAPL ${this} OpenMP::OpenMP_Fortran)
ecbuild_add_executable (TARGET mkCubeFVRaster.x SOURCES mkCubeFVRaster.F90 LIBS MAPL ${this})
ecbuild_add_executable (TARGET mkLandRaster.x SOURCES mkLandRaster.F90 LIBS MAPL ${this})
ecbuild_add_executable (TARGET mkLatLonRaster.x SOURCES mkLatLonRaster.F90 LIBS MAPL ${this})
ecbuild_add_executable (TARGET mkMITAquaRaster.x SOURCES mkMITAquaRaster.F90 LIBS MAPL ${this})
ecbuild_add_executable (TARGET mkMOMAquaRaster.x SOURCES mkMOMAquaRaster.F90 LIBS MAPL ${this})
ecbuild_add_executable (TARGET FillMomGrid.x SOURCES FillMomGrid.F90 LIBS MAPL ${this})
ecbuild_add_executable (TARGET mk_runofftbl.x SOURCES mk_runofftbl.F90 LIBS MAPL ${this})
ecbuild_add_executable (TARGET mkEASETilesParam.x SOURCES mkEASETilesParam.F90 LIBS MAPL ${this})

install(PROGRAMS clsm_plots.pro create_README.csh  DESTINATION bin)
file(GLOB MAKE_BCS_PYTHON CONFIGURE_DEPENDS "./make_bcs*.py")
list(FILTER MAKE_BCS_PYTHON EXCLUDE REGEX "make_bcs_shared.py")
install(PROGRAMS  ${MAKE_BCS_PYTHON} DESTINATION bin)

set(file ./make_bcs_shared.py)
configure_file(${file} ${file} @ONLY)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${file} DESTINATION bin)
