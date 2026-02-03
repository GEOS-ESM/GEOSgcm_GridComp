esma_set_this ()

set(srcs
  routing_model_constants.f90
  river_read.f90
  k_module_cali.f90
)

set (exe_srcs
   get_Pfaf_file.f90
   get_num_sub_catchment.f90
   get_lonlat_bond.f90
   get_isub.f90
   get_area.f90
   get_Qr_clmt.f90
   get_river_length.f90
   get_K_model_calik.f90
   read_input_TopoCat.f90
)

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES MAPL GEOS_LandShared ESMF::ESMF NetCDF::NetCDF_Fortran)  

foreach (src ${exe_srcs})
  string (REGEX REPLACE ".f90" ".x" exe ${src})
  ecbuild_add_executable (
    TARGET ${exe}
    SOURCES ${src}
    NOINSTALL
    LIBS MAPL GFTL_SHARED::gftl-shared GEOSroute_GridComp GEOS_LandShared
    ${this}
    LIBS    ${this} NetCDF::NetCDF_Fortran MPI::MPI_Fortran)
endforeach ()

# copy to the build directory
set(SCRIPT_FILES
  create_river_input.py
  get_latloni_cellarea.py
  process_lake_data.py
  get_dam_data.py
  get_lonlati_maptile.py
  routing_model_constants.py
  run_routing_preproc.py
)

foreach (SCRIPT ${SCRIPT_FILES})
  configure_file(${SCRIPT} ${SCRIPT})
endforeach()
