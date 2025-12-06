esma_set_this ()

set(srcs
  constant.f90
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
)
foreach (src ${exe_srcs})
  string (REGEX REPLACE ".f90" ".x" exe ${src})
  ecbuild_add_executable (
    TARGET ${exe}
    SOURCES ${src}
    NOINSTALL
    LIBS    ${this} NetCDF::NetCDF_Fortran MPI::MPI_Fortran)
endforeach ()

# copy to the build directory
configure_file(run_routing_preproc.py run_routing_preproc.py)