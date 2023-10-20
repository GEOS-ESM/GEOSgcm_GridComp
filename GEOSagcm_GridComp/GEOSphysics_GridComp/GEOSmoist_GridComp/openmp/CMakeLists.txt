esma_set_this(OVERRIDE MOIST_OPENMP_OFFLOAD)

message(STATUS "Building GPU (OpenMP) ports of Moist routines")

set (srcs
  ../aer_cloud.F90
  ../aer_actv_single_moment.F90
  Process_Library.F90
  GEOS_GFDL_1M_InterfaceMod.F90
  gfdl_cloud_microphys.F90)

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES MAPL GEOS_Shared GMAO_mpeu)
  # DEPENDENCIES GEOS_Shared GMAO_mpeu MAPL Chem_Shared Chem_Base esmf)

get_target_property (extra_incs fms_r4 INCLUDE_DIRECTORIES)
target_include_directories(${this} PRIVATE
   $<BUILD_INTERFACE:${extra_incs}>
   )

# We need the SHELL prefix to get around the error where the flags are quoted as in
# gfortran: error: unrecognized command-line option '-fopenacc -foffload=nvptx-none’
target_compile_options(${this} PRIVATE SHELL:${GEOS_Fortran_OpenMP_Offload_Flags})
