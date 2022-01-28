esma_set_this (OVERRIDE RRTMGP)

# These are from RRTMGP
set (RRTMGP_SRCS
  rte/mo_rte_kind.F90
  rte/mo_rte_config.F90
  rte/mo_rte_util_array.F90
  rte/kernels/mo_rte_solver_kernels.F90
  rte/kernels/mo_optical_props_kernels.F90
  rte/kernels/mo_fluxes_broadband_kernels.F90
  rte/mo_optical_props.F90
  rte/mo_source_functions.F90
  rte/mo_fluxes.F90
  rrtmgp/mo_rrtmgp_constants.F90
  rrtmgp/mo_rrtmgp_util_reorder.F90
  rrtmgp/mo_rrtmgp_util_string.F90
  rrtmgp/mo_gas_optics.F90
  rrtmgp/mo_gas_optics_rrtmgp.F90
  rrtmgp/mo_gas_concentrations.F90
  rrtmgp/kernels/mo_gas_optics_kernels.F90
  rrtmgp/kernels/mo_rrtmgp_util_reorder_kernels.F90
  extensions/mo_compute_bc.F90
  extensions/mo_fluxes_byband.F90
  extensions/mo_fluxes_bygpoint.F90
  extensions/mo_heating_rates.F90
  extensions/mo_rrtmgp_clr_all_sky.F90
  extensions/cloud_optics/mo_cloud_optics.F90
  extensions/cloud_optics/mo_cloud_sampling.F90
  examples/mo_simple_netcdf.F90
  examples/mo_load_coefficients.F90
  examples/all-sky/mo_load_cloud_coefficients.F90
  rte/mo_rte_lw.F90
  rte/mo_rte_sw.F90
  )

esma_mepo_style(RRTMGP RRTMGP_rel_path REL_PATH ..)
set (RRTMGP_path ${CMAKE_CURRENT_SOURCE_DIR}/${RRTMGP_rel_path})

set (SRCS)
foreach (file ${RRTMGP_SRCS})
  list (APPEND SRCS ${RRTMGP_path}/${file})
endforeach ()

# Now we append some sources not in RRTMGP proper
list (APPEND SRC
  rng/mo_rng.F90
  )

# These files need MKL
if (APPLE)
   set(MKL_Fortran TRUE)
endif ()
find_package(MKL)
if (MKL_FOUND)
  set(BLA_VENDOR Intel10_64lp_seq)

  list (APPEND SRCS
    rng/mo_rng_mklvsl.F90
    rng/mo_rng_mklvsl_plus.F90
    )

else ()
  ecbuild_warn("MKL was not found. The mklvsl rng will not be compiled. RRTMGP should not be used.")
endif()

esma_add_library(${this}
  SRCS ${SRCS}
  DEPENDENCIES NetCDF::NetCDF_Fortran)

string(REPLACE " " ";" tmp ${FREAL8})
foreach(flag ${tmp})
  target_compile_options (${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${flag}>)
endforeach()

if (MKL_FOUND)
  target_include_directories (${this} PRIVATE ${MKL_INCLUDE_DIRS})
endif()

set_target_properties (${this} PROPERTIES Fortran_MODULE_DIRECTORY ${include_${this}})
