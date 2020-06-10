esma_set_this ()

set (srcs
   rte/mo_rte_kind.F90
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
   extensions/rng/mo_rng.F90
   extensions/mo_compute_bc.F90
   extensions/mo_fluxes_byband.F90
   extensions/mo_fluxes_byband_kernels.F90
   extensions/mo_fluxes_bygpoint.F90
   extensions/mo_heating_rates.F90
   extensions/mo_rrtmgp_clr_all_sky.F90
   extensions/cloud_optics/mo_cloud_optics.F90
   extensions/cloud_optics/mo_cloud_sampling.F90
   examples/mo_simple_netcdf.F90
   examples/mo_load_coefficients.F90
   examples/mo_load_cloud_coefficients.F90
   rte/mo_rte_lw.F90
   rte/mo_rte_sw.F90
   )

if (MKL_FOUND)
   list (APPEND srcs
      extensions/rng/mo_rng_mklvsl.F90
      extensions/rng/mo_rng_mklvsl_plus.F90
      )
else ()
   ecbuild_warn("MKL was not found. The mklvsl rng will not be compiled. RRTMGP should not be used.")
endif()

esma_add_library(${this} 
   SRCS ${srcs})

string(REPLACE " " ";" tmp ${FREAL8})
foreach(flag ${tmp})
  target_compile_options (${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${flag}>)
endforeach()

target_include_directories (${this} PRIVATE ${MKL_INCLUDE_DIRS} ${INC_NETCDF})

set_target_properties (${this} PROPERTIES Fortran_MODULE_DIRECTORY ${include_${this}})
target_link_libraries (${this} ${NETCDF_LIBRARIES})

