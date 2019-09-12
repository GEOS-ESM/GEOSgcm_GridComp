esma_set_this ()

set (srcs
  tools/fv_treat_da_inc.F90
  model/lin_cloud_microphys.F90
  model/fv_cmp.F90
  geos_utils/fill_corner.F90
  geos_utils/flow_proj.F90
  geos_utils/ghost_cubsph.F90
  geos_utils/grid_utils_mod.F90
  geos_utils/cub2latlon.F90
  geos_utils/cub2cub.F90
  model/fv_arrays.F90
  tools/fv_mp_mod.F90
  tools/fv_eta.F90
  tools/fv_timing.F90
  model/fv_grid_utils.F90
  tools/sorted_index.F90
  tools/fv_grid_tools.F90
  model/a2b_edge.F90
  model/fv_fill.F90
  tools/init_hydro.F90
  model/fv_mapz.F90
  tools/fv_surf_map.F90
  model/fv_update_phys.F90
  model/fv_sg.F90
  tools/test_cases.F90
  tools/fv_diagnostics.F90
  tools/external_ic.F90
  tools/fv_restart.F90
  model/fv_control.F90
  model/tp_core.F90
  model/sw_core.F90
  model/fastexp.c
  model/nh_utils.F90
  model/nh_core.F90
  model/dyn_core.F90
  model/fv_tracer2d.F90
  model/fv_dynamics.F90
  tools/fv_io.F90
  tools/sim_nc_mod.F90
  model/boundary.F90
  tools/external_sst.F90
  tools/fv_nudge.F90
  model/fv_nesting.F90
  geos_utils/fv_regrid_c2c.F90
  geos_utils/fv_regridding_utils.F90
  geos_utils/rs_scaleMod.F90
  )

if (FV_PRECISION STREQUAL R4)
elseif (FV_PRECISION STREQUAL R4R8) # FV is R4 but FMS is R8
   get_target_property (extra_incs GFDL_fms_r4 INCLUDE_DIRECTORIES)
elseif (FV_PRECISION STREQUAL R8)
endif ()

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES MAPL_Base 
  INCLUDES ${extra_incs} ${INC_ESMF} ${INC_NETCDF})

if (FV_PRECISION STREQUAL R4)
   target_compile_definitions (${this} PRIVATE -DSINGLE_FV -DOVERLOAD_R4)
   target_link_libraries (${this} PUBLIC GFDL_fms_r4)
elseif (FV_PRECISION STREQUAL R4R8) # FV is R4 but FMS is R8
   target_compile_definitions (${this} PRIVATE -DSINGLE_FV -DOVERLOAD_R4)
   target_link_libraries (${this} PUBLIC GFDL_fms_r8)
elseif (FV_PRECISION STREQUAL R8)
   target_link_libraries (${this} PUBLIC GFDL_fms_r8)
   string(REPLACE " " ";" tmp ${FREAL8})
   foreach(flag ${tmp})
     target_compile_options (${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${flag}>)
   endforeach()
endif ()

set (CMAKE_Fortran_FLAGS_RELEASE "${GEOS_Fortran_FLAGS_VECT}")

if (CRAY_POINTER)
  set_target_properties (${this} PROPERTIES COMPILE_FLAGS ${CRAY_POINTER})
endif()

target_compile_definitions (${this} PRIVATE MAPL_MODE SPMD TIMING)




