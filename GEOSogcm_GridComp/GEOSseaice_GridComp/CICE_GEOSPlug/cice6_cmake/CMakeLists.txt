esma_set_this (OVERRIDE cice6)


# got below list from following two commands:

# 1. ls -1 icepack/columnphysics/*F90
# 2. ls -1 cicecore/shared/*0
# 3. find cicecore/cicedynB/ -not -path "*cicedynB/infrastructure*" -print | grep .F90
# 4. find cicecore/cicedynB/infrastructure/ -not -path "*comm/serial*" -not -path "*io/io_pio2*" -not -path "*io/io_binary*" -print | grep .F90
#
list( APPEND CICE6_SRCS
   icepack/columnphysics/icepack_aerosol.F90
   icepack/columnphysics/icepack_age.F90
   icepack/columnphysics/icepack_algae.F90
   icepack/columnphysics/icepack_atmo.F90
   icepack/columnphysics/icepack_brine.F90
   icepack/columnphysics/icepack_firstyear.F90
   icepack/columnphysics/icepack_flux.F90
   icepack/columnphysics/icepack_fsd.F90
   icepack/columnphysics/icepack_intfc.F90
   icepack/columnphysics/icepack_isotope.F90
   icepack/columnphysics/icepack_itd.F90
   icepack/columnphysics/icepack_kinds.F90
   icepack/columnphysics/icepack_mechred.F90
   icepack/columnphysics/icepack_meltpond_lvl.F90
   icepack/columnphysics/icepack_meltpond_topo.F90
   icepack/columnphysics/icepack_mushy_physics.F90
   icepack/columnphysics/icepack_ocean.F90
   icepack/columnphysics/icepack_orbital.F90
   icepack/columnphysics/icepack_parameters.F90
   icepack/columnphysics/icepack_shortwave.F90
   icepack/columnphysics/icepack_snow.F90
   icepack/columnphysics/icepack_therm_bl99.F90
   icepack/columnphysics/icepack_therm_itd.F90
   icepack/columnphysics/icepack_therm_mushy.F90
   icepack/columnphysics/icepack_therm_shared.F90
   icepack/columnphysics/icepack_therm_vertical.F90
   icepack/columnphysics/icepack_tracers.F90
   icepack/columnphysics/icepack_warnings.F90
   icepack/columnphysics/icepack_wavefracspec.F90
   icepack/columnphysics/icepack_zbgc.F90
   icepack/columnphysics/icepack_zbgc_shared.F90
   icepack/columnphysics/icepack_zsalinity.F90
   cicecore/shared/ice_arrays_column.F90
   cicecore/shared/ice_calendar.F90
   cicecore/shared/ice_constants.F90
   cicecore/shared/ice_distribution.F90
   cicecore/shared/ice_domain_size.F90
   cicecore/shared/ice_fileunits.F90
   cicecore/shared/ice_init_column.F90
   cicecore/shared/ice_kinds_mod.F90
   cicecore/shared/ice_restart_column.F90
   cicecore/shared/ice_restart_shared.F90
   cicecore/shared/ice_spacecurve.F90
   cicecore/cicedynB/dynamics/ice_dyn_eap.F90
   cicecore/cicedynB/dynamics/ice_dyn_shared.F90
   cicecore/cicedynB/dynamics/ice_transport_driver.F90
   cicecore/cicedynB/dynamics/ice_dyn_evp_1d.F90
   cicecore/cicedynB/dynamics/ice_transport_remap.F90
   cicecore/cicedynB/dynamics/ice_dyn_vp.F90
   cicecore/cicedynB/dynamics/ice_dyn_evp.F90
   cicecore/cicedynB/general/ice_init.F90
   cicecore/cicedynB/general/ice_state.F90
   cicecore/cicedynB/general/ice_flux_bgc.F90
   cicecore/cicedynB/general/ice_step_mod.F90
   cicecore/cicedynB/general/ice_forcing_bgc.F90
   cicecore/cicedynB/general/ice_flux.F90
   cicecore/cicedynB/general/ice_forcing.F90
   cicecore/cicedynB/analysis/ice_diagnostics_bgc.F90
   cicecore/cicedynB/analysis/ice_history_drag.F90
   cicecore/cicedynB/analysis/ice_history_shared.F90
   cicecore/cicedynB/analysis/ice_history_mechred.F90
   cicecore/cicedynB/analysis/ice_history_bgc.F90
   cicecore/cicedynB/analysis/ice_history.F90
   cicecore/cicedynB/analysis/ice_history_fsd.F90
   cicecore/cicedynB/analysis/ice_diagnostics.F90
   cicecore/cicedynB/analysis/ice_history_snow.F90
   cicecore/cicedynB/analysis/ice_history_pond.F90
   cicecore/cicedynB/infrastructure/ice_domain.F90
   cicecore/cicedynB/infrastructure/ice_grid.F90
   cicecore/cicedynB/infrastructure/ice_memusage.F90
   cicecore/cicedynB/infrastructure/ice_restoring.F90
   cicecore/cicedynB/infrastructure/ice_restart_driver.F90
   cicecore/cicedynB/infrastructure/comm/mpi/ice_global_reductions.F90
   cicecore/cicedynB/infrastructure/comm/mpi/ice_reprosum.F90
   cicecore/cicedynB/infrastructure/comm/mpi/ice_exit.F90
   cicecore/cicedynB/infrastructure/comm/mpi/ice_timers.F90
   cicecore/cicedynB/infrastructure/comm/mpi/ice_communicate.F90
   cicecore/cicedynB/infrastructure/comm/mpi/ice_broadcast.F90
   cicecore/cicedynB/infrastructure/comm/mpi/ice_gather_scatter.F90
   cicecore/cicedynB/infrastructure/comm/mpi/ice_boundary.F90
   cicecore/cicedynB/infrastructure/io/io_netcdf/ice_restart.F90
   cicecore/cicedynB/infrastructure/io/io_netcdf/ice_history_write.F90
   cicecore/cicedynB/infrastructure/ice_blocks.F90
   cicecore/cicedynB/infrastructure/ice_read_write.F90
)

list( APPEND CICE6_SRCS
   # a single c source file
   cicecore/cicedynB/infrastructure/ice_shr_reprosum86.c
   cicecore/cicedynB/infrastructure/ice_memusage_gptl.c
   # drivers for mapl: ls -1 cicecore/drivers/mapl/geos/*0
   cicecore/drivers/mapl/geos/CICE_FinalMod.F90
   cicecore/drivers/mapl/geos/CICE_InitMod.F90
   cicecore/drivers/mapl/geos/CICE_RunMod.F90
   cicecore/drivers/mapl/geos/ice_import_export.F90
   cicecore/drivers/mapl/geos/ice_shr_methods.F90
   cicecore/drivers/mapl/geos/ice_prescribed_mod.F90
)


esma_mepo_style(cice6 CICE6_rel_path REL_PATH ..)
set (CICE6_path ${CMAKE_CURRENT_SOURCE_DIR}/${CICE6_rel_path})

set (SRCS)
foreach (file ${CICE6_SRCS})
  list (APPEND SRCS ${CICE6_path}/${file})
endforeach ()

esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES ESMF NetCDF::NetCDF_Fortran
  TYPE SHARED
)

target_compile_definitions (${this} PRIVATE USE_NETCDF GEOSCOUPLED FORTRANUNDERSCORE)

