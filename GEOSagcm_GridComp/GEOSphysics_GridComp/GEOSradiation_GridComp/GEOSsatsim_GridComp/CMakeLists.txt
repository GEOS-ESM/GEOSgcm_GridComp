esma_set_this ()

set (srcs
  actsim/lidar_simulator.F90
  actsim/lmd_ipsl_stats.F90
  llnl/cosp_radar.F90
  llnl/llnl_stats.F90
  llnl/pf_to_mr.f
  llnl/prec_scops.f
  quickbeam/array_lib.f90
  quickbeam/atmos_lib.f90
  quickbeam/dsd.f90
  quickbeam/format_input.f90
  quickbeam/gases.f90
  quickbeam/load_hydrometeor_classes.f90
#  quickbeam/load_mie_table.f90
  quickbeam/math_lib.f90
  quickbeam/mrgrnk.f90
  quickbeam/optics_lib.f90
  quickbeam/radar_simulator.f90
  quickbeam/radar_simulator_types.f90
  quickbeam/zeff.f90
  MISR_simulator.f
#  congvec.f
  cosp.F90
  cosp_constants.F90
  cosp_isccp_simulator.F90
  cosp_lidar.F90
  cosp_misr_simulator.F90
  cosp_modis_simulator.F90
  cosp_simulator.F90
  cosp_stats.F90
  cosp_types.F90
  cosp_utils.F90
  geos5_cosp_const.F90
  geos5_read.f90
  icarus.f
  isccp_cloud_types.f
  modis_simulator.F90
  scops.f
  GEOS_SatsimGridComp.F90
  )

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES GEOS_Shared MAPL GEOS_RadiationShared)

