esma_set_this ()

set (srcs
  shr_kind_mod.F90 FVperf_module.F90 dynamics_vars.F90 diag_module.F90 FVdycore_wrapper.F90
  FVdycore_GridCompMod.F90 fill_module.F90 pft_module.F90 mapz_module.F90 fft99.F90
  tp_core.F90 par_vecsum.F90 par_xsum.F90 pmaxmin.F90 sw_core.F90
  benergy.F90 cd_core.F90 geopk.F90 pkez.F90 te_map.F90 trac2d.F90
  epvd.F90 mfz_comp.F90 remap.F90 gmap.F90 glosum.F90 BlendingMod.F90 G3_AVRX.F90
  )

add_definitions(-DSPMD -DMAPL_MODE -DTWOD_YZ -DrFV=r4)


esma_add_library (${this} SRCS ${srcs} DEPENDENCIES GEOS_Shared GMAO_pilgrim MAPL)
target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})

if (CRAY_POINTER)
  set_target_properties (${this} PROPERTIES COMPILE_FLAGS ${CRAY_POINTER})
endif()

file (GLOB rc_files *.rc)
install(FILES ${rc_files} DESTINATION etc)
