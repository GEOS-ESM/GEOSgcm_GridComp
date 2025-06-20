esma_set_this ()

set (srcs
  g3_dynamics_lattice_module.F g3_dynamics_state_module.F
  g3_mpi_util.F
  g3_mymalloc.F g3_fft.F
  g3_wrapper.F g3_dycore.F g3_gcmutil.F g3_grids.F
  shr_kind_mod.F90 dynamics_vars.F90
  pft_module.F90 par_xsum.F90 fft99.F90 remap.F90 gmap.F90
  ${this}Mod.F90
  )

set (files_that_need_extended_source
  g3_dycore.F g3_dynamics_lattice_module.F
  g3_dynamics_state_module.F g3_gcmutil.F
  g3_grids.F
  g3_mpi_util.F
  g3_mymalloc.F g3_fft.F
  g3_wrapper.F
  )
if (EXTENDED_SOURCE)
  foreach (file ${files_that_need_extended_source})
    set_source_files_properties(${file} APPEND_STRING PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
    set_source_files_properties(${file} APPEND_STRING PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
  endforeach ()
endif ()

add_definitions(-DSPMD -DMAPL_MODE -DTWOD_YZ -DrFV=r4)


esma_add_library (${this} SRCS ${srcs}
   DEPENDENCIES GEOS_Shared GMAO_pilgrim MAPL GMAO_transf)
target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})
if (CRAY_POINTER)
  set_target_properties (${this} PROPERTIES COMPILE_FLAGS ${CRAY_POINTER})
endif()

