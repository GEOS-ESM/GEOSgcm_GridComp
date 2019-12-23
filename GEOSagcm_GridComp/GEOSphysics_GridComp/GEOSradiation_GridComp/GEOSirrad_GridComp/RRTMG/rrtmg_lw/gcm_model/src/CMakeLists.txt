esma_set_this(OVERRIDE RRTMG)

set (srcs
  mcica_random_numbers.F90 rrtmg_lw_init.F90
  rrtmg_lw_k_g.F90 rrtmg_lw_rad.F90
  rrtmg_lw_rtrn.F90 rrtmg_lw_rtrnmr.F90
  rrtmg_lw_setcoef.F90 WaterDistributionMod.F90
  gmcica_subcol_gen_lw.F90 gpu_rrtmg_lw_rtrnmc.F90
  gpu_rrtmg_lw_setcoef.F90 gpu_rrtmg_lw_taumol.F90
  gpu_rrtmg_lw_cldprmc.F90
  )

# MAT With GCC -O1 -g with these files causes an issue. But disabling var tracking
#     seems to make compilation bearable on release
if (CMAKE_BUILD_TYPE MATCHES Release)
   if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
      set_source_files_properties( rrtmg_lw_k_g.F90 PROPERTIES COMPILE_FLAGS "-O1 -fno-var-tracking")
   else ()
      set_source_files_properties( rrtmg_lw_k_g.F90 PROPERTIES COMPILE_FLAGS "-O1")
   endif ()
endif ()

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES GEOS_Shared MAPL_Base RRTMG_mods)
