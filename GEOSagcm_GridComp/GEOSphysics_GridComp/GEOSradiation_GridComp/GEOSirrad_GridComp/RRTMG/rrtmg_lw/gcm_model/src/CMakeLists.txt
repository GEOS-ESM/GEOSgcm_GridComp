esma_set_this(OVERRIDE RRTMG)

set (srcs
  rrtmg_lw_init.F90
  rrtmg_lw_cldprmc.F90
  rrtmg_lw_setcoef.F90
  rrtmg_lw_taumol.F90
  rrtmg_lw_rtrnmc.F90
  rrtmg_lw_rad.F90
  )

set (k_g_srcs
   rrtmg_lw_k_g_01.F90
   rrtmg_lw_k_g_02.F90
   rrtmg_lw_k_g_03.F90
   rrtmg_lw_k_g_04.F90
   rrtmg_lw_k_g_05.F90
   rrtmg_lw_k_g_06.F90
   rrtmg_lw_k_g_07.F90
   rrtmg_lw_k_g_08.F90
   rrtmg_lw_k_g_09.F90
   rrtmg_lw_k_g_10.F90
   rrtmg_lw_k_g_11.F90
   rrtmg_lw_k_g_12.F90
   rrtmg_lw_k_g_13.F90
   rrtmg_lw_k_g_14.F90
   rrtmg_lw_k_g_15.F90
   rrtmg_lw_k_g_16.F90
   )

# MAT With GCC -O1 -g with these files causes an issue. But disabling var tracking
#     seems to make compilation bearable on release
foreach (file ${k_g_srcs})
   if (CMAKE_BUILD_TYPE MATCHES Release OR CMAKE_BUILD_TYPE MATCHES Aggressive)
      if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
         set_source_files_properties( ${file} PROPERTIES COMPILE_FLAGS "-O1 -fno-var-tracking")
      else ()
         set_source_files_properties( ${file} PROPERTIES COMPILE_FLAGS "-O1")
      endif ()
   endif ()
endforeach ()

esma_add_library (${this}
  SRCS ${srcs} ${k_g_srcs}
  DEPENDENCIES GEOS_Shared MAPL RRTMG_mods GEOS_RadiationShared)
