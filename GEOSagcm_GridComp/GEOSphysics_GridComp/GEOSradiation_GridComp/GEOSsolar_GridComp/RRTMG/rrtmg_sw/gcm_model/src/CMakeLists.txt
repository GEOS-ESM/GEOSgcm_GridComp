# Non default name
esma_set_this (OVERRIDE RRTMG_SW)
set (srcs
  NRLSSI2.F90
  rrtmg_sw_init.F90
  rrtmg_sw_cldprmc.F90
  rrtmg_sw_setcoef.F90
  rrtmg_sw_taumol.F90
  rrtmg_sw_spcvmc.F90
  rrtmg_sw_rad.F90
  )

set (k_g_srcs
   rrtmg_sw_k_g_16.F90
   rrtmg_sw_k_g_17.F90
   rrtmg_sw_k_g_18.F90
   rrtmg_sw_k_g_19.F90
   rrtmg_sw_k_g_20.F90
   rrtmg_sw_k_g_21.F90
   rrtmg_sw_k_g_22.F90
   rrtmg_sw_k_g_23.F90
   rrtmg_sw_k_g_24.F90
   rrtmg_sw_k_g_25.F90
   rrtmg_sw_k_g_26.F90
   rrtmg_sw_k_g_27.F90
   rrtmg_sw_k_g_28.F90
   rrtmg_sw_k_g_29.F90
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

esma_add_library (${this} SRCS ${srcs} ${k_g_srcs}
   DEPENDENCIES GEOS_Shared MAPL RRTMG_SW_mods GEOS_RadiationShared)
target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})


