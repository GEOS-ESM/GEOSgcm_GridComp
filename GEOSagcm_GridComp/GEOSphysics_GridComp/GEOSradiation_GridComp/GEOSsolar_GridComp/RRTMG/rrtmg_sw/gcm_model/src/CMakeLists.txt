# Non default name
esma_set_this (OVERRIDE RRTMG_SW)
set (srcs
  mcica_random_numbers.F90 mcica_subcol_gen_sw.F90 rrtmg_sw_cldprmc.F90
  rrtmg_sw_cldprop.F90 rrtmg_sw_init.F90 rrtmg_sw_k_g.F90
  rrtmg_sw_rad.F90 rrtmg_sw_setcoef.F90 rrtmg_sw_spcvmc.F90
  rrtmg_sw_taumol.F90 tab_xcw.F90
  )

# MAT With GCC -O1 -g with these files causes an issue. But disabling var tracking
#     seems to make compilation bearable on release
if (CMAKE_BUILD_TYPE MATCHES Release)
   if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
      set_source_files_properties( rrtmg_sw_k_g.F90 PROPERTIES COMPILE_FLAGS "-O1 -fno-var-tracking")
   else ()
      set_source_files_properties( rrtmg_sw_k_g.F90 PROPERTIES COMPILE_FLAGS "-O1")
   endif ()
endif ()

esma_add_library (${this} SRCS ${srcs}
   DEPENDENCIES GEOS_Shared MAPL RRTMG_SW_mods)
target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})


