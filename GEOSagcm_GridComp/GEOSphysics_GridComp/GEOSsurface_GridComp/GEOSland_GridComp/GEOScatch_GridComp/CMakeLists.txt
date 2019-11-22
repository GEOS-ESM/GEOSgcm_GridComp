esma_set_this ()
string (REPLACE GEOScatch_GridComp_ "" is_openmp ${this})


set (srcs
  GEOS_CatchGridComp.F90
  catchment.F90
  catch_incr.F90  
  )

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES MAPL_Base GEOS_LandShared)
target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})
#target_compile_definitions (${this} PRIVATE DBG_CATCH_INPUTS)
ecbuild_add_executable (TARGET dbg_clsm_offline.x SOURCES dbg_clsm_offline.F90  LIBS MAPL_Base GEOS_LandShared ${this})

add_subdirectory(mk_restarts)

if (is_openmp)
  target_compile_options(${this} PRIVATE ${OpenMP_Fortran_FLAGS})
endif ()

