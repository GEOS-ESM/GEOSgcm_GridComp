esma_set_this()

set (srcs
  held_suarez_oacc.F90
  GEOS_HsGridComp.F90
  )

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES GEOS_Shared MAPL)
target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})

target_compile_options(${this} PRIVATE -fopenacc -foffload=nvptx-none)
