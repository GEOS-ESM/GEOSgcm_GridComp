esma_set_this()


set (dependencies MAPL esmf)

esma_add_library (${this}
  SRCS GigaTraj_GridComp.F90
  DEPENDENCIES ${dependencies})

add_subdirectory (@GigaTraj)
