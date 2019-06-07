# Note that directory and library name do not correspond!
esma_set_this (OVERRIDE GEOS_SurfaceShared)

set (srcs
  StieglitzSnow.F90
  SurfParams.F90)

esma_add_library (${this} SRCS ${srcs} DEPENDENCIES GEOS_Shared MAPL_Base)
#target_compile_options (${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${FREAL4}>)

