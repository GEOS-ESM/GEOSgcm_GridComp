# Note that directory and library name do not correspond!
esma_set_this (OVERRIDE GEOS_SurfaceShared)

set (resource_files
    GEOS_SurfaceGridComp.rc
   ) 

install( FILES ${resource_files} 
   DESTINATION etc
   )

set (srcs
  getids.F90
  StieglitzSnow.F90
  SurfParams.F90)

esma_add_library (${this} SRCS ${srcs} DEPENDENCIES GEOS_Shared MAPL
TYPE SHARED)
#target_compile_options (${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${FREAL4}>)
