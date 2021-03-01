# Note that directory and library name do not correspond!
esma_set_this (OVERRIDE GEOS_SurfaceShared)
esma_add_subdirectories (Raster)

set (resource_files
    GEOS_SurfaceGridComp.rc
   ) 

install( FILES ${resource_files} 
   DESTINATION etc
   )

set (srcs
  StieglitzSnow.F90
  SurfParams.F90)

esma_add_library (${this} SRCS ${srcs} DEPENDENCIES GEOS_Shared MAPL)
#target_compile_options (${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${FREAL4}>)
