esma_set_this()

set (srcs
  GEOS_OpenWaterGridComp.F90
  GEOS_SaltWaterGridComp.F90 
  GEOS_SimpleSeaiceGridComp.F90
  GEOS_ObioGridComp.F90)

esma_add_library (${this}
  SRCS ${srcs}
  SUBCOMPONENTS GEOS_CICE4ColumnPhysGridComp GEOS_SeaiceInterfaceGridComp 
  DEPENDENCIES GEOS_Shared MAPL esmf NetCDF::NetCDF_Fortran)

