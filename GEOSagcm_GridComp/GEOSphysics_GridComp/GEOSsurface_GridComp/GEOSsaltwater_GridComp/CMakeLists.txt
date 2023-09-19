esma_set_this()

set (srcs
  GEOS_CICE4ColumnPhysGridComp.F90 GEOS_OpenWaterGridComp.F90
  GEOS_SaltWaterGridComp.F90 GEOS_SimpleSeaiceGridComp.F90
  GEOS_SeaiceInterfaceGridComp.F90
  GEOS_ObioImportsGridComp.F90)

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES GEOS_Shared MAPL CICE4 esmf NetCDF::NetCDF_Fortran)

