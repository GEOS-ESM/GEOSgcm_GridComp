esma_set_this ()

esma_add_library (${this}
   SRCS GEOS_MITDynaGridComp.F90
   DEPENDENCIES MAPL CICE4 ESMF::ESMF NetCDF::NetCDF_Fortran
)



