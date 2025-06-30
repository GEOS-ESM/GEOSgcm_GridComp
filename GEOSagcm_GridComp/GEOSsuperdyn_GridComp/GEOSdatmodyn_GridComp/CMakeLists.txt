esma_set_this ()

set (srcs
  ppm.F90 cptread.F90 reader.F90 gate_reader.F90 bomex_reader.F90
  stratus_ic.F90 cfmip_ic.F90 NeuralNet.F90 GEOS_DatmoDynGridComp.F90
  )

esma_add_library (${this} SRCS ${srcs} DEPENDENCIES GEOS_Shared MAPL ESMF::ESMF NetCDF::NetCDF_Fortran)

