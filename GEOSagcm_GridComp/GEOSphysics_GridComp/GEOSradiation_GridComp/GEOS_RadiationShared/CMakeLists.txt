esma_set_this ()

set (srcs
  radconstants.F90
  gettau.F90
  cloud_condensate_inhomogeneity.F90
  cloud_subcol_gen.F90
  )

esma_add_subdirectories(
  RRTMGP_cmake
  )

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES GEOS_Shared MAPL RRTMGP)
