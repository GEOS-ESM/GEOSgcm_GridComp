esma_set_this ()

set (srcs
  getco2.F90
  radconstants.F90
  gettau.F90
  )

esma_add_library (${this}
  SRCS ${srcs}
  SUBDIRS RRTMGP
  DEPENDENCIES GEOS_Shared MAPL_Base)
