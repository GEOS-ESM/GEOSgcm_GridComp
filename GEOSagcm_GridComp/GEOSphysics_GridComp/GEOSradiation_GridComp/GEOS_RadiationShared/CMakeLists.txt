esma_set_this ()

set (srcs
  radconstants.F90
  gettau.F90
  )

esma_add_library (${this}
  SRCS ${srcs}
  SUBDIRS RRTMGP
  DEPENDENCIES GEOS_Shared MAPL_Base)
