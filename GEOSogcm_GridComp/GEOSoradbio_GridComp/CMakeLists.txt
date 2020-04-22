esma_set_this ()

set (SRCS
  slingo.F90
  sfcirr.F90
  setsfclte.F90
  setlte.F90
  ocalbedo.F90
  lidatatm.F90
  lidata.F90
  glight.F90
  edeu.F90
  daysetrad.F90
  clrtrans.F90
  aasack.F90
  GEOS_OradBioGridComp.F90
)


esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES MAPL
  INCLUDES ${INC_ESMF})

