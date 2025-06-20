esma_set_this ()

set (SRCS
  setlte.F90
  sfcrfl.F90
  ocalbedo.F90
  lidata.F90
  glight.F90
  edeu.F90
  daysetrad.F90
  aasack.F90
  GEOS_OradBioGridComp.F90
)


esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES MAPL ESMF::ESMF)

