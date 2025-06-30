esma_set_this()

set (SRCS
  GEOS_OceanbiogeochemGridComp.F90
  carbon.F90
  co2calc_SWS.F90
  daysetbio.F90
  drtsafe.F90
  kloop.F90
  ppco2.F90
  ptend.F90
  setbio.F90
  sink.F90
  ta_iter_SWS.F90
  )

esma_add_library(${this}
  SRCS ${SRCS}
  DEPENDENCIES MAPL ESMF::ESMF)






