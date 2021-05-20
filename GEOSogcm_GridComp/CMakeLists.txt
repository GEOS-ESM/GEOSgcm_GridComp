esma_set_this ()


set (alldirs
  GEOS_OceanBioGeoChemGridComp
  GEOS_OradGridComp
  GEOS_OradBioGridComp
  )

set (dependencies GEOSdatasea_GridComp  GEOSdataseaice_GridComp GuestOcean_GridComp GEOSCICEDyna_GridComp MAPL CICE4)

esma_add_library (${this}
  SRCS GEOS_OgcmGridComp.F90
  SUBCOMPONENTS ${alldirs}
  SUBDIRS GEOSseaice_GridComp GEOSocean_GridComp 
  DEPENDENCIES ${dependencies}
  INCLUDES ${INC_ESMF})

