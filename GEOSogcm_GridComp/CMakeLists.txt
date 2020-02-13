esma_set_this ()


set (alldirs
  GEOSoceanbiogeochem_GridComp
  GEOSoceanbiosimple_GridComp
  GEOSorad_GridComp
  GEOSoradbio_GridComp
  )

set (rename_GEOSoceanbiogeochem_GridComp GEOS_OceanBioGeoChemGridComp)
set (rename_GEOSoceanbiosimple_GridComp GEOS_OceanBioGridComp)
set (rename_GEOSoradbio_GridComp GEOS_OradBioGridComp)

set (dependencies GEOSdatasea_GridComp  GEOSdataseaice_GridComp GuestOcean_GridComp GEOSCICEDyna_GridComp MAPL LANL_cice)

esma_add_library (${this}
  SRCS GEOS_OgcmGridComp.F90
  SUBCOMPONENTS ${alldirs}
  SUBDIRS GEOSseaice_GridComp GEOSocean_GridComp 
  DEPENDENCIES ${dependencies}
  INCLUDES ${INC_ESMF})

