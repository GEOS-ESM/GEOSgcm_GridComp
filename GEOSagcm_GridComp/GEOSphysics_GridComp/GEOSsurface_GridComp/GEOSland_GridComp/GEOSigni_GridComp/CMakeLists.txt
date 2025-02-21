esma_set_this ()


set (resource_files
    GEOS_IgniGridComp.rc
   )

install( FILES ${resource_files}
   DESTINATION etc
   )

set (srcs
  GEOS_IgniGridComp.F90
  cffwi.F90
  )


esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES MAPL GEOS_Shared ESMF::ESMF
  )

