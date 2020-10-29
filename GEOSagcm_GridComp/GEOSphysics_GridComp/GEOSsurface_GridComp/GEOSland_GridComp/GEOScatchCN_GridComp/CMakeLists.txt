esma_set_this ()
#
# cmake command line
# -DCLM4.0=ON -DCLM4.5=OFF
#
option(CLM4.0 "running clm4.0" ON)
option(CLM4.5 "running clm4.5" OFF)

if (CLM4.0) 
  if (CLM4.5)
    message( FATAL_ERROR "CLM 4.0 and 4.5 should be exclusive")
  endif()
endif()

if (CLM4.0)
  set (CATCHCN_COMPONENT CLM40)
elseif (CLM4.5)
  set (CATCHCN_COMPONENT CLM45)
else()
  message( FATAL_ERROR "One of CLM 4.0 and 4.5 should be ON")
endif

esma_add_subdirectories (${CATCHCN_COMPONENT})

set (srcs
  ${CATCHCN_COMPONENT}/GEOS_CatchCNGridComp.F90) 

esma_add_library (${this}
    SRCS ${srcs}
    DEPENDENCIES MAPL GEOS_Shared GEOS_LandShared GEOS_CatchCNShared)

