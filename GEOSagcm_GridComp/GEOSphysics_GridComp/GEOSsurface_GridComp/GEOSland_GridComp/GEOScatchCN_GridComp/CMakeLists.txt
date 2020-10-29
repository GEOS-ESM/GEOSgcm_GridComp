esma_set_this ()

#
# cmake command line
# -DCLM4.0=ON -DCLM4.5=OFF
option(CLM4.0 "running clm4.0" ON)
option(CLM4.5 "running clm4.5" OFF)

if (CLM4.0) 
  if (CLM4.5)
    message( FATAL_ERROR "CLM 4.0 and 4.5 should be exclusive")
  endif()
endif()

if (CLM4.0)
  set (srcs
  GEOS_CatchCNclm40GridComp.F90) 
  esma_add_subdirectories (CLM40)
  esma_add_library (${this}
    SRCS ${srcs}
    DEPENDENCIES MAPL GEOS_Shared GEOS_LandShared GEOS_CatchCNShared)
endif()

if (CLM4.5)
  set (srcs
  GEOS_CatchCNclm45GridComp.F90) 
  esma_add_subdirectories (CLM45)

  esma_add_library (${this}
    SRCS ${srcs}
    DEPENDENCIES MAPL GEOS_Shared GEOS_LandShared CLM45 GEOS_CatchCNShared)
endif()

