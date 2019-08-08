esma_set_this ()

set (srcs
  sw.f90 jw.f90 testcases_3_4_5_6_stand_alone.f90
  GetWeightsC2C.F90
  GetWeights.F90
  CubeHalo.F90
  Cube2LatLon.F90 LatLon2Cube.F90 AppGridCreate.F90 FV_StateMod.F90
  AdvCore_GridCompMod.F90
  DynCore_GridCompMod.F90 CreateInterpWeights_GridCompMod.F90
  StandAlone_DynAdvCore_GridCompMod.F90
  CubeToLatLonRegridder.F90 
  LatLonToCubeRegridder.F90 
  CubeToCubeRegridder.F90
  CubeToLatLon.F90 
  CubeGridPrototype.F90
  )

set(dependencies MAPL_Base GMAO_hermes GEOS_Shared)

if (FV_PRECISION STREQUAL R4)
   set (GFDL GFDL_fms_r4)
elseif (FV_PRECISION STREQUAL R4R8) # FV is R4 but FMS is R8
   get_target_property (extra_incs GFDL_fms_r4 INCLUDE_DIRECTORIES)
   set (GFDL GFDL_fms_r8)
elseif (FV_PRECISION STREQUAL R8)
   set (GFDL GFDL_fms_r8)
endif ()

esma_add_library (${this}
  SRCS ${srcs}
  SUBCOMPONENTS fvdycore
  DEPENDENCIES ${dependencies}
  DEPENDENCIES ${GFDL}
  INCLUDES ${INC_ESMF} ${extra_incs})

if (FV_PRECISION MATCHES R4)
  target_compile_definitions (${this} PRIVATE -DSINGLE_FV -DOVERLOAD_R4)
elseif (FV_PRECISION MATCHES R4R8) # FV is R4 but FMS is R8
  target_compile_definitions (${this} PRIVATE -DSINGLE_FV -DOVERLOAD_R4)
endif ()


set (CMAKE_Fortran_FLAGS_RELEASE "${GEOS_Fortran_FLAGS_VECT}")

if (CRAY_POINTER)
  set_target_properties (${this} PROPERTIES COMPILE_FLAGS ${CRAY_POINTER})
endif()

add_definitions (-DSPMD -DMAPL_MODE)

foreach(flag ${tmp})
   target_compile_options (${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${flag}>)
endforeach()

ecbuild_add_executable (
  TARGET StandAlone_FV3_Dycore.x
  SOURCES StandAlone_FV3_Dycore.F90
  LIBS ${this})
set_target_properties (StandAlone_FV3_Dycore.x PROPERTIES LINK_FLAGS "${OpenMP_Fortran_FLAGS}")

ecbuild_add_executable (
  TARGET rs_scale.x
  SOURCES rs_scale.F90
  LIBS ${this} GMAO_pFIO)
set_target_properties (rs_scale.x PROPERTIES LINK_FLAGS "${OpenMP_Fortran_FLAGS}")

ecbuild_add_executable (
   TARGET StandAlone_AdvCore.x
   SOURCES StandAlone_AdvCore.F90
   LIBS ${this})
set_target_properties (StandAlone_AdvCore.x PROPERTIES LINK_FLAGS "${OpenMP_Fortran_FLAGS}")

ecbuild_add_executable (
   TARGET StandAlone_DynAdvCore.x
   SOURCES StandAlone_DynAdvCore.F90
   LIBS ${this})
set_target_properties (StandAlone_DynAdvCore.x PROPERTIES LINK_FLAGS "${OpenMP_Fortran_FLAGS}")

ecbuild_add_executable (
   TARGET c2c.x
   SOURCES c2c.F90
   LIBS ${this})
set_target_properties (c2c.x PROPERTIES LINK_FLAGS "${OpenMP_Fortran_FLAGS}")

ecbuild_add_executable (
   TARGET interp_restarts.x
   SOURCES interp_restarts.F90
   LIBS ${this})
set_target_properties (interp_restarts.x PROPERTIES LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
