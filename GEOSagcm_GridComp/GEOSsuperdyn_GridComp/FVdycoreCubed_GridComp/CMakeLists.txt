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
  CubedSphereGridFactory.F90
  CubeToLatLonRegridder.F90 
  LatLonToCubeRegridder.F90 
  CubeToCubeRegridder.F90
  CubeToLatLon.F90
#***DAS***#  CubeGridPrototype.F90
  )

set(dependencies MAPL_Base GMAO_hermes GEOS_Shared)

if (FV_PRECISION MATCHES R4)
  get_target_property (extra_incs GFDL_fms_r4 INCLUDE_DIRECTORIES)
else ()
  get_target_property (extra_incs GFDL_fms_r8 INCLUDE_DIRECTORIES)
endif ()


esma_add_library (${this}
  SRCS ${srcs}
  SUBCOMPONENTS fvdycore
  DEPENDENCIES ${dependencies}
  INCLUDES ${INC_ESMF} ${extra_incs})

if (FV_PRECISION MATCHES R4)
   target_compile_definitions (${this} PRIVATE -DSINGLE_FV -DOVERLOAD_R4)
#  list(APPEND dependencies GFDL_fms_r4)
else ()
#  list(APPEND dependencies GFDL_fms_r8)

#  string(REPLACE " " ";" tmp ${FREAL8})
#  foreach(flag ${tmp})
#    target_compile_options (${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${flag}>)
#  endforeach()
endif ()


set (CMAKE_Fortran_FLAGS_RELEASE "${GEOS_Fortran_FLAGS_VECT}")

if (CRAY_POINTER)
  set_target_properties (${this} PROPERTIES COMPILE_FLAGS ${CRAY_POINTER})
endif()

add_definitions (-DSPMD -DMAPL_MODE)

foreach(flag ${tmp})
   target_compile_options (${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${flag}>)
endforeach()

add_executable (StandAlone_FV3_Dycore.x StandAlone_FV3_Dycore.F90)
set_target_properties (StandAlone_FV3_Dycore.x PROPERTIES LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
target_link_libraries(StandAlone_FV3_Dycore.x ${this})

add_executable (rs_scale.x rs_scale.F90)
set_target_properties (rs_scale.x PROPERTIES LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
target_link_libraries(rs_scale.x ${this} GMAO_pFIO)


