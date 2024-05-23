esma_set_this()

option(BUILD_SERIALBOX_SER "Build with Serialbox serialization statements" OFF)

set (srcs
  GEOS_RAS_InterfaceMod.F90 RASPARAMS.F90 ras.F90 ddf.F90 rascnvv2_v.F90
  GEOS_BACM_1M_InterfaceMod.F90 cloudnew.F90
  GEOS_MGB2_2M_InterfaceMod.F90 micro_mg3_0.F90 micro_mg_utils.F90
  cldwat2m_micro.F90 wv_saturation.F90 aer_cloud.F90
  wv_sat_methods.F90
  GEOS_GFDL_1M_InterfaceMod.F90 gfdl_cloud_microphys.F90
  GEOS_THOM_1M_InterfaceMod.F90 module_mp_thompson.F90 module_mp_radar.F90 machine.F
  GEOS_GF_InterfaceMod.F90 ConvPar_GF_GEOS5.F90 ConvPar_GF2020.F90 ConvPar_GF_Shared.F90 module_gate.F90
  GEOS_UW_InterfaceMod.F90 uwshcu.F90
  # aer_actv_single_moment.F90
  # GEOS_MoistGridComp.F90
  # files that are sometimes present?
  # ras00.F90 cloud.F90
  )
if (GPU_PHYSICS)
  list (APPEND srcs
    openacc/Process_Library.F90
    openacc/gfdl_cloud_microphys.F90)
else ()
  list (APPEND srcs
    Process_Library.F90
    gfdl_cloud_microphys.F90)
endif ()

if (BUILD_SERIALBOX_SER)
  list (APPEND srcs
    sb/aer_actv_single_moment.F90
    sb/GEOS_MoistGridComp.F90)
else ()
  list (APPEND srcs
    aer_actv_single_moment.F90
    GEOS_MoistGridComp.F90
    )
endif()

if (CMAKE_Fortran_COMPILER_ID MATCHES Intel AND CMAKE_BUILD_TYPE MATCHES Aggressive)
   set (CMAKE_Fortran_FLAGS_AGGRESSIVE "${GEOS_Fortran_FLAGS_VECT}")
endif ()

if (CMAKE_Fortran_COMPILER_ID MATCHES GNU AND CMAKE_BUILD_TYPE MATCHES Release)
   string (REPLACE "${FOPT3}" "${FOPT2}" CMAKE_Fortran_FLAGS_RELEASE ${CMAKE_Fortran_FLAGS_RELEASE})
endif ()

# Note For unknown reasons, BACM_1M_Interface takes 20 minutes to compile at O3
#      and 10 minutes at O2. But only 7 seconds with O1. So we compile at O1
if (CMAKE_Fortran_COMPILER_ID MATCHES Intel AND CMAKE_BUILD_TYPE MATCHES Release)
  set_source_files_properties(GEOS_BACM_1M_InterfaceMod.F90 PROPERTIES COMPILE_OPTIONS ${FOPT1})
  set_source_files_properties(GEOS_MGB2_2M_InterfaceMod.F90 PROPERTIES COMPILE_OPTIONS ${FOPT1})
endif ()

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES GEOS_Shared GMAO_mpeu MAPL Chem_Shared Chem_Base esmf)

if(BUILD_SERIALBOX_SER)
  message(STATUS "Building with Serialbox serialization statements")
  message(STATUS "${SERIALBOX_ROOT}/lib/")
  add_definitions(-DSERIALIZE)
  target_include_directories(${this} PRIVATE ${SERIALBOX_ROOT}/include)
  find_library(SERIALBOX_LIB
      NAMES SerialboxFortran
      HINTS ${SERIALBOX_ROOT}/lib/
      REQUIRED
  )
  target_link_libraries(${this} PRIVATE ${SERIALBOX_ROOT}/lib/libSerialboxFortran.a ${SERIALBOX_ROOT}/lib/libSerialboxC.a ${SERIALBOX_ROOT}/lib/libSerialboxCore.a)
endif()

get_target_property (extra_incs fms_r4 INCLUDE_DIRECTORIES)
target_include_directories(${this} PRIVATE
   $<BUILD_INTERFACE:${extra_incs}>
   )
if (GPU_PHYSICS)
  message(STATUS "Building GPU ports of Moist routines")
  set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${GEOS_Fortran_OpenACC_Flags}")
  set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${GEOS_Fortran_OpenACC_Flags}")
endif ()

file (GLOB_RECURSE rc_files CONFIGURE_DEPENDS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.rc *.yaml)
foreach ( file ${rc_files} )
   get_filename_component( dir ${file} DIRECTORY )
   install( FILES ${file} DESTINATION etc/${dir} )
endforeach()



