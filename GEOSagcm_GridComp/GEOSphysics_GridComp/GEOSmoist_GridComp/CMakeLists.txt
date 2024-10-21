esma_set_this()

set (srcs
  Process_Library.F90
  GEOS_RAS_InterfaceMod.F90 RASPARAMS.F90 ras.F90 ddf.F90 rascnvv2_v.F90
  GEOS_BACM_1M_InterfaceMod.F90 cloudnew.F90
  GEOS_MGB2_2M_InterfaceMod.F90 micro_mg3_0.F90 micro_mg_utils.F90
  cldwat2m_micro.F90 wv_saturation.F90 aer_cloud.F90
  wv_sat_methods.F90
  GEOS_GFDL_1M_InterfaceMod.F90 gfdl_cloud_microphys.F90
  GEOS_THOM_1M_InterfaceMod.F90 module_mp_thompson.F90 module_mp_radar.F90 machine.F
  GEOS_GF_InterfaceMod.F90 ConvPar_GF_GEOS5.F90 ConvPar_GF2020.F90 ConvPar_GF_Shared.F90 module_gate.F90
  GEOS_UW_InterfaceMod.F90 uwshcu.F90
  aer_actv_single_moment.F90
  GEOS_MoistGridComp.F90
  # files that are sometimes present?
  # ras00.F90 cloud.F90
  )

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
  DEPENDENCIES GEOS_Shared GMAO_mpeu MAPL Chem_Shared Chem_Base ESMF::ESMF)

get_target_property (extra_incs fms_r4 INCLUDE_DIRECTORIES)
target_include_directories(${this} PRIVATE
   $<BUILD_INTERFACE:${extra_incs}>
   )

file (GLOB_RECURSE rc_files CONFIGURE_DEPENDS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.rc *.yaml)
foreach ( file ${rc_files} )
   get_filename_component( dir ${file} DIRECTORY )
   install( FILES ${file} DESTINATION etc/${dir} )
endforeach()
