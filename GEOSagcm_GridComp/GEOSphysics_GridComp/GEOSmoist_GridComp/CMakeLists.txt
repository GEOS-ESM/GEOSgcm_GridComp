esma_set_this()

set (srcs
  CLDPARAMS.F90 RASPARAMS.F90 ras.F90 cloudnew.F90 ddf.F90 rascnvv2_v.F90
  cldmacro.F90 cldwat2m_micro.F90 wv_saturation.F90 aer_cloud.F90
  GEOS_MoistGridComp.F90 partition_pdf.F90
  ConvPar_GF_GEOS5.F90 module_gate.F90 uwshcu.F90 SHLWPARAMS.F90
  aer_actv_single_moment.F90 gfdl_cloud_microphys.F90
  micro_mg3_0.F90 micro_mg_utils.F90 wv_sat_methods.F90
  # files that are sometimes present?
  # ras00.F90 cloud.F90
  )

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES GEOS_Shared GMAO_mpeu MAPL Chem_Shared Chem_Base
  INCLUDES ${INC_ESMF})

get_target_property (extra_incs fms_r4 INCLUDE_DIRECTORIES)
target_include_directories(${this} PRIVATE
   $<BUILD_INTERFACE:${extra_incs}>
   )

file (GLOB_RECURSE rc_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.rc)
foreach ( file ${rc_files} )
   get_filename_component( dir ${file} DIRECTORY )
   install( FILES ${file} DESTINATION etc/${dir} )
endforeach()
