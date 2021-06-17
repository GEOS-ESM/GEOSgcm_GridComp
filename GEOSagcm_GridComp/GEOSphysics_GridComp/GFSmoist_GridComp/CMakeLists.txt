esma_set_this()

set (srcs
  GFS_MoistGridComp.F90
  @physics/physics/machine.F
  @physics/physics/gfdl_cloud_microphys.F90
  @physics/physics/module_gfdl_cloud_microphys.F90
  @physics/physics/GFS_suite_interstitial.F90
  @physics/physics/GFS_DCNV_generic.F90
  @physics/physics/samfdeepcnv.f
  @physics/physics/samfshalcnv.f
  @physics/physics/GFS_SCNV_generic.F90
  @physics/physics/cnvc90.f
  @physics/physics/GFS_MP_generic.F90
  @physics/physics/maximum_hourly_diagnostics.F90
  @physics/physics/phys_tend.F90
  @physics/physics/calpreciptype.f90
  @physics/physics/module_mp_thompson_make_number_concentrations.F90
  @physics/physics/physcons.F90
  @physics/physics/module_mp_radar.F90
  @physics/physics/samfaerosols.F
  @physics/physics/funcphys.f90
  @physics/physics/get_prs_fv3.F90
  )

#if (CMAKE_Fortran_COMPILER_ID MATCHES Intel AND CMAKE_BUILD_TYPE MATCHES Aggressive)
#   set (CMAKE_Fortran_FLAGS_AGGRESSIVE "${GEOS_Fortran_FLAGS_VECT}")
#endif ()

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES GEOS_Shared GMAO_mpeu MAPL Chem_Shared Chem_Base
  INCLUDES ${INC_ESMF})

get_target_property (extra_incs fms_r4 INCLUDE_DIRECTORIES)
target_include_directories(${this} PRIVATE
   $<BUILD_INTERFACE:${extra_incs}>
   )

file (GLOB_RECURSE rc_files CONFIGURE_DEPENDS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.rc)
foreach ( file ${rc_files} )
   get_filename_component( dir ${file} DIRECTORY )
   install( FILES ${file} DESTINATION etc/${dir} )
endforeach()

target_compile_definitions (${this} PRIVATE -DSINGLE_PREC)
