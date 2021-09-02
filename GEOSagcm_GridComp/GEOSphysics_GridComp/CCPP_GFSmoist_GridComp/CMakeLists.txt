esma_set_this()

set (srcs
  CCPP_GFS_MoistGridComp.F90
  GFS_typedefs_stubs.F90
  @ccpp-physics/physics/machine.F
  @ccpp-physics/physics/gfdl_cloud_microphys.F90
  @ccpp-physics/physics/module_gfdl_cloud_microphys.F90
  @ccpp-physics/physics/GFS_suite_interstitial.F90
  @ccpp-physics/physics/GFS_DCNV_generic.F90
  @ccpp-physics/physics/samfdeepcnv.f
  @ccpp-physics/physics/samfshalcnv.f
  @ccpp-physics/physics/GFS_SCNV_generic.F90
  @ccpp-physics/physics/cnvc90.f
  @ccpp-physics/physics/GFS_MP_generic.F90
  @ccpp-physics/physics/maximum_hourly_diagnostics.F90
  @ccpp-physics/physics/phys_tend.F90
  @ccpp-physics/physics/calpreciptype.f90
  @ccpp-physics/physics/module_mp_thompson_make_number_concentrations.F90
  @ccpp-physics/physics/physcons.F90
  @ccpp-physics/physics/module_mp_radar.F90
  @ccpp-physics/physics/samfaerosols.F
  @ccpp-physics/physics/funcphys.f90
  @ccpp-physics/physics/get_prs_fv3.F90
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
