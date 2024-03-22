esma_set_this (OVERRIDE umwm)

list (APPEND UMWM_SRCS
   src/umwm_mpi.F90
   src/umwm_io.F90
   src/umwm_util.F90
   src/umwm_module.F90
   src/umwm_constants.F90
   src/umwm_init.F90
   src/umwm_sheltering.F90
   src/umwm_source_functions.F90
   src/umwm_physics.F90
   src/umwm_advection.F90
   src/umwm_stress.F90
   src/umwm_stokes.F90
)

esma_mepo_style(umwm UMWM_rel_path REL_PATH ..)
set (UMWM_path ${CMAKE_CURRENT_SOURCE_DIR}/${UMWM_rel_path})

set (SRCS)
foreach (file ${UMWM_SRCS})
   list (APPEND SRCS ${UMWM_path}/${file})
endforeach ()

esma_add_library (${this}
   SRCS ${SRCS}
   DEPENDENCIES esmf NetCDF::NetCDF_Fortran
)

target_compile_definitions(${this} PRIVATE MPI ESMF GEOS)

if (CMAKE_Fortran_COMPILER_ID MATCHES Intel AND CMAKE_BUILD_TYPE MATCHES Release)
  target_compile_options (${this} PRIVATE "-ftrapuv")
endif ()

