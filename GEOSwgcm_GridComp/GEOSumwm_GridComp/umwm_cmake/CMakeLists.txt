esma_set_this (OVERRIDE umwm)

add_compile_definitions(MPI ESMF GEOS)

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

set (UMWM_path ${CMAKE_CURRENT_SOURCE_DIR}/../@umwm)
set (SRCS)
foreach (file ${UMWM_SRCS})
   list (APPEND SRCS ${UMWM_path}/${file})
endforeach ()

include_directories (${INC_ESMF} ${INC_NETCDF})

esma_add_library (${this} 
   SRCS ${SRCS}
#  DEPENDENCIES MAPL
)

target_include_directories (${this} 
   PUBLIC ${INC_ESMF} ${INC_NETCDF}
)

target_compile_options (${this} PRIVATE "-ftrapuv")
#target_link_libraries (${this} PRIVATE ${NETCDF_LIBRARIES}  ${MPI_Fortran_LIBRARIES} MPI:MPI_LIBS)

