esma_set_this()

option(BUILD_SERIALBOX_SER "Build with Serialbox serialization statements" OFF)
option(BUILD_PYMOIST_INTERFACE "Build pyMoist interface" OFF)

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
  aer_actv_single_moment.F90
  GEOS_MoistGridComp.F90
  # files that are sometimes present?
  # ras00.F90 cloud.F90
  )

if (BUILD_PYMOIST_INTERFACE)
  list (APPEND srcs
  pyMoist/pyMoist/interface/cffi_lib/interface.f90
  pyMoist/pyMoist/interface/cffi_lib/interface.c)
endif ()

if (BUILD_PYMOIST_INTERFACE)

  message(STATUS "Building pyMoist interface")

  add_definitions(-DPYMOIST_INTEGRATION)
  
  # The Python library creation requires mpiexec/mpirun to run on a
  # compute node. Probably a weird SLURM thing?
  find_package(MPI REQUIRED)
  find_package(Python3 COMPONENTS Interpreter REQUIRED)

  # Set up some variables in case names change
  set(PYMOIST_INTERFACE_LIBRARY     ${CMAKE_CURRENT_BINARY_DIR}/libpyMoist_interface_py.so)
  set(PYMOIST_INTERFACE_HEADER_FILE ${CMAKE_CURRENT_BINARY_DIR}/pymoist_interface_py.h)
  set(PYMOIST_INTERFACE_FLAG_HEADER_FILE ${CMAKE_CURRENT_SOURCE_DIR}/pyMoist/pyMoist/interface/cffi_lib/moist.h)
  set(PYMOIST_INTERFACE_SRCS        ${CMAKE_CURRENT_SOURCE_DIR}/pyMoist/pyMoist/interface/cffi_lib/interface.py)

  # This command creates the shared object library from Python
  add_custom_command(
    OUTPUT ${PYMOIST_INTERFACE_LIBRARY}
    # Note below is essentially:
    #  mpirun -np 1 python file
    # but we use the CMake options as much as we can for flexibility
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PYMOIST_INTERFACE_FLAG_HEADER_FILE} ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${Python3_EXECUTABLE} ${PYMOIST_INTERFACE_SRCS}
    BYPRODUCTS ${PYMOIST_INTERFACE_HEADER_FILE}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    MAIN_DEPENDENCY ${PYMOIST_INTERFACE_SRCS}
    COMMENT "Building pyMoist interface library with Python"
    VERBATIM
  )

  # This creates a target we can use for dependencies and post build
  add_custom_target(generate_pymoist_interface_library DEPENDS ${PYMOIST_INTERFACE_LIBRARY})

  # Because of the weird hacking of INTERFACE libraries below, we cannot
  # use the "usual" CMake calls to install() the .so. I think it's because
  # INTERFACE libraries don't actually produce any artifacts as far as
  # CMake is concerned. So we add a POST_BUILD custom command to "install"
  # the library into install/lib
  add_custom_command(TARGET generate_pymoist_interface_library
    POST_BUILD
    # We first need to make a lib dir if it doesn't exist. If not, then
    # the next command can copy the script into a *file* called lib because
    # of a race condition (if install/lib/ isn't mkdir'd first)
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_INSTALL_PREFIX}/lib
    # Now we copy the file (if different...though not sure if this is useful)
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${PYMOIST_INTERFACE_LIBRARY}" ${CMAKE_INSTALL_PREFIX}/lib
    )

  # We use INTERFACE libraries to create a sort of "fake" target library we can use
  # to make libFVdycoreCubed_GridComp.a depend on. It seems to work!
  add_library(pymoist_interface_py INTERFACE)

  # The target_include_directories bits were essentially stolen from the esma_add_library
  # code...
  target_include_directories(pymoist_interface_py INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}> # stubs
    # modules and copied *.h, *.inc
    $<BUILD_INTERFACE:${esma_include}/${this}>
    $<INSTALL_INTERFACE:include/${this}>
    )
  target_link_libraries(pymoist_interface_py INTERFACE ${PYMOIST_INTERFACE_LIBRARY})

  # This makes sure the library is built first
  add_dependencies(pymoist_interface_py generate_pymoist_interface_library)

  # This bit is to resolve an issue and Google told me to do this. I'm not
  # sure that the LIBRARY DESTINATION bit actually does anything since
  # this is using INTERFACE
  install(TARGETS pymoist_interface_py
    EXPORT ${PROJECT_NAME}-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
    )

endif ()


if (GPU_PHYSICS)
  list (APPEND srcs
    openacc/Process_Library.F90
    openacc/gfdl_cloud_microphys.F90)
else ()
  list (APPEND srcs
    Process_Library.F90
    gfdl_cloud_microphys.F90)
endif ()

#if (BUILD_SERIALBOX_SER)
#  list (APPEND srcs
#    #sb/GEOS_GFDL_1M_InterfaceMod.F90
#    sb/aer_actv_single_moment.F90
#    sb/GEOS_MoistGridComp.F90
#    )
#else ()
#  list (APPEND srcs
#    # GEOS_GFDL_1M_InterfaceMod.F90
#    aer_actv_single_moment.F90
#    GEOS_MoistGridComp.F90
#    )
#endif()

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

if (BUILD_PYMOIST_INTERFACE)
  esma_add_library (${this}
    SRCS ${srcs}
    DEPENDENCIES pymoist_interface_py GEOS_Shared GMAO_mpeu MAPL Chem_Shared Chem_Base esmf)
else ()
  esma_add_library (${this}
    SRCS ${srcs}
    DEPENDENCIES GEOS_Shared GMAO_mpeu MAPL Chem_Shared Chem_Base esmf)
endif ()

if(BUILD_SERIALBOX_SER)
  message(STATUS "Building with Serialbox serialization statements")
  message(STATUS "  Serialbox root ${SERIALBOX_ROOT}")
  add_definitions(-DSERIALIZE)
  target_include_directories(${this} PRIVATE ${SERIALBOX_ROOT}/include)
  target_link_directories(${this} PRIVATE ${SERIALBOX_ROOT}/lib/)
  find_library(SERIALBOX_FORTRAN_LIB
      NAMES SerialboxFortran
      HINTS ${SERIALBOX_ROOT}/lib/
      REQUIRED
  )
  find_library(SERIALBOX_C_LIB
      NAMES SerialboxC
      HINTS ${SERIALBOX_ROOT}/lib/
      REQUIRED
  )
  message(STATUS "  Serialbox Fortran library found: ${SERIALBOX_FORTRAN_LIB}")
  message(STATUS "  Serialbox C library found: ${SERIALBOX_C_LIB}")
  target_link_libraries(${this} PRIVATE ${SERIALBOX_FORTRAN_LIB} ${SERIALBOX_C_LIB})
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



