esma_set_this()

option(BUILD_PYMKIAU_INTERFACE "Build pyMKIAU interface" OFF)

set (srcs
  IAU_GridCompMod.F90
  GEOS_mkiauGridComp.F90
  DFI_GridComp.F90
  mkiau_specmod.F90
  DynVec_GridComp.F90
  )

if (BUILD_PYMKIAU_INTERFACE)
  list (APPEND srcs
    pyMKIAU/interface/interface.f90
    pyMKIAU/interface/interface.c)

  message(STATUS "Building pyMKIAU interface")

  add_definitions(-DPYMKIAU_INTEGRATION)

  # The Python library creation requires mpiexec/mpirun to run on a
  # compute node. Probably a weird SLURM thing?
  find_package(Python3 COMPONENTS Interpreter REQUIRED)

  # Set up some variables in case names change
  set(PYMKIAU_INTERFACE_LIBRARY     ${CMAKE_CURRENT_BINARY_DIR}/libpyMKIAU_interface_py.so)
  set(PYMKIAU_INTERFACE_HEADER_FILE ${CMAKE_CURRENT_BINARY_DIR}/pyMKIAU_interface_py.h)
  set(PYMKIAU_INTERFACE_FLAG_HEADER_FILE ${CMAKE_CURRENT_SOURCE_DIR}/pyMKIAU/interface/interface.h)
  set(PYMKIAU_INTERFACE_SRCS        ${CMAKE_CURRENT_SOURCE_DIR}/pyMKIAU/interface/interface.py)

  # This command creates the shared object library from Python
  add_custom_command(
    OUTPUT ${PYMKIAU_INTERFACE_LIBRARY}
    # Note below is essentially:
    #  mpirun -np 1 python file
    # but we use the CMake options as much as we can for flexibility
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PYMKIAU_INTERFACE_FLAG_HEADER_FILE} ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${Python3_EXECUTABLE} ${PYMKIAU_INTERFACE_SRCS}
    BYPRODUCTS ${PYMKIAU_INTERFACE_HEADER_FILE}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    MAIN_DEPENDENCY ${PYMKIAU_INTERFACE_SRCS}
    COMMENT "Building pyMKIAU interface library with Python"
    VERBATIM
  )

  # This creates a target we can use for dependencies and post build
  add_custom_target(generate_pyMKIAU_interface_library DEPENDS ${PYMKIAU_INTERFACE_LIBRARY})

  # Because of the weird hacking of INTERFACE libraries below, we cannot
  # use the "usual" CMake calls to install() the .so. I think it's because
  # INTERFACE libraries don't actually produce any artifacts as far as
  # CMake is concerned. So we add a POST_BUILD custom command to "install"
  # the library into install/lib
  add_custom_command(TARGET generate_pyMKIAU_interface_library
    POST_BUILD
    # We first need to make a lib dir if it doesn't exist. If not, then
    # the next command can copy the script into a *file* called lib because
    # of a race condition (if install/lib/ isn't mkdir'd first)
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_INSTALL_PREFIX}/lib
    # Now we copy the file (if different...though not sure if this is useful)
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${PYMKIAU_INTERFACE_LIBRARY}" ${CMAKE_INSTALL_PREFIX}/lib
    )

  # We use INTERFACE libraries to create a sort of "fake" target library we can use
  # to make libFVdycoreCubed_GridComp.a depend on. It seems to work!
  add_library(pyMKIAU_interface_py INTERFACE)

  # The target_include_directories bits were essentially stolen from the esma_add_library
  # code...
  target_include_directories(pyMKIAU_interface_py INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}> # stubs
    # modules and copied *.h, *.inc
    $<BUILD_INTERFACE:${esma_include}/${this}>
    $<INSTALL_INTERFACE:include/${this}>
    )
  target_link_libraries(pyMKIAU_interface_py INTERFACE ${PYMKIAU_INTERFACE_LIBRARY})

  # This makes sure the library is built first
  add_dependencies(pyMKIAU_interface_py generate_pyMKIAU_interface_library)

  # This bit is to resolve an issue and Google told me to do this. I'm not
  # sure that the LIBRARY DESTINATION bit actually does anything since
  # this is using INTERFACE
  install(TARGETS pyMKIAU_interface_py
    EXPORT ${PROJECT_NAME}-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
    )

endif ()

if (BUILD_PYMKIAU_INTERFACE)
  set(dependencies pyMKIAU_interface_py MAPL_cfio_r4 NCEP_sp_r4i4 GEOS_Shared GMAO_mpeu MAPL FVdycoreCubed_GridComp ESMF::ESMF NetCDF::NetCDF_Fortran)
else ()
  set(dependencies MAPL_cfio_r4 NCEP_sp_r4i4 GEOS_Shared GMAO_mpeu MAPL FVdycoreCubed_GridComp ESMF::ESMF NetCDF::NetCDF_Fortran)
endif ()

esma_add_library (${this} SRCS ${srcs} DEPENDENCIES ${dependencies})
