file(STRINGS "cases.txt" TEST_CASES)

if(APPLE)
  set(LD_PATH "DYLD_LIBRARY_PATH")
else()
  set(LD_PATH "LD_LIBRARY_PATH")
  # OpenMPI flags
  string(REPLACE " " ";" MPI_Fortran_LIBRARY_VERSION_LIST ${MPI_Fortran_LIBRARY_VERSION_STRING})
  list(GET MPI_Fortran_LIBRARY_VERSION_LIST 0 MPI_Fortran_LIBRARY_VERSION_FIRSTWORD)
  if(MPI_Fortran_LIBRARY_VERSION_FIRSTWORD MATCHES "Open")
    set(MPIEXEC_PREFLAGS "--oversubscribe;--bind-to;core;--map-by;core" CACHE STRING "Flags for mpiexec" FORCE)
  endif()
endif()

# Sync regression test data from S3
set(HS_PATH GEOSagcm_GridComp/GEOShs_GridComp)
set(REGRESSION_DATA_DIR ${CMAKE_BINARY_DIR}/regression-data/${HS_PATH})
if (NOT EXISTS ${REGRESSION_DATA_DIR})
  message(STATUS "Syncing regression data from S3: ${HS_PATH}")
  esma_sync_data(
    TARGET sync_regression_data_hs
    URI "s3://gmao-usw2-si-team/regression-data/${HS_PATH}"
    DESTINATION "${REGRESSION_DATA_DIR}"
    CREDENTIALS_URL "https://llvsm4u7ij.execute-api.us-east-1.amazonaws.com/credentials"
    ALL
  )
endif()

# Regression tests
foreach(TEST_CASE ${TEST_CASES})
  add_test(
    NAME "${TEST_CASE}"
    COMMAND ${CMAKE_COMMAND}
      -DTEST_CASE=${TEST_CASE}
      -DMPIEXEC_EXECUTABLE=${MPIEXEC_EXECUTABLE}
      -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
      -DMY_BINARY_DIR=${CMAKE_BINARY_DIR}/bin
      "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
      -DREGRESSION_DATA_DIR=${REGRESSION_DATA_DIR}
      -P ${CMAKE_CURRENT_SOURCE_DIR}/run_case.cmake
  )
  set_tests_properties("${TEST_CASE}" PROPERTIES
    LABELS "REGRESSION"
    ENVIRONMENT "${LD_PATH}=${CMAKE_BINARY_DIR}/lib:$ENV{${LD_PATH}}"
  )
endforeach()
