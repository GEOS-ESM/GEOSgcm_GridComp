# Detect if we are using Open MPI and add oversubscribe/bind-to/map-by flags
string(REPLACE " " ";" MPI_Fortran_LIBRARY_VERSION_LIST ${MPI_Fortran_LIBRARY_VERSION_STRING})
list(GET MPI_Fortran_LIBRARY_VERSION_LIST 0 MPI_Fortran_LIBRARY_VERSION_FIRSTWORD)
if(MPI_Fortran_LIBRARY_VERSION_FIRSTWORD MATCHES "Open")
  set(MPIEXEC_PREFLAGS "--oversubscribe;--bind-to;core;--map-by;core" CACHE STRING "Flags for mpiexec" FORCE)
endif()

file(STRINGS "cases.txt" TEST_CASES)

if(APPLE)
  set(LD_PATH "DYLD_LIBRARY_PATH")
else()
  set(LD_PATH "LD_LIBRARY_PATH")
endif()

string(ASCII 27 ESC)
set(ORANGE "${ESC}[1;38;5;214m")
set(RESET  "${ESC}[0m")

# # Sync regression test data from S3 (skipped if aws CLI, API key, or internet is not available)
# find_program(AWS_CLI aws)

# set(_api_key_exists FALSE)
# if(DEFINED ENV{AWS_API_KEY_GMAO_SITEAM_S3} AND NOT "$ENV{AWS_API_KEY_GMAO_SITEAM_S3}" STREQUAL "")
#   set(_api_key_exists TRUE)
# endif()

# set(_internet_is_available FALSE)
# file(DOWNLOAD https://www.nasa.gov/ ${CMAKE_BINARY_DIR}/.nasa_connectivity_test
#   TIMEOUT 10
#   STATUS _download_status
# )
# list(GET _download_status 0 _download_result)
# if(_download_result EQUAL 0)
#   set(_internet_is_available TRUE)
# endif()
# file(REMOVE ${CMAKE_BINARY_DIR}/.nasa_connectivity_test)

# if(AWS_CLI AND _api_key_exists AND _internet_is_available)
#   set(S3_URI "s3://gmao-usw2-si-team/regression-data")
#   set(LOCAL_DIR "${CMAKE_BINARY_DIR}/regression-data")
#   add_custom_target(sync_regression_data ALL
#     COMMAND ${CMAKE_COMMAND} -E make_directory "${LOCAL_DIR}"
#     COMMAND ${CMAKE_COMMAND} -DAWS_CLI=${AWS_CLI} -DS3_URI=${S3_URI} -DLOCAL_DIR=${LOCAL_DIR}
#       -P ${CMAKE_CURRENT_SOURCE_DIR}/sync_data.cmake
#     BYPRODUCTS "${LOCAL_DIR}"
#     COMMENT "Syncing regression test data from ${S3_URI}"
#     VERBATIM
#   )
# else()
#   message(WARNING "${ORANGE}AWS CLI or AWS_API_KEY_GMAO_SITEAM_S3 not available - regression test data sync and comparisons will be skipped${RESET}")
# endif()

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
      -P ${CMAKE_CURRENT_SOURCE_DIR}/run_case.cmake
  )
  set_tests_properties("${TEST_CASE}" PROPERTIES
    LABELS "REGRESSION"
    ENVIRONMENT "${LD_PATH}=${CMAKE_BINARY_DIR}/lib:$ENV{${LD_PATH}}"
  )
endforeach()
