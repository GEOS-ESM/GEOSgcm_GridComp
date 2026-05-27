macro(run_case CASE)
  string(RANDOM LENGTH 24 tempdir)
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E make_directory ${tempdir}
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_LIST_DIR}/${CASE} ${tempdir}
  )

  # Run the case
  set(num_procs "6")
  execute_process(
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${num_procs} -prepend-rank ${MPIEXEC_PREFLAGS} ${MY_BINARY_DIR}/GEOS.x cap.yaml
    RESULT_VARIABLE CMD_RESULT
    WORKING_DIRECTORY ${tempdir}
    COMMAND_ECHO STDOUT
  )
  if(EXISTS ${tempdir}/PET0.ESMF_LogFile)
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E cat ${tempdir}/PET0.ESMF_LogFile
    )
  endif()
  if(CMD_RESULT)
    set(MSG "Error running ${CASE}")
    message(FATAL_ERROR "${MSG}")
  endif()

  # Compare output against baseline
  set(BASELINE_DIR ${MY_BINARY_DIR}/../regression-data/GEOSagcm_GridComp/${CASE})
  if(NOT EXISTS ${BASELINE_DIR})
    string(ASCII 27 ESC)
    set(ORANGE "${ESC}[1;38;5;214m")
    set(RESET  "${ESC}[0m")
    message(WARNING "${ORANGE}Baseline directory not found, skipping comparison: ${BASELINE_DIR}${RESET}")
  else()
    file(GLOB baseline_files ${BASELINE_DIR}/*.nc)
    foreach(baseline_file IN LISTS baseline_files)
      get_filename_component(fname ${baseline_file} NAME)
      message(STATUS "Comparing ${fname}")
      execute_process(
        COMMAND cmp ${baseline_file} ${tempdir}/checkpoints/last/${fname}
        RESULT_VARIABLE CMP_RESULT
        OUTPUT_VARIABLE CMP_OUTPUT
        ERROR_VARIABLE CMP_OUTPUT
      )
      if(CMP_RESULT)
        message(FATAL_ERROR "Files differ: ${fname}\n${CMP_OUTPUT}")
      endif()
    endforeach()
  endif()

  # execute_process(
  #   COMMAND ${CMAKE_COMMAND} -E rm -rf ${tempdir}
  # )
endmacro()
run_case(${TEST_CASE})
