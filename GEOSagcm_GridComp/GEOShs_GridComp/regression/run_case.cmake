include("${ESMA_REGRESSION_HELPERS}")

function(run_case case_name regression_data_dir)
  string(RANDOM LENGTH 24 expdir)
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E make_directory ${expdir}
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_LIST_DIR}/${case_name} ${expdir}
  )

  set(num_procs "6")
  set(checkpoints_dir ${regression_data_dir}/${case_name}/checkpoints/last)

  run_geos(${num_procs} ${case_name} ${expdir})
  if (NOT EXISTS ${checkpoints_dir})
    message(WARNING "Baseline directory [${checkpoints_dir}] not found. Skipping comparison.")
  else()
    compare_results(${checkpoints_dir} ${expdir}/checkpoints/last)
  endif()

  file(REMOVE_RECURSE ${expdir})
endfunction()

run_case(${TEST_CASE} ${REGRESSION_DATA_DIR})
