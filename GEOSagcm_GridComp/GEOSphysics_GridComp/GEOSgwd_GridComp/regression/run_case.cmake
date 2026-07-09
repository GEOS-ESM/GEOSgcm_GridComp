include("${ESMA_REGRESSION_HELPERS}")

function(run_case case_name regression_data_dir)
  string(RANDOM LENGTH 24 expdir)
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E make_directory ${expdir}
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_LIST_DIR}/${case_name} ${expdir}
  )

  set(root_dir ${regression_data_dir}/${case_name})
  set(num_procs "6")
  set(start_date_time "1891-03-01T00:00:00")
  set(restart_dir ${root_dir}/checkpoints/${start_date_time})
  set(checkpoints_dir ${root_dir}/checkpoints/last)

  if(NOT EXISTS ${root_dir})
    message(STATUS "Regression data not found for ${case_name}: ${root_dir} -- skipping")
    return()
  endif()

  copy_directory(${restart_dir} ${expdir}/checkpoints/${start_date_time})
  copy_file(${regression_data_dir}/newmfspectra40_dc25.nc ${expdir})
  run_geos(${num_procs} ${case_name} ${expdir})
  compare_results(${checkpoints_dir} ${expdir}/checkpoints/last)

  # execute_process(
  #   COMMAND ${CMAKE_COMMAND} -E rm -rf ${expdir}
  # )
endfunction()

run_case(${TEST_CASE} ${REGRESSION_DATA_DIR})
