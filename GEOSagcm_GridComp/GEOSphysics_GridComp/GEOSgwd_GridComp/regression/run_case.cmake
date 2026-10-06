include("${ESMA_REGRESSION_HELPERS}")

function(run_case case_name regression_data_dir)
  string(RANDOM LENGTH 24 expdir)
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E make_directory ${expdir}
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_LIST_DIR}/${case_name} ${expdir}
  )

  set(root_dir ${regression_data_dir}/${case_name})
  set(num_procs "6")
  set(checkpoints_dir ${root_dir}/checkpoints/last)

  if(NOT EXISTS ${root_dir})
    message(STATUS "Regression data not found for ${case_name}: ${root_dir} -- skipping")
    return()
  endif()

  copy_restarts(${root_dir} ${expdir})
  copy_file(${regression_data_dir}/newmfspectra40_dc25.nc ${expdir})
  run_geos(${num_procs} ${case_name} ${expdir})
  if(FORTRAN_COMPILER_ID STREQUAL "IntelLLVM") # only compare against IntelLLVM baselines
    # compare_results() globs the baseline directory, so stage links to just the
    # files we want compared. GWD_export.nc is left out until the baselines are
    # regenerated for the NCAR GWD science update.
    set(staged_baseline_dir ${expdir}/baseline)
    execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${staged_baseline_dir})
    foreach(fname IN ITEMS GWD_import.nc GWD_internal.nc)
      execute_process(
        COMMAND ${CMAKE_COMMAND} -E create_symlink ${checkpoints_dir}/${fname} ${staged_baseline_dir}/${fname}
      )
    endforeach()

    # The NCAR GWD update renamed import DTDT_DC to HT_dc, so the two files
    # disagree on the name of a field that is otherwise identical. Exclude both
    # spellings until the baselines are regenerated.
    compare_results(
      ${staged_baseline_dir} ${expdir}/checkpoints/last
      NANS_ARE_EQUAL
      EXCLUDE_VARS lons lats corner_lons corner_lats DQIDT DQLDT HT_mi WSPD_STABLE300M DTDT_DC HT_dc
    )
  endif()

  file(REMOVE_RECURSE ${expdir})
endfunction()

run_case(${TEST_CASE} ${REGRESSION_DATA_DIR})
