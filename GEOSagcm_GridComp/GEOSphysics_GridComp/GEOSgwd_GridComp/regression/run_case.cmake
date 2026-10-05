include("${ESMA_REGRESSION_HELPERS}")

# The NCAR GWD update renamed import DTDT_DC to HT_dc and reworded its
# standard_name. Put the generated file back into the baseline's naming so the
# field is compared value-for-value rather than excluded. Rewriting the output
# rather than the baseline keeps the read-only reference data untouched.
function(align_import_to_baseline import_file)
  find_program(NCRENAME_EXECUTABLE ncrename)
  find_program(NCATTED_EXECUTABLE ncatted)
  if(NOT NCRENAME_EXECUTABLE OR NOT NCATTED_EXECUTABLE)
    message(FATAL_ERROR "ncrename/ncatted (NCO) are required to compare the GWD import state")
  endif()

  execute_process(
    COMMAND ${NCRENAME_EXECUTABLE} -O -h -v HT_dc,DTDT_DC ${import_file}
    COMMAND_ERROR_IS_FATAL ANY
  )
  # ncrename stamps an NCO global attribute that the baseline does not carry.
  execute_process(
    COMMAND ${NCATTED_EXECUTABLE} -O -h
            -a standard_name,DTDT_DC,o,c,T\ tendency\ due\ to\ deep\ convection
            -a NCO,global,d,,
            ${import_file}
    COMMAND_ERROR_IS_FATAL ANY
  )
endfunction()

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
    align_import_to_baseline(${expdir}/checkpoints/last/GWD_import.nc)

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

    compare_results(
      ${staged_baseline_dir} ${expdir}/checkpoints/last
      NANS_ARE_EQUAL
      EXCLUDE_VARS lons lats corner_lons corner_lats DQIDT DQLDT HT_mi WSPD_STABLE300M
    )
  endif()

  file(REMOVE_RECURSE ${expdir})
endfunction()

run_case(${TEST_CASE} ${REGRESSION_DATA_DIR})
