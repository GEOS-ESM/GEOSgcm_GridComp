file(STRINGS "cases.txt" TEST_CASES)

include(esma_regression_tests)
esma_add_regression_tests(
  NAME gwd
  DATA_PATH GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSgwd_GridComp
  TEST_CASES ${TEST_CASES}
  # gwd-ncar-sa runs GWD threaded (mapl: misc: num_threads in gwd.yaml); the
  # OpenMP worker threads need more than the default stack for the NCAR GWD
  # column physics (same setting as gcm_run.j)
  EXTRA_ENV OMP_NUM_THREADS=1 OMP_STACKSIZE=16M
)