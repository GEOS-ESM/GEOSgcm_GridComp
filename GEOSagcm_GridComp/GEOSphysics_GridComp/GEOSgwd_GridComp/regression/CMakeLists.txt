file(STRINGS "cases.txt" TEST_CASES)

include(esma_regression_tests)
esma_add_regression_tests(
  NAME gwd
  DATA_PATH GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSgwd_GridComp
  TEST_CASES ${TEST_CASES}
  EXTRA_ENV OMP_NUM_THREADS=1
)