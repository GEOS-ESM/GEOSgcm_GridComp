# Cases that couple components across the AGCM's children (here FV3's ADV with
# GOCART2G under Physics). Each child's own cases live in its regression/.
file(STRINGS "cases.txt" TEST_CASES)

include(esma_regression_tests)
esma_add_regression_tests(
  NAME agcm
  DATA_PATH GEOSagcm_GridComp
  TEST_CASES ${TEST_CASES}
  EXTRA_ENV OMP_NUM_THREADS=1
)
