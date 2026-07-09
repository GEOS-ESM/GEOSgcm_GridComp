include(esma_regression_tests)

file(STRINGS "cases.txt" TEST_CASES)

esma_add_regression_tests(
  NAME hs
  DATA_PATH GEOSagcm_GridComp/GEOShs_GridComp
  TEST_CASES ${TEST_CASES}
)
