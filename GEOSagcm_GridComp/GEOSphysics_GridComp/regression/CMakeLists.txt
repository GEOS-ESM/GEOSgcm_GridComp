file(STRINGS "cases.txt" TEST_CASES)

# moist-gfdl1m-gf-uw-sa uses FakeGOCART_GridComp from @GEOSradiation_GridComp/regression
include(esma_regression_tests)
esma_add_regression_tests(
  NAME moist
  DATA_PATH GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp
  TEST_CASES ${TEST_CASES}
  EXTRA_ENV OMP_NUM_THREADS=1
)
