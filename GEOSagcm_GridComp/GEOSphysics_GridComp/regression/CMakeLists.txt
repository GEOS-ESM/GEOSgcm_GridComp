# Cases that couple several Physics children (here SOLAR and IRRAD with GOCART2G
# as their aerosol provider). Each child's own standalone cases live in its regression/.
file(STRINGS "cases.txt" TEST_CASES)

include(esma_regression_tests)
esma_add_regression_tests(
  NAME physics
  DATA_PATH GEOSagcm_GridComp/GEOSphysics_GridComp
  TEST_CASES ${TEST_CASES}
  EXTRA_ENV OMP_NUM_THREADS=1
)
