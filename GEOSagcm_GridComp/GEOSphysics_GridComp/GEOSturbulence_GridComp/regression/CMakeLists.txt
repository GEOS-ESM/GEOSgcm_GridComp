# Root of the run2 case: runs TURBULENCE's Run1 and Run2 with Surface's import
# changes in between (README.md). Test only, so not in the component library.
add_subdirectory(TurbulencePhasesDriver_GridComp)

file(STRINGS "cases.txt" TEST_CASES)

include(esma_regression_tests)
esma_add_regression_tests(
  NAME turbulence
  DATA_PATH GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSturbulence_GridComp
  TEST_CASES ${TEST_CASES}
  EXTRA_ENV OMP_NUM_THREADS=1
)
