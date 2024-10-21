#!/bin/bash

set -e -x

TEST_DATA_PATH="../../test_data/" #11.5.2/Moist/
mkdir -p $TEST_DATA_PATH
cd $TEST_DATA_PATH
wget -r -nH --cut-dir=5 -np -R "index.html*" https://portal.nccs.nasa.gov/datashare/astg/smt/geos-fp/translate/11.5.2/x86_GNU/Moist/TBC_C24_L72_Debug/
