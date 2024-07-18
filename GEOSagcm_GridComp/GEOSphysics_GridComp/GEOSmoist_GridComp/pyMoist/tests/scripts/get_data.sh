#!/bin/bash

set -e -x

TEST_DATA_PATH="../../test_data"
mkdir -p $TEST_DATA_PATH
wget https://portal.nccs.nasa.gov/datashare/astg/smt/geos-fp/translate/11.5.2/Moist/AerActivation_nc.tar.gz
mv AerActivation_nc.tar.gz $TEST_DATA_PATH
cd $TEST_DATA_PATH
tar -xzvf AerActivation_nc.tar.gz
rm AerActivation_nc.tar.gz
