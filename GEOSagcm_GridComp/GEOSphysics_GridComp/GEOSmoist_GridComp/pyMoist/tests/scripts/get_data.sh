#!/bin/bash

set -e -x

TEST_DATA_PATH="../../test_data"
mkdir -p $TEST_DATA_PATH
wget https://portal.nccs.nasa.gov/datashare/astg/smt/geos-fp/translate/11.5.2/Moist/RadCouple.tar.gz
mv RadCouple.tar.gz $TEST_DATA_PATH
cd $TEST_DATA_PATH
tar -xzvf RadCouple.tar.gz
rm RadCouple.tar.gz
