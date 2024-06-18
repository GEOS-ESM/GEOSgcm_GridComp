#!/bin/bash

set -e -x

TEST_DATA_PATH="../../test_data"
mkdir -p $TEST_DATA_PATH
wget https://portal.nccs.nasa.gov/datashare/astg/smt/geos-fp/translate/11.5.2/Moist/aer_activation_data_nc.tar.gz
mv aer_activation_data_nc.tar.gz $TEST_DATA_PATH
cd $TEST_DATA_PATH
tar -xzvf aer_activation_data_nc.tar.gz
rm aer_activation_data_nc.tar.gz
