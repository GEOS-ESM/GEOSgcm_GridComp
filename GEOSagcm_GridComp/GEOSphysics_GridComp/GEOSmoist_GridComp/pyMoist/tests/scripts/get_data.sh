#!/bin/bash

set -e -x

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}"; )" &> /dev/null && pwd; )"

TEST_DATA_PATH="$SCRIPT_DIR/../../test_data/11.5.2/"
mkdir -p $TEST_DATA_PATH
cd $TEST_DATA_PATH
wget https://portal.nccs.nasa.gov/datashare/astg/smt/geos-fp/translate/11.5.2/x86_GNU/Moist/TBC_C24L72.tar.gz
tar -xzvf TBC_C24L72.tar.gz
rm TBC_C24L72.tar.gz
