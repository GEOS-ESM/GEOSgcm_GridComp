#!/bin/bash
#
# simple script for generating issm meshes at different resolutions
# 

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd $SCRIPT_DIR
cd ../

h_max=$1
h_min=$2


find . -mindepth 1 -type d ! -name utils_issm -exec bash -c '
h_max="$1"
h_min="$2"
shift 2

for dir do
(
    hdir=$(pwd)
    cd "$dir" || exit
    name=$(basename "$dir")
    
    for f in "${name}_meshgen.py" "${name}_parameterize.py" "${name}_control.py" "${name}_finalize.py"; do
        [[ -f "$f" ]] || exit
    done

    source "$hdir/issm_env"
    rm -f ISSM_${name}.bin ISSM_${name}.outbin ISSM_${name}.errlog 
    rm -rf netcdfs
   
    (export LD_LIBRARY_PATH=$PYTHON_LIB:$LD_LIBRARY_PATH; python ${name}_meshgen.py "$h_max" "$h_min")
)
done
' bash "$h_max" "$h_min" {} +

source issm_env
domain_name=$(
  LD_LIBRARY_PATH="$PYTHON_LIB:$LD_LIBRARY_PATH" \
  python ./utils_issm/domain_name.py
)


echo ""
echo "================================================================================================================"
echo "Candidate ISSM mesh:"
echo ""
echo "Domain name: $domain_name"
echo "(ME=mean edge length [meters], N = total nodes)"
echo ""
echo "================================================================================================================"
echo ""
