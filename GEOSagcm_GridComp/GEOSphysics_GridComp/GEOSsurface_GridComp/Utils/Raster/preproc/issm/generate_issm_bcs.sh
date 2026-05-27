#!/bin/bash
#SBATCH --job-name=preproc_issm
#SBATCH --time=1:00:00
#SBATCH --ntasks=126
#SBATCH --qos=debug
#SBATCH --constraint=mil

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
    (export LD_LIBRARY_PATH=$PYTHON_LIB:$LD_LIBRARY_PATH; python ${name}_control.py)
    
    if [ -n "$SLURM_JOB_ID" ]; then
        mpirun -np $SLURM_NTASKS ${ISSM_DIR}/bin/issm.exe StressbalanceSolution $(pwd) ISSM_${name} 2>> ISSM_${name}.errlog
    else
        ${ISSM_DIR}/bin/issm.exe StressbalanceSolution $(pwd) ISSM_${name} 2>> ISSM_${name}.errlog
    fi
    
    (export LD_LIBRARY_PATH=$PYTHON_LIB:$LD_LIBRARY_PATH; python ${name}_finalize.py)
)
done
' bash "$h_max" "$h_min" {} +

source issm_env
domain_name=$(
  LD_LIBRARY_PATH="$PYTHON_LIB:$LD_LIBRARY_PATH" \
  python ./utils_issm/domain_name.py
)

rm -rf "$domain_name" && mkdir "$domain_name"

find . -type f -name "ISSM*.bin" -not -path "./ISSM_ME*/*" -exec cp -t "$domain_name" {} +
find . -type f -name "ISSM*.toolkits" -not -path "./ISSM_ME*/*" -exec cp -t "$domain_name" {} +

echo ""
echo "================================================================================================================"
echo "Created ISSM BCs!"
echo ""
echo "Domain name: $domain_name"
echo "(ME=mean edge length [meters], N = total nodes)"
echo ""
echo "ISSM BCs copied to $(pwd)/$domain_name"
echo "================================================================================================================"
echo ""
