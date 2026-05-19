find . -mindepth 1 -type d ! -name utils_issm -exec bash -c '
for dir do
(
    hdir=$(pwd)
    cd "$dir" || exit
    name=$(basename "$dir")
    
    for f in "${name}_meshgen.py" "${name}_parameterize.py" "${name}_control.py" "${name}_finalize.py"; do
    [[ -f "$f" ]] || exit
    done

    source "$hdir/issm_env"
    rm -f ISSM_${name}.bin ISSM_${name}.outbin ISSM_${name}.outlog ISSM_${name}.errlog
    rm -rf netcdfs
    (export LD_LIBRARY_PATH=$PYTHON_LIB:$LD_LIBRARY_PATH; python ${name}_meshgen.py)
    (export LD_LIBRARY_PATH=$PYTHON_LIB:$LD_LIBRARY_PATH; python ${name}_control.py)
    ${ISSM_DIR}/bin/issm.exe StressbalanceSolution $(pwd) ISSM_${name} 2>> ISSM_${name}.errlog
    (export LD_LIBRARY_PATH=$PYTHON_LIB:$LD_LIBRARY_PATH; python ${name}_finalize.py)
)
done
' bash {} +

source issm_env
domain_name=$(
  LD_LIBRARY_PATH="$PYTHON_LIB:$LD_LIBRARY_PATH" \
  python ./utils_issm/domain_name.py
)

mkdir "$domain_name"

find . -type f -name "*.bin" -exec cp -t "$domain_name" {} +
find . -type f -name "*.toolkits" -exec cp -t "$domain_name" {} +

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
