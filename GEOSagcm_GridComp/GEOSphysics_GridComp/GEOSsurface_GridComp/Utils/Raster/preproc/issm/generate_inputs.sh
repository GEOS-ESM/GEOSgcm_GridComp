find . -mindepth 1 -type d ! -name utils_issm -exec bash -c '
for dir do
(
    cd "$dir" || exit

    name=$(basename "$dir")
    script="generate_${name}_input.sh"

    if [[ -x "$script" ]]; then
        "./$script"
    elif [[ -f "$script" ]]; then
        bash "$script"
    fi
)
done
' bash {} +

#source issm_env
#(export LD_LIBRARY_PATH=$PYTHON_LIB:$LD_LIBRARY_PATH; python ./utils_issm/domain_name.py)
