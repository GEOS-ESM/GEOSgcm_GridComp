#!/bin/bash
# ==============================================================================
# Script: prepare_upscale_dirs.sh
# Description: Creates year-specific directories (2003-2025) inside 'upscaled/'
#              and copies all files from the GMAO shared template directory into them.
# ==============================================================================

echo "Initializing directory structure for upscale processing (2003-2025)..."

# Define the source template directory on Discover
SRC_DIR="/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/land/veg/lai_grn/v2/MODIS_8-DayClim"

# Check if the template directory exists
if [ ! -d "${SRC_DIR}" ]; then
    echo "Error: Template directory not found at ${SRC_DIR}!"
    exit 1
fi

# Check if the template directory is empty
if [ -z "\((ls -A "\){SRC_DIR}")" ]; then
    echo "Warning: Template directory is empty. No files will be copied."
fi

# Loop through years to create directories and copy templates
for year in {2003..2025}; do
    target_dir="upscaled/${year}"
    
    # Create the target directory (if it does not exist)
    mkdir -p "${target_dir}"
    
    # Copy all files from the source directory to the target directory
    # Use cp -p to preserve original timestamps and file attributes
    cp -p ${SRC_DIR}/* ${target_dir}
    
    echo "Created template for ${year}"
done

echo "Directory preparation complete! You can now run upscale_global_array.py."
