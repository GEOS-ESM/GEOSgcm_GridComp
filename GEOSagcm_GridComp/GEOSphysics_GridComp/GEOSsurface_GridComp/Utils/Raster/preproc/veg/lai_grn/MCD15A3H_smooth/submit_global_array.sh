#!/bin/bash
#SBATCH --job-name=MODIS_Global_10x10
#SBATCH --output=logs/global_%A_%a.out       # %A is master job ID, %a is array ID
#SBATCH --error=logs/global_%A_%a.err
#SBATCH --time=04:00:00                     # Max time per node
#SBATCH --qos=allnccs
#SBATCH --nodes=1                           # 1 node per array task
#SBATCH --ntasks=16                         # 16 cores per node
#SBATCH --account=s1583

##Please run the script by:
##mkdir -p logs temporary output
##sbatch --array=0-347%20 submit_global_array.sh
##when all the jobs are done, then:
##sbatch --array=348-647%20 submit_global_array.sh

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

MODIS_DATA_DIR="/discover/nobackup/projects/lis/LS_PARAMETERS/MODIS/MCD15A3H.061"

if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    echo "FATAL ERROR: Missing array parameter!"
    echo "You MUST submit this script with an array range, for example:"
    echo "sbatch --array=0-347%20 submit_global_array.sh"
    echo "sbatch --array=348-647%20 submit_global_array.sh"
    exit 1
fi

mkdir -p logs
mkdir -p temporary
mkdir -p output

# ==============================================================================
# [GLOBAL GRID MAPPING] Convert Array ID (0-647) to Lat & Lon
# Grid size: 10x10. Lat range: -90 to 90 (18 rows). Lon range: -180 to 180 (36 cols).
# ==============================================================================
LAT_IDX=$((SLURM_ARRAY_TASK_ID / 36))
LON_IDX=$((SLURM_ARRAY_TASK_ID % 36))

# Calculate numerical values using bc (Starts from Bottom-Left: Lat -90, Lon -180)
CALC_LAT=$(echo "-90 + $LAT_IDX * 10" | bc)
CALC_LON=$(echo "-180 + $LON_IDX * 10" | bc)

# Format to ensure .0 suffix (e.g., -90 -> -90.0)
BASE_LAT=$(printf "%.1f" $CALC_LAT)
BASE_LON=$(printf "%.1f" $CALC_LON)

echo "=========================================================="
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Node assigned: $SLURM_NODELIST"
echo "Processing Global Grid at Base Lat: $BASE_LAT, Lon: $BASE_LON"
echo "=========================================================="

# ==============================================================================
# BATCH 1: Execute sub-tasks 0 to 7
# ==============================================================================
echo "--- Launching Batch 1 (Tasks 0 to 7) ---"
for task_id in {0..7}
do
    python3 process_tile_array.py --task_id $task_id --base_lat $BASE_LAT --base_lon $BASE_LON --data_dir "$MODIS_DATA_DIR" &
done
wait
echo "Batch 1 completed successfully!"

# ==============================================================================
# BATCH 2: Execute sub-tasks 8 to 15
# ==============================================================================
echo "--- Launching Batch 2 (Tasks 8 to 15) ---"
for task_id in {8..15}
do
    python3 process_tile_array.py --task_id $task_id --base_lat $BASE_LAT --base_lon $BASE_LON --data_dir "$MODIS_DATA_DIR" &
done
wait
echo "Batch 2 completed successfully!"

# ==============================================================================
# STITCHING & COMPRESSION
# ==============================================================================
echo "--- Initiating Stitching and Compression Process ---"
python3 stitch_tiles.py --base_lat $BASE_LAT --base_lon $BASE_LON

if [ $? -eq 0 ]; then
    echo "Stitching script executed successfully!"
    echo "Cleaning up 16 uncompressed temporary sub-tiles..."
    
    for row in {0..3}; do
        for col in {0..3}; do
            LAT_MIN=$(echo "$BASE_LAT + $row * 2.5" | bc)
            LON_MIN=$(echo "$BASE_LON + $col * 2.5" | bc)
            LAT_MIN_FMT=$(printf "%.1f" $LAT_MIN)
            LON_MIN_FMT=$(printf "%.1f" $LON_MIN)
            rm -f "temporary/MCD15A3H_2003_2025_lat_${LAT_MIN_FMT}_lon_${LON_MIN_FMT}.nc"
        done
    done
    echo "Cleanup complete!"
else
    echo "ERROR: Stitching failed! Temporary files preserved for debugging."
    exit 1
fi


echo "=========================================================="
echo "Global Task ID $SLURM_ARRAY_TASK_ID (Lat: $BASE_LAT, Lon: $BASE_LON) finished!"
echo "=========================================================="

