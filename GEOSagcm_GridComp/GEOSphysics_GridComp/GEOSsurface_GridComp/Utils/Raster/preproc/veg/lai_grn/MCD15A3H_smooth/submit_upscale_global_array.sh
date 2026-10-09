#!/bin/bash
#SBATCH --job-name=LAI_Global         # Job name
#SBATCH --array=1-36                  # Launch 36 array tasks (one for each H tile)
#SBATCH --time=02:00:00               # Time limit per task
#SBATCH --nodes=1                     # 1 node per array task
#SBATCH --ntasks=1                    # 1 Python script per node
#SBATCH --cpus-per-task=18            # 18 parallel threads per node for V01-V18
#SBATCH --account=s1583               # REPLACE THIS with your actual GMAO project code (e.g. s0827, j1010)
#SBATCH --constraint="[mil|cas]"
#SBATCH --qos=allnccs
#SBATCH --output=slurm_LAI_%A_%a.out  # %A is the array master Job ID, %a is the specific task ID (1-36)
#SBATCH --error=slurm_LAI_%A_%a.err   # Error log

echo "=========================================================="
echo "Starting LAI Array Task ID: $SLURM_ARRAY_TASK_ID at $(date)"
echo "Running on node: $SLURM_NODELIST"
echo "=========================================================="

# Pass the Array Task ID (1 to 36) directly to the Python script
python3 upscale_global_array.py $SLURM_ARRAY_TASK_ID

echo "Array Task $SLURM_ARRAY_TASK_ID finished at $(date)"
