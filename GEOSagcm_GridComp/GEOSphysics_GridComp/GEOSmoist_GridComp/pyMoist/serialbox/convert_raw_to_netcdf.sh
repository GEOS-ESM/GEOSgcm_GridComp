EXPERIMENT_SCRATCH_DIR=/home/mad/work/fp/experiments/TBC_C24_L72/scratch
OUTPUT_NETCDF_DIR=/home/mad/work/fp/savepoints
LOG_RANK0_FILE=$EXPERIMENT_SCRATCH_DIR/logs/1/rank.0/stdout

python extract_nml_from_log.py $LOG_RANK0_FILE
mv -f input.nml $EXPERIMENT_SCRATCH_DIR/serialized_data
ndsl-serialbox_to_netcdf $EXPERIMENT_SCRATCH_DIR/serialized_data $OUTPUT_NETCDF_DIR
