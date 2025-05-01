EXPERIMENT_SCRATCH_DIR=/Users/kfandric/experiments/TBC_C24_L72/scratch
OUTPUT_NETCDF_DIR=/Users/kfandric/netcdf
LOG_RANK0_FILE=$EXPERIMENT_SCRATCH_DIR/logs/1/rank.0/stdout

#python extract_nml_from_log.py $LOG_RANK0_FILE
cp /Users/kfandric/misc/input.nml $EXPERIMENT_SCRATCH_DIR/moist_serialbox_data
ndsl-serialbox_to_netcdf $EXPERIMENT_SCRATCH_DIR/moist_serialbox_data $OUTPUT_NETCDF_DIR