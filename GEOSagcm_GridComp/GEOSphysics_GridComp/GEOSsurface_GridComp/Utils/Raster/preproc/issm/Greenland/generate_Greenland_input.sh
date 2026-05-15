source ../issm_env
rm -f GreenlandGEOS.bin GreenlandGEOS.outbin GreenlandGEOS.outlog GreenlandGEOS.errlog
rm -rf netcdfs
(export LD_LIBRARY_PATH=$PYTHON_LIB:$LD_LIBRARY_PATH; python Greenland_input_control.py)
${ISSM_DIR}/bin/issm.exe StressbalanceSolution $(pwd) GreenlandGEOS 2>> GreenlandGEOS.errlog
(export LD_LIBRARY_PATH=$PYTHON_LIB:$LD_LIBRARY_PATH; python Greenland_input_transient.py)
