source ../issm_env
rm -f AntarcticaGEOS.bin AntarcticaGEOS.outbin AntarcticaGEOS.outlog AntarcticaGEOS.errlog
rm -rf netcdfs
(export LD_LIBRARY_PATH=$PYTHON_LIB:$LD_LIBRARY_PATH; python Antarctica_input_control.py)
${ISSM_DIR}/bin/issm.exe StressbalanceSolution $(pwd) AntarcticaGEOS 2>> AntarcticaGEOS.errlog
(export LD_LIBRARY_PATH=$PYTHON_LIB:$LD_LIBRARY_PATH; python Antarctica_input_transient.py)
