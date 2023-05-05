programName = RUN_EDMF_MINIAPP

FC = nvfortran
#OPT = -O3 -Mflushz -Mfunc32 -Minfo #NVIDIA compiler options
OPT = -O3 -Mfunc32 -acc=gpu -cudalib=curand -gpu=flushz -Kieee #-Minfo=acc #NVIDIA compiler options
#OPT = -O3 -Mfunc32 -Mflushz -acc=multicore -Kieee
#OPT = -O3 # Gfortran compiler options
#OPT = -O3 -g -march=core-avx2 -fma -qopt-report0 -ftz -align all -fno-alias -align array32byte -traceback -assume realloc_lhs -fpe3 -fp-model consistent -assume noold_maxminloc -align dcommons # Ifort compiler options
OBJ = MAPL_Constants.o GEOS_Utilities.o edmfparams.o edmf.o standalone.o
INC = -I/usr/include
LIB = -L/usr/lib/x86_64-linux-gnu -lnetcdff

%.o : %.F90
	$(FC) -c -o $@ $^ $(OPT) $(INC)

$(programName) : $(OBJ)
	$(FC) -o $@ $^ $(OPT) $(INC) $(LIB)

clean :
	rm -rf *.o *.mod $(programName) *~ *.out output_data s_*
