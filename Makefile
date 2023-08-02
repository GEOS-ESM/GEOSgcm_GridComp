programName = run_UWSC
sourceName = test_UWSC.F90
sourceName2 =test_UWSC.o

FC = nvfortran
OPT = -O3 -Mfunc32 -Kieee #-acc=gpu -gpu=flushz #NVIDIA compiler options
#OPT = -O3 -fopenmp # Gfortran compiler options
#OPT = -O3 -qopenmp -fp-model source -fimf-archi-consistency=true# Ifort compiler options
OBJ = SHLWPARAMS.o MAPL_Constants.o GEOS_Utilities.o uwshcu.o $(sourceName2)

%.o : %.F90
	$(FC) -c -o $@ $^ $(OPT) $(INC)

$(programName) : $(OBJ)
	$(FC) -o $@ $^ $(OPT) $(INC) $(LIB)

clean :
	rm -rf *.o *.mod $(programName) *~ *.out output_data s_*
