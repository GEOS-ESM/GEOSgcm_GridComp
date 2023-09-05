programName = run_UWSC

FC = gfortran
#OPT = -O3 -Mfunc32 -Kieee #-acc=gpu -gpu=flushz #NVIDIA compiler options
OPT = -g -O3  -ffree-line-length-0#-fopenmp # Gfortran compiler options
#OPT = -O3 -qopenmp -fp-model source -fimf-archi-consistency=true# Ifort compiler options
OBJ = MAPL_Constants.o GEOS_Utilities.o moist_subroutines_process_library.o read_tracer_data.o uwshcu.o test_UWSC.o 

%.o : %.F90
	$(FC) -c -o $@ $^ $(OPT) $(INC)

$(programName) : $(OBJ)
	$(FC) -o $@ $^ $(OPT) $(INC) $(LIB)

clean :
	rm -rf *.o *.mod $(programName) *~ *.out output_data s_*
