programName = TEST_MOIST

FC = gfortran
#OPT = -O3 -Mflushz -Mfunc32 -Kieee #NVIDIA compiler options
#OPT = -O3 -Mfunc32 -Kieee -acc=gpu -gpu=flushz -Minfo=acc #NVIDIA compiler options
OPT = -g -O3 -fPIC -ffree-line-length-0 -fopenacc -foffload="-lgfortran -lgomp -lm" -foffload=nvptx-none # Gfortran compiler options
#OPT = -O3 -g -march=core-avx2 -fma -qopt-report0 -ftz -align all -fno-alias -align array32byte -traceback -assume realloc_lhs -fpe3 -fp-model consistent -assume noold_maxminloc -align dcommons #-prof-gen=srcpos# Ifort compiler options
OBJ = MAPL_Constants.o \
	moist_subroutines_aer_activation.o \
	test_moist_subroutines_aer_activation.o \
	standalone.o

%.o : %.F90
	$(FC) -c -o $@ $^ $(OPT) $(INC)

$(programName) : $(OBJ)
	$(FC) -o $@ $^ $(OPT) $(INC) $(LIB)

clean :
	rm -rf *.o *.mod $(programName) *~ *.out output_data s_* coverage* *.gcda *.gcno *.spi *.spl *.dyn
