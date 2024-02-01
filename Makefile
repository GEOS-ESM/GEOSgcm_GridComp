programName = TEST_MOIST

FC = gfortran
#OPT = -O3 -Mflushz -Mfunc32 -Kieee #NVIDIA compiler options
#OPT = -O3 -Mfunc32 -Kieee -acc=gpu -gpu=flushz -Minfo=acc #NVIDIA compiler options
OPT = -g -O3 -fPIC -ffree-line-length-0 -fopenacc -foffload="-lgfortran -lgomp -lm" -foffload=nvptx-none # Gfortran compiler options
#OPT = -O3 -g -march=core-avx2 -fma -qopt-report0 -ftz -align all -fno-alias -align array32byte -traceback -assume realloc_lhs -fpe3 -fp-model consistent -assume noold_maxminloc -align dcommons #-prof-gen=srcpos# Ifort compiler options
OBJ = MAPL_SatVapor.o MAPL_Constants.o GEOS_Utilities.o \
      moist_subroutines_process_library.o \
	  timing_module.o \
	  test_moist_subroutines_process_library.o \
	  standalone.o

FFLAGS ?=$(OPT)

%.o : %.F90
	$(FC) -c -o $@ $^ $(FFLAGS) $(INC)

$(programName) : $(OBJ)
	$(FC) -o $@ $^ $(FFLAGS) $(INC) $(LIB)

clean :
	rm -rf *.o *.mod $(programName) *~ *.out output_data s_* coverage* *.gcda *.gcno *.spi *.spl *.dyn
