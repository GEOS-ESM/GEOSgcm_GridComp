programName = TEST_MOIST

FC = nvfortran
#OPT = -O3 -Mflushz -Mfunc32 -Kieee #NVIDIA compiler options
OPT = -O3 -Mfunc32 -Kieee -acc=gpu -gpu=flushz -Minfo=acc #NVIDIA compiler options
#OPT = -g -O3 -fPIC -ffree-line-length-0 #--coverage #-fopenacc -foffload="-lgfortran -lgomp -lm" -foffload=nvptx-none # Gfortran compiler options
#OPT = -O3 -g -march=core-avx2 -fma -qopt-report0 -ftz -align all -fno-alias -align array32byte -traceback -assume realloc_lhs -fpe3 -fp-model consistent -assume noold_maxminloc -align dcommons #-prof-gen=srcpos# Ifort compiler options
OBJ = MAPL_SatVapor.o MAPL_Constants.o GEOS_Utilities.o module_gate.o ConvPar_GF_Shared.o moist_subroutines_GEOS5.o \
      moist_subroutines_GF2020.o moist_subroutines_process_library.o moist_subroutines_cloud_microphys.o \
	  moist_subroutines_aer_activation.o moist_subroutines_evap_subl_pdf_loop.o \
	  moist_subroutines_radcoup_loop.o \
	  timing_module.o \
	  test_moist_subroutines_process_library.o test_moist_subroutines_aer_activation.o \
	  test_moist_subroutines_GEOS5.o test_moist_subroutines_GF2020.o \
	  test_moist_subroutines_cloud_microphys.o test_moist_subroutines_evap_subl_pdf_loop.o \
	  test_moist_subroutines_radcoup_loop.o \
	  standalone.o

%.o : %.F90
	$(FC) -c -o $@ $^ $(OPT) $(INC)

$(programName) : $(OBJ)
	$(FC) -o $@ $^ $(OPT) $(INC) $(LIB)

clean :
	rm -rf *.o *.mod $(programName) *~ *.out output_data s_* coverage* *.gcda *.gcno *.spi *.spl *.dyn
