programName = TEST_MOIST 

FC = ifort
#OPT = -O3 -Mfunc32 -Kieee #-acc=gpu -gpu=flushz #NVIDIA compiler options
#OPT = -g -O3  -ffree-line-length-0#-fopenmp # Gfortran compiler options
OPT =  -O3 -g -march=core-avx2 -fma -qopt-report0 -ftz -align all -fno-alias -align array32byte -traceback -assume realloc_lhs -fpe3 -fp-model consistent -assume noold_maxminloc -align dcommons# Ifort compiler options
OBJ = InternalConstants.o MathConstants.o PhysicalConstants.o Constants.o GEOS_Utilities.o m_fpe.o aer_cloud.o MAPL_SatVapor.o aer_actv_single_moment.o Process_Library.o ConvPar_GF_Shared.o module_gate.o ConvPar_GF2020.o GEOS_GF_InterfaceMod.o standalone.o

%.o : %.F90
	$(FC) -c -o $@ $^ $(OPT) $(INC)

$(programName) : $(OBJ)
	$(FC) -o $@ $^ $(OPT) $(INC) $(LIB)

clean :
	rm -rf *.o *.mod $(programName) *~ *.out output_data s_*
