programName = GEOS_ncar_gwd

FC = nvfortran
#OPT = -O3 -Kieee -Mflushz -Mfunc32 -Minfo #NVIDIA compiler options
OPT = -O3 -Mfunc32 -Kieee -acc=gpu -gpu=flushz -Minfo=acc #NVIDIA compiler options
#OPT = -O3 -fPIC -fopenacc -foffload="-lgfortran -lgomp -lm" -foffload=nvptx-none #Gfortran compiler options
# OPT = -O3 -g -march=core-avx2 -fma -qopt-report0 -ftz -align all -fno-alias -align array32byte -traceback -assume realloc_lhs -fpe3 -fp-model consistent -assume noold_maxminloc -align dcommons # Ifort compiler options
OBJ = MAPL_Constants.o cesm_const_mod.F90 coords_1d.o linear_1d_operators.o interpolate_data.o vdiff_lu_solver.o gw_utils.o gw_diffusion.o gw_common.o gw_convect.o gw_oro.o gw_drag.o standalone.o
INC = -I/usr/include
LIB = -L/usr/lib/x86_64-linux-gnu -lnetcdff

# *** Uncomment the below line if compiling using gfortran with OpenACC offload ***
#gw_oro.o : OPT = -O0 -fPIC -fopenacc -foffload=nvptx-none -foffload="-lgfortran -lgomp -lm"

%.o : %.F90
	$(FC) -c -o $@ $^ $(OPT) $(INC)

$(programName) : $(OBJ)
	$(FC) -o $@ $^ $(OPT) $(INC) $(LIB)

clean :
	rm -rf *.o *.mod $(programName) *~ *.out output_data s_*
