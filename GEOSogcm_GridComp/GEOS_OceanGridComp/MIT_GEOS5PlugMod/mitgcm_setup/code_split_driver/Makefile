#
# Makefile for ESMA components.
#
# REVISION HISTORY:
#
# 09Jun2003  da Silva  First crack.
#
.NOTPARALLEL: # Force disabling of -j flag

# Make sure ESMADIR is defined
# ----------------------------
ifndef ESMADIR
       ESMADIR := $(PWD)/..
endif

# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies
  include $(ESMADIR)/Config/GMAO_base.mk  # Generic stuff

#                  ---------------------
#                  Standard ESMA Targets
#                  ---------------------

esma_help help:
	@echo "Standard ESMA targets:"
	@echo "% make esma_install    (builds and install under ESMADIR)"
	@echo "% make esma_clean      (removes deliverables: *.[aox], etc)"
	@echo "% make esma_distclean  (leaves in the same state as cvs co)"
	@echo "% make esma_doc        (generates PDF, installs under ESMADIR)"
	@echo "% make esma_help       (this message)"
	@echo "Environment:"
	@echo "      ESMADIR = $(ESMADIR)"
	@echo "      BASEDIR = $(BASEDIR)"
	@echo "         ARCH = $(ARCH)"
	@echo "         SITE = $(SITE)"

THIS = mitgcmuv_driver
LIB  = lib$(THIS).a
esma_install install pinstall: links
	@$(MAKE) --no-print-directory target

target: $(LIB)
	$(MKDIR) $(ESMALIB) $(ESMAINC)/$(THIS)
	$(CP) -p *.a            $(ESMALIB)
	$(CP) -p *.mod          $(ESMAINC)/$(THIS)

esma_clean clean:
	-$(RM) *~ *.[aox] *.[Mm][Oo][Dd]

esma_distclean distclean:
	-$(RM) *~ *.[aoxd] *.F?? *.[Mm][Oo][Dd]


#                  --------------------
#                  User Defined Targets
#                  --------------------

# ESMF bits
ESMF_F90COMPILEOPTS=#-threads
ESMF_F90COMPILEPATHS=-I$(BASEDIR)/Linux/include/esmf

#ESMF_F90LINKLIBS=-lcprts -limf -lm -lcxa -lunwind -lrt -ldl
#ESMF_F90ESMFLINKLIBS=-lesmf   -lmpich -lpmpich++ -lcprts -limf -lm -lcxa -lunwind -lrt -ldl
#ESMF_F90LINKPATHS=-L$(BASEDIR)/Linux/lib

MAPL_F90COMPILEPATHS=-I$(ESMADIR)/Linux/include/MAPL_Base

#MAPL_F90LINKPATHS=-L$(ESMADIR)/Linux/lib -lMAPL_Base -lesmf   -lmpich -lpmpich++ -lcprts -limf -lm -lcxa -lunwind -lrt -ldl

F_MAIN=main.F
F_SRCS=
C_SRCS= setdir.c

F90_FIXED_SRCS =driver/driver_init_mod.FOR         \
                driver/driver_run_mod.FOR          \
                driver/driver_get_export_state.FOR \
                driver/driver_set_import_state.FOR
F90_FREE_SRCS  =driver/mitgcm_driver_mod.F90

F90_FIXED_SRCS+=state/dynvars_h/dynvars_h_state_save_restore_mod.FOR
F90_FREE_SRCS +=state/dynvars_h/dynvars_h_state_mod.F90                \
                state/dynvars_h/dynvars_h_state_types_mod.F90          \
                state/dynvars_h/dynvars_h_state_create_destroy_mod.F90

F90_FIXED_SRCS+=state/export/export_state_fill_mod.FOR
F90_FREE_SRCS +=state/export/export_state_create_destroy_mod.F90 \
                state/export/export_state_mod.F90                \
                state/export/export_state_types_mod.F90

F90_FREE_SRCS+= state/ffields_h/ffields_h_state_mod.F90

F90_FIXED_SRCS+=state/import/import_state_fill_mod.FOR
F90_FREE_SRCS +=state/import/import_state_create_destroy_mod.F90 \
                state/import/import_state_mod.F90                \
                state/import/import_state_types_mod.F90

F90_FREE_SRCS +=state/mitgcm_state/mitgcm_state_mod.F90                \
                state/mitgcm_state/mitgcm_state_types_mod.F90          \
                state/mitgcm_state/mitgcm_state_create_destroy_mod.F90 \
                state/mitgcm_state/mitgcm_state_save_restore_mod.F90   \
                state/mitgcm_state/mitgcm_state_getdp_mod.F90

F90_FREE_SRCS +=state/stackvars/stackvars_state_mod.F90

F90_FIXED_SRCS+=state/timevars/timevars_state_save_restore_mod.FOR
F90_FREE_SRCS +=state/timevars/timevars_state_mod.F90                \
                state/timevars/timevars_state_types_mod.F90          \
                state/timevars/timevars_state_create_destroy_mod.F90

F90_FREE_SRCS +=utils/str4c.F90

F_OBJS=$(F_SRCS:.F=.o)
F90_OBJS1=$(F90_FREE_SRCS:.F90=.o)
F90_OBJS2=$(F90_FIXED_SRCS:.FOR=.o)
F90_OBJS=$(notdir $(F90_OBJS1)) \
         $(notdir $(F90_OBJS2))
C_OBJS=$(C_SRCS:.c=.o)
OBJS=$(F_OBJS) \
     $(C_OBJS) \
     $(F90_OBJS)
MAIN_OBJS=$(F_MAIN:.F=.o)
CPP_FLAGS=-DALLOW_USE_MPI -DALWAYS_USE_MPI
#FC=mpif90
FREAL8 =

FC_FLAGS=-I../inc -I../build $(FIXED_SOURCE) $(EXTENDED_SOURCE) $(FREAL8) $(BIG_ENDIAN) $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(MAPL_F90COMPILEPATHS) #-W0 -WB -assume byterecl -fpconstant 
FC_OPT=-O3 $(ALIGNCOM)
#CC=mpicc
CC_FLAGS=-I../inc -I../build
CC_OPT=-g
F90C=$(FC)
FOR_FLAGS=$(FC_FLAGS)
F90_FLAGS=-I../inc -I../build $(FREAL8) $(BIG_ENDIAN) $(F90FLAGS) $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILERATHS) $(MAPL_F90COMPILEPATHS) #-W0 -WB -assume byterecl -fpconstant 
F90_OPT=$(FC_OPT) -g #-check -ftrapuv
LIBS=-L../lib -lmitgcmuv $(ESMF_F90LINKPATHS) $(ESMF_F90ESMFLINKLIBS) $(MAPL_F90LINKPATHS)

.SUFFIXES:
.SUFFIXES: .o .c .FOR .F .f .c .F90

a.out: $(OBJS) $(MAIN_OBJS)
	@ echo Linking \"$?\"
	$(FC) $(OBJS) $(LIBS)
	cp a.out ../run_split

lib: libmitgcmuv_driver.a

libmitgcmuv_driver.a: $(OBJS)
	@ echo Creating \"$@\"
	ar rcv $@ $(OBJS)

%.o: %.F
	@ echo Compiling \"$<\"
	$(FC) $(FC_FLAGS) $(FC_OPT) $(CPP_FLAGS) -c $<

%.o: %.FOR
	@ echo Compiling \"$<\"
	$(F90C) $(FOR_FLAGS) $(F90_OPT) $(CPP_FLAGS) -c $<

%.o: %.F90
	@ echo Compiling \"$<\"
	$(F90C) $(F90_FLAGS) $(F90_OPT) $(CPP_FLAGS) -c $<

%.o: %.c
	@ echo Compiling \"$<\"
	$(CC) $(CC_FLAGS) $(CC_OPT) $(CPP_FLAGS) -c $<

links:
	@ln -sf $(F90_FIXED_SRCS) .
	@ln -sf $(F90_FREE_SRCS) .

CLEAN:
	@rm -f *.o *.mod STD* a.out
	@find . -type l -exec \rm -f {} \;

# Dependencies due to modules and headers (should be auto generated)
main.o: mitgcm_driver_mod.o
main.o: mitgcm_state_mod.o

mitgcm_driver_mod.o: driver_init_mod.o
mitgcm_driver_mod.o: driver_run_mod.o

driver_init_mod.o:   mitgcm_state_mod.o

driver_run_mod.o:    mitgcm_state_mod.o

mitgcm_state_mod.o:  mitgcm_state_create_destroy_mod.o
mitgcm_state_mod.o:  mitgcm_state_save_restore_mod.o
mitgcm_state_mod.o:  mitgcm_state_getdp_mod.o
mitgcm_state_mod.o:  mitgcm_state_types_mod.o

mitgcm_state_types_mod.o:  dynvars_h_state_mod.o
mitgcm_state_types_mod.o:  timevars_state_mod.o
mitgcm_state_types_mod.o:  stackvars_state_mod.o
mitgcm_state_types_mod.o:  ffields_h_state_mod.o
mitgcm_state_types_mod.o:  export_state_mod.o
mitgcm_state_types_mod.o:  import_state_mod.o

mitgcm_state_create_destroy_mod.o:  mitgcm_state_types_mod.o

mitgcm_state_save_restore_mod.o:    mitgcm_state_types_mod.o

mitgcm_state_getdp_mod.o:           mitgcm_state_types_mod.o

dynvars_h_state_mod.o:  dynvars_h_state_create_destroy_mod.o
dynvars_h_state_mod.o:  dynvars_h_state_save_restore_mod.o
dynvars_h_state_mod.o:  dynvars_h_state_types_mod.o

dynvars_h_state_create_destroy_mod.o: dynvars_h_state_types_mod.o

dynvars_h_state_save_restore_mod.o: dynvars_h_state_types_mod.o

timevars_state_mod.o:  timevars_state_create_destroy_mod.o
timevars_state_mod.o:  timevars_state_save_restore_mod.o
timevars_state_mod.o:  timevars_state_types_mod.o

timevars_state_create_destroy_mod.o: timevars_state_types_mod.o

timevars_state_save_restore_mod.o: timevars_state_types_mod.o

export_state_mod.o:  export_state_types_mod.o export_state_create_destroy_mod.o export_state_fill_mod.o
export_state_create_destroy_mod.o: export_state_types_mod.o

import_state_mod.o:  import_state_types_mod.o import_state_create_destroy_mod.o import_state_fill_mod.o
import_state_create_destroy_mod.o: import_state_types_mod.o
