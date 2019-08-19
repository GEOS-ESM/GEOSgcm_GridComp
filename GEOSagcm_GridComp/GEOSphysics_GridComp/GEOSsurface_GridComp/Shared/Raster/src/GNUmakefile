
# $Id: GNUmakefile,v 1.7.4.2 2018/04/03 14:45:47 mathomp4 Exp $

#.SILENT:


# Make sure ESMADIR is defined
# ----------------------------
ifndef ESMADIR
       ESMADIR := ../..
endif

# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies
  include $(ESMADIR)/Config/GMAO_base.mk  # Generic stuff

  FREAL = $(FREAL4)

ifdef CEREBUS
	NETCDFDIR?=/home/suarez/lib/netcdf/
	F90=f90
	LIBS:=$(LIBS) -lU77
#	MODARG=-p
else

	NETCDFDIR=$(BASEDIR)/$(ARCH)
endif

LIBS =  -lz -L../lib -lraster $(LIB_NETCDF) $(LIB_SYS)
F90 = $(FC)
MODARG = $(M)

#ALT this has no effect on compilation
#OMP_NUM_THREADS:= 12 

#ALT: this is first attempt to pass -openmp in a more portable way
#     without changing files inside Config. Once we enable OMPFLAG,
#     this should be removed
RASTER_OMPFLAG =
ifeq ($(ESMA_FC), ifort)
#	RASTER_OMPFLAG = $(OMPFLAG)
endif
OPENMP_FLAG = $(OMPFLAG)

DIRS= ../mod ../bin ../lib ../til ../rst ../clsm ../data

SRDDIR=
MODDIR=../mod/
BINDIR=../bin/
LIBDIR=../lib/
INCDIR=../include/

INCS=$(INCDIR)/Raster.h

SRCS=$(wildcard mk*.F90 CombineRasters.F90 rmTinyTiles.F90 chk_clsm_params.F90)
BINS=$(addprefix $(BINDIR),$(basename $(notdir $(SRCS))))

SRCS_OPENMP= mkCatchParam.F90  
BINS_OPENMP=$(addsuffix _openmp, $(addprefix $(BINDIR),$(basename $(notdir $(SRCS_OPENMP)))))

F90SCRS= MAPL_Constants.F90 MAPL_Sort.F90 MAPL_Hash.F90 CubedSphere_GridMod.F90 rasterize.F90 leap_year.F90 date_time_util.F90\
 rmTinyCatchParaMod.F90 mod_process_hres_data.F90 easeV1_conv.F90 easeV2_conv.F90 comp_CATCHCN_AlbScale_parameters.F90

F90OBJS=$(addprefix $(BINDIR),$(addsuffix .o, $(basename $(F90SCRS))))
F90OBJS_OPENMP=$(addprefix $(BINDIR),$(addsuffix _openmp.o, $(basename $(F90SCRS))))

COBJS=$(addprefix $(BINDIR),zip.o util.o sort.o hash.o)

#f90OBJS=$(addprefix $(BINDIR),zenith.o)

fOBJS=$(addprefix $(BINDIR),sibalb.o zenith.o)

OBJS=$(COBJS) $(F90OBJS) $(f90OBJS) $(fOBJS)
OBJS_OPENMP=$(COBJS) $(F90OBJS_OPENMP) $(f90OBJS) $(fOBJS)

SRCH=-I$(NETCDFDIR)/include/netcdf -L$(NETCDFDIR)/lib -L$(LIBDIR) -I$(INCDIR)


MAPLDIR=$(ESMADIR)/src/GMAO_Shared/MAPL_Base

VPATH=$(MAPLDIR)

OPTS=$(F90FLAGS) $(RASTER_OMPFLAG) -DEIGHT_BYTE $(MODARG)$(MODDIR)

OPTS_OPENMP=$(F90FLAGS) $(OPENMP_FLAG) -DEIGHT_BYTE $(MODARG)$(MODDIR)

all: $(DIRS) $(LIBDIR)libraster.a $(INCS) $(BINS) $(BINS_OPENMP)

$(BINS):$(BINDIR)%: %.F90 $(LIBDIR)libraster.a
	$(F90) -o $@ $(OPTS) $(SRCH) $^ $(LIBS)

$(BINS_OPENMP):$(BINDIR)%_openmp: %.F90 $(LIBDIR)libraster_openmp.a
	$(F90) -o $@ $(OPTS_OPENMP) $(SRCH) $^ $(LIBS)

$(LIBDIR)libraster.a: $(OBJS)
	$(AR) -r $@ $? 2> /dev/null

$(LIBDIR)libraster_openmp.a: $(OBJS_OPENMP)
	$(AR) -r $@ $? 2> /dev/null

hash.o: $(MAPLDIR)/hash.c
	$(CC) -o $@ -c $(COPT3) $<
sort.o: $(MAPLDIR)/sort.c
	$(CC) -o $@ -c $(COPT3) $<
MAPL_Hash.o: $(MAPLDIR)/MAPL_Hash.F90
	$(F90) -o $@ -c $(OPTS) $<
MAPL_Sort.o: $(MAPLDIR)/MAPL_Sort.F90
	$(F90) -o $@ -c $(OPTS) $<

$(COBJS):$(BINDIR)%.o: $(SRCDIR)%.c $(INCS)
	$(CC) -o $@ -c $(COPT3) $<

$(F90OBJS):$(BINDIR)%.o: $(SRCDIR)%.F90 $(INCS)
	$(F90) -o $@ -c $(OPTS) $(SRCH) $<
	/bin/mv *.mod $(MODDIR)

$(F90OBJS_OPENMP):$(BINDIR)%_openmp.o: $(SRCDIR)%.F90 $(INCS)
	$(F90) -o $@ -c $(OPTS_OPENMP) $(SRCH) $<
	/bin/mv *.mod $(MODDIR)

$(f90OBJS):$(BINDIR)%.o: $(SRCDIR)%.f90 $(INCS)
	$(F90) -o $@ -c $(OPTS)  $(SRCH) $<

$(fOBJS):$(BINDIR)%.o: $(SRCDIR)%.f $(INCS)
	$(F90) -o $@ -c $(OPTS)  $(SRCH) $<

$(DIRS):
	$(MKDIR) $@

clean:
	-@$(RM) $(BINDIR)* $(SRCDIR)*~ $(LIBDIR)* $(MODDIR)* 2>/dev/null

