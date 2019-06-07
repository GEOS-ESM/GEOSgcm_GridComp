#
# Makefile for ESMA components.
#
# REVISION HISTORY:
#
# 09Jun2003  da Silva  First crack.
#

# Make sure ESMADIR is defined
# ----------------------------
ifndef ESMADIR
       ESMADIR := $(PWD)/../..
endif

# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies
  include $(ESMADIR)/Config/GMAO_base.mk  # System dependencies

  include $(ESMADIR)/Config/vectorize.mk  # Vectorize

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

THIS := $(shell basename `pwd`)
LIB   = lib$(THIS).a
#BINS  = StandAlone_FV3_Dycore.x interp_restarts.x c2c.x c2l_CFIO_offline.x 
BINS  = StandAlone_DynAdvCore.x StandAlone_AdvCore.x StandAlone_FV3_Dycore.x interp_restarts.x c2c.x rs_scale.x
#BINS  = StandAlone_DynAdvCore.x StandAlone_AdvCore.x StandAlone_FV3_Dycore.x c2l_CFIO_offline.x CreateTopography.x CreateInterpWeights.x interp_restarts.x c2c.x

ALLDIRS = fvdycore

SUBDIRS = $(wildcard $(ALLDIRS))

TARGETS = esma_install esma_clean esma_distclean esma_doc \
          install clean distclean doc 

export ESMADIR BASEDIR ARCH SITE

PROTEX_FLAGS = -bM

export ESMADIR BASEDIR ARCH SITE

$(TARGETS): 
	@ t=$@; argv="$(SUBDIRS)" ;\
	  for d in $$argv; do			 \
	    ( cd $$d				;\
	      echo ""; echo Making $$t in `pwd`          ;\
	      $(MAKE) -e $$t ) \
	  done
	$(MAKE) local_$@

local_esma_install local_install: $(LIB) $(BINS)
	$(MKDIR) $(ESMALIB) $(ESMAETC) $(ESMAINC)/$(THIS)
	$(CP) -p *.a                   $(ESMALIB)
	$(CP) -p *.[Mm][Oo][Dd]        $(ESMAINC)/$(THIS)
	$(CP) -p $(BINS)               $(ESMABIN)

local_esma_clean local_clean:
	-$(RM) *~ *.[aox] *.[Mm][Oo][Dd] 

local_esma_distclean local_distclean:
	-$(RM) *~ *.[aoxd] *.[Mm][Oo][Dd]

local_esma_doc local_doc:
	@$(PROTEX) $(PROTEX_FLAGS) *GridComp*.[fF]* > $(ESMADOC)/$(THIS).tex


#                  --------------------
#                  User Defined Targets
#                  --------------------


SRCS := sw.f90 jw.f90 testcases_3_4_5_6_stand_alone.f90 \
	GetWeightsC2C.F90 \
	GetWeights.F90 \
	CubeHalo.F90 \
        Cube2LatLon.F90 LatLon2Cube.F90 AppGridCreate.F90 FV_StateMod.F90 \
        AdvCore_GridCompMod.F90 \
        DynCore_GridCompMod.F90 CreateInterpWeights_GridCompMod.F90 \
        StandAlone_DynAdvCore_GridCompMod.F90\
	CubedSphereGridFactory.F90\
        CubeToLatLonRegridder.F90 \
        LatLonToCubeRegridder.F90 \
        CubeToCubeRegridder.F90 \
        CubeToLatLon.F90

OBJS := $(addsuffix .o, $(basename $(SRCS))) 
DEPS := $(addsuffix .d, $(basename $(SRCS))) \
        $(addsuffix .d, $(basename $(BINS)))

INC_DIRS = . $(INC_ESMF) $(INC_HERMES) $(INC_GMAO_SHARED) $(INC_MPI) $(INC_GEOS_FV1) $(INC_GMAO_pFIO) $(INC_FLAP)
MOD_DIRS = . $(INC_DIRS)

USER_FINCS  = $(foreach dir,$(INC_DIRS),$(I)$(dir))
USER_FINCS += $(I). $(I)$(INC_MPI)
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 
vpath % $(INC_DIRS) $(MOD_DIRS)

#USER_FFLAGS += $(BIG_ENDIAN) 
FDEFS       += $(D)SPMD $(D)MAPL_MODE

USER_LDFLAGS += $(OMPFLAG)

# MAT This option is used to enable single-precision
#     FV code in the dycore
#export FV_PRECISION = R4

ifeq ("$(FV_PRECISION)","R4")
USER_FDEFS += $(D)SINGLE_FV $(D)OVERLOAD_R4
THIS_GFDL_FMS = GFDL_fms_r4
FREAL= $(FREAL4)
else
THIS_GFDL_FMS = GFDL_fms_r8
FREAL= $(FREAL8)
endif

LIB_SHARED = $(wildcard \
             $(ESMALIB)/libGEOS_Shared.a \
             $(ESMALIB)/libMAPL_Base.a \
             $(ESMALIB)/libGMAO_pFIO.a \
             $(ESMALIB)/libMAPL_cfio_r4.a \
             $(ESMALIB)/libGMAO_gfio_r4.a \
             $(ESMALIB)/libGMAO_mpeu.a )

LIB_FV     = $(ESMALIB)/libfvdycore.a
LIB_SHARED += $(LIB) $(LIB_FV)

FREAL = #$(FREAL8)  # avoid the nasty -r4 and -r8 flags at all costs

$(LIB) lib : $(DEPS) $(OBJS)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)
	$(RANLIB) $(RANLIB_FLAGS) $(LIB)

StandAlone_DynAdvCore.x: StandAlone_DynAdvCore.o AppGridCreate.o LatLon2Cube.o Cube2LatLon.o $(LIB)
	$(LD) -o $@ $+ $(LDFLAGS) $(LIB_FV) \
				   $(LIB_HERMES) $(LIB_SHARED) $(LIB_ESMF) $(LIB_FLAP) \
				   $(LIB_GFDL_FMS)  \
				   $(LIB_SDF) $(LIB_MPI)  \
				   $(LIB_SYS)

StandAlone_AdvCore.x: StandAlone_AdvCore.o AppGridCreate.o LatLon2Cube.o Cube2LatLon.o $(LIB)
	$(LD) -o $@ $+ $(LDFLAGS) $(LIB_FV) \
				   $(LIB_HERMES) $(LIB_SHARED) $(LIB_ESMF) $(LIB_FLAP) \
				   $(LIB_GFDL_FMS)  \
				   $(LIB_SDF) $(LIB_MPI)  \
				   $(LIB_SYS)

StandAlone_FV3_Dycore.x: StandAlone_FV3_Dycore.o AppGridCreate.o LatLon2Cube.o Cube2LatLon.o $(LIB)
	$(LD) -o $@ $+ $(LDFLAGS) $(LIB_FV) \
				   $(LIB_HERMES) $(LIB_SHARED) $(LIB_ESMF) $(LIB_FLAP) \
				   $(LIB_GFDL_FMS)  \
				   $(LIB_SDF) $(LIB_MPI)  \
				   $(LIB_SYS)

CreateTopography.x: CreateTopography.o AppGridCreate.o $(LIB)
	$(LD) -o $@ $+ $(LDFLAGS) $(LIB_FV) \
				   $(LIB_HERMES) $(LIB_SHARED) $(LIB_ESMF) $(LIB_FLAP) \
				   $(LIB_GFDL_FMS)  \
				   $(LIB_SDF) $(LIB_MPI) \
				   $(LIB_SYS)

CreateInterpWeights.x: CreateInterpWeights.o CreateInterpWeights_GridCompMod.o AppGridCreate.o $(LIB)
	$(LD) -o $@ $+ $(LDFLAGS) $(LIB_FV) \
				   $(LIB_HERMES) $(LIB_SHARED) $(LIB_ESMF) $(LIB_FLAP) \
				   $(LIB_GFDL_FMS)  \
				   $(LIB_SDF) $(LIB_MPI)  \
				   $(LIB_SYS)

c2l_CFIO_offline.x: c2l_CFIO_offline.o AppGridCreate.o LatLon2Cube.o Cube2LatLon.o $(LIB)
	$(LD) -o $@ $+ $(LDFLAGS) $(LIB_FV) \
				   $(LIB_HERMES) $(LIB_SHARED) $(LIB_ESMF) $(LIB_FLAP) \
				   $(LIB_GFDL_FMS)  \
				   $(LIB_SDF) $(LIB_MPI)  \
				   $(LIB_SYS)

c2c.x: c2c.o AppGridCreate.o $(LIB)
	$(LD) -o $@ $+ $(LDFLAGS) $(LIB_FV) \
				   $(LIB_SHARED) $(LIB_ESMF) $(LIB_FLAP) \
				   $(LIB_GFDL_FMS)  \
				   $(LIB_SDF) $(LIB_MPI)  \
				   $(LIB_SYS)

rs_scale.x: rs_scale.o AppGridCreate.o $(LIB)
	$(LD) -o $@ $+ $(LDFLAGS) $(LIB_FV) \
				   $(LIB_SHARED) $(LIB_ESMF) $(LIB_FLAP) \
				   $(LIB_GFDL_FMS)  \
				   $(LIB_SDF) $(LIB_MPI)  \
				   $(LIB_SYS)

c2c_rst.x: c2c_rst.o AppGridCreate.o $(LIB)
	$(LD) -o $@ $+ $(LDFLAGS) $(LIB_FV) \
				   $(LIB_SHARED) $(LIB_ESMF) $(LIB_FLAP) \
				   $(LIB_GFDL_FMS)  \
				   $(LIB_SDF) $(LIB_MPI)  \
				   $(LIB_SYS)

interp_restarts.x: interp_restarts.o AppGridCreate.o LatLon2Cube.o Cube2LatLon.o $(LIB)
	$(LD) -o $@ $+ $(LDFLAGS) $(LIB_FV) \
				   $(LIB_HERMES) $(LIB_SHARED) $(LIB_ESMF) $(LIB_FLAP) \
				   $(LIB_GFDL_FMS)  \
				   $(LIB_SDF) $(LIB_MPI) \
				   $(LIB_SYS)

dyn_interp_rst.x: dyn_interp_rst.o $(LIB)
	$(LD) -o $@ $+ $(LDFLAGS)  $(LIB_FV) \
                                   $(LIB_SDF) $(LIB_MPI) $(LIB_SCI) \
                                   $(LIB_SYS) 

moist_interp_rst.x: moist_interp_rst.o 
	$(LD) -o $@ $+ $(LDFLAGS)  $(LIB_FV) \
                                   $(LIB_SDF) $(LIB_MPI) $(LIB_SCI) \
                                   $(LIB_SYS) 

#                  --------------------
#                      Dependencies
#                  --------------------

# Make sure dep files are not remade during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

#.

  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

