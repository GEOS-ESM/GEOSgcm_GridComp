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
  include $(ESMADIR)/Config/GMAO_base.mk  # Generic stuff

# include $(ESMADIR)/Config/vectorize.mk  # Vectorize

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

esma_install install: $(LIB)
	$(MKDIR) $(ESMALIB) $(ESMAINC)/$(THIS)
	$(CP) -p *.a            $(ESMALIB)
	$(CP) -p *.[Mm][Oo][Dd] $(ESMAINC)/$(THIS)
#	$(CP) -p GEOS_ErrLog.h  $(ESMAINC)/$(THIS)

esma_clean clean: 
	-$(RM) *~ *.[aox] *.[Mm][Oo][Dd] 

esma_distclean distclean:
	-$(RM) *~ *.[aoxd] *.[Mm][Oo][Dd] 

esma_doc doc:
	@$(PROTEX) $(PROTEX_FLAGS) *GridComp*.[fF]* > $(ESMADOC)/$(THIS).tex


#                  --------------------
#                  User Defined Targets
#                  --------------------


SRCS := $(wildcard gw_drag.F90 gw_drag_v2.F90 \
                   machine.F90 gwdps.f gwdc.f \
                   GEOS_GwdGridComp.F90 )
OBJS := $(addsuffix .o, $(basename $(SRCS))) 
DEPS := $(addsuffix .d, $(basename $(SRCS))) 

FREAL = $(FREAL4) # for now, require 32 bit reals (R4)

INC_DIRS = . $(INC_GMAO_SHARED) $(INC_ESMF)
MOD_DIRS = . $(INC_DIRS)

ifdef GPU_TARGET
   USER_FDEFS += $(D)MAXPGWV=4  # Maximum number of waves allowed set at compile time for GPU
   INC_DIRS   += $(INC_PGI)     # Needed to pick up CUDA libraries
endif

USER_FINCS  = $(foreach dir,$(INC_DIRS),$(I)$(dir)) 
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 
USER_FFLAGS = $(BIG_ENDIAN)

vpath % $(MOD_DIRS)

ifeq ($(ESMA_FC),ifort)
.f.o:
	$(FC) $(fFLAGS) -fixed -noextend-source -c $<
#	$(FC) -c $(F90FLAGS) -nofixed $<
else
ifeq ($(ESMA_FC),gfortran)
.f.o:
	$(FC) $(fFLAGS) -ffixed-line-length-72 -c $<
else
ifeq ($(ESMA_FC),pgfortran)
.f.o:
	$(FC) $(patsubst $(EXTENDED_SOURCE),$(FIXED_SOURCE),$(fFLAGS)) -c $<
endif
endif
endif

$(LIB) lib : $(OBJS)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)
	$(RANLIB) $(RANLIB_FLAGS) $(LIB)

# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

#.


  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

