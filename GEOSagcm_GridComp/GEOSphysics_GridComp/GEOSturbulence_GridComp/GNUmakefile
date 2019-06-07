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
	@echo "         SRCS = $(SRCS)"

THIS = $(shell basename `pwd`)
LIB  = lib$(THIS).a
NAME = $(shell basename $(THIS) _GridComp)


esma_install install: $(LIB)
	$(MKDIR) $(ESMALIB) $(ESMAINC)/$(THIS)
	$(CP) -p *.a            $(ESMALIB)
	$(CP) -p *.mod          $(ESMAINC)/$(THIS)
#	$(CP) -p GEOS_ErrLog.h  $(ESMAINC)/$(THIS)

esma_clean clean:
	-$(RM) *~ *.[aox] *.mod *.x

esma_distclean distclean:
	-$(RM) *~ *.[aoxd] *.mod *.x

#                  --------------------
#                  User Defined Targets
#                  --------------------

SRCS = GEOS_TurbulenceGridComp.F90 LockEntrain.F90 shoc.F90
OBJS = $(SRCS:.F90=.o) 
DEPS = $(OBJS:.o=.d) 

INC_DIRS = . $(INC_GMAO_SHARED) $(INC_ESMF)
MOD_DIRS = . $(INC_DIRS)

ifdef GPU_TARGET
   INC_DIRS += $(INC_PGI)                             # Needed to pick up CUDA libraries
   INC_DIRS += $(ESMADIR)/src/GMAO_Shared/GEOS_Shared # Needed for internal DQSAT
endif

USER_FFLAGS = $(BIG_ENDIAN) 
USER_FINCS  = $(foreach dir,$(INC_DIRS),$(I)$(dir))
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 

vpath % $(INC_DIRS) 

FREAL = $(FREAL4) # for now, require 32 bit reals (R4)

$(LIB) lib : $(DEPS) $(OBJS)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)
	$(RANLIB) $(RANLIB_FLAGS) $(LIB)

GGUG_SRCS = $(SRCS)
esma_doc doc:
	@$(PROTEX) $(PROTEX_FLAGS) *GridComp*.[fF]* > $(ESMADOC)/$(THIS).tex

# Make sure dep files are not remade during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

#.












  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

