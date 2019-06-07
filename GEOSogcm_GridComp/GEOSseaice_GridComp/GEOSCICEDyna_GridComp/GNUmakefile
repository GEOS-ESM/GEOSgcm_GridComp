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
	$(CP) -p *.mod          $(ESMAINC)/$(THIS)
#	$(CP) -p GEOS_ErrLog.h  $(ESMAINC)/$(THIS)

esma_clean clean:
	-$(RM) *~ *.[aox] *.mod *.x

esma_distclean distclean:
	-$(RM) *~ *.[aoxd] *.mod *.x

esma_doc doc:
	@$(PROTEX) $(PROTEX_FLAGS) *GridComp*.[fF]* > $(ESMADOC)/$(THIS).tex


#                  --------------------
#                  User Defined Targets
#                  --------------------


SRCS = $(shell echo *.F90)
OBJS = $(SRCS:.F90=.o)
DEPS := $(addsuffix .d, $(basename $(SRCS)))

FREAL = $(FREAL4) # for now, require 32 bit reals (R4)

INC_DIRS = . $(INC_GMAO_SHARED) $(INC_ESMF) $(ESMAINC)/LANL_cice $(ESMAINC)/GMAO_gems $(ESMAINC)/GMAOodas_GridComp
MOD_DIRS =   $(INC_DIRS)

USER_FINCS  = $(foreach dir,$(INC_DIRS),$(I)$(dir)) 
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 
USER_FFLAGS = $(BIG_ENDIAN) -DDIAGOUT -DMODIFY_TOPOGRAPHY -DUSE_R8

ifeq ($(findstring iodas,$(options)),iodas)
   USER_FDEFS += -DUSE_ODAS
   FOPT += -fpe0   
endif

vpath % $(MOD_DIRS)

$(LIB) lib : $(OBJS) $(DEPS)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)


# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

#.


  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

