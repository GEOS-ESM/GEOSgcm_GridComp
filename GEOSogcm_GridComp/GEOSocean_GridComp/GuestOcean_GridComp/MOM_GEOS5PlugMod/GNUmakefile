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
	@echo "       MOMDIR = $(MOMDIR)"
	@echo "         ARCH = $(ARCH)"
	@echo "         SITE = $(SITE)"

THIS := $(shell basename `pwd`)
LIB   = lib$(THIS).a

SUBDIRS = mom

TARGETS = esma_install esma_clean esma_distclean esma_doc \
          install clean distclean doc

export ESMADIR BASEDIR MOMDIR ARCH SITE

$(TARGETS): mom-stamp
	@ t=$@; argv="$(SUBDIRS)" ;\
	  for d in $$argv; do			 \
	    ( cd $$d				;\
	      echo ""; echo Making $$t in `pwd`          ;\
	      $(MAKE) -e $$t ) \
	  done
	$(MAKE) local_$@

esma_install install:

local_esma_install local_install: $(LIB)
	$(MKDIR) $(ESMALIB) $(ESMAINC)/$(THIS)
	$(CP) -p *.a            $(ESMALIB)
	$(CP) -p *.mod          $(ESMAINC)/$(THIS)
#	$(CP) -p GEOS_ErrLog.h  $(ESMAINC)/$(THIS)

local_esma_clean local_clean:
	-$(RM) *~ *.[aox] *.mod *.x

local_esma_distclean local_distclean:
	-$(RM) *~ *.[aoxd] *.mod *.x

local_esma_doc local_doc:
	@$(PROTEX) $(PROTEX_FLAGS) *GridComp*.[fF]* > $(ESMADOC)/$(THIS).tex


#                  --------------------
#                  User Defined Targets
#                  --------------------


SRCS := MOM_GEOS5PlugMod.F90
USER_FDEFS += $(M)$(ESMAINC)/GMAO_gems $(M)$(ESMAINC)/GMAOodas_GridComp

ifneq ( $(wildcard seawater.F90), $(null) ) 
   SRCS += seawater.F90
endif


GITBRANCH = geos5
GITREPOSITORY = $(MOMDIR)/mom

OBJS := $(addsuffix .o, $(basename $(SRCS))) 
DEPS := $(addsuffix .d, $(basename $(SRCS))) 

FREAL = $(FREAL4) # for now, require 32 bit reals (R4)

INC_DIRS = . $(ESMAINC)/MAPL $(ESMAINC)/MAPL_Base 
MOD_DIRS = . $(INC_ESMF) \
             $(INC_GMAO_SHARED) \
             $(INC_GFDL_FMS) \
             $(ESMAINC)/GEOS_Shared $(ESMAINC)/MAPL_Base \
             $(foreach dir,$(SUBDIRS),$(ESMAINC)/$(dir))

USER_FINCS  = $(foreach dir,$(INC_DIRS),$(I)$(dir)) 
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 
#USER_FFLAGS = $(BIG_ENDIAN) -i4 -r8 
USER_FFLAGS = $(BIG_ENDIAN) $(FREAL8) 

vpath % $(MOD_DIRS) $(ESMAINC)/GMAO_gems $(ESMAINC)/GMAOodas_GridComp

$(LIB) lib : $(OBJS)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)


mom-stamp:
	git clone -b $(GITBRANCH) $(GITREPOSITORY)
	touch $@

# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

#.


#
# Earth System Modeling Applications (ESMA) post processing makefile fragment.
# Optionally, include this at the end of your makefile.
#
# Notes:
# 1. Building the "pinstall" target in a directory does the following:
#    a) If the directory has SUBDIRS, then
#       - each SUBDIR is built with "gmake pinstall"
#       - the directory is then built with the "local_install" target
#    b) If the directory has no SUBDIRS, then 
#       - the directory is built with the "install" target
# 2. Use $(subdir)_install dependency statements in the GNUmakefile
#    to indicate subdirectory dependencies.
#
# REVISION HISTORY:
#
# 15Dec2006  da Silva  First Crack
# 02Feb2010  Stassi    Added notes and made small changes to the code structure
#                      for clarity (i.e. to help me remember how it works).
#--------------------------------------------------------------------------

#                     -----------------
#                     Parallel Install
#                     -----------------

ifdef SUBDIRS
	PINSTALL_DIRS = $(foreach dir,$(SUBDIRS),$(dir)_install) 
	LOCAL_INSTALL = local_install
	PINSTALL_TARGET = $(LOCAL_INSTALL)
else
	PINSTALL_TARGET = install
endif

pinstall: $(PINSTALL_DIRS) 
	$(MAKE) -e $(PINSTALL_TARGET)

%_install: mom-stamp
	@$(ESMA_TIMER_BEG)
	$(MAKE) -e -C $* pinstall
	@$(ESMA_TIMER_END)

pinstall_skip:
	@echo "Skipping local_install in `pwd`"

