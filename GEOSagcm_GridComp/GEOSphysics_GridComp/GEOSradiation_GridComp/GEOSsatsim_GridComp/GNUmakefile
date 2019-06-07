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

#ALLDIRS = actsim
#SUBDIRS = $(wildcard $(ALLDIRS))

TARGETS = esma_install esma_clean esma_distclean esma_doc \
          install clean distclean doc 

export ESMADIR BASEDIR ARCH SITE

#$(TARGETS): 
#	@ t=$@; argv="$(SUBDIRS)" ;\
#	  for d in $$argv; do			 \
#	    ( cd $$d				;\
#	      echo ""; echo Making $$t in `pwd`          ;\
#	      $(MAKE) -e $$t ) \
#	  done
#	$(MAKE) local_$@

local_esma_install local_install install: $(LIB)
	$(MKDIR) $(ESMALIB) $(ESMAINC)/$(THIS)
	$(CP) -p *.a            $(ESMALIB)
	$(CP) -p *.mod          $(ESMAINC)/$(THIS)
#	$(CP) -p GEOS_ErrLog.h  $(ESMAINC)/$(THIS)

local_esma_clean local_clean clean esma_clean:
	-$(RM) *~ *.[aox] *.[Mm][Oo][Dd]

local_esma_distclean local_distclean distclean esma_distclean:
	-$(RM) *~ *.[aoxd] *.[Mm][Oo][Dd]

esma_doc doc:
	@$(PROTEX) $(PROTEX_FLAGS) *GridComp*.[fF]* > $(ESMADOC)/$(THIS).tex


#                  --------------------
#                  User Defined Targets
#                  --------------------

SRC_DIRS := . \
            actsim \
            llnl \
            quickbeam

UNSRC = ./congvec.f quickbeam/load_mie_table.f90

SRCS := $(filter-out $(UNSRC), \
        $(foreach dir,$(SRC_DIRS), \
        $(wildcard $(dir)/*.[fFc]) $(wildcard $(dir)/*.[fF]90) ) )

# notdir seems to be the trick to moving the object files to the gridded comp directory
OBJS := $(notdir $(addsuffix .o, $(basename $(SRCS))))
DEPS := $(notdir $(addsuffix .d, $(basename $(SRCS))))


FREAL = $(FREAL4) # for now, require 32 bit reals (R4)

INC_DIRS = . $(INC_GMAO_SHARED) $(INC_GEOS_RADSHR) $(INC_ESMF) $(SRC_DIRS)

MOD_DIRS = . $(INC_DIRS) \
             $(foreach dir,$(SUBDIRS),$(ESMAINC)/$(dir))

USER_FFLAGS = $(BIG_ENDIAN) $(PP)
USER_FINCS  = $(foreach dir,$(INC_DIRS),$(I)$(dir)) 
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 

vpath % $(MOD_DIRS)

$(LIB) lib : $(OBJS)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)
	$(RANLIB) $(RANLIB_FLAGS) $(LIB)

#              ------------------------------------------
#              Package Dependencies for Parallel Install
#              ------------------------------------------
#  GEOSsolar_GridComp_install : GEOS_RadiationShared_install 


# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros
#.
