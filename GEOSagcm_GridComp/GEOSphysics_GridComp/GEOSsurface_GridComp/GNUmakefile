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

ALLDIRS = Shared \
          GEOSlandice_GridComp \
          GEOSlake_GridComp \
          GEOSland_GridComp  \
          GEOSsaltwater_GridComp

SUBDIRS = $(wildcard $(ALLDIRS))

TARGETS = esma_install esma_clean esma_distclean esma_doc \
          install clean distclean doc 

export ESMADIR BASEDIR ARCH SITE

$(TARGETS): 
	@ t=$@; argv="$(SUBDIRS)" ;\
	  for d in $$argv; do			 \
	    ( cd $$d				;\
	      echo ""; echo Making $$t in `pwd`          ;\
	      $(MAKE) -e $$t ) \
	  done
	$(MAKE) local_$@

ifneq ("$(wildcard GEOS_SurfaceGridComp.F90)","")
   local_esma_install local_install: $(LIB)
	$(MKDIR) $(ESMALIB) $(ESMAINC)/$(THIS)
	$(CP) -p *.a            $(ESMALIB)
	$(CP) -p *.mod          $(ESMAINC)/$(THIS)
#	$(CP) -p GEOS_ErrLog.h  $(ESMAINC)/$(THIS)
else
  local_esma_install local_install: ;
endif

local_esma_clean local_clean:
	-$(RM) *~ *.[aox] *.mod *.x *.ps

local_esma_distclean local_distclean:
	-$(RM) *~ *.[aoxd] *.mod *.x *.ps

local_esma_doc local_doc:
	echo ""; echo Making *GridComp*.[fF]* in `pwd`          ;\
	$(PROTEX) $(PROTEX_FLAGS) *GridComp*.[fF]* > $(ESMADOC)/$(THIS).tex


#                  --------------------
#                  User Defined Targets
#                  --------------------


SRCS := GEOS_SurfaceGridComp.F90
OBJS := $(addsuffix .o, $(basename $(SRCS)))
DEPS := $(addsuffix .d, $(basename $(SRCS)))

FREAL = $(FREAL4) # for now, require 32 bit reals (R4)

INC_DIRS = . $(INC_GMAO_SHARED) $(INC_ESMF) $(INC_GEOS_SURF) $(INC_MPI)
MOD_DIRS =   $(INC_DIRS)

USER_FINCS  = $(foreach dir,$(INC_DIRS),$(I)$(dir)) 
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 
USER_FFLAGS = $(BIG_ENDIAN) 

ifeq ($(findstring obio,$(options)),obio)
   USER_FDEFS += -DUSE_OBIO
endif

vpath % $(MOD_DIRS)

$(LIB) lib : $(OBJS) $(DEPS)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)
	$(RANLIB) $(RANLIB_FLAGS) $(LIB)


# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

#              ------------------------------------------
#              Package Dependencies for Parallel Install
#              ------------------------------------------

     GEOSland_GridComp_install : Shared_install
  GEOSlandice_GridComp_install : Shared_install


  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros
