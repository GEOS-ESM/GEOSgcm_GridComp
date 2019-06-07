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

ALLDIRS = GEOSorad_GridComp \
          GEOSoradbio_GridComp \
          GEOSoceanbiogeochem_GridComp \
	  GEOSoceanbiosimple_GridComp \
          GEOSseaice_GridComp \
          GEOSocean_GridComp

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


SRCS := GEOS_OgcmGridComp.F90 
ifeq ($(wildcard GEOSoceanbiogeochem_GridComp),$(null))
    SRCS += GEOS_OceanBioGeoChemGridCompMod___.F90
endif
ifeq ($(wildcard GEOSoceanbiosimple_GridComp),$(null))
    SRCS += GEOS_OceanBioGridCompMod___.F90
endif
ifeq ($(wildcard GEOSoradbio_GridComp),$(null))
    SRCS += GEOS_OradBioGridCompMod___.F90
endif
ifeq ($(wildcard GEOSseaice_GridComp/GEOSCICEDyna_GridComp),$(null))
    SRCS += GEOS_CICEDynaGridCompMod___.F90
endif

OBJS := $(addsuffix .o, $(basename $(SRCS))) 
DEPS := $(addsuffix .d, $(basename $(SRCS))) 

FREAL = $(FREAL4) # for now, require 32 bit reals (R4)

THIS_CFIO = MAPL_cfio_r4

INC_DIRS = . $(INC_ESMF) $(INC_GMAO_SHARED) $(ESMAINC)/GMAO_gems $(INC_ODS) $(INC_GEOS_OGCM)

MOD_DIRS = . $(INC_DIRS) \
             $(ESMAINC)/GEOSdatasea_GridComp \
             $(ESMAINC)/GEOSdataseaice_GridComp \
             $(ESMAINC)/GEOSCICEDyna_GridComp \
             $(ESMAINC)/GuestOcean_GridComp \
             $(foreach dir,$(SUBDIRS),$(ESMAINC)/$(dir))

USER_FINCS  = $(foreach dir,$(INC_DIRS),$(I)$(dir)) 
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 
USER_FFLAGS = $(BIG_ENDIAN) 

vpath % $(MOD_DIRS)

$(LIB) lib : $(OBJS)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)
	$(RANLIB) $(RANLIB_FLAGS) $(LIB)

# Create stubs
# ------------
%___.F90: $(STUB) 
	$(STUB) $(patsubst %___,%Mod,$*) > $@ 

# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

#.


  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

