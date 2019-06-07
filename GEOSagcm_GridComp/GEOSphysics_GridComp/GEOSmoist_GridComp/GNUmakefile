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

THIS = GEOSmoist_GridComp
LIB  = lib$(THIS).a

esma_install install: $(LIB)
	$(MKDIR) $(ESMALIB) $(ESMAINC)/$(THIS)
	$(CP) -p *.a            $(ESMALIB)
	$(CP) -p *.mod          $(ESMAINC)/$(THIS)
	$(CP) -p *.rc           $(ESMAETC)
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


SRCS := $(wildcard CLDPARAMS.F90 RASPARAMS.F90 SHLWPARAMS.F90 ras.F90 ras00.F90 ConvPar_GF_GEOS5.F90 module_gate.F90 uwshcu.F90 cloud.F90 cloudnew.F90 ddf.F90 rascnvv2_v.F90 \
          cldmacro.F90 cldwat2m_micro.F90 wv_saturation.F90 aer_cloud.F90 aer_actv_single_moment.F90 gfdl_cloud_microphys.F90 \
          GEOS_MoistGridComp.F90 )

OBJS := $(addsuffix .o, $(basename $(SRCS))) 
DEPS := $(addsuffix .d, $(basename $(SRCS))) 

INC_DIRS = . $(INC_GMAO_SHARED)
MOD_DIRS = . $(INC_DIRS) $(INC_ESMF) \
             $(foreach dir,$(SUBDIRS),$(ESMAINC)/$(dir))

ifdef GPU_TARGET
   INC_DIRS += $(INC_PGI)                             # Needed to pick up CUDA libraries
   INC_DIRS += $(ESMADIR)/src/GMAO_Shared/GEOS_Shared # Needed for internal DQSAT
endif

USER_FFLAGS = $(BIG_ENDIAN) 
USER_FINCS  = $(foreach dir,$(INC_DIRS),$(I)$(dir)) 
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir))

FREAL = $(FREAL4) # for now, require 32 bit reals (R4)

ifeq ($(ESMA_FC), pgfortran)
   ifneq ("$(BOPT)", "g")

      # MAT PGI 16.5 is pathologically slow with Moist GC. More than 10 minutes.

      GEOS_MoistGridComp.o: FOPT = -O -Kieee -g $(GPU_TARGET)

   endif
endif


vpath % $(MOD_DIRS)

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

