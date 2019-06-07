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

THIS := $(shell basename `pwd` )

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

LIB  = lib$(THIS).a
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
	@$(PROTEX) $(PROTEX_FLAGS) GEOS_CICE4ColumnPhysGridComp.F90 > $(ESMADOC)/GEOS_CICE4ColumnPhysGridComp.tex
	@$(PROTEX) $(PROTEX_FLAGS) GEOS_OpenWaterGridComp.F90 > $(ESMADOC)/GEOS_OpenWaterGridComp.tex
	@$(PROTEX) $(PROTEX_FLAGS) GEOS_SaltWaterGridComp.F90 > $(ESMADOC)/GEOS_SaltWaterGridComp.tex
	@$(PROTEX) $(PROTEX_FLAGS) GEOS_SimpleSeaiceGridComp.F90 > $(ESMADOC)/GEOS_SimpleSeaiceGridComp.tex



#                  --------------------
#                  User Defined Targets
#                  --------------------


SRCS = $(shell echo *.F90)
OBJS := $(addsuffix .o, $(basename $(SRCS)))
DEPS := $(addsuffix .d, $(basename $(SRCS)))

INC_DIRS = . $(INC_GMAO_SHARED) $(INC_LANL_CICE) $(INC_ESMF)
MOD_DIRS = . $(INC_DIRS) $(INC_ESMF) \
             $(foreach dir,$(SUBDIRS),$(ESMAINC)/$(dir))

USER_FINCS  = $(foreach dir,$(INC_DIRS),$(I)$(dir))
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir))
USER_FFLAGS = $(BIG_ENDIAN)

FREAL = $(FREAL4) # for now, require 32 bit reals (R4)

GEOS_SaltWaterGridComp.o: GEOS_OpenWaterGridComp.o GEOS_SimpleSeaiceGridComp.o GEOS_CICE4ColumnPhysGridComp.o 

ifeq ($(ESMA_FC), pgfortran)
   ifneq ("$(BOPT)", "g")

      # MAT PGI 16.5 is much too slow with the Salt GC. Reduce optimiation.

      GEOS_SaltWaterGridComp.o: FOPT = -O -Kieee -g $(GPU_TARGET)

   endif
endif
vpath % $(MOD_DIRS)

$(LIB) lib : $(OBJS)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)
	$(RANLIB) $(RANLIB_FLAGS) $(LIB)


#.


  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

