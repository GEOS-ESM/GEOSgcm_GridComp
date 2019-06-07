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

THIS = $(shell basename `pwd`)
LIB  = lib$(THIS).a

ALLDIRS = RRTMG

SUBDIRS = $(wildcard $(ALLDIRS))

TARGETS = esma_install esma_clean esma_distclean esma_doc \
			 install clean distclean doc

$(TARGETS): 
	@ t=$@; argv="$(SUBDIRS)" ;\
     for d in $$argv; do          \
       ( cd $$d            ;\
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


#SRCS = GEOS_IrradGridComp.F90
#OBJS = $(SRCS:.F90=.o) irradconstants.o irrad.o
#DEPS = $(OBJS:.o=.d)

SRCS = GEOS_IrradGridComp.F90 irradconstants.F90 irrad.F90
OBJS = $(addsuffix .o, $(basename $(SRCS)))
DEPS = $(addsuffix .d, $(basename $(SRCS)))

INC_DIRS = . $(INC_GMAO_SHARED) $(INC_GEOS_RADSHR) $(INC_ESMF) $(INC_GEOS_RRTMG)

ifdef GPU_TARGET
   USER_FDEFS += $(D)MAXNS=1    # Number of subgrids set at compile time for GPU
   INC_DIRS   += $(INC_PGI)     # Needed to pick up CUDA libraries
   INC_DIRS   += ../GEOS_RadiationShared # Needed for tau code
endif

# Uncomment this USER_FDEF if you wish to run the simpler, overcast cloud
# IRRAD code
#USER_FDEFS += $(D)OVERCAST

USER_FFLAGS = $(BIG_ENDIAN)
USER_FINCS  = $(foreach dir,$(INC_DIRS),$(I)$(dir))
USER_FMODS  = $(foreach dir,$(INC_DIRS),$(M)$(dir))

ifeq ($(ESMA_FC), pgfortran)
   ifneq ("$(BOPT)", "g")

      # MAT There seems to be a bug with PGI 16.5 at -O2 and higher with
      #     Irrad GC. This seems to express itself on discover only. For now
      #     just compile -O1 until a fix can be determined. The GPU seems to
      #     -O0. This is a regression back to PGI 14

      ifeq ("$(BOPT)","GPU")
         GEOS_IrradGridComp.o: FOPT = -O0 -Kieee -g $(GPU_TARGET)
         irradconstants.o: FOPT = -fast -Kieee -g
      else
         GEOS_IrradGridComp.o: FOPT = -O1 -Kieee -g $(GPU_TARGET)
      endif

   endif
endif

vpath % $(INC_DIRS)

FREAL = $(FREAL4) # for now, require 32 bit reals (R4)

$(LIB) lib : $(DEPS) $(OBJS)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)
	$(RANLIB) $(RANLIB_FLAGS) $(LIB)


# Make sure dep files are not remade during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

#.


  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

