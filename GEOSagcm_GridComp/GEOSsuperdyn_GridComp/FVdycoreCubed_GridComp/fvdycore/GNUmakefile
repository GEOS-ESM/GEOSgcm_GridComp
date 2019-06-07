#
# Makefile for ESMA components.
#
#
# REVISION HISTORY:
#
# 09Jun2003  da Silva  First crack.
# 07Aug2003  Sawyer    Modified for FVdycore
# 27Aug2003  Sawyer    Added a2d3d and d2a3d
# 22Sep2004  Sawyer    Modifications for merge with current code
# 06Oct2004  Sawyer    Removed spmd_dyn
# 17Feb2005  Sawyer    Added epvd and ppme
# 17May2005  Sawyer    Added FVdycore_wrapper
# 18oct2005  da Silva  Moved ALTIX specific flags to FVdycore.mk
# 18Jan2006  Putman    Added mfz_comp
# 17May2006  Sawyer    Added diag_module.F90
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
  include               FVdycore_arch.mk  # arch dependent flags 

  include $(ESMADIR)/Config/vectorize.mk  # Vectorize

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
LIB  = lib$(THIS).a

esma_install install: $(LIB)
	$(MKDIR) $(ESMALIB) $(ESMAETC) $(ESMAINC)/$(THIS)
	$(CP) -p *.a            $(ESMALIB)
	$(CP) -p *.[Mm][Oo][Dd] $(ESMAINC)/$(THIS)
#	$(CP) -p *.rc   $(ESMAETC)

esma_clean clean:
	-$(RM) *~ *.[aox] *.[Mm][Oo][Dd] 

esma_distclean distclean:
	-$(RM) *~ *.[aoxd] *.[Mm][Oo][Dd]

esma_doc doc:
	@ echo "nothing to do here"
	#@$(PROTEX) $(PROTEX_FLAGS) *GridComp*.[fF]* > $(ESMADOC)/$(THIS).tex


#                  --------------------
#                  User Defined Targets
#                  --------------------

SRCS := $(shell cat path_names)
OBJS := $(notdir $(addsuffix .o, $(basename $(SRCS)))) 
DEPS := $(notdir $(addsuffix .d, $(basename $(SRCS))))

DEPS := $(filter-out nsclock.d threadloc.d, $(DEPS)) # filter these out

INC_DIRS = . $(INC_GFDL_FMS) \
            $(INC_GMAO_SHARED) $(INC_SDF) $(INC_ESMF) $(INC_MPI)
MOD_DIRS = . $(INC_DIRS) 
SRC_DIRS = $(sort $(dir $(SRCS)))

USER_FINCS  = $(foreach dir,$(INC_DIRS),$(I)$(dir))
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 

vpath % $(SRC_DIRS) $(INC_DIRS) $(MOD_DIRS)

USER_FDEFS   = $(D)MAPL_MODE $(D)SPMD $(D)TIMING

# If Non-Hydrostatic FV3 uncomment this line
#USER_FDEFS  += $(D)USE_COND $(D)MOIST_CAPPA

USER_CDEFS   = $(USER_FDEFS)

# The FV_PRECISION variable is exported from the parent
# SINGLE PRECISION FV_CUBED
ifeq ($(FV_PRECISION),R4)
   USER_FDEFS += $(D)SINGLE_FV $(D)OVERLOAD_R4
   THIS_GFDL_FMS = GFDL_fms_r4
   FREAL= $(FREAL4)
else
   THIS_GFDL_FMS = GFDL_fms_r8
   FREAL= $(FREAL8)
endif

#ifeq ($(IFORT_MAJOR),15)
#   USER_FFLAGS += -nowarn -sox -align array64byte -xCORE-AVX2 -qno-opt-dynamic-align -qoverride-limits -qopt-prefetch=3
#endif

#ifeq ($(HYDRO),Y)
#CPPDEFS += -Duse_libMPI -Duse_netCDF -DSPMD -DUSE_LOG_DIAG_FIELD_INFO -Duse_LARGEFILE -DUSE_GFSL63 -DGFS_PHYS
#else
#CPPDEFS += -Duse_libMPI -Duse_netCDF -DSPMD -DUSE_LOG_DIAG_FIELD_INFO -Duse_LARGEFILE -DUSE_GFSL63 -DGFS_PHYS -DMOIST_CAPPA -DUSE_COND
#endif

$(LIB) lib : $(DEPS) $(OBJS)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)
	$(RANLIB) $(RANLIB_FLAGS) $(LIB)

#                  --------------------
#                      Dependencies
#                  --------------------

# Make sure dep files are not remade during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

#.

  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros



