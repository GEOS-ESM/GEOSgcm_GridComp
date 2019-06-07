#
# Makefile for ESMA components.
#
#
# REVISION HISTORY:
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

THIS = ARIESg3_GridComp
LIB  = lib$(THIS).a

esma_install install: $(LIB)
	$(MKDIR) $(ESMALIB) $(ESMAETC) $(ESMAINC)/$(THIS)
	$(CP) -p *.a    $(ESMALIB)
	$(CP) -p *.mod  $(ESMAINC)/$(THIS)
#	$(CP) -p *.rc   $(ESMAETC)

esma_clean esma_distclean clean distclean:
	-$(RM) *~ *.[aoxd] *.[Mm][Oo][Dd] *.x

esma_doc doc:
	@$(PROTEX) $(PROTEX_FLAGS) *GridComp*.[fF]* > $(ESMADOC)/$(THIS).tex


#                  --------------------
#                  User Defined Targets
#                  --------------------

SRCS = g3_dynamics_lattice_module.F g3_dynamics_state_module.F \
       g3_mpi_util.F g3_mymalloc.F g3_fft.F \
       g3_wrapper.F g3_dycore.F g3_gcmutil.F g3_grids.F \
       shr_kind_mod.F90 dynamics_vars.F90 \
       pft_module.F90 par_xsum.F90 fft99.F90 remap.F90 gmap.F90 \
       $(THIS)Mod.F90

OBJS := $(addsuffix .o, $(basename $(SRCS)))
DEPS := $(addsuffix .d, $(basename $(SRCS)))

#OBJS = $(SRCS:.F90=.o) 
#DEPS = $(SRCS:.F90=.d) 

INC_DIRS = . $(INC_GMAO_SHARED) $(INC_ESMF) $(INC_MPI)
MOD_DIRS = . $(INC_DIRS)

USER_FINCS  = $(foreach dir,$(INC_DIRS),$(I)$(dir))
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 
vpath % $(INC_DIRS) $(MOD_DIRS)

USER_FFLAGS += $(BIG_ENDIAN) 
USER_FDEFS  = -DSPMD -DMAPL_MODE -DTWOD_YZ -DrFV=r4

PPAWK=./pp.awk
FREAL = # $(FREAL8)  # avoid the nasty -r4 and -r8 flags at all costs

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

