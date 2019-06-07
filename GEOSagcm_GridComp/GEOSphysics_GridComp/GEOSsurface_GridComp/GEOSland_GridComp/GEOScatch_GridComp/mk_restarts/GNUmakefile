SHELL=/bin/sh

#.......................................................................
#
# Makefile for ESMA components
#
# REVISION HISTORY:
#
#  13Mar2012  jstassi  Created with gmm script.
#
#.......................................................................

#---------------------------------------------------
# Get GNUmakefile name
# note: MAKEFILE_LIST defined for gmake versions >= 3.80
#---------------------------------------------------
ifdef MAKEFILE_LIST
   GNUMAKEFILE := $(word $(words $(MAKEFILE_LIST)),$(MAKEFILE_LIST))
else
   GNUMAKEFILE := GNUmakefile
endif

#------------------------------
# Make sure ESMADIR is defined
# -----------------------------
ifndef ESMADIR
   ESMADIR := $(PWD)/../../../../../../../../..
endif

#-------------------------------
# Compilation rules, flags, etc
#-------------------------------
include $(ESMADIR)/Config/ESMA_base.mk   # Generic stuff
include $(ESMADIR)/Config/ESMA_arch.mk   # System dependencies
include $(ESMADIR)/Config/GMAO_base.mk   # System dependencies

########################################################################

LIB_GEOS_SURFSHR = $(ESMALIB)/libGEOS_SurfaceShared.a

INC_LIST = $(INC_MAPL_BASE) $(INC_NETCDF) $(INC_GEOS_ROUTE) $(INC_GEOS_LANDSHR) \
           $(INC_GEOS_SURFSHR)
LIB_LIST = $(LIB_MAPL_BASE) \
           $(LIB_CFIO) $(LIB_GFIO) $(LIB_MPEU) $(LIB_GEOS_ROUTE) $(LIB_GEOS_SURFSHR) \
           $(LIB_ESMF) $(LIB_SDF) $(LIB_MPI)
OMP_NUM_THREADS:= $(shell /usr/bin/lscpu | grep 'CPU(s)' | cut -d ':' -f2)
export OMP_NUM_THREADS
FOPT += $(OMPFLAG)

########################################################################

BINS = Scale_Catch \
       Scale_CatchCN \
       cv_SaltRestart\
       SaltIntSplitter \
       SaltImpConverter\
       mk_CICERestart \
       mk_CatchCNRestarts\
       mk_CatchRestarts \
       mk_LDASsaRestarts\
       mk_LakeLandiceSaltRestarts \
       mk_RouteRestarts \
       newcatch \
       strip_vegdyn

UT =

########################################################################

INCFILES = $(wildcard *.[Mm][Oo][Dd]) $(wildcard *.[Hh])

NMLFILES =

PERLFILES = $(wildcard *.pl) $(wildcard *.pm) mk_Restarts

RCFILES =

SHELLSCRIPTS = $(wildcard *.csh)

SEDCMD = s^\@DASPERL^\$(PERL)^;s^\@DASSED^\$(SED)^;s^\@DASHELL^\$(SHELL)^


########################################################################

ALLDIRS =

SUBDIRS = $(wildcard $(ALLDIRS))

########################################################################

THIS = $(shell basename `pwd`)
ALLSRCS =
SRCS = $(wildcard $(ALLSRCS))

OBJS := $(addsuffix .o, $(basename $(SRCS)))
DEPS := $(addsuffix .d, $(basename $(SRCS)))

ALLDEPS := $(DEPS) $(addsuffix .d, $(basename $(BINS)))

ifneq ("$(SRCS)", "")
  LIB = lib$(THIS).a
else
  LIB = 
endif

########################################################################

#-----------------------
# Standard ESMA Targets
#-----------------------
esma_doc doc:
	@echo no documentation here

esma_install install:
	@echo
	@if [ "$(SUBDIRS)" != "" ]; then ( \
	   SUBTARGET=$@; export SUBTARGET; $(MAKE) -f $(GNUMAKEFILE) subdirs \
	) fi

	$(MAKE) -f $(GNUMAKEFILE) local

local_install local:
	$(MKDIR) $(ESMABIN) $(ESMAETC) $(ESMALIB)
	$(MAKE) -f $(GNUMAKEFILE) bins
	$(MAKE) -f $(GNUMAKEFILE) etc
	$(MAKE) -f $(GNUMAKEFILE) inc
	$(MAKE) -f $(GNUMAKEFILE) library
	$(MAKE) -f $(GNUMAKEFILE) scripts

esma_clean clean:
	$(RM) *~ $(BINS) *.[aox] *.[Mm][Oo][Dd]

	@if [ "$(SUBDIRS)" != "" ]; then ( \
	   SUBTARGET=$@; export SUBTARGET; $(MAKE) -f $(GNUMAKEFILE) subdirs \
	) fi

esma_distclean distclean:
	$(RM) *~ $(BINS) *.[aoxd] *.[Mm][Oo][Dd]

	@if [ "$(SUBDIRS)" != "" ]; then ( \
	   SUBTARGET=$@; export SUBTARGET; $(MAKE) -f $(GNUMAKEFILE) subdirs \
	) fi

esma_help help:
	@echo "Standard ESMA targets:"
	@echo "% gmake install       (builds and install under ESMADIR)"
	@echo "% gmake clean         (removes deliverables: *.[aox], etc)"
	@echo "% gmake distclean     (leaves in the same state as cvs co)"
	@echo "% gmake esma_doc      (generates PDF, installs under ESMADIR)"
	@echo "% gmake help          (this message)"
	@echo
	@echo "Subtargets:"
	@echo "% gmake bins          (compiles and installs binary files)"
	@echo "% gmake csh           (installs shell scripts)"
	@echo "% gmake etc           (installs rc and nml files)"
	@echo "% gmake inc           (installs include files)"
	@echo "% gmake library       (builds and installs library files)"
	@echo "% gmake nml           (installs nml files)"
	@echo "% gmake perl          (installs perl scripts and packages)"
	@echo "% gmake rc            (installs rc files)"
	@echo "% gmake scripts       (installs shell scripts and perl files)"
	@echo "% gmake subdirs       (builds subdirectories)"
	@echo "% gmake ut            (compiles utility programs)"
	@echo
	@echo "Environment:"
	@echo "      ESMADIR = $(ESMADIR)"
	@echo "      BASEDIR = $(BASEDIR)"
	@echo "         ARCH = $(ARCH)"
	@echo "         SITE = $(SITE)"

########################################################################

bins: $(BINS)
	@if [ "$(BINS)" != "" ]; then ( \
	   echo "................................."; \
	   for file in $(BINS); do \
	      ( echo "copy to \$$(ESMABIN): $$file"; \
	        $(CP) $$file $(ESMABIN) ) \
	   done; \
	   echo "\$$(ESMABIN) = $(ESMABIN)"; \
	   echo "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"; \
	) fi

csh:
	@if [ "$(SHELLSCRIPTS)" != "" ]; then ( \
	   echo "................................."; \
	   for scr in $(SHELLSCRIPTS); do \
	      ( $(SED) -e "$(SEDCMD)"< $$scr > $(ESMABIN)/$$scr; \
	        echo "script to \$$(ESMABIN): $$scr"; \
	        chmod 755 $(ESMABIN)/$$scr ) \
	   done; \
	   echo "\$$(ESMABIN) = $(ESMABIN)"; \
	   echo "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"; \
	) fi

etc: nml rc

inc:
	@if [ "$(INCFILES)" != " " ]; then ( \
	   $(MKDIR) $(ESMAINC)/$(THIS); \
	   echo "................................."; \
	   for file in $(INCFILES); do \
	      ( echo "copy to \$$(ESMAINC)/\$$(THIS): $$file"; \
	        $(CP) $$file $(ESMAINC)/$(THIS) ) \
	   done; \
	   echo "\$$(ESMAINC)/\$$(THIS) = $(ESMAINC)/$(THIS)"; \
	   echo "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"; \
	) fi

library: $(LIB)
	@if [ "$(LIB)" != "" ]; then ( \
	   echo "................................."; \
	   for file in $(LIB); do \
	      ( echo "library to \$$(ESMALIB): $$file"; \
	        $(CP) $$file $(ESMALIB) ) \
	   done; \
	   echo "\$$(ESMALIB) = $(ESMALIB)"; \
	   echo "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"; \
	) fi

nml:
	@if [ "$(NMLFILES)" != "" ]; then ( \
	   echo "................................."; \
	   for file in $(NMLFILES); do \
	      ( echo "copy to \$$(ESMAETC): $$file"; \
	        $(CP) $$file $(ESMAETC) ) \
	   done; \
	   echo "\$$(ESMAETC) = $(ESMAETC)"; \
	   echo "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"; \
	) fi

perl:
	@if [ "$(PERLFILES)" != " " ]; then ( \
	   echo "................................."; \
	   for scr in $(PERLFILES); do \
	      ( $(SED) -e "$(SEDCMD)"< $$scr > $(ESMABIN)/$$scr; \
	        echo "script to \$$(ESMABIN): $$scr"; \
	        chmod 755 $(ESMABIN)/$$scr ) \
	   done; \
	   echo "\$$(ESMABIN) = $(ESMABIN)"; \
	   echo "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"; \
	) fi

rc:
	@if [ "$(RCFILES)" != "" ]; then ( \
	   echo "................................."; \
	   for file in $(RCFILES); do \
	      ( echo "copy to \$$(ESMAETC): $$file"; \
	        $(CP) $$file $(ESMAETC) ) \
	   done; \
	   echo "\$$(ESMAETC) = $(ESMAETC)"; \
	   echo "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"; \
	) fi

scripts: csh perl

subdirs:
	@for subdir in $(SUBDIRS); do \
	   ( echo "> cd $$subdir; $(MAKE) -e $(SUBTARGET)"; \
	     cd $$subdir; $(MAKE) -e $(SUBTARGET) ) \
	done; \

ut UT:
	@if [ "$(UT)" != "" ]; then ( \
	   for prog in $(UT); do \
	      ( echo ""; echo "making $$prog in local directory"; \
	        $(MAKE) -f $(GNUMAKEFILE) $$prog ) \
	   done \
	) fi

########################################################################

#----------------------
# User Defined Targets
#----------------------

vpath % $(INC_LIST)

USER_FFLAGS = $(DGEOS5)
USER_FMODS  = $(foreach dir,$(INC_LIST),$(M)$(dir)) 
USER_FINCS  = $(foreach dir,$(INC_LIST),$(I)$(dir)) 
FREAL = $(FREAL4)

$(LIB) lib : $(DEPS) $(OBJS)
	@if [ "$(OBJS)" != "" ]; then ( \
	   $(RM) $(LIB); \
	   $(AR) $(AR_FLAGS) $(LIB) $(OBJS) \
	   $(RANLIB) $(RANLIB_FLAGS) $(LIB) \
	) fi

% : %.o getids.H
	$(FC) $(LDFLAGS) -o $@ $*.o $(LIB) $(LIB_LIST) $(LIB_SYS) $(FOPT)

#----------------------------------------------------
# Hack to prevent remaking dep files during cleaning
#----------------------------------------------------
ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
ifneq ($(findstring clean,$(SUBTARGET)),clean)
   -include $(ALLDEPS)
endif
endif

  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

#.
