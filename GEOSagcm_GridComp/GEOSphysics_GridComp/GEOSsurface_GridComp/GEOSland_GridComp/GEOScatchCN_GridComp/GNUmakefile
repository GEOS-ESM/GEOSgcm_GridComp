#
# Makefile for ESMA components.
#
# REVISION HISTORY:
#
# 20oct2004  da Silva  First ESMA version.
# 03Mar2008  Stassi    Modified version from GMAO_Shared/MAPL_cfio
# 27Aug2008  Reichle   Renamed (mislabeled) "mpi/nompi" to "openmp/noopenmp"
#

# Make sure ESMADIR is defined
# ----------------------------
ifndef ESMADIR
       ESMADIR := $(shell pwd)/../../../../../../../..
endif

THIS := $(shell basename `pwd`)

# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies
  include $(ESMADIR)/Config/GMAO_base.mk

#                  ---------------------
#                  Standard ESMA Targets
#                  ---------------------

esma_help :
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
	@echo "         SITE = $(SITE) "
	@echo "         THIS = $(THIS) "


#                  --------------------------------
#                   Recurse Make in Sub-directories
#                  --------------------------------

ifeq ($(strip $(OMPFLAG)),)
  SUBDIRS_OPENMP = noopenmp
else
  SUBDIRS_OPENMP = openmp noopenmp
endif

#SUBDIRS = mk_restarts # not ready commented out
SUBDIRS = 

TARGETS = esma_install esma_clean esma_distclean esma_doc \
          help install doc test benchmark clean distclean

$(TARGETS): 
	@ t=$@; argv="$(SUBDIRS_OPENMP)" ;\
	for d in $$argv; do			 \
	  ( $(MKDIR) $$d            		;\
	    $(MAKE) -C $$d -f ../GNUmakefile_openmp -e $$t NTHIS=$(THIS) OPENMP=$$d ) \
	done

	@ t=$@; argv="$(SUBDIRS)" ;\
	for d in $$argv; do			 \
	  ( cd $$d            		        ;\
	    echo ""; echo Making $$t in `pwd`   ;\
	    $(MAKE) -e $$t ) \
	done

	$(MAKE) local_$@

local_esma_install local_install: $(LIB)
	@echo No local install in here

#                  --------------------
#                  User Defined Targets
#                  --------------------

openmp_install:
	$(MKDIR) openmp
	@echo;echo Making install for openmp in $(PWD)
	$(MAKE) -C openmp -f ../GNUmakefile_openmp -e install NTHIS=$(THIS) OPENMP=openmp

noopenmp_install: openmp_install   # let openmp_install finish before competing for memory
	$(MKDIR) noopenmp
	@echo;echo Making install for noopenmp in $(PWD)
	$(MAKE) -C noopenmp -f ../GNUmakefile_openmp -e install NTHIS=$(THIS) OPENMP=noopenmp

mk_restarts_install: openmp_install noopenmp_install

local_esma_distclean local_distclean:
	-$(RM) -r *~ *.[aoxd] *.[Mm][Oo][Dd] openmp noopenmp

local_esma_clean local_clean:
	-$(RM) -r *~ *.[aox] *.[Mm][Oo][Dd]

local_doc local_help local_test local_benchmark:
	@$(PROTEX) $(PROTEX_FLAGS) *GridComp*.[fF]* > $(ESMADOC)/$(THIS).tex

  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

#.
