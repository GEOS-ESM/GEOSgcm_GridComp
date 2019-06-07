#
# Makefile for ESMA components.
#

ifeq ("$(USER)","suarez")
.SILENT:
endif

# Make sure ESMADIR is defined
# ----------------------------
ifndef ESMADIR
       ESMADIR := $(PWD)/..
endif


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


#                  --------------------------------
#                   Recurse Make in Sub-directories
#                  --------------------------------

ALLDIRS = GEOSdataseaice_GridComp\
          GEOSCICEDyna_GridComp\

SUBDIRS = $(wildcard $(ALLDIRS) )

TARGETS = esma_install esma_clean esma_distclean esma_doc \
          help install clean doc

FREAL = $(FREAL4)

.PHONY: install FORCE

export ESMADIR BASEDIR ARCH SITE FREAL OCEAN_FDEFS AGCM_FDEFS

$(TARGETS):
	  for d in $(SUBDIRS); \
	  do			 \
	      echo Making $@ in $$d	;\
	      $(MAKE) -C $$d -e $@	;\
	      [ $$? -ne 0 ] && exit	;\
	  done || exit 0;

#                  --------------------
#                  User Defined Targets
#                  --------------------

local_esma_install local_install: $(LIB)
	@echo No local install in here

dist:
	$(MAKE) distclean
	( cd ..  ;\
          $(TAR)  cvf evac-`date '+%d%b%Y'`.tar src ;\
          $(GZIP) evac-`date '+%d%b%Y'`.tar )

distclean:
	-$(RM) *~ *.[aox] *.mod *.x
	$(MAKE) esma_distclean

binclean:
	-$(RM) -r $(ESMADIR)/Config 
	-$(RM) -r $(ESMADIR)/$(ARCH)


realclean:
	$(MAKE) distclean
	$(MAKE) binclean


#                  --------------------
#                      Configuration
#                  --------------------

$(ESMADIR)/Config/ESMA_base.mk:
	gmake -C Config esma_install

$(ESMADIR)/Config/ESMA_arch.mk: ;

#.











  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

