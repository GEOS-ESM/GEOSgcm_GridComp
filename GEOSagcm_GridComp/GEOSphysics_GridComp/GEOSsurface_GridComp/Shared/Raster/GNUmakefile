# $Id$
#!/bin/sh -x

#SILENT:

XY = -x 8640 -y 4320
 T = -t 235000000

#LAND_SURFACE_DATA = /discover/nobackup/ltakacs/bcs/LandSurfaceData_for_Catchment_Parameters
LAND_SURFACE_DATA  = /discover/nobackup/projects/gmao/ssd/land/l_data/LandBCs_files_for_mkCatchParam/V001/
GLOBAL_CATCH_DATA ?= $(LAND_SURFACE_DATA)/global.cat_id.catch.DL

OMP_NUM_THREADS ?= 1
export OMP_NUM_THREADS
export MASKFILE

ifdef NC
	ATMOSGRID =CUBE
	ATMNAME=CF$(NC)x6C
	ifeq ($(NC), 0768)
		XY:=-x 17280 -y 8640
	endif
	ifeq ($(NC), 1440)
		XY:=-x 17280 -y 8640
	endif
	ifeq ($(NC), 1536)
		XY:=-x 43200 -y 21600
	endif	
	ifeq ($(NC), 2000)
		XY:=-x 21600 -y 10800
	endif
	ifeq ($(NC), 2880)
		XY:=-x 21600 -y 10800
	endif
	ifeq ($(NC), 3072)
		XY:=-x 43200 -y 21600
	endif
	ifeq ($(NC), 5760)
		XY:=-x 43200 -y 21600
	endif
else
	ATMOSGRID =LATLON
	IM ?= 0288
	JM ?= 0181
	ATMNAME=DC$(IM)xPC$(JM)
endif

ifeq ($(MASKFILE),$(filter $(MASKFILE), GEOS5_10arcsec_mask.nc GEOS5_10arcsec_mask GEOS5_10arcsec_mask_freshwater-lakes.nc))
	XY:=-x 43200 -y 21600 
endif

ifdef IMM
	OCNNAME=TM$(IMM)xTM$(JMM)
	LENX:=$(shell expr $(IMM) : '0*')
	LENX:=$(shell expr $(LENX) + 1)
	IMO:=$(shell expr substr $(IMM) $(LENX) 4)
	LENY:=$(shell expr $(JMM) : '0*')
	LENY:=$(shell expr $(LENY) + 1)
	JMO:=$(shell expr substr $(JMM) $(LENY) 4)
	OCEANGRID =data/MOM/$(IMO)x$(JMO)/grid_spec.nc
else
ifdef IMO
	OCEANGRID =LATLON
	IMO ?= 0360
	JMO ?= 0180
	OCNNAME=DE$(IMO)xPE$(JMO)
else
	OCEANGRID =CUBEDO
	OCNNAME=$(ATMNAME)
endif
endif

LNDNAME=Pfafstetter
SFCNAME=$(OCNNAME)-$(LNDNAME)
TILNAME=$(ATMNAME)_$(OCNNAME)-$(LNDNAME)

CATCH=data/CATCH/global.cat_id.catch.DL
LAND =til/$(LNDNAME).til rst/$(LNDNAME).rst
OCEAN=til/$(OCNNAME).til rst/$(OCNNAME).rst
ATMOS=til/$(ATMNAME).til rst/$(ATMNAME).rst 
SURFC=til/$(SFCNAME).til rst/$(SFCNAME).rst 
TILEF=til/$(TILNAME).til rst/$(TILNAME).rst 

DIRS= mod bin lib til rst clsm data

##all: $(DIRS) BUILD UNTAR $(TILEF) CLSM

all: $(DIRS) BUILD $(TILEF) CLSM

$(TILEF): $(ATMOS) $(SURFC)
	bin/CombineRasters $(T) $(ATMNAME) $(OCNNAME)-$(LNDNAME)
ifeq ($(OCEANGRID),data/MOM/$(IMO)x$(JMO)/grid_spec.nc)
	bin/mk_runofftbl $(TILNAME)
endif

$(ATMOS):
ifeq ($(ATMOSGRID),LATLON)
	bin/mkLatLonRaster $(XY) -t -1  $(IM) $(JM) >/dev/null
else
	bin/mkCubeFVRaster $(XY) $(NC) >/dev/null 
endif

$(SURFC): $(OCEAN) $(LAND)
	bin/CombineRasters -f 0 $(T) $(OCNNAME) $(LNDNAME) >/dev/null

$(OCEAN):
ifeq ($(OCEANGRID),LATLON)
	bin/mkLatLonRaster $(XY) -b DE -p PE -t 0 $(IMO) $(JMO) >/dev/null
else
ifneq  ($(OCEANGRID),CUBEDO)	
	bin/mkMOMAquaRaster $(XY) $(OCEANGRID) > /dev/null
endif
endif

$(LAND): $(CATCH)
	bin/mkLandRaster $(XY) -v $(T)

# The variable:  GLOBAL_CATCH_DATA 
# allows different catchment definitions to be used.
# The DEFAULT points to those defined by the LAND group
# -----------------------------------------------------
$(CATCH):
	mkdir -p data/CATCH
	/bin/ln -s  $(LAND_SURFACE_DATA)/* data/CATCH
	/bin/rm -rf                        data/CATCH/$(MASKFILE)
	/bin/ln -s  $(GLOBAL_CATCH_DATA)   data/CATCH/$(MASKFILE)
	/bin/ls -lat data/CATCH

$(DIRS):
	mkdir -p $@

BUILD: $(DIRS)
	gmake -C src

#UNTAR: $(DIRS)
#	cd data/CATCH; tar -xvzf ../CatchData.tar.gz; cd ../../

clean:
	/bin/rm -rf $(DIRS)
#	cvs upd -d > /dev/null

CLSM: $(TILEF)
	/bin/rm -rf data/CATCH
	/bin/ln -s  $(LAND_SURFACE_DATA)/ data/CATCH
	bin/mkCatchParam $(XY) -g $(TILNAME)
	bin/mkCatchParam_openmp $(XY) -g $(TILNAME)
	/bin/rm -rf data/CATCH
