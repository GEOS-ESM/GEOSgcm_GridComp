#!/bin/bash

ml SMTStack/2024.04.00

MOIST_COMP=/home/mad/work/fp/geos/src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp

# find ../../ -iname *F90.SER

PPSER_VERB="$SERIALBOX_ROOT/python/pp_ser/pp_ser.py --verbose --ignore-identical -m utils_ppser_kbuff"


FILE=$MOIST_COMP/GEOS_MoistGridComp.F90
python3 $PPSER_VERB $FILE.SER > $FILE
tail -n +2 $FILE > $FILE.swp && mv $FILE.swp $FILE

FILE=$MOIST_COMP/GEOS_GFDL_1M_InterfaceMod.F90
python3 $PPSER_VERB $FILE.SER > $FILE
tail -n +2 $FILE > $FILE.swp && mv $FILE.swp $FILE

