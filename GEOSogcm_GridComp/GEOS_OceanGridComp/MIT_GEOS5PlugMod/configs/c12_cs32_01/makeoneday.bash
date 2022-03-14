#!/bin/bash

# -----
# Usage
# -----

usage ()
{
   echo "Usage: $0 [SIX] [TEN] [1MO] [RRMTG] [REPLAY] [BENCH] [BRO]  [HISTORY] [G40] [TINY] [RRTMG_SW] [RRTMG_LW] [POLICE] [LINK]  [NC4] [BIN] [WW3] [MIC] [HYPERQ] [DAS]"
   echo ""
   echo " Common Well-Tested Options "
   echo " ========================== "
   echo "        SIX: Make a six-hour experiment"
   echo "        TEN: Make a ten-day experiment"
   echo "        1MO: Use single-moment moist"
   echo "        MEM: Add memory stats"
   echo "      RRTMG: Enable both shortwave and longwave RRTMG code"
   echo "     REPLAY: Turn on regular replay"
   echo "      BENCH: Use benchmark qos"
   echo "        BRO: Turn on Broadwell at NAS"
   echo ""
   echo " Less Common Options "
   echo " =================== "
   echo "    HISTORY: Enables smaller HISTORY.rc"
   echo "        G40: Use Ganymed-4_0 directories"
   echo "       TINY: Use minimal BCs"
   echo "   RRTMG_SW: Enable shortwave RRTMG code"
   echo "   RRTMG_LW: Enable longwave RRTMG code"
   echo "     POLICE: Enables PoliceMe functionality"
   echo "       LINK: Link restarts in scratch"
   echo ""
   echo " Rare or Obsolete Options (Use at your own risk) "
   echo " =============================================== "
   echo "        NC4: Convert to use NC4 restarts"
   echo "        BIN: Convert to use binary restarts"
   echo "        WW3: Turn on WAVEWATCH III capabilities"
   echo "        MIC: Enable MIC code"
   echo "     HYPERQ: Enable Hyper-Q"
   echo "        DAS: Add DAS SBATCH pragmas"
   echo ""
   echo " Note: Lowercase also allowed"
}

# ------------------------------
# Process command line arguments
# ------------------------------

# Set defaults
# ------------

USEG40=FALSE
SIX=FALSE
TEN=FALSE
ONEMO=FALSE
NC4=FALSE
BIN=FALSE
WW3=FALSE
MEM=FALSE
TINY=FALSE
MIC=FALSE
RRTMG=FALSE
RRTMG_SW=FALSE
RRTMG_LW=FALSE
HYPERQ=FALSE
HISTORY=FALSE
POLICE=FALSE
LINK=FALSE
DAS=FALSE
BENCH=FALSE
REPLAY=FALSE
BRO=FALSE

while [ "${1+defined}" ]
do
   case "$1" in
      "G40" | "g40")
         USEG40=TRUE
         shift
         ;;
      "SIX" | "six")
         SIX=TRUE
         shift
         ;;
      "TEN" | "ten")
         TEN=TRUE
         shift
         ;;
      "1MO" | "1mo")
         ONEMO=TRUE
         shift
         ;;
      "NC4" | "nc4")
         NC4=TRUE
         shift
         ;;
      "BIN" | "bin")
         BIN=TRUE
         shift
         ;;
      "MEM" | "mem")
         MEM=TRUE
         shift
         ;;
      "TINY" | "tiny" | "MIN" | "min")
         TINY=TRUE
         shift
         ;;
      "WW3" | "ww3")
         WW3=TRUE
         shift
         ;;
      "MIC" | "mic")
         MIC=TRUE
         shift
         ;;
      "RRTMG" | "rrtmg")
         RRTMG=TRUE
         shift
         ;;
      "RRTMG_SW" | "rrtmg_sw")
         RRTMG_SW=TRUE
         shift
         ;;
      "RRTMG_LW" | "rrtmg_lw")
         RRTMG_LW=TRUE
         shift
         ;;
      "HYPERQ" | "hyperq" | "hq")
         HYPERQ=TRUE
         shift
         ;;
      "HIS" | "his" | "HISTORY" | "history" )
         HISTORY=TRUE
         shift
         ;;
      "POL" | "pol" | "POLICE" | "police" )
         POLICE=TRUE
         shift
         ;;
      "LINK" | "link" )
         LINK=TRUE
         shift
         ;;
      "DAS" | "das" )
         DAS=TRUE
         shift
         ;;
      "BENCH" | "bench" )
         BENCH=TRUE
         shift
         ;;
      "REPLAY" | "replay" )
         REPLAY=TRUE
         shift
         ;;
      "BRO" | "bro" )
         BRO=TRUE
         shift
         ;;
      -h | --help)
         usage
         exit 0
         ;;
      *)
         echo "Unknown option: $1"
         echo ""
         usage
         exit 1
         ;;
   esac
done

if [[ $USEG40 == TRUE ]] 
then
   echo "Using Ganymed-4_0 directories"
   BCDIRNAME="G40"
else
   echo "Using Heracles-5_0 directories"
   BCDIRNAME="H50"
fi

if [[ $SIX == TRUE ]] 
then
   echo "Making six-hour experiment with "
elif [[ $TEN == TRUE ]]
then
   echo "Making ten-day experiment with "
else
   echo "Making one-day experiment with "
fi

if [[ $NC4 == TRUE && $BIN == TRUE ]]
then
   echo "You can't have both NC4 and BIN set to true"
   exit 9
fi

echo "     USEG40: $USEG40"
echo "        SIX: $SIX"
echo "        TEN: $TEN"
echo "        1MO: $ONEMO"
echo "        NC4: $NC4"
echo "        BIN: $BIN"
echo "        MEM: $MEM"
echo "       TINY: $TINY"
echo "        WW3: $WW3"
echo "        MIC: $MIC"
echo "      RRTMG: $RRTMG"
echo "   RRTMG_SW: $RRTMG_SW"
echo "   RRTMG_LW: $RRTMG_LW"
echo "     HYPERQ: $HYPERQ"
echo "    HISTORY: $HISTORY"
echo "     POLICE: $POLICE"
echo "       LINK: $LINK"
echo "        DAS: $DAS"
echo "      BENCH: $BENCH"
echo "     REPLAY: $REPLAY"
echo "        BRO: $BRO"
echo ""

# -------------------
# Locate where we are
# -------------------

NODENAME=$(uname -n)

if [[ $NODENAME == discover* || $NODENAME == dali* || $NODENAME == warp* || $NODENAME == borg* ]]
then
   SITE=NCCS

   COLORDIFF=/home/mathomp4/bin/colordiff
   UTIL_DIR=$WorkingDir/TinyFromMatt
   PBZIP2=/home/mathomp4/bin/pbzip2

elif [[ $NODENAME == pfe* || $NODENAME == r[0-9]*i[0-9]*n[0-9]* || $NODENAME == bridge* || $NODENAME == maia* ]]
then
   SITE=NAS

   COLORDIFF=/nobackup/gmao_SIteam/Utilities/bin/colordiff
   UTIL_DIR=$WorkingDir/TinyFromMatt
   PBZIP2=/usr/bin/pbzip2

elif [[ $NODENAME == jibb* || $NODENAME == jcc* ]]
then
   SITE=JIBB

   COLORDIFF=/jibb/nobackup/gmaosi/Utilities/bin/colordiff
   UTIL_DIR=/jibb/nobackup/gmaosi/ModelData
   PBZIP2=/jibb/nobackup/gmaosi/Utilities/bin/pbzip2

else
   SITE=DESKTOP

   COLORDIFF=colordiff
   UTIL_DIR=$HOME/geos5/TinyFromMatt
   PBZIP2=/ford1/share/gmao_SIteam/Utilities/bin/pbzip2

fi

MIN_DIR=$UTIL_DIR/TinyBCs-$BCDIRNAME
MIN_BCS_DIR=$MIN_DIR/bcs
MIN_SST_DIR=$MIN_DIR/sst
MIN_CHM_DIR=$MIN_DIR/chem
MIN_AERO_DIR=$MIN_CHM_DIR/g5chem/L72/aero_clm
MIN_RESTART_DIR=$MIN_DIR/rs

RESTARTS_H10_DIR=Restarts-H10

# ----------------------------------------
# Replay only works at NCCS. Die otherwise
# ----------------------------------------

if [[ $REPLAY == TRUE && $SITE != NCCS ]]
then
   echo "Detected site: $SITE and REPLAY: $REPLAY"
   echo "REPLAY only works at NCCS"
   exit 400
fi

# ------------------------------------------
# Broadwell only works at NAS. Die otherwise
# ------------------------------------------

if [[ $BRO == TRUE && $SITE != NAS ]]
then
   echo "Detected site: $SITE and BRO: $BRO"
   echo "BRO only works at NAS"
   exit 401
fi

# ---------------
# Local Functions
# ---------------

restore_save ()
{
   if [ -e $1.save ]
   then
      echo "Restoring $1.save to $1..."
      mv $1.save $1
   fi
}

copy_save ()
{
   if [ ! -e $1.save ]
   then
      echo "Copying $1 to $1.save..."
      cp $1 $1.save
   fi
}

print_changes ()
{
   DIFF=$(diff "$1.save" "$1")
   if [ $? -ne 0 ]
   then
      echo "Changes made to $1:"
      $COLORDIFF $1.save $1
   fi
   echo
}

convert_rrtmg_lw ()
{
   AERODIR=$UTIL_DIR/AerosolTables-$BCDIRNAME/ChouS-RRTMGI
   sed -i -r -e "/^DU_OPTICS:/ s#ExtData/.*/x/opticsBands_DU.v14_2.nc#$AERODIR/opticsBands_DU.ChouS-RRTMGI.v14_2.nc#" \
             -e "/^DU_OPTICS:/ s#ExtData/.*/x/opticsBands_DU.v15_3.nc#$AERODIR/opticsBands_DU.ChouS-RRTMGI.v15_3.nc#" \
             -e "/^SS_OPTICS:/ s#ExtData/.*/x/opticsBands_SS.v3_3.nc#$AERODIR/opticsBands_SS.ChouS-RRTMGI.v3_3.nc#" \
             -e "/^SU_OPTICS:/ s#ExtData/.*/x/opticsBands_SU.v1_3.nc#$AERODIR/opticsBands_SU.ChouS-RRTMGI.v1_3.nc#" \
             -e "/^OC_OPTICS:/ s#ExtData/.*/x/opticsBands_OC.v1_3.nc#$AERODIR/opticsBands_OC.ChouS-RRTMGI.v1_3.nc#" \
             -e "/^BC_OPTICS:/ s#ExtData/.*/x/opticsBands_BC.v1_3.nc#$AERODIR/opticsBands_BC.ChouS-RRTMGI.v1_3.nc#" \
             -e "/^NI_OPTICS:/ s#ExtData/.*/x/opticsBands_NI.v1_5.nc#$AERODIR/opticsBands_NI.ChouS-RRTMGI.v1_5.nc#" \
             -e "/^NI_OPTICS:/ s#ExtData/.*/x/opticsBands_NI.v2_5.nc#$AERODIR/opticsBands_NI.ChouS-RRTMGI.v2_5.nc#" \
             -e "/^NUM_BANDS:/ s/18/24/" AGCM.rc

   echo "USE_RRTMG_IRRAD: 1.0" >> AGCM.rc
}

convert_rrtmg_sw ()
{
   AERODIR=$UTIL_DIR/AerosolTables-$BCDIRNAME/RRTMGS-ChouI
   sed -i -r -e "/^DU_OPTICS:/ s#ExtData/.*/x/opticsBands_DU.v14_2.nc#$AERODIR/opticsBands_DU.RRTMGS-ChouI.v14_2.nc#" \
             -e "/^DU_OPTICS:/ s#ExtData/.*/x/opticsBands_DU.v15_3.nc#$AERODIR/opticsBands_DU.RRTMGS-ChouI.v15_3.nc#" \
             -e "/^SS_OPTICS:/ s#ExtData/.*/x/opticsBands_SS.v3_3.nc#$AERODIR/opticsBands_SS.RRTMGS-ChouI.v3_3.nc#" \
             -e "/^SU_OPTICS:/ s#ExtData/.*/x/opticsBands_SU.v1_3.nc#$AERODIR/opticsBands_SU.RRTMGS-ChouI.v1_3.nc#" \
             -e "/^OC_OPTICS:/ s#ExtData/.*/x/opticsBands_OC.v1_3.nc#$AERODIR/opticsBands_OC.RRTMGS-ChouI.v1_3.nc#" \
             -e "/^BC_OPTICS:/ s#ExtData/.*/x/opticsBands_BC.v1_3.nc#$AERODIR/opticsBands_BC.RRTMGS-ChouI.v1_3.nc#" \
             -e "/^NI_OPTICS:/ s#ExtData/.*/x/opticsBands_NI.v1_5.nc#$AERODIR/opticsBands_NI.RRTMGS-ChouI.v1_5.nc#" \
             -e "/^NI_OPTICS:/ s#ExtData/.*/x/opticsBands_NI.v2_5.nc#$AERODIR/opticsBands_NI.RRTMGS-ChouI.v2_5.nc#" \
             -e "/^NUM_BANDS:/ s/18/24/" AGCM.rc

   echo "USE_RRTMG_SORAD: 1.0" >> AGCM.rc
}

convert_rrtmg_swlw ()
{
   AERODIR=$UTIL_DIR/AerosolTables-$BCDIRNAME/RRTMGS-RRTMGI
   sed -i -r -e "/^DU_OPTICS:/ s#ExtData/.*/x/opticsBands_DU.v14_2.nc#$AERODIR/opticsBands_DU.RRTMGS-RRTMGI.v14_2.nc#" \
             -e "/^DU_OPTICS:/ s#ExtData/.*/x/opticsBands_DU.v15_3.nc#$AERODIR/opticsBands_DU.RRTMGS-RRTMGI.v15_3.nc#" \
             -e "/^SS_OPTICS:/ s#ExtData/.*/x/opticsBands_SS.v3_3.nc#$AERODIR/opticsBands_SS.RRTMGS-RRTMGI.v3_3.nc#" \
             -e "/^SU_OPTICS:/ s#ExtData/.*/x/opticsBands_SU.v1_3.nc#$AERODIR/opticsBands_SU.RRTMGS-RRTMGI.v1_3.nc#" \
             -e "/^OC_OPTICS:/ s#ExtData/.*/x/opticsBands_OC.v1_3.nc#$AERODIR/opticsBands_OC.RRTMGS-RRTMGI.v1_3.nc#" \
             -e "/^BC_OPTICS:/ s#ExtData/.*/x/opticsBands_BC.v1_3.nc#$AERODIR/opticsBands_BC.RRTMGS-RRTMGI.v1_3.nc#" \
             -e "/^NI_OPTICS:/ s#ExtData/.*/x/opticsBands_NI.v1_5.nc#$AERODIR/opticsBands_NI.RRTMGS-RRTMGI.v1_5.nc#" \
             -e "/^NI_OPTICS:/ s#ExtData/.*/x/opticsBands_NI.v2_5.nc#$AERODIR/opticsBands_NI.RRTMGS-RRTMGI.v2_5.nc#" \
             -e "/^NUM_BANDS:/ s/18/30/" AGCM.rc

   echo "USE_RRTMG_IRRAD: 1.0" >> AGCM.rc
   echo "USE_RRTMG_SORAD: 1.0" >> AGCM.rc
}

# ------------------------------------
# Restore to original all edited files
# ------------------------------------

restore_save "AGCM.rc"
copy_save "AGCM.rc"

restore_save "CAP.rc"
copy_save "CAP.rc"

restore_save "HISTORY.rc"
copy_save "HISTORY.rc"

restore_save "gcm_run.j"
copy_save "gcm_run.j"

restore_save "regress/gcm_regress.j"
copy_save "regress/gcm_regress.j"

#restore_save "RC/GOCARTdata_ExtData.rc"
#copy_save "RC/GOCARTdata_ExtData.rc"

restore_save "RC/GEOS_ChemGridComp.rc"
copy_save "RC/GEOS_ChemGridComp.rc"

# -----------------------------
# Detect Atmospheric Resolution
# -----------------------------

AGCM_IM=$(grep "^ \+AGCM_IM:" AGCM.rc | awk '{print $2}')

case $AGCM_IM in
   12)
   ATMOS_RES="c12"
   ;;
   24)
   ATMOS_RES="c24"
   ;;
   48)
   ATMOS_RES="c48"
   ;;
   90)
   ATMOS_RES="c90"
   ;;
   180)
   ATMOS_RES="c180"
   ;;
   360)
   ATMOS_RES="c360"
   ;;
   720)
   ATMOS_RES="c720"
   ;;
   1000)
   ATMOS_RES="c1000"
   ;;
   1440)
   ATMOS_RES="c1440"
   ;;
   2880)
   ATMOS_RES="c2880"
   ;;
   72)
   ATMOS_RES="72x46"
   ;;
   144)
   ATMOS_RES="144x91"
   ;;
   288)
   ATMOS_RES="288x181"
   ;;
   576)
   ATMOS_RES="576x361"
   ;;
   1152)
   ATMOS_RES="1152x721"
   ;;
   *)
   ATMOS_RES="UNKNOWN"
   echo "$ATMOS_RES atmospheric resolution found!"
   ;;
esac

# --------------------------------
# Detect Number of Pressure Levels
# --------------------------------

AGCM_LM=$(grep "^ \+AGCM_LM:" AGCM.rc | awk '{print $2}')

case $AGCM_LM in
   72)
   ATMOS_LEVS="72"
   ;;
   132)
   ATMOS_LEVS="132"
   ;;
   137)
   ATMOS_LEVS="137"
   ;;
   144)
   ATMOS_LEVS="144"
   ;;
   *)
   ATMOS_LEVS="UNKNOWN"
   echo "Unknown number of atmospheric levels detected: $AGCM_LM"
   exit 338
   ;;
esac

# For now, append L<levs) if we are running other than 72
# -------------------------------------------------------

if [[ $ATMOS_LEVS == 72 ]]
then
   ATMOS_RES_LEVS=$ATMOS_RES
else
   ATMOS_RES_LEVS=$ATMOS_RES-L$ATMOS_LEVS
fi

# -----------------------
# Detect Ocean Resolution
# -----------------------

OGCM_IM=$(grep "^ \+OGCM_IM:" AGCM.rc | awk '{print $2}')

case $OGCM_IM in
   360)
   OCEAN_RES="Reynolds"
   ;;
   1440)
   OCEAN_RES="MERRA-2"
   ;;
   2880)
   OCEAN_RES="Ostia"
   ;;
   192)
   OCEAN_RES="Reynolds"
   ;;
   *)
   OCEAN_RES="UNKNOWN"
   echo "$OCEAN_RES ocean resolution found!"
   ;;
esac

# --------------------
# Detect Restart Types
# --------------------

DYN_RESTART_TYPE=$(grep "^DYN_INTERNAL_RESTART_TYPE:" AGCM.rc | awk '{print $2}')

case $DYN_RESTART_TYPE in
   binary)
   RESTART_TYPE="binary"
   ;;
   pbinary)
   RESTART_TYPE="binary"
   ;;
   pnc4)
   RESTART_TYPE="nc4"
   ;;
   *)
   echo "Unknown DYN_INTERNAL_RESTART_TYPE detected: $DYN_RESTART_TYPE"
   exit 1
   ;;
esac

# -----------------------------------------
# Do we want to convert to NC4 or to binary
# -----------------------------------------

if [[ $NC4 == TRUE && $RESTART_TYPE == "binary" ]]
then
   echo "NC4 requested and binary restart found. Converting to NC4"

   #sed -i -e "/_rst$/ s/_rst/_rst.nc4/" \
          #-e "/_checkpoint$/ s/_checkpoint/_checkpoint.nc4/" AGCM.rc

   sed -i -e "/.*_TYPE: \+pbinary$/ s/pbinary/pnc4/" \
          -e "/.*_TYPE: \+binary$/ s/binary/pnc4/" AGCM.rc

   # The above clobbers VEGDYN. Undo.
   sed -i -e "/VEGDYN/ s/pnc4/binary/" AGCM.rc

   RESTART_TYPE="nc4"

elif [[ $BIN == TRUE && $RESTART_TYPE == "nc4" ]]
then
   echo "Binary requested and NC4 restart found. Converting to binary"

   #sed -i -e "/_rst$/ s/_rst/_rst.nc4/" \
          #-e "/_checkpoint$/ s/_checkpoint/_checkpoint.nc4/" AGCM.rc

   sed -i -e "/.*_TYPE: \+pnc4$/ s/pnc4/binary/" AGCM.rc

   # The above clobbers VEGDYN. Undo.
   #sed -i -e "/VEGDYN/ s/pnc4/binary/" AGCM.rc

   RESTART_TYPE="binary"
fi


# ----------------------------------
# Do our restart files have suffixes
# ----------------------------------

FVNAME=$(grep "^DYN_INTERNAL_RESTART_FILE:" AGCM.rc | awk '{print $2}')
EXTENSION="${FVNAME##*.}"

# --------------------------
# Detect Boundary Conditions
# --------------------------

BCSDIR=$(grep "^setenv BCSDIR" gcm_run.j | awk '{print $3}')

#echo $BCSDIR

USING_4_0_BCS=$(echo $BCSDIR | grep -ow "Ganymed-4_0" )
USING_HNL_BCS=$(echo $BCSDIR | grep -ow "Heracles-NL" )

#echo $USING_4_0_BCS

if [[ $USING_4_0_BCS == "Ganymed-4_0" ]]
then
   if [[ $TINY == TRUE ]]
   then
      RESTARTS_DIR=$MIN_RESTART_DIR
   else
      RESTARTS_DIR=$UTIL_DIR/$RESTARTS_H10_DIR
   fi
elif [[ $USING_HNL_BCS == "Heracles-NL" ]]
then
   RESTARTS_DIR=$UTIL_DIR/$RESTARTS_H10_DIR
else
   echo "You seem to be using an unknown BCSDIR:"
   echo "   $BCSDIR"
   echo
   echo "This script can handle Ganymed-4_0."
   exit
fi

# -------------
# Link Restarts
# -------------

if [[ ! -e fvcore_internal_rst ]]
then
   if [[ ! -d $RESTARTS_DIR/$RESTART_TYPE/$OCEAN_RES/$ATMOS_RES_LEVS ]]
   then
      echo "FAILURE!"
      echo "Restarts of type $RESTART_TYPE for resolution $ATMOS_RES using $ATMOS_LEVS levels on ocean $OCEAN_RES do not exist in $RESTARTS_DIR"
      exit 2
   else
      echo "Linking $RESTART_TYPE restarts for resolution $ATMOS_RES using $ATMOS_LEVS levels on ocean $OCEAN_RES from $RESTARTS_DIR..."
      if [[ ${EXTENSION} == ${FVNAME} ]]
      then
         echo "Restarts have no EXTENSION..."
         ln -sv $RESTARTS_DIR/$RESTART_TYPE/$OCEAN_RES/$ATMOS_RES_LEVS/*_rst .
      else
         echo "Restarts have EXTENSION: $EXTENSION"
         if [[ -e $RESTARTS_DIR/$RESTART_TYPE/$OCEAN_RES/$ATMOS_RES_LEVS/fvcore_internal_rst.$EXTENSION ]]
         then
            ln -sv $RESTARTS_DIR/$RESTART_TYPE/$OCEAN_RES/$ATMOS_RES_LEVS/*_rst.$EXTENSION .
         else
            for file in $RESTARTS_DIR/$RESTART_TYPE/$OCEAN_RES/$ATMOS_RES_LEVS/*_rst
            do
               filename=$(basename $file)
               ln -sv $file ./$filename.$EXTENSION
            done
         fi
      fi
      if [[ ! -e cap_restart ]]
      then
         cp $RESTARTS_DIR/$RESTART_TYPE/$OCEAN_RES/$ATMOS_RES_LEVS/cap_restart .
         if [[ ! -w cap_restart ]]
         then
            echo "cap_restart seems to be read only. For safety's sake, we make it writable"
            chmod -v u+w cap_restart
         fi
      else
         echo "cap_restart already exists. Not copying."
      fi
   fi
else
   echo "Found fvcore_internal_rst. Assuming you have needed restarts!"
fi

if [[ $WW3 == TRUE ]]
then
   echo -n "Linking WW3 mod_def.ww3 file for $OCEAN_RES"
   if [[ $OCEAN_RES == "Reynolds" ]] 
   then
      WW3_DIR=$UTIL_DIR/GridGen-v2.1/geos_1deg_latlon_grid_dateline_edge_poleline_edge
      echo " from $WW3_DIR"
      MOD_FILE=mod_def.ww3.CFL_for_30m
      ln -s $WW3_DIR/$MOD_FILE .
      ln -s $MOD_FILE mod_def.ww3
   elif [[ $OCEAN_RES == "Ostia" ]]
   then
      WW3_DIR=$UTIL_DIR/GridGen-v2.1/ostia_eighth_latlon_grid_dateline_edge_poleline_edge
      echo " from $WW3_DIR"
      MOD_FILE=mod_def.ww3.150_50_75_15
      ln -s $WW3_DIR/$MOD_FILE .
      ln -s $MOD_FILE mod_def.ww3
   elif [[ $OCEAN_RES == "MERRA-2" ]]
   then
      WW3_DIR=$UTIL_DIR/GridGen-v2.1/merra2_quart_deg_latlon_grid_dateline_edge_poleline_edge
      echo " from $WW3_DIR"
      MOD_FILE=mod_def.ww3.300_100_150_30
      ln -s $WW3_DIR/$MOD_FILE .
      ln -s $MOD_FILE mod_def.ww3
   else
      echo ""
      echo "ERROR: No WW3 mod_def.ww3 available for $OCEAN_RES"
      exit 2
   fi
fi


# ------
# CAP.rc
# ------

if [[ $SIX == TRUE ]]
then
   sed -i -e '/^JOB_SGMT:/ s/000000[0-9][0-9] 000000/00000000 060000/' CAP.rc
elif [[ $TEN == TRUE ]]
then
   sed -i -e '/^JOB_SGMT:/ s/000000[0-9][0-9]/00000010/' CAP.rc
else
   sed -i -e '/^JOB_SGMT:/ s/000000[0-9][0-9]/00000001/' CAP.rc
fi

sed -i -e '/^NUM_SGMT:/ s/[0-9][0-9]*/1/' \
       -e '/^MAPL_ENABLE_MEMUTILS:/ s/NO/NO/' \
       -e '/^MAPL_ENABLE_TIMERS:/ s/NO/YES/' CAP.rc

if [[ $MEM == TRUE ]]
then
   sed -i -e '/^MAPL_ENABLE_MEMUTILS:/ s/NO/YES/' CAP.rc
fi

sed -i -e '/^CoresPerNode:/ d' CAP.rc

print_changes "CAP.rc"

# ----------
# HISTORY.rc
# ----------

if [[ $HISTORY == TRUE ]]
then

   # This command should add a # to all lines that have geosgcm
   # or tavg with spaces at the beginning between COLLECTIONS and ::

   sed -i -e "/^COLLECTIONS:/,/.*  ::/ {
                                       /^ .*geosgcm/ s/ /#/
                                       /^ .*tavg/    s/ /#/}" HISTORY.rc

fi

if [[ $TINY == TRUE ]]
then

   echo "Since we turn off TR, tracer collection cannot run"

   sed -i -e "/ *'geosgcm_tracer'/ s/ /#/" HISTORY.rc
fi

print_changes "HISTORY.rc"

# ---------
# gcm_run.j
# ---------

sed -i -e '/^echo GEOSgcm/ a exit' \
       -e '/^#SBATCH --account/ a #SBATCH --mail-type=ALL' \
       -e '/^#SBATCH -A / a #SBATCH --mail-type=ALL' \
       -e '/^#PBS -W group_list/ a #SBATCH --mail-type=ALL' \
       -e '/numrs == 0/ s/== 0/== 1/ ' gcm_run.j

if [[ $TEN != TRUE ]]
then

   if (( $AGCM_IM < 180 ))
   then
      sed -i -e '/^#PBS -l walltime/ s/=12:00:/=0:15:/ ' \
            -e '/^#PBS -l walltime/ s/=8:00:/=0:15:/ ' \
            -e '/^#SBATCH --time/ s/=12:00:/=0:15:/ ' gcm_run.j
   elif (( $AGCM_IM < 720 ))
   then
      sed -i -e '/^#PBS -l walltime/ s/=12:/=1:/ ' \
            -e '/^#PBS -l walltime/ s/=8:/=1:/ ' \
            -e '/^#SBATCH --time/ s/=12:/=1:/ ' gcm_run.j
   else
      sed -i -e '/^#PBS -l walltime/ s/=12:/=2:/ ' \
            -e '/^#PBS -l walltime/ s/=8:/=2:/ ' \
            -e '/^#SBATCH --time/ s/=12:/=2:/ ' gcm_run.j
   fi
fi

if [[ $LINK == TRUE ]]
then
   sed -i -e '/^  if(-e $EXPDIR\/$rst ) \/bin/ s/cp/ln -s/' gcm_run.j
fi

if [[ $SITE == NAS ]] 
then
   sed -i -e '/^#PBS -l walltime/ a #PBS -W umask=0022' \
          -e '/^#PBS -l walltime/ a #PBS -m abe' gcm_run.j
fi

if [[ $DAS == TRUE ]]
then
   sed -i -e '/^#SBATCH --account/ a #SBATCH --reservation=das' \
          -e '/^#SBATCH --account/ a #SBATCH --qos=das2015' \
          -e '/^#SBATCH -A / a #SBATCH --reservation=das' \
          -e '/^#SBATCH -A/ a #SBATCH --qos=das2015' gcm_run.j
fi

if [[ $BENCH == TRUE ]]
then
   sed -i -e '/^#SBATCH --account/ a #SBATCH --partition=preops' \
          -e '/^#SBATCH --account/ a #SBATCH --qos=benchmark' \
          -e '/^#SBATCH -A / a #SBATCH --partition=preops' \
          -e '/^#SBATCH -A/ a #SBATCH --qos=benchmark' gcm_run.j
fi

if [[ $BRO == TRUE ]]
then
   sed -i -e '/^#PBS -l select/ s/ncpus=24:mpiprocs=24:model=has/ncpus=28:mpiprocs=28:model=bro/' gcm_run.j
fi

if [[ $WW3 == TRUE ]]
then
   echo "Enabling WW3 in gcm_run.j"
   sed -i -e '/                             \/bin\/cp \-f/ a\
                             /bin/cp -f $EXPDIR/*.ww3 .' gcm_run.j
fi 

if [[ $POLICE == TRUE ]]
then
   echo "Enabling PoliceMe in gcm_run.j"
   sed -i -e '/limit stacksize/ a\
\
if { /home/mathomp4/bin/amibatch } then \
   echo "I am running in batch mode. Executing PoliceMe" \
   mkdir -p policeme \
   /usr/local/other/policeme/policeme.exe -d policeme \
endif ' gcm_run.j
fi

if [[ $MIC == TRUE ]]
then
   echo "Enabling for MIC"
   sed -i -e '/limit stacksize/ a \
\
   setenv PATH /opt/intel/mic/bin:${PATH} \
   setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/opt/intel/mic/coi/host-linux-release/lib \
 \
   unsetenv MIC_ENV_PREFIX \
   unsetenv MIC_OMP_NUM_THREADS \
   unsetenv OMP_NUM_THREADS \
   unsetenv OFFLOAD_INIT \
   unsetenv MIC_USE_2MB_BUFFERS \
   unsetenv MIC_KMP_AFFINITY \
   unsetenv MKL_MIC_ENABLE \
   unsetenv OFFLOAD_REPORT \
   unsetenv I_MPI_MIC  \
   unsetenv I_MPI_FABRIC \
   unsetenv I_MPI_DEBUG \
 \
   setenv MIC_ENV_PREFIX MIC \
   setenv MIC_OMP_NUM_THREADS 236 #multiple of 59 4 threads per core \
   setenv OMP_NUM_THREADS 16 \
   setenv OFFLOAD_INIT on_start #good for timing prepares MiC \
   setenv MIC_USE_2MB_BUFFERS 2K \
   setenv MIC_KMP_AFFINITY balanced,granularity=fine \
   setenv MKL_MIC_ENABLE 1 \
   setenv OFFLOAD_REPORT 1 \
   setenv I_MPI_MIC 1  \
   setenv I_MPI_FABRIC shm:ofa \
   setenv I_MPI_DEBUG 5 ' \
   -e '/mpirun_rsh/ s/mpirun_rsh/mpirun_rsh -export/ ' gcm_run.j
fi

if [[ $HYPERQ == TRUE ]]
then
   echo "Enabling for Hyper-Q"
   if [[ $SITE == NCCS ]]
   then

   cat > HYPERQSEDFILE1 << _EOF_

####################################################################### 
#               Set up the MPS Server on Each GPU Node          
####################################################################### 
 
setenv NODELIST \`scontrol show hostnames\` 
 
setenv CUDA_TMPDIR /tmp/nvidia-hyperq-\$SLURM_JOBID 
_EOF_
   elif [[ $SITE == NAS ]]
   then

   cat > HYPERQSEDFILE1 << _EOF_

####################################################################### 
#               Set up the MPS Server on Each GPU Node          
####################################################################### 
 
setenv NODELIST \`cat \$PBS_NODEFILE | uniq\` 
 
setenv CUDA_TMPDIR /tmp/nvidia-hyperq-\$PBS_JOBID 
_EOF_
   fi

   cat >> HYPERQSEDFILE1 << _EOF_
setenv CUDA_VISIBLE_DEVICES 0 
setenv CUDA_MPS_CLIENT 1 
setenv CUDA_MPS_PIPE_DIRECTORY \$CUDA_TMPDIR/mps_0 
setenv CUDA_MPS_LOG_DIRECTORY  \$CUDA_TMPDIR/mps_log_0 
 
env | grep CUDA 
 
foreach node (\$NODELIST) 
   ssh -f \$node "setenv CUDA_TMPDIR \$CUDA_TMPDIR; $HOME/bin/kill_mps_server.bash" 
end 
 
sleep 3 
 
foreach node (\$NODELIST) 
   ssh -f \$node "setenv CUDA_TMPDIR \$CUDA_TMPDIR; $HOME/bin/kill_mps_server.bash" 
end 
 
sleep 3 
 
foreach node (\$NODELIST) 
   echo "Running MPS Server on \$node" 
   ssh -f \$node "setenv CUDA_TMPDIR \$CUDA_TMPDIR; $HOME/bin/run_mps_server.bash" 
end 
 
sleep 3 
_EOF_

   sed -i -e '/limit stacksize/r HYPERQSEDFILE1' gcm_run.j

   cat > HYPERQSEDFILE2 << _EOF2_

foreach node (\$NODELIST) 
   ssh -f \$node "setenv CUDA_TMPDIR \$CUDA_TMPDIR; $HOME/bin/kill_mps_server.bash" 
end

_EOF2_

   sed -i -e '/echo GEOSgcm Run/r HYPERQSEDFILE2' gcm_run.j

   /bin/rm HYPERQSEDFILE1
   /bin/rm HYPERQSEDFILE2

fi

if [[ $TINY == TRUE ]]
then
   echo "Using minimal boundary datasets"
   sed -i -e "/^setenv BCSDIR/ s#\(setenv BCSDIR\)\(.*\)#\1   ${MIN_BCS_DIR}#" \
          -e "/^setenv SSTDIR/ s#\(setenv SSTDIR\)\(.*\)#\1   ${MIN_SST_DIR}#" \
          -e "/^setenv CHMDIR/ s#\(setenv CHMDIR\)\(.*\)#\1   ${MIN_CHM_DIR}#" \
          -e "s/^#\!\/bin\/ksh$/#\!\/bin\/bash/"                               \
          -e "/pchem.species/  s#1870-2097#1999-2000#"                         \
          -e "/sst_1971/       s#1971-current#2000#"                           \
          -e "/fraci_1971/     s#1971-current#2000#"                               gcm_run.j
fi

print_changes "gcm_run.j"

# -------
# AGCM.rc
# -------

DOING_GOCART=$(grep AERO_PROVIDER AGCM.rc | awk '{print $2}')

FOUND_BOOTSTRAP=$(grep MAPL_ENABLE_BOOTSTRAP AGCM.rc | awk '{print $1}')

sed -i -e "/^NUM_WRITERS:/ s/4/1/" AGCM.rc

if [[ "$FOUND_BOOTSTRAP" == "MAPL_ENABLE_BOOTSTRAP:" ]]
then
   sed -i -e "/^MAPL_ENABLE_BOOTSTRAP:/ s/NO/YES/" AGCM.rc
else
   if [[ "$DOING_GOCART" == "GOCART" && ! -e 'gocart_internal_rst' ]]
   then
      echo "Didn't see MAPL_ENABLE_BOOTSTRAP"
      echo "For safety's sake, we bootstrap gocart"
      sed -i -e '/GOCART_INTERNAL_RESTART_FILE:/ s/ gocart_internal_rst/-gocart_internal_rst/' AGCM.rc
   fi
fi


if [[ $RRTMG_LW == TRUE ]]
then
   convert_rrtmg_lw 
fi

if [[ $RRTMG_SW == TRUE ]]
then
   convert_rrtmg_sw 
fi

if [[ $RRTMG == TRUE ]]
then
   convert_rrtmg_swlw 
fi

if [[ $ONEMO == TRUE ]]
then
   echo "Enabling single-moment in AGCM.rc"
   echo "CLDMICRO: 1MOMENT" >> AGCM.rc
fi

if [[ $WW3 == TRUE ]]
then
   echo "Enabling WW3 in AGCM.rc"
   echo "USE_WW3: 1" >> AGCM.rc
   echo "WRITE_WW3_RESTART: 0" >> AGCM.rc
fi


if [[ $MIC == TRUE ]]
then
   echo "Enabling for MIC in AGCM.rc"
   echo "SOLAR_LOAD_BALANCE: 0" >> AGCM.rc
fi

if [[ $REPLAY == TRUE ]]
then
   echo "Turning on regular replay"
   sed -i -e "/^#   REPLAY_MODE: Regular/ s/#/ /" \
          -e "/REPLAY_MODE: Regular/{n;s/#/ /}" AGCM.rc
fi

if [[ $TINY == TRUE ]]
then

   # Before we turned of GOCART.data. We now provide a
   # tiny version of it in Tiny-BCs

   #sed -i -e "/AERO_PROVIDER:/ s/GOCART.data  /None  /" \
          #-e "/pchem_clim_years:/ s/228/2/" AGCM.rc

   sed -i -e "/pchem_clim_years:/ s/228/2/" AGCM.rc

fi

print_changes "AGCM.rc"

# ---------------------
# regress/gcm_regress.j
# ---------------------

if [[ $REPLAY == TRUE ]]
then
   echo "Altering regression test for replay"
   sed -i -e "/test_duration = 21/ s/21/18/" \
          -e "/test_duration = 03/ s/03/06/" regress/gcm_regress.j
fi

print_changes "regress/gcm_regress.j"

# -----------
# cap_restart
# -----------

if [[ "$DOING_GOCART" == "GOCART" ]]
then
   USING_NR=$(grep DU_OPTICS AGCM.rc | grep NR)
   if [[ $? -eq 0 ]]
   then
      echo "     You seem to be using the Nature Run GOCART."
      echo "     Setting cap_restart to be in 2005"
      echo ""

      restore_save "cap_restart"
      copy_save "cap_restart"

      sed -i -e "/^2000/ s/2000/2005/" cap_restart

      print_changes "cap_restart"
   fi

   USING_OPS=$(grep DU_OPTICS AGCM.rc | grep g5chem)
   if [[ $? -eq 0 ]]
   then
      echo "     You seem to be using the OPS GOCART."
      echo "     Setting cap_restart to be in 2015"
      echo ""

      restore_save "cap_restart"
      copy_save "cap_restart"

      sed -i -e "/^2000/ s/2000/2015/" cap_restart

      print_changes "cap_restart"
   fi
fi

# ------------------------
# RC/GEOS_ChemGridComp.rc
# ------------------------

if [[ $TINY == TRUE ]]
then

   # Before we turned of GOCART.data. We now provide a
   # tiny version of it in Tiny-BCs

   #echo "Turning off GOCART.data and TR in GEOS_ChemGridComp"
   #sed -i -e "/ENABLE_GOCART_DATA:/ s/TRUE/FALSE/" \
          #-e "/ENABLE_TR:/ s/TRUE/FALSE/" RC/GEOS_ChemGridComp.rc

   echo "Turning off TR in GEOS_ChemGridComp"
   sed -i -e "/ENABLE_TR:/ s/TRUE/FALSE/" RC/GEOS_ChemGridComp.rc
   
   print_changes "RC/GEOS_ChemGridComp.rc"
fi

# -------------
# src directory
# -------------

if [[ -d src ]]
then
   echo "src directory found. tarring to save inodes"
   tar cf src.tar src
   if [[ -x $PBZIP2 ]]
   then
      echo "pbzip2 found: $PBZIP2, compressing"
      $PBZIP2 -l src.tar
   fi
   echo "removing src directory"
   rm -rf src
fi
