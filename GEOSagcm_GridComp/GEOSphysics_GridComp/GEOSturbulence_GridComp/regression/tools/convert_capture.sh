#!/bin/bash
# Converts a MAPL2 capture (TURBULENCE_{import,internal,export}_{before,after}_runPhase{1,2},
# old cubed-sphere layout, no .nc suffix) to the new layout MAPL3 reads.
#
# usage: convert_capture.sh OLD_DIR NEW_DIR [convertOldCStoNewCS.py]
#   the converter defaults to $GEOS_RESTART_UTILITIES/convertOldCStoNewCS.py
#   (GEOS_Restart_Utilities); no example file is given, so the output has no lons/lats.
set -e
old=$1; new=$2
conv=${3:-${GEOS_RESTART_UTILITIES:?set GEOS_RESTART_UTILITIES or pass the converter}/convertOldCStoNewCS.py}
mkdir -p "$new"
for f in "$old"/TURBULENCE_*; do
  python3 "$conv" -i "$f" -o "$new/$(basename "$f").nc"
done
ls "$new"
