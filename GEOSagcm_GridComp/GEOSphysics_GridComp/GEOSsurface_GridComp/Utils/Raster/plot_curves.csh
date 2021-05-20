#!/usr/bin/csh -x

/bin/cp ../src/clsm_plots.pro .

module load tool/idl-8.1

idl  <<EOB

.compile clsm_plots

check_satparam

exit
EOB

xview img.0000001.jpg &

