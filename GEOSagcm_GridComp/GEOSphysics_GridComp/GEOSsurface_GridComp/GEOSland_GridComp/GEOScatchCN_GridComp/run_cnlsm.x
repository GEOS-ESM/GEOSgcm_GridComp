#!/bin/csh -f

cat << EOF > /dev/null 
WORK_PATH: /discover/nobackup/smahanam/ROUTING_MODEL/RESEARCH/FCST/may31/
# MET_PATH: /discover/nobackup/smahanam/GEOS5_s2s/data/v2.1/runoff/may31.geosgcm_vis2d.2020.nc4
MET_PATH: /discover/nobackup/projects/gmao/merra2/
# MET_PATH: /css/smapl4/public/L4_SM/OL4001/gph/
# MET_PATH: /css/smapl4/public/L4_SM/Vv4030/gph/
BEG_DATE: 20200601
END_DATE: 20210301
# Specify extremities of lat/lon rectangle (DEFAULT: Global):
#MIN_LON:  -129.9
#MAX_LON:  -60.
#MIN_LAT: 20.
#MAX_LAT: 55.
EOF

#goto PROCESS_S2S

#######################################################################
#                     Batch Parameters and Run Job
#######################################################################

set EXPDIR = `head -n 25 run_rrm.x | grep -m1 WORK_PATH | cut -d':' -f2`

#SBATCH --output=${EXPDIR}/GEOSroute_log.txt
#SBATCH --error=${EXPDIR}/GEOSroute_err_txt
#SBATCH --account=s1583
#SBATCH --time=12:00:00
#SBATCH --ntasks=2
#SBATCH --job-name=run_rrm
#SBATCH --constraint=cssro

#######################################################################
#                     Set up Experiment  Directory
#######################################################################

setenv    MYNAME         `finger $USER | cut -d: -f3 | head -1`
set numprocs = `head -n 35 run_rrm.x | grep ntasks | cut -d'=' -f2`
set  YEARB = `head -n 25 run_rrm.x | grep -m1 BEG_DATE | cut -d':' -f2 | cut -c1-5`
set  YEARE = `head -n 25 run_rrm.x | grep -m1 END_DATE | cut -d':' -f2 | cut -c1-5`
set HOMDIR = `pwd`

mkdir -p $EXPDIR
cd $EXPDIR
mkdir rs output

while ($YEARB <= $YEARE)
mkdir -p rs/$YEARB
mkdir -p output/$YEARB
@ YEARB++
end

cd $PWD

#######################################################################
#                            Submit the job
#######################################################################

umask 022
limit stacksize unlimited
setenv ARCH `uname`

setenv GEOSBIN    $HOMDIR/bin/
source $GEOSBIN/g5_modules

#setenv RUN_CMD "$GEOSBIN/esma_mpirun -np "
setenv RUN_CMD "mpirun -np "
$RUN_CMD $numprocs $GEOSBIN/GEOSroute.x

exit

#######################################################################
#               EXTRACT S2S RUNOFF AND CREATE RRM INPUT FILE
#######################################################################

PROCESS_S2S:

set MMDDM = (dec27 jan31 feb25 mar27 apr26 may31 jun30 jul30 aug29 sep28 oct28 nov27)
set month = 1 
set   PWD = `pwd`
module use -a /discover/swdev/gmao_SIteam/modulefiles-SLES12
module load nco/4.8.1   

while ($month <= 12)

    set yearB = 1981
    set yearE = 2019
    if ($month <= 6) set yearE = 2020
    echo $month $yearB $yearE
    while ($yearB <= $yearE)
        set this_list = /gpfsm/dnb04/projects/p71/aogcm/g5fcst/forecast/production/geos-s2s/runx/$yearB/$MMDDM[$month]/ens1/geosgcm_vis2d/*.geosgcm_vis2d.daily.*.nc4.tar

        mkdir -p tar_dir
        cd tar_dir
        set N_TAR = `echo $#this_list`
        set f = 2
        while ($f <= $N_TAR)
            tar -xvf $this_list[$f]
            @ f++
        end
        exit
        set daily_files = `ls`
        ncrcat -h -v RUNOFF $daily_files $PWD/${yearB}.nc4
        cd $PWD
        nccopy -d6 ${yearB}.nc4 s2s_runoff/$MMDDM[$month].geosgcm_vis2d.${yearB}.nc4
        /bin/rm ${yearB}.nc4
        /bin/rm -r tar_dir
        @ yearB++
    end
    @ month++
end
