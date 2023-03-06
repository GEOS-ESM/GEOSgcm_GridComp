#!/bin/csh -fx
#=============
#PBS -N lossper_day
#PBS -S /bin/csh
#PBS -W group_list=s1178
#PBS -l select=1:ncpus=16:mpiprocs=16
#PBS -o /discover/nobackup/smahanam//Ganymed-4_0_UNSTABLE/GEOSagcm/src/GEOSgcs_GridComp/GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/Shared/Raster/src/loss_during_day.txt1
#PBS -j oe
#PBS -q general_small
#PBS -l walltime=12:00:00
#=============
cd /discover/nobackup/smahanam//Ganymed-4_0_UNSTABLE/GEOSagcm/src/GEOSgcs_GridComp/GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/Shared/Raster/
setenv OMP_NUM_THREADS 16
setenv MPI_USE_XPMEM
setenv MPI_BUFFER_MAX 2000
setenv MPI_TYPE_MAX   6553600
setenv MPI_MSGS_MAX   10485760

source /discover/nobackup/smahanam//Ganymed-4_0_UNSTABLE/GEOSagcm/src/g5_modules; 
./src/loss_during_day -js 1
