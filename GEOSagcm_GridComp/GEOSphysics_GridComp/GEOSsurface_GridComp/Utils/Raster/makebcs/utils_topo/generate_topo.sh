#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1 --ntasks-per-node=126
#SBATCH --job-name=topo_intel
#SBATCH --constraint=mil
#SBATCH --qos=benchmark
#SBATCH --partition=preops
#SBATCH --account=g0620
#SBATCH --mail-type=ALL
#SBATCH --no-requeue

# put path to path to bin here
export TOPODIR=
source ${TOPODIR}/g5_modules.sh

echo "STARTING"
# generate intermediate cube
raw_latlon_data="/discover/nobackup/bmauer/gmted_topo/gmted_fix_superior/gmted_fixed_anartica_superior_caspian.nc4"
intermediate_cube="c3000.gmted_fixedanarticasuperior.nc"
source_topo="gmted_intel"
 
cat << _EOF_ > bin_to_cube.nl
&binparams
  raw_latlon_data_file='$raw_latlon_data'
  output_file='$intermediate_cube'
  ncube=3000
/
_EOF_

if [[ ! -e landm_coslat.nc ]]; then
   ln -s $TOPODIR/landm_coslat.nc landm_coslat.nc 
fi

${TOPODIR}/bin_to_cube.x

res=("12" "24" "48" "90" "180" "360" "720" "1120" "1440" "2880")
cutoff=25
smoothmap[12]="773.91"
smoothmap[24]="386.52"
smoothmap[48]="193.07"
smoothmap[90]="102.91"
smoothmap[180]="51.44"
smoothmap[360]="25.71"
smoothmap[720]="12.86"
smoothmap[1120]="8.26"
smoothmap[1440]="6.43"
smoothmap[2880]="3.21"

for n in "${res[@]}";
do

   let jm=$n*6
   echo $n
   echo ${smoothmap[$n]}

   export output_dir=output_${n}
   if [[ ! -e $output_dir ]]; then
      mkdir $output_dir
   fi

   config_file=GenScrip.rc
   let jm=$n*6
   echo $n
   echo $jm
   scripfile=PE${n}x${jm}-CF.nc4
   echo $scripfile
   echo $config_file
   cat << _EOF_ > ${config_file}
CUBE_DIM: $n
output_scrip: ${scripfile}
output_geos: c${n}_coords.nc4
_EOF_
   mpirun -np 6 ${TOPODIR}/generate_scrip_cube_topo.x
   rm GenScrip.rc
   output_grid=PE${im}x${jm}-CF
   if (( $res < $cutoff )); then
      ${TOPODIR}/cube_to_target.x --grid_descriptor_file="PE${n}x${jm}-CF.nc4" --intermediate_cs_name=${intermediate_cube} --output_grid="PE${n}x${jm}-CF.nc4" --output_data_directory=${output_dir} --smoothing_scale=${smoothmap[$n]} --name_email_of_creator='gmao' --fine_radius=0 --output_grid=${output_grid} --source_data_identifier=${source_topo} --jmax_segments=100000
   else
      ${TOPODIR}/cube_to_target.x --grid_descriptor_file="PE${n}x${jm}-CF.nc4" --intermediate_cs_name=${intermediate_cube} --output_grid="PE${n}x${jm}-CF.nc4" --output_data_directory=${output_dir} --smoothing_scale=${smoothmap[$n]} --name_email_of_creator='gmao' --fine_radius=0 --output_grid=${output_grid} --source_data_identifier=${source_topo}
   fi
   rm $scripfile

   #convert to gmao
   cd ${output_dir}
   arr=(*.nc)
   echo ${arr[0]}
   ${TOPODIR}/scrip_to_restart_topo.py -i ${arr[0]} -o gwd_internal_rst
   ${TOPODIR}/convert_to_gmao_output_topo.x -i ${arr[0]} --im $n --jm$n 
   cd ..
done
