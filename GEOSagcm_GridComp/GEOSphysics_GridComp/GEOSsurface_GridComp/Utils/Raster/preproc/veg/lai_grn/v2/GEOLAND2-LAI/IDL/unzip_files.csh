#!/bin/csh 

setenv numB 679

setenv numE 701

set PWD=`pwd`

while ($numB <= $numE)

cd catftp.vgt.vito.be/M0000${numB}/

echo Processing M0000${numB}
  
ls *.zip > zip_list

foreach line ("`cat zip_list`")

### getting date_stamp and file tail of the xml file
set date_stamp=`echo $line | cut -d'_' -f4 | cut -c1-8`
set ftail=`echo $line | cut -d'.' -f1 |  cut -c15-200`

unzip $line "*.h5" -q -d /l_data/GEOLAND2-LAI/download/extracts/
unzip $line "*.xml" -q -d /l_data/GEOLAND2-LAI/download/extracts/

## WestBound Longitude
set ln = `grep -n "<westBoundLongitude>" /l_data/GEOLAND2-LAI/download/extracts/${date_stamp}/g2_BIOPAR_LAI_PROD-DESC_${ftail}.3.xml`
set lnp1=`echo $ln | cut -d':' -f1`
@ lnp1 = $lnp1 + 1 
set m1=`sed -n ${lnp1}p /l_data/GEOLAND2-LAI/download/extracts/${date_stamp}/g2_BIOPAR_LAI_PROD-DESC_${ftail}.3.xml | cut -d'>' -f2 | cut -d'<' -f1`

## EastBound Longitude
set ln = `grep -n "<eastBoundLongitude>" /l_data/GEOLAND2-LAI/download/extracts/${date_stamp}/g2_BIOPAR_LAI_PROD-DESC_${ftail}.3.xml`
set lnp1=`echo $ln | cut -d':' -f1`
@ lnp1 = $lnp1 + 1 
set m2=`sed -n ${lnp1}p /l_data/GEOLAND2-LAI/download/extracts/${date_stamp}/g2_BIOPAR_LAI_PROD-DESC_${ftail}.3.xml | cut -d'>' -f2 | cut -d'<' -f1`

## SouthBound Latitude
set ln = `grep -n "<southBoundLatitude>" /l_data/GEOLAND2-LAI/download/extracts/${date_stamp}/g2_BIOPAR_LAI_PROD-DESC_${ftail}.3.xml`
set lnp1=`echo $ln | cut -d':' -f1`
@ lnp1 = $lnp1 + 1 
set p1=`sed -n ${lnp1}p /l_data/GEOLAND2-LAI/download/extracts/${date_stamp}/g2_BIOPAR_LAI_PROD-DESC_${ftail}.3.xml | cut -d'>' -f2 | cut -d'<' -f1`

## NorthBound Latitude
set ln = `grep -n "<northBoundLatitude>" /l_data/GEOLAND2-LAI/download/extracts/${date_stamp}/g2_BIOPAR_LAI_PROD-DESC_${ftail}.3.xml`
set lnp1=`echo $ln | cut -d':' -f1`
@ lnp1 = $lnp1 + 1 
set p2=`sed -n ${lnp1}p /l_data/GEOLAND2-LAI/download/extracts/${date_stamp}/g2_BIOPAR_LAI_PROD-DESC_${ftail}.3.xml | cut -d'>' -f2 | cut -d'<' -f1`

##Time Slice
set T1 = `grep "beginPosition" /l_data/GEOLAND2-LAI/download/extracts/${date_stamp}/g2_BIOPAR_LAI_PROD-DESC_${ftail}.3.xml | cut -d'>' -f2 | cut -d'T' -f1`
set T2 = `grep "endPosition" /l_data/GEOLAND2-LAI/download/extracts/${date_stamp}/g2_BIOPAR_LAI_PROD-DESC_${ftail}.3.xml | cut -d'>' -f2 | cut -d'T' -f1`

##Angle
set angle = `grep "Angle" /l_data/GEOLAND2-LAI/download/extracts/${date_stamp}/g2_BIOPAR_LAI_PROD-DESC_${ftail}.3.xml | cut -d'>' -f2 | cut -d'<' -f1`

##Dimesions
set ln = `grep -n "<dimensionSize>" /l_data/GEOLAND2-LAI/download/extracts/${date_stamp}/g2_BIOPAR_LAI_PROD-DESC_${ftail}.3.xml`
set ln1=`echo $ln | cut -d':' -f1`
set ln2=`echo $ln | cut -d':' -f2 | cut -d'>' -f2 | cut -d' ' -f2`
@ ln1 = $ln1 + 1
@ ln2 = $ln2 + 1
set rows=`sed -n ${ln1}p /l_data/GEOLAND2-LAI/download/extracts/${date_stamp}/g2_BIOPAR_LAI_PROD-DESC_${ftail}.3.xml | cut -d'>' -f2 | cut -d'<' -f1`
set cols=`sed -n ${ln2}p /l_data/GEOLAND2-LAI/download/extracts/${date_stamp}/g2_BIOPAR_LAI_PROD-DESC_${ftail}.3.xml | cut -d'>' -f2 | cut -d'<' -f1`

cat << _EOF_ > /l_data/GEOLAND2-LAI/download/extracts/${date_stamp}/${line}.txt
${m1}                : WestBound Longitude
${m2}                : EastBound Longitude
${p1}                : SouthBound Latitude
${p2}                : NorthBound Latitude
${T1}                : Start Date
${T2}                : End Date
${angle}             : Resolution
${rows}              : Number of rows
${cols}              : Number of columns
_EOF_

end

cd $PWD

@ numB = $numB    +  1

end 


