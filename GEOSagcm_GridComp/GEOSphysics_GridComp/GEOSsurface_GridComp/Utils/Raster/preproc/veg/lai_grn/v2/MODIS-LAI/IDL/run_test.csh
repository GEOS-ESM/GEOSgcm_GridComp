#!/bin/csh
#$ -cwd
#$ -m e


setenv filename /l_data/MODIS/LAI/2000.02.18/MOD15A2.A2000049.h35v10.005.2006268212134.hdf

echo $filename

/usr/local/bin/idl <<EOB

filename=GETENV('filename')

print,filename 

test_geos5,filename

exit
EOB
