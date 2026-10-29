#!/bin/csh 

cd /l_data/GEOLAND2-LAI/download/extracts

set PWD=`pwd`

##ls > directory_list2

foreach dir ("`cat directory_list4`")

cd $dir

ls *h5 > h5_files

set nfiles = `wc h5_files | cut -d' ' -f3`

echo $nfiles > temp.file
cat h5_files >> temp.file
sed -e "s/.h5/ /g" temp.file > h5_files

cd $PWD

end
