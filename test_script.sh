#!/bin/bash

EXE=TEST_MOIST

#for C in c24_data c90_data c180_data;
for C in c180_data;
do
	for test in aer_activation buoyancy cup_gf_GEOS5 cup_gf_GF2020 cup_gf_sh fillq2zero gfdl_cloud_microphys_driver hystpdf evap_subl_pdf_loop radcoup_loop
	do
		for index in {0..5}
		do
			echo "******************************"
			echo ./$EXE ./$C/$test/ $index
			./$EXE ./$C/$test/ $index
		done
	done
done
