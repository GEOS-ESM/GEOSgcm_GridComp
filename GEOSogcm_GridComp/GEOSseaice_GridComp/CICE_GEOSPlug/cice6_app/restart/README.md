# How to invoke
python3 regrid.py -i iced_tx1_v5.nc -ig tx1_bathy.nc -o iced_mom1_v5.nc -og MAPL_Tripolar.nc -fs True


python3 make_import_internal.py -i iced_mom6_x5.nc -g MAPL_Tripolar.nc -im seaicethermo_import_rst -in seaicethermo_internal_rst -t CF0012x6C_TM0072xTM0036-Pfafstetter.til -d yyymmdd 
