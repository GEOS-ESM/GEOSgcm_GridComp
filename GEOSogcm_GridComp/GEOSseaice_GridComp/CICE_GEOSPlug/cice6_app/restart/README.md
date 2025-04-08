# How to invoke

```python
python3 regrid.py -i input/iced_tx1_v5.nc -ig input/MAPL_Tripolar.nc -o output/iced_mom1_v5.nc -og output/MAPL_Tripolar.nc
```

Above command will generate restart with prognostic salinity if the original file has them.

If a fixed salinity profile (like BL99) is desired while the original restart has prognostic ones, add ```-fs``` switch, e.g.

```python
python3 regrid.py -i input/iced_tx1_v5.nc -ig input/MAPL_Tripolar.nc -o output/iced_mom1_v5.nc -og output/MAPL_Tripolar.nc -fs
```    
 
```python
python3 make_import_internal.py -i iced_mom6_x5.nc -g MAPL_Tripolar.nc -im seaicethermo_import_rst -in seaicethermo_internal_rst -t CF0012x6C_TM0072xTM0036-Pfafstetter.til -d yyymmdd 
```

```python
python3 cice4_to_cice6.py -i seaicethermo_internal_rst  -ig tilefile -o iced.nc -ot iced_template.nc -og MAPL_Tripolar.nc
```
