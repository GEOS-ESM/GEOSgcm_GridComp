begin

ns=22612
nlat=21600
nlon=43200

lat=asciiread("outputs/outlet_sinklat.txt",ns,"double")
lon=asciiread("outputs/outlet_sinklon.txt",ns,"double")

lat30s=asciiread("inputs/lat_30s.txt",nlat,"double")
lon30s=asciiread("inputs/lon_30s.txt",nlon,"double")

lati=ind_nearest_coord(lat,lat30s,0)
loni=ind_nearest_coord(lon,lon30s,0)

asciiwrite("outputs/outlet_sinky.txt",lati+1)
asciiwrite("outputs/outlet_sinkx.txt",loni+1)

end
