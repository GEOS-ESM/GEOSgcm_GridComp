function grid,data,lat,lon,region=rgn,nlat=nlat,nlon=nlon, $
                   npts=npts,mis_val=mval

; function to [re]grid data given lat and lon position of the elements. 
; Output is a 2D array on a lat/lon grid, holding means of the elements that fell 
; into each gridbox. Default regridded data resolution is 1deg (both lat/lon)

; Optional arguments include: 
; - number of points along lat-direction (nlat) and lon-direction (nlon) - controles the resolution of the output
; - a user-defined missing value (mval) - defines missing value
; - a user-defined region to [re]grid over. [min_lat,mat_lat, min_lon, max_lon]
; - flag '/npts'. Will add the number of elements per grid box to the output 

if (n_elements(rgn)  EQ 0) then rgn=[-90.0,90.0,-180.0,180.0]
if (n_elements(nlat) EQ 0) then nlat=180
if (n_elements(nlon) EQ 0) then nlon=360
if (n_elements(mval) EQ 0) then mval=-999.0
if (n_elements(eps)  EQ 0) then eps=0.1
if max(lon) gt 180. then lon[where(lon gt 180.,/null)]=lon[where(lon gt 180.,/null)]-360.
minlat = rgn[0]
maxlat = rgn[1]
minlon = rgn[2]
maxlon = rgn[3]
xinc=(maxlon-minlon)/float(NLON)
yinc=(maxlat-minlat)/float(NLAT)
ngrd=lonarr(NLON,NLAT)
rgrd=fltarr(NLON,NLAT)
sgrd=fltarr(NLON,NLAT)

; Eliminate missing data

lat2=lat
lon2=lon
data2=data
m=where(abs(data2-mval) GT eps AND abs(lat2-mval) GT eps AND abs(lon2-mval) GT eps,cnt)
if (cnt GT 0) then begin
  lat2=lat[m]
  lon2=lon[m]
  data2=data[m]
endif

; Determine if data is inside specified grid region

xind=long((lon2 - minlon)/double(xinc))
yind=long((lat2 - minlat)/double(yinc))
m=where(xind GE 0 AND xind LT NLON AND yind GE 0 AND yind LT NLAT,cnt)

lat2=lat2[m]
lon2=lon2[m]
data2=data2[m]
xind=xind[m]
yind=yind[m]

ind=xind+NLON*yind
h=histogram(ind,MIN=0,MAX=long(NLON)*long(NLAT),locations=x)
m=where(h GT 0,mcnt)
for i=0L,mcnt-1 do begin
  ind2=where(ind EQ x[m[i]],ncnt)
  ngrd[m[i]] = ngrd[m[i]] + ncnt
  rgrd[m[i]] = rgrd[m[i]] + total(data2[ind2])
endfor

; Compute grid averages

ind=where(ngrd EQ 0,cnt)
if (cnt GT 0) then rgrd[ind] = mval
ind=where(ngrd GT 0,cnt)
if (cnt GT 0) then rgrd[ind] = rgrd[ind] / float(ngrd[ind])

n=1
if (keyword_set(npts))  then n=n+1
a=fltarr(NLON,NLAT,n)
a[*,*,0] = rgrd
n=1
if (keyword_set(npts)) then begin
  a[*,*,n] = float(ngrd)
  n = n + 1
endif
return,a
end
