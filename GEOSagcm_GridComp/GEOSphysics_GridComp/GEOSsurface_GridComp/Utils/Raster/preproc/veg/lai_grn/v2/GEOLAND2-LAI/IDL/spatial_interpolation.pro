;pro spatial_interpolation
;sp,1,1,/f,/h,/color,/land
;UL=fltarr(36,18,2)
;LR=fltarr(36,18,2)
loadct2,0
ndirs = 25
ncols = 1120
nrows = 1120
resol = 0.00892857142857

nc_global = 40320l
nr_global = 20160l

path = '/l_data/GEOLAND2-LAI/download/extracts/'

openr,1,path+'directory_list5'
wlon =0.
elon =0.
slat =0.
nlat =0.
sdate ='    '
edate ='    '
rdate1='    '
rdate2='    '
dx    =0.
dy    =0.
nck   =0
nrk   =0

for dir_no = 1,ndirs do begin

    lai_global =  fltarr (nc_global,nr_global)
    lai_global(*,*) = !VALUES.F_NAN
    dir_name='      '
    nfiles=0
    
    readf,1,dir_name
    print, 'Processing :',dir_name
    openr,2,path+dir_name+'h5_files'
    readf,2,nfiles
    fname ='    '
    for n =1,nfiles do begin
        readf,2,fname
        openr,3,path+dir_name+strtrim(fname)+'.zip.txt'
        readf,3,wlon
        readf,3,elon
        readf,3,slat
        readf,3,nlat
        readf,3,format='(a10)',rdate1
        readf,3,format='(a10)',rdate2
        readf,3,dx,dy
        readf,3,nrk
        readf,3,nck
        
        if (n eq 1) then begin
            sdate = rdate1
            edate = rdate2
        endif
        
        if ((dx ne resol) or (dy ne resol)) then begin
            print, 'DX or DY /= RESOL :',dir_name+fname, dx,dy
            stop
        endif

        if ((nrk ne nrows) or (nck ne ncols)) then begin
            print, 'NR or NC mismatch :',dir_name+fname, nck,nrk
;            stop
        endif

        if (n gt 1) then begin
            if ((sdate ne rdate1) or (edate ne rdate2)) then begin
                print, 'Date mismatch :',dir_name+fname, rdate1,rdate2
                stop
            endif
        endif

        ; reading HDF5 file

        h5_name = path+dir_name+strtrim(fname)+'.h5'
        file_id=H5F_OPEN(h5_name)
        data_id=H5D_OPEN(file_id,'LAI')
        lai=H5D_READ(data_id)
        latitude_id=H5A_OPEN_NAME(file_id, 'LAT')
        lat_nw = H5A_READ(latitude_id)
        H5A_CLOSE,latitude_id

        longitude_id=H5A_OPEN_NAME(file_id, 'LONG')
        lon_nw = H5A_READ(longitude_id)
        H5A_CLOSE,longitude_id

        if ((lat_nw ne nlat) or (lon_nw ne wlon)) then begin
            print, 'Bunding box mismatch :',h5_name 
            stop 
        endif

        H5D_CLOSE, data_id
        H5F_CLOSE, file_id
        lai_size=size(lai)

        icol = ncols*(wlon + 180)/10 - 1
        irow = nrows*(90 - nlat)/10 - 1

        if(icol eq -1) then icol =0
        if(irow eq -1) then irow =0

        for j =0,lai_size(2) - 1 do begin
            for i=0,lai_size(1) -1 do begin
                if((lai(i,j) ge 0) and (lai(i,j) le 250)) then begin
                    lai_global(icol +i,irow + j) = lai(i,j)/30.
                endif
            endfor
        endfor
        ;
        close,3
    endfor
    close,2
    lai_global = reverse(lai_global,2)
    save,lai_global,file=path+dir_name+'lai.idl'
endfor
close,1
goto, jump

nc_low=2880
nr_low=1440

k1 = nc_global/nc_low
k2 = k1 -1 

pgrid = fltarr (nc_low,nr_low)
for j = 0,nr_low -1 do begin
    for i =0,nc_low-1 do begin
        pgrid (i,j) = mean (lai_global(i*long(k1):i*long(k1)+k2,j*long(k1):j*long(k1)+k2),/nan)
    endfor
endfor
sp,1,3,/f,/h,/color,landscape=0

loadct2,0
levp=indgen(21)*0.5
colors=indgen(21) +80

contour,pgrid,title='LAI',c_colors=colors,levels=levp,xstyle=4,ystyle=4,/CELL_FILL

closeps,'lai.ps'
jump:

end
