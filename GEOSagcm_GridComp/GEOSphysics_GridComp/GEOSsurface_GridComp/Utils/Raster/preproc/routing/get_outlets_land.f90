program main

  use routing_constant,only : nc,nl,ng,nlon=>nlon1m,nlat=>nlat1m
  implicit none
  include 'netcdf.inc'
  
  real*8,allocatable  :: lon(:),lat(:),long(:),latg(:),lons(:),lats(:)
  integer,allocatable :: catchind(:,:)
  real,allocatable    :: acah(:,:)
  integer,allocatable :: down(:),sx(:),sy(:),msk(:)
  real,allocatable    :: acas(:)
  
  integer :: id,xi,yi,i,k,xis,yis,ntot,ncid,ret,varid
  
  ntot=nl+ng
  allocate(catchind(nlon,nlat),acah(nlon,nlat))
  allocate(lon(nlon),lat(nlat))
  allocate(sx(nc),sy(nc),acas(nc),down(nc),msk(nc))
  allocate(long(ng),latg(ng),lons(ntot),lats(ntot))
  
  
  ret=nf_open("inputs/CatchIndex.nc",0,ncid)
  ret=nf_inq_varid(ncid,"lon",varid)
  ret=nf_get_var_double(ncid,varid,lon)
  ret=nf_close(ncid)
  ret=nf_open("inputs/CatchIndex.nc",0,ncid)
  ret=nf_inq_varid(ncid,"lat",varid)
  ret=nf_get_var_double(ncid,varid,lat)
  ret=nf_close(ncid)
  
  ret=nf_open("inputs/CatchIndex.nc",0,ncid)
  ret=nf_inq_varid(ncid,"data",varid)
  ret=nf_get_var_int(ncid,varid,catchind)
  ret=nf_close(ncid)
  
  ret=nf_open("inputs/HydroSHEDS_drainage_area.nc",0,ncid)
  ret=nf_inq_varid(ncid,"data",varid)
  ret=nf_get_var_real(ncid,varid,acah)
  ret=nf_close(ncid)
  
  open(77,file="inputs/downstream_1D_new_noadj.txt")
  read(77,*)down
  open(77,file="inputs/Pfaf_msk.txt")
  read(77,*)msk
  
  acas=-9999.
  sx=0
  sy=0
  do xi=1,nlon
     do yi=1,nlat
        if(catchind(xi,yi)>=1)then
           id=catchind(xi,yi)
           if(down(id)==-1.and.acah(xi,yi)>=acas(id))then
              acas(id)=acah(xi,yi)
              sx(id)=xi
              sy(id)=yi
           endif
        endif
     enddo
  enddo
  
  where(down/=-1)sx=-1
  where(down/=-1)sy=-1
  k=0
  do i=1,nc
     if(msk(i)==2)then
        k=k+1
        lons(k)=lon(sx(i))
        lats(k)=lat(sy(i))
     endif
  enddo
  
  open(77,file="inputs/Greenland_outlets_lat.txt")
  read(77,*)latg
  open(77,file="inputs/Greenland_outlets_lon.txt")
  read(77,*)long
  
  lons(k+1:ntot)=long
  lats(k+1:ntot)=latg
  
  
  open(88,file="outputs/outlet_sinklat.txt")
  do i=1,ntot
     write(88,*)lats(i)
  enddo
  open(88,file="outputs/outlet_sinklon.txt")
  do i=1,ntot
     write(88,*)lons(i)
  enddo
  
end program main
