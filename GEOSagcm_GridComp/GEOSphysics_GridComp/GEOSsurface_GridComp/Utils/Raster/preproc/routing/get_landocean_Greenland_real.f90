program main
  
  use routing_constant, only : nlon,nlat,nlon_G,nlat_G,loni_min,loni_max,lati_min,lati_max,id_glac,id_lake,id_landend
  
  implicit none
  include 'netcdf.inc'
  
  real*8,allocatable,dimension(:)    :: lon,lat,lon_G,lat_G
  integer,allocatable,dimension(:,:) :: landocean,Greenland
  integer,allocatable,dimension(:)   :: Pfaf_real, countc
  
  integer :: i,j,ret,ncid,varid
  
  allocate(landocean(nlon,nlat))
  allocate(lon(nlon),lat(nlat))
  
  ret=nf_open("inputs/Pfafstetter.nc",0,ncid)
  ret=nf_inq_varid(ncid,"lon",varid)
  ret=nf_get_var_double(ncid,varid,lon)
  ret=nf_close(ncid)
  ret=nf_open("inputs/Pfafstetter.nc",0,ncid)
  ret=nf_inq_varid(ncid,"lat",varid)
  ret=nf_get_var_double(ncid,varid,lat)
  ret=nf_close(ncid)
  ret=nf_open("inputs/Pfafstetter.nc",0,ncid)
  ret=nf_inq_varid(ncid,"data",varid)
  ret=nf_get_var_int(ncid,varid,landocean)
  ret=nf_close(ncid)
  
  
  allocate(Greenland(nlon_G,nlat_G))
  allocate(lon_G(nlon_G),lat_G(nlat_G))
  ret=nf_open("inputs/GreenlandID_30s.nc",0,ncid)
  ret=nf_inq_varid(ncid,"lon",varid)
  ret=nf_get_var_double(ncid,varid,lon_G)
  ret=nf_close(ncid)
  ret=nf_open("inputs/GreenlandID_30s.nc",0,ncid)
  ret=nf_inq_varid(ncid,"lat",varid)
  ret=nf_get_var_double(ncid,varid,lat_G)
  ret=nf_close(ncid)
  ret=nf_open("inputs/GreenlandID_30s.nc",0,ncid)
  ret=nf_inq_varid(ncid,"data",varid)
  ret=nf_get_var_int(ncid,varid,Greenland)
  ret=nf_close(ncid)
  
  
  where(Greenland/=-9999.and.(landocean(loni_min:loni_max,lati_min:lati_max)<=id_landend.or.&
       landocean(loni_min:loni_max,lati_min:lati_max)==id_glac ))&
       landocean(loni_min:loni_max,lati_min:lati_max)=Greenland
  
  
  where(landocean>id_landend.and.landocean<id_lake) landocean=-9999
  where(landocean==id_lake.or.landocean==id_glac) landocean=0
  
  allocate(Pfaf_real(id_landend))
  open(77,file="inputs/Pfaf_real.txt")
  read(77,*)Pfaf_real
  
  do i=1,nlon
     do j=1,nlat
        if(landocean(i,j)<=id_landend.and.landocean(i,j)>=1)then
           landocean(i,j)=Pfaf_real(landocean(i,j))
        else if(landocean(i,j)>=700000000)then
           landocean(i,j)=landocean(i,j)-700000000+291284
        endif
     enddo
  enddo
  
  open(30,file="outputs/Pfafstetter_Greenland_real",form="unformatted")
  do j = 1,nlat
     write (30) landocean(:,j)
  end do
  
end program main
