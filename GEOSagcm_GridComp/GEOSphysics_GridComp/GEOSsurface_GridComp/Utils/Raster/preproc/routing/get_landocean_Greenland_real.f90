program main
  
  use routing_constant, only : nc,nlon,nlat,nlon_G,nlat_G,loni_min,loni_max,lati_min,lati_max,&
                               nlon1m,nlat1m
  
  implicit none
  include 'netcdf.inc'
  
  real*8,allocatable,dimension(:)    :: lon_G,lat_G
  integer,allocatable,dimension(:,:) :: catchind,landocean,Greenland
  integer,allocatable,dimension(:)   :: Pfaf_real, countc
  
  integer :: i,j,ret,ncid,varid,ntile
  real    :: val(4)
  
  allocate(catchind(nlon1m,nlat1m))
  ret=nf_open("/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/land/topo/v1/SRTM-TopoData/SRTM_PfafData.nc",0,ncid)
  ret=nf_inq_varid(ncid,"CatchIndex",varid)
  ret=nf_get_var_int(ncid,varid,catchind)
  ret=nf_close(ncid)
  allocate(landocean(nlon,nlat))
  landocean=-9999
  do i=1,nlon1m
    do j=1,nlat1m
      if(catchind(i,j)/=-9999)then
        landocean(2*i-1:2*i,2*j-1:2*j)=catchind(i,j)
      endif
    enddo
  enddo


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
  
  
  where(Greenland/=-9999) landocean(loni_min:loni_max,lati_min:lati_max)=Greenland
  
  
  do i=1,nlon
     do j=1,nlat
        if(landocean(i,j)>=700000000)then
           landocean(i,j)=landocean(i,j)-700000000+nc
        endif
     enddo
  enddo
  
  open(30,file="outputs/Pfafstetter_Greenland_real",form="unformatted")
  do j = 1,nlat
     write (30) landocean(:,j)
  end do
  
end program main
