program main
  
  use routing_constant, only : nc,nlon,nlat,nlon_G,nlat_G,loni_min,loni_max,lati_min,lati_max
  
  implicit none
  include 'netcdf.inc'
  
  real*8,allocatable,dimension(:)    :: lon_G,lat_G
  integer,allocatable,dimension(:,:) :: landocean,Greenland
  integer,allocatable,dimension(:)   :: Pfaf_real, countc
  
  integer :: i,j,ret,ncid,varid,ntile,id_glac,id_lake,id_landend
  real    :: val(4)
  
  allocate(landocean(nlon,nlat))
  open(77,file="/discover/nobackup/projects/gmao/bcs_shared/legacy_bcs/Icarus-NLv5/Icarus-NLv5_EASE/SMAP_EASEv2_M09/rst/Pfafstetter.rst",form="unformatted",status="old")
  do j=1,nlat
    read(77) landocean(:,j)
  enddo
  close(77)   

  open(77,file="/discover/nobackup/projects/gmao/bcs_shared/legacy_bcs/Icarus-NLv5/Icarus-NLv5_EASE/SMAP_EASEv2_M09/til/Pfafstetter.til",form="formatted",status="old")
  read(77,*) ntile
  id_glac=ntile
  id_lake=ntile-1
  id_landend=ntile-3
  allocate(Pfaf_real(id_landend))  
  do i=1,4
    read(77,*)
  enddo
  do i=1,id_landend
    read(77,*) val(1:4),Pfaf_real(i)
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
  
  
  where(Greenland/=-9999.and.(landocean(loni_min:loni_max,lati_min:lati_max)<=id_landend.or.&
       landocean(loni_min:loni_max,lati_min:lati_max)==id_glac ))&
       landocean(loni_min:loni_max,lati_min:lati_max)=Greenland
  
  
  where(landocean>id_landend.and.landocean<id_lake) landocean=-9999
  where(landocean==id_lake.or.landocean==id_glac) landocean=0
  
  
  do i=1,nlon
     do j=1,nlat
        if(landocean(i,j)<=id_landend.and.landocean(i,j)>=1)then
           landocean(i,j)=Pfaf_real(landocean(i,j))
        else if(landocean(i,j)>=700000000)then
           landocean(i,j)=landocean(i,j)-700000000+nc
        endif
     enddo
  enddo
  
  open(30,file="outputs/Pfafstetter_Greenland_real",form="unformatted")
  do j = 1,nlat
     write (30) landocean(:,j)
  end do
  
end program main
