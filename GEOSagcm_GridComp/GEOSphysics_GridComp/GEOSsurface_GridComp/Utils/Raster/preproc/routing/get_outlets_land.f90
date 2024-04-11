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
  
  character(len=100) :: file_path1 !/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/land/topo/v1/SRTM-TopoData/SRTM_PfafData.nc
  character(len=100) :: file_path2 !/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/routing/HydroSHEDS_drainage_area.nc
  character(len=100) :: file_path3 !/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/routing/Greenland_outlets_lat.txt
  character(len=100) :: file_path4 !/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/routing/Greenland_outlets_lon.txt  

  if (command_argument_count() /= 4) then
      print *, "no <file_path1> <file_path2> <file_path3> <file_path4> found"
      stop
  endif
  call get_command_argument(1, file_path1)
  call get_command_argument(2, file_path2)
  call get_command_argument(3, file_path3)
  call get_command_argument(4, file_path4)

  ntot=nl+ng
  allocate(catchind(nlon,nlat),acah(nlon,nlat))
  allocate(lon(nlon),lat(nlat))
  allocate(sx(nc),sy(nc),acas(nc),down(nc),msk(nc))
  allocate(long(ng),latg(ng),lons(ntot),lats(ntot))
  
  
  ret=nf_open(file_path1,0,ncid)
  ret=nf_inq_varid(ncid,"longitude",varid)
  ret=nf_get_var_double(ncid,varid,lon)
  ret=nf_close(ncid)
  ret=nf_open(file_path1,0,ncid)
  ret=nf_inq_varid(ncid,"latitude",varid)
  ret=nf_get_var_double(ncid,varid,lat)
  ret=nf_close(ncid)
  
  ret=nf_open(file_path1,0,ncid)
  ret=nf_inq_varid(ncid,"CatchIndex",varid)
  ret=nf_get_var_int(ncid,varid,catchind)
  ret=nf_close(ncid)
  
  ret=nf_open(file_path2,0,ncid)
  ret=nf_inq_varid(ncid,"data",varid)
  ret=nf_get_var_real(ncid,varid,acah)
  ret=nf_close(ncid)
  
  open(77,file="outputs/Pfaf_downid.txt")
  read(77,*)down
  open(77,file="outputs/Pfaf_msk.txt")
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
  
  open(77,file=file_path3)
  read(77,*)latg
  open(77,file=file_path4)
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
