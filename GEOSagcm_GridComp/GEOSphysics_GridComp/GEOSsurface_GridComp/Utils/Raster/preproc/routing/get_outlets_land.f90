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

! Get sink points on land or in Greenland (from Lauren Andrews) by picking the point (i.e., 1 minute grid cell) within each sink catchment that has the largest drainage area per the HydroSHEDS (https://www.hydrosheds.org/) dataset.
! acah = map of drainage area with a resolution of 1m from HydroSHEDS.
! down = catchment index of the downstream catchment.
! sx = the lon index (on a 1m map) of the outlet point (>0 only for sink catchments).
! sy = the lat index (on a 1m map) of the outlet point (>0 only for sink catchments). 
! msk: 1 = has downstream catchment; 2 = drains to ocean; 3 = drains to inland lake
! acas = maximum drainage area, defined (>0) only for sink catchments 
! ntot = number of the total outlets (including Greenland) 
! nl = number of outlets to ocean in land (not including Greenland) 
! ng = number of outlets to ocean in Greenland

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

! For each long/lat location, determine if the catchment holding it is an outlet catchment to the ocean, 
! and if so, determine if this point has the maximum drainage area  
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

! Construct arrays of longitudes and latitudes of sink points.  
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

! Append Greenland values to the longitude and latitude arrays.  
  open(77,file=file_path3)
  read(77,*)latg
  open(77,file=file_path4)
  read(77,*)long
  
  lons(k+1:ntot)=long
  lats(k+1:ntot)=latg
  
! Write out arrays of sink point longitudes and latitudes
  open(88,file="outputs/outlet_sinklat.txt")
  do i=1,ntot
     write(88,*)lats(i)
  enddo
  open(88,file="outputs/outlet_sinklon.txt")
  do i=1,ntot
     write(88,*)lons(i)
  enddo
  
end program main
