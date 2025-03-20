program main

use river_read
implicit none

integer :: nc=291284
real :: cur_avg=1.4
real :: cur_min=0.5
real :: cur_max=5.

integer,parameter :: nlon=21600
integer,parameter :: nlat=10800
real*8,allocatable :: lon(:),lat(:)
real,allocatable :: ldn1m(:,:),elev1m(:,:)
integer,allocatable :: catid(:,:),flag_slp(:)

integer,parameter :: nlonh=86400
integer,parameter :: nlath=33600
real*8,allocatable :: lonh(:),lath(:)
real,allocatable :: ldnh(:,:),elev_15s(:,:)

real,allocatable,dimension(:) :: lon_dn,lat_dn,lon_up,lat_up,dist_ref,dist_ref2,ldn_min,ldn_max,riv_len,str_len,slp
real,allocatable,dimension(:) :: lon_min,lon_max,lat_min,lat_max,area,elevdiff_ref,elevdiff
integer,allocatable,dimension(:) :: xi_min,yi_min,xi_max,yi_max
integer,allocatable,dimension(:) :: downid



integer xi,yi
integer :: num,i,j,cid,did,k
integer :: data1,data12
real*8  :: data2
real    :: data7,data9,data10
real    :: elev_temp

!-----------------------------------------------------------------------
!Regrid LDN from HydroSHEDS

allocate(ldn1m(nlon,nlat),catid(nlon,nlat))
allocate(lon(nlon),lat(nlat))
call read_ncfile_double1d("input/SRTM_PfafData.nc","longitude",lon,nlon)
call read_ncfile_double1d("input/SRTM_PfafData.nc","latitude",lat,nlat)
call read_ncfile_int2d("input/SRTM_PfafData.nc","CatchIndex",catid,nlon,nlat)
ldn1m=-1.
where(catid==-9999) ldn1m=-9999.

allocate(ldnh(nlonh,nlath))
call read_ncfile_real2d("input/hyd_glo_ldn_15s.nc","Band1",ldnh,nlonh,nlath)
where(ldnh.lt.4.e9) ldnh=ldnh/1.e3 !m -> km

do xi=1,nlon
  do yi=2041,10440
    if(ldn1m(xi,yi).ne.-9999.)then
      ldn1m(xi,yi)=minval(ldnh(4*xi-3:4*xi,4*yi-3-8160:4*yi-8160))
      if(ldn1m(xi,yi).gt.4.e9)ldn1m(xi,yi)=-1. 
    end if
  enddo
enddo
print *,maxval(ldn1m)

allocate(ldn_min(nc),ldn_max(nc),xi_min(nc),yi_min(nc),xi_max(nc),yi_max(nc))
ldn_min=1.e20
ldn_max=-9999.
xi_min=-9999;yi_min=-9999;xi_max=-9999;yi_max=-9999
do i=1,nlon
  do j=1,nlat
    if(catid(i,j)>=1)then
      cid=catid(i,j)
      if(ldn1m(i,j)>0. .and. ldn1m(i,j)<ldn_min(cid))then
        ldn_min(cid)=ldn1m(i,j)
        xi_min(cid)=i
        yi_min(cid)=j
      endif
      if(ldn1m(i,j)>0. .and. ldn1m(i,j)>ldn_max(cid))then
        ldn_max(cid)=ldn1m(i,j)
        xi_max(cid)=i
        yi_max(cid)=j        
      endif      
    endif
  enddo
enddo
where(ldn_min==1.e20)ldn_min=-9999

!open(88,file="xi_yi_min.txt")
!do i=1,nc
!  write(88,*)xi_min(i),yi_min(i)
!enddo


allocate(elev_15s(nlonh,nlath),elev1m(nlon,nlat))
call read_ncfile_real2d("input/hyd_glo_dem_15s.nc","Band1",elev_15s,nlonh,nlath)
where(elev_15s>30000.)elev_15s=0.
elev1m=0.
do xi=1,nlon
  do yi=2041,10440
      elev1m(xi,yi)=sum(elev_15s(4*xi-3:4*xi,4*yi-3-8160:4*yi-8160))/16.
  enddo
enddo

!call create_ncfile_real2d("elev_1m.nc","data",elev1m,lon,lat,nlon,nlat)

deallocate(ldnh,elev_15s)
!-----------------------------------------------------------------------
!Get reference distance

open(77,file="input/Pfafcatch-routing.dat", form="formatted", status="old")
read(77,*)num
allocate(lon_dn(nc),lat_dn(nc),lon_up(nc),lat_up(nc),dist_ref(nc),dist_ref2(nc))
allocate(lon_min(nc),lon_max(nc),lat_min(nc),lat_max(nc),area(nc),elevdiff_ref(nc),elevdiff(nc))

do i=1,nc
  read(77,*)data1,data2,lon_min(i),lon_max(i),lat_min(i),lat_max(i),data7,area(i),data9,data10,elevdiff_ref(i),data12,lon_dn(i),lat_dn(i),lon_up(i),lat_up(i)
enddo

do i=1,nc
  dist_ref(i)=spherical_distance(lon_dn(i), lat_dn(i), lon_up(i), lat_up(i))
  dist_ref2(i)=spherical_distance(lon_min(i), lat_min(i), lon_max(i), lat_max(i))
enddo
where(dist_ref>dist_ref2.or.dist_ref==0.)dist_ref=0.5*dist_ref2


!--------------------------------------------------------------------
! Get intial guess of river length
allocate(riv_len(nc),downid(nc),flag_slp(nc))
open(77,file="output/downstream_1D_new_noadj.txt")
read(77,*)downid

!open(88,file="temp/xi_yi_min.txt")
!do i=1,nc
!  write(88,*)xi_min(i),yi_min(i)
!enddo
!open(88,file="temp/xi_yi_max.txt")
!do i=1,nc
!  write(88,*)xi_max(i),yi_max(i)
!enddo


flag_slp=1

riv_len=-9999.
elevdiff=-9999.
do i=1,nc
  if(downid(i)>=1)then
    did=downid(i)
    if(.not. (riv_len(did)>=cur_min*dist_ref(did).and.riv_len(did)<=cur_max*dist_ref(did)) )then
      riv_len(did)=ldn_min(i)-ldn_min(did)
      if(xi_min(i)>0.and.xi_min(did)>0)then
        elevdiff(did)=max(0.,elev1m(xi_min(i),yi_min(i)) - elev1m(xi_min(did),yi_min(did)))
        flag_slp(did)=1
      else
        elevdiff(did)=elevdiff_ref(did)
        flag_slp(did)=0
      endif
    else if(flag_slp(did)==0.or.elevdiff(did)==0.)then
      riv_len(did)=ldn_min(i)-ldn_min(did)
      if(xi_min(i)>0.and.xi_min(did)>0)then
        elevdiff(did)=max(0.,elev1m(xi_min(i),yi_min(i)) - elev1m(xi_min(did),yi_min(did)))
        flag_slp(did)=1
      else
        elevdiff(did)=elevdiff_ref(did)
        flag_slp(did)=0
      endif
    endif
  endif
enddo

do i=1,nc
  if(riv_len(i)==-9999.)then
    riv_len(i)=(ldn_max(i)-ldn_min(i))*0.5
    if(xi_min(i)>0)then
      elevdiff(i)=max(0.,0.5*elev1m(xi_max(i),yi_max(i)) - 0.5*elev1m(xi_min(i),yi_min(i)) )
    else
      elevdiff(i)=elevdiff_ref(i)
      flag_slp(did)=0
    endif
  endif
enddo

k=0
do i=1,nc
  if(.not. (riv_len(i)>=cur_min*dist_ref(i).and.riv_len(i)<=cur_max*dist_ref(i)) )then
    riv_len(i)=cur_avg*dist_ref(i)
    elevdiff(i)=elevdiff_ref(i)
    flag_slp(i)=0
    k=k+1
  endif
enddo
open(88,file="output/Pfaf_lriv_PR.txt")
do i=1,nc
  write(88,*) riv_len(i)
enddo



!--------------------------------------------------------------------
! Calculate the length scale of local streams
allocate(str_len(nc))
str_len=area/riv_len/4.*cur_avg
open(88,file="output/Pfaf_lstr_PR.txt")
do i=1,nc
  write(88,*) str_len(i)
enddo
!--------------------------------------------------------------------
! Calculate the Catchment slope
allocate(slp(nc))
slp=elevdiff*1.e-3/riv_len
where(slp.lt.1.e-5) flag_slp=0
where(slp.lt.1.e-5) slp=1.e-5
print *,sum(flag_slp)
open(88,file="temp/Pfaf_slope.txt")
do i=1,nc
  write(88,*) slp(i)
enddo
print *,minval(slp)
open(88,file="temp/Pfaf_slope_flag.txt")
do i=1,nc
  write(88,*)flag_slp(i)
enddo

!--------------------------------------------------------------------
contains

function spherical_distance(lon_dn, lat_dn, lon_up, lat_up) result(distance)
    implicit none
    ! Declare variables
    real, intent(in) :: lon_dn, lat_dn   ! Input coordinates (downstream point)
    real, intent(in) :: lon_up, lat_up   ! Input coordinates (upstream point)
    real :: distance                     ! Output distance (in kilometers)
    real :: R, dlon, dlat, a, c          ! Intermediate variables

    ! Radius of the Earth in kilometers
    R = 6371.0

    ! Convert degrees to radians
    dlon = (lon_up - lon_dn) * (acos(-1.0) / 180.0)
    dlat = (lat_up - lat_dn) * (acos(-1.0) / 180.0)

    ! Haversine formula
    a = sin(dlat / 2.0)**2 + cos(lat_dn * (acos(-1.0) / 180.0)) * &
        cos(lat_up * (acos(-1.0) / 180.0)) * sin(dlon / 2.0)**2
    c = 2.0 * atan2(sqrt(a), sqrt(1.0 - a))

    ! Distance calculation
    distance = R * c

end function spherical_distance

end program