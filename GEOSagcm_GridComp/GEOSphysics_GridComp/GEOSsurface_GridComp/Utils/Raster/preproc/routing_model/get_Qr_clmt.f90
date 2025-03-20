program main

  use omp_lib
  use river_read
  implicit none

  integer,parameter   :: nlat=1624,nlon=3856,nc=291284,nupmax=34
  character(len=500)  :: filename="input/SMAPL4_OL7000_runoff_mean_2016_2023.nc"
  real,allocatable    :: runoff(:,:),qrunf(:),temp(:,:),qri(:),qin(:)
  integer,allocatable :: nts(:),downid(:),upstream(:,:)

  integer :: i,j,nmax,did
  
  allocate(runoff(nlon,nlat),qrunf(nc),temp(nlon,nlat))
  call read_ncfile_real2d(trim(filename),"mean_runoff_flux",runoff,nlon,nlat)
  where(runoff==-9999.)runoff=0.
  temp=runoff(:,nlat:1:-1)
  runoff=temp
  runoff=runoff*86400. !mm/d
  qrunf=M09_to_cat(runoff,nlon,nlat,nc) !kg/s
  qrunf=qrunf/1.e3 !m3/s
  open(88,file="output/Pfaf_qstr.txt")
  do i=1,nc
    write(88,*)qrunf(i)
  enddo

  allocate(nts(nc),downid(nc),qri(nc))
  open(77,file="output/Pfaf_tosink.txt")
  read(77,*)nts
  open(77,file="output/downstream_1D_new_noadj.txt")  
  read(77,*)downid

  nmax=maxval(nts)
  qri=qrunf
  do j=nmax,1,-1
    do i=1,nc  
    if(nts(i)==j)then
      did=downid(i)
      qri(did)=qri(did)+qri(i)
    endif
    enddo
  enddo
  open(88,file="output/Pfaf_qri.txt")
  do i=1,nc
    write(88,*)qri(i)
  enddo  

  allocate(upstream(nupmax,nc),qin(nc))
  open(77,file="output/upstream_1D.txt")
  read(77,*)upstream
  qin=-9999.
  where(upstream(1,:)/=-1)qin=qri-qrunf
  where(upstream(1,:)==-1)qin=qrunf/2.
  open(88,file="output/Pfaf_qin.txt")
  do i=1,nc
    write(88,*)qin(i)
  enddo   

contains
!------------------------------------------------------------------------------
! This function maps runoff data from M09 resolution to catchments (cat)
function M09_to_cat(runoff,nlon,nlat,ncat) result(Qrunf)

  integer,intent(in) :: nlon,nlat,ncat    ! Input: number of longitude, latitude, and catchments
  real,intent(in) :: runoff(nlon,nlat)  ! Input: runoff array of size (nlon, nlat)
  real  :: Qrunf(ncat)                  ! Output: runoff mapped to catchments

  real,parameter :: small=1.e-12        ! Small value to avoid division by zero

  integer,parameter :: nmax=458           ! Maximum number of sub-areas per catchment
  integer,parameter :: nc=291284          ! Total number of catchments

  real,allocatable,dimension(:,:) :: subarea,frac  ! Arrays for sub-area and fractions
  integer,allocatable,dimension(:,:) :: subx,suby    ! Arrays for x and y coordinates of sub-areas
  real,allocatable,dimension(:) :: tot,runfC,fracA ! Arrays for total area, calculated runoff, and fraction
  integer,allocatable,dimension(:) :: nsub           ! Array for number of sub-areas per catchment

  integer :: i,j,sx,sy                   ! Loop variables and coordinates for sub-areas

  ! Allocate memory for arrays
  allocate(nsub(nc),subarea(nmax,nc),subx(nmax,nc),suby(nmax,nc),tot(nc))

  ! Read sub-area data from text files
  open(77,file="output/Pfaf_nsub_M09.txt"); read(77,*)nsub
  open(77,file="output/Pfaf_asub_M09.txt"); read(77,*)subarea
  open(77,file="output/Pfaf_xsub_M09.txt"); read(77,*)subx
  open(77,file="output/Pfaf_ysub_M09.txt"); read(77,*)suby
  open(77,file="output/Pfaf_area.txt"); read(77,*)tot

  ! Allocate memory for fraction array
  allocate(frac(nmax,nc))

  ! Compute fraction of each sub-area relative to the total catchment area
  do i=1,nc
    frac(:,i)=subarea(:,i)/tot(i)
  enddo

  ! Allocate memory for runoff and fraction arrays
  allocate(runfC(nc),fracA(nc))
  runfC=0.                     ! Initialize runoff array to zero
  fracA=0.                     ! Initialize fraction array to zero

  !$OMP PARALLEL default(shared) private(i,j,sx,sy)   ! Start OpenMP parallel region
  !$OMP DO
  ! Loop over all catchments and sub-areas
  do i=1,nc
    do j=1,nsub(i)
      sy=suby(j,i)  ! Get y-coordinate of the sub-area
      sx=subx(j,i)  ! Get x-coordinate of the sub-area
      ! Check for valid fraction and runoff values
      if(frac(j,i)>0..and.runoff(sx,sy)<1.e14.and.runoff(sx,sy)>=0.)then
        runfC(i)=runfC(i)+frac(j,i)*runoff(sx,sy)     ! Accumulate runoff for the catchment
        fracA(i)=fracA(i)+frac(j,i)                   ! Accumulate fraction
      endif
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL    ! End OpenMP parallel region

  ! Convert to kg/s by multiplying by area (in m²) and dividing by time (in seconds)
  Qrunf=runfC*(tot*1.e6)/86400.  

  ! Deallocate arrays to free memory
  deallocate(subarea,subx,suby,tot,frac,&
             runfC,fracA,nsub)

end function M09_to_cat
!------------------------------------------------------------------------------

end program main
