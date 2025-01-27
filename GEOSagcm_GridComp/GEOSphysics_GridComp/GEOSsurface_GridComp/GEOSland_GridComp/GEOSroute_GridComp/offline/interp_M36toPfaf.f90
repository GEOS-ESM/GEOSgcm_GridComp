module interp

use omp_lib          ! Use OpenMP library for parallel processing
use rwncfile         ! Use custom module for reading NetCDF files
implicit none

private
public :: M36_to_cat   ! Make the M36_to_cat function public
public :: M09_to_cat   ! Make the M09_to_cat function public

contains

!------------------------------------------------------------------------------
! This function maps runoff data from M36 resolution to catchments (cat)
function M36_to_cat(runoff,nlon,nlat,ncat,inputdir) result(Qrunf)

  integer,intent(in) :: nlon,nlat,ncat    ! Input: number of longitude, latitude, and catchments
  real*8,intent(in) :: runoff(nlon,nlat)  ! Input: runoff array of size (nlon, nlat)
  character(len=500),intent(in) :: inputdir  ! Input: directory path for input files
  real*8  :: Qrunf(ncat)                  ! Output: runoff mapped to catchments

  real*8,parameter :: small=1.D-12        ! Small value to avoid division by zero

  integer,parameter :: nmax=150           ! Maximum number of sub-areas per catchment
  integer,parameter :: nc=291284          ! Total number of catchments

  real*8,allocatable,dimension(:,:) :: subarea,frac  ! Arrays for sub-area and fractions
  integer,allocatable,dimension(:,:) :: subx,suby    ! Arrays for x and y coordinates of sub-areas
  real*8,allocatable,dimension(:) :: tot,runfC,fracA ! Arrays for total area, calculated runoff, and fraction
  integer,allocatable,dimension(:) :: nsub           ! Array for number of sub-areas per catchment

  integer :: i,j,sx,sy                   ! Loop variables and coordinates for sub-areas

  ! Allocate memory for arrays
  allocate(nsub(nc),subarea(nmax,nc),subx(nmax,nc),suby(nmax,nc),tot(nc))

  ! Read sub-area data from text files
  open(77,file=trim(inputdir)//"/Pfaf_nsub_M36.txt"); read(77,*)nsub
  open(77,file=trim(inputdir)//"/Pfaf_asub_M36.txt"); read(77,*)subarea
  open(77,file=trim(inputdir)//"/Pfaf_xsub_M36.txt"); read(77,*)subx
  open(77,file=trim(inputdir)//"/Pfaf_ysub_M36.txt"); read(77,*)suby
  open(77,file=trim(inputdir)//"/Pfaf_area.txt"); read(77,*)tot

  ! Allocate memory for fraction array
  allocate(frac(nmax,nc))

  ! Compute fraction of each sub-area relative to the total catchment area
  do i=1,nc
    frac(:,i)=subarea(:,i)/tot(i)
  enddo

  ! Allocate memory for runoff and fraction arrays
  allocate(runfC(nc),fracA(nc))
  runfC=0.D0                     ! Initialize runoff array to zero
  fracA=0.D0                     ! Initialize fraction array to zero

  !$OMP PARALLEL default(shared) private(i,j,sx,sy)   ! Start OpenMP parallel region
  !$OMP DO
  ! Loop over all catchments and sub-areas
  do i=1,nc
    if(nsub(i)>=1)then
      do j=1,nsub(i)
        sy=suby(j,i)  ! Get y-coordinate of the sub-area
        sx=subx(j,i)  ! Get x-coordinate of the sub-area
        ! Check for valid fraction and runoff values
        if(frac(j,i)>0.D0.and.runoff(sx,sy)<1.D14)then
          runfC(i)=runfC(i)+frac(j,i)*runoff(sx,sy)     ! Accumulate runoff for the catchment
          fracA(i)=fracA(i)+frac(j,i)                   ! Accumulate fraction
        endif
      enddo
    endif
  enddo
  !$OMP END DO
  !$OMP END PARALLEL    ! End OpenMP parallel region

  ! Convert to kg/s by multiplying by area (in m²) and dividing by time (in seconds)
  Qrunf=runfC*(tot*1.D6)/86400.D0  

  ! Deallocate arrays to free memory
  deallocate(subarea,subx,suby,tot,frac,&  
             runfC,fracA,nsub)

end function M36_to_cat
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! This function maps runoff data from M09 resolution to catchments (cat)
function M09_to_cat(runoff,nlon,nlat,ncat,inputdir) result(Qrunf)

  integer,intent(in) :: nlon,nlat,ncat    ! Input: number of longitude, latitude, and catchments
  real*8,intent(in) :: runoff(nlon,nlat)  ! Input: runoff array of size (nlon, nlat)
  character(len=500),intent(in) :: inputdir  ! Input: directory path for input files
  real*8  :: Qrunf(ncat)                  ! Output: runoff mapped to catchments

  real*8,parameter :: small=1.D-12        ! Small value to avoid division by zero

  integer,parameter :: nmax=458           ! Maximum number of sub-areas per catchment
  integer,parameter :: nc=291284          ! Total number of catchments

  real*8,allocatable,dimension(:,:) :: subarea,frac  ! Arrays for sub-area and fractions
  integer,allocatable,dimension(:,:) :: subx,suby    ! Arrays for x and y coordinates of sub-areas
  real*8,allocatable,dimension(:) :: tot,runfC,fracA ! Arrays for total area, calculated runoff, and fraction
  integer,allocatable,dimension(:) :: nsub           ! Array for number of sub-areas per catchment

  integer :: i,j,sx,sy                   ! Loop variables and coordinates for sub-areas

  ! Allocate memory for arrays
  allocate(nsub(nc),subarea(nmax,nc),subx(nmax,nc),suby(nmax,nc),tot(nc))

  ! Read sub-area data from text files
  open(77,file=trim(inputdir)//"/Pfaf_nsub_M09.txt"); read(77,*)nsub
  open(77,file=trim(inputdir)//"/Pfaf_asub_M09.txt"); read(77,*)subarea
  open(77,file=trim(inputdir)//"/Pfaf_xsub_M09.txt"); read(77,*)subx
  open(77,file=trim(inputdir)//"/Pfaf_ysub_M09.txt"); read(77,*)suby
  open(77,file=trim(inputdir)//"/Pfaf_area.txt"); read(77,*)tot

  ! Allocate memory for fraction array
  allocate(frac(nmax,nc))

  ! Compute fraction of each sub-area relative to the total catchment area
  do i=1,nc
    frac(:,i)=subarea(:,i)/tot(i)
  enddo

  ! Allocate memory for runoff and fraction arrays
  allocate(runfC(nc),fracA(nc))
  runfC=0.D0                     ! Initialize runoff array to zero
  fracA=0.D0                     ! Initialize fraction array to zero

  !$OMP PARALLEL default(shared) private(i,j,sx,sy)   ! Start OpenMP parallel region
  !$OMP DO
  ! Loop over all catchments and sub-areas
  do i=1,nc
    do j=1,nsub(i)
      sy=suby(j,i)  ! Get y-coordinate of the sub-area
      sx=subx(j,i)  ! Get x-coordinate of the sub-area
      ! Check for valid fraction and runoff values
      if(frac(j,i)>0.D0.and.runoff(sx,sy)<1.D14.and.runoff(sx,sy)>=0.D0)then
        runfC(i)=runfC(i)+frac(j,i)*runoff(sx,sy)     ! Accumulate runoff for the catchment
        fracA(i)=fracA(i)+frac(j,i)                   ! Accumulate fraction
      endif
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL    ! End OpenMP parallel region

  ! Convert to kg/s by multiplying by area (in m²) and dividing by time (in seconds)
  Qrunf=runfC*(tot*1.D6)/86400.D0  

  ! Deallocate arrays to free memory
  deallocate(subarea,subx,suby,tot,frac,&
             runfC,fracA,nsub)

end function M09_to_cat
!------------------------------------------------------------------------------

end module interp