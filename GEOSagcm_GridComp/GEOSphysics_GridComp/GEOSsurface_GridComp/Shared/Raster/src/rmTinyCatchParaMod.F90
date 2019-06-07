
#include "Raster.h"

!
! A Collection subroutines needed by mkCatchParam.F90
!   Contact: Sarith Mahanama  sarith.p.mahanama@nasa.gov
!   Email  : sarith.p.mahanama@nasa.gov

module rmTinyCatchParaMod

  use date_time_util  
  use leap_year
  use MAPL_ConstantsMod
  

  implicit none
  logical, parameter :: error_file=.true.
  integer, parameter :: n_SoilClasses = 253
  real, parameter :: gnu=1.0, zks = 2.0
  integer, parameter :: i_raster = 8640, j_raster=4320
  integer, parameter :: ncat_gswp2 = 15238
  REAL, PARAMETER :: undef = 1.e+20
  integer, parameter :: arr_len = 1734915,ip1 =0 
  real, parameter :: dx_gswp2 =1.,dy_gswp2=1. 
  integer, parameter :: MAX_NOF_GRID = ncat_gswp2
  integer, PARAMETER :: nbdep=150, NAR=1000,nwt=81,nrz=41
  real, parameter :: slice=0.1, lim =5.,grzdep =1.1
  logical, parameter ::  bug =.false.
  include 'netcdf.inc'	
  logical :: preserve_soiltype = .false.
  character*100 :: c_data = 'data/CATCH/'

  private

  public remove_tiny_tiles,modis_alb_on_tiles,modis_scale_para  
  public make_75,catchment_def,soil_para_high
  public create_soil_types_files,compute_mosaic_veg_types
  public cti_stat_file, create_model_para_woesten
  public create_model_para, modis_lai,regridraster,regridrasterreal
  public i_raster, j_raster,regridraster1,regridraster2,n_SoilClasses,gnu,zks
  public mineral_perc, process_gswp2_veg,center_pix, soil_class
  public tgen, sat_param,REFORMAT_VEGFILES,base_param,ts_param
  public :: Get_MidTime, Time_Interp_Fac, compute_stats, c_data	
  public :: ascat_r0, jpl_canoph

type :: mineral_perc
   real :: clay_perc
   real :: silt_perc
   real :: sand_perc
end type mineral_perc

contains
! _____________________________________________________________________________________________
!

  SUBROUTINE Get_MidTime (      &
                          yr1,mn1,dy1,yr2,mn2,dy2, &
                          MIDT)
  
    implicit none 
    real, intent (in)                  :: yr1,mn1,dy1,yr2,mn2,dy2
    type(date_time_type), intent(out ) :: MIDT
    type(date_time_type)               :: TIME1, TIME2
    integer                            :: TIMEDIF

    TIME1%year  = NINT(yr1) + 2001
    TIME1%month = NINT(mn1)
    TIME1%day   = NINT(dy1)
    TIME1%hour  = 0
    TIME1%min   = 0
    TIME1%sec   = 0

    call get_dofyr_pentad(TIME1)
    MIDT = TIME1

    TIME2%year  = NINT(yr2) + 2001
    TIME2%month = NINT(mn2)
    TIME2%day   = NINT(dy2)
    TIME2%hour  = 23
    TIME2%min   = 59
    TIME2%sec   = 59
    call get_dofyr_pentad(TIME2)

    TIMEDIF = datetime2_minus_datetime1(TIME1,TIME2)
    TIMEDIF = TIMEDIF/2

    call augment_date_time(TIMEDIF, MIDT) 
!    print *,'MIDTIME'
!    print *,'TIME1:', time1
!    print *,'MIDT :', midt
!    print *,'TIME2:', time2

  END SUBROUTINE Get_MidTime

!
! ---------------------------------------------------------------------------------------------
!

SUBROUTINE Time_Interp_Fac (TIME0, TIME1, TIME2, FAC1, FAC2)

implicit none
!  PURPOSE:
!  ========
!
!    Compute interpolation factors, fac, to be used 
!    in the calculation of the instantaneous boundary 
!    conditions, ie:
!
!     q(i,j) = fac1*q1(i,j) + (1.-fac1)*q2(i,j)
!
!    where:
!     q(i,j)  => Boundary Data valid    at time0
!     q1(i,j) => Boundary Data centered at time1
!     q2(i,j) => Boundary Data centered at time2

!  INPUT:
!  ======
!    time0    : Time of current timestep
!    time1    : Time of boundary data 1 
!    time2    : Time of boundary data 2 

!  OUTPUT:
!  =======
!     fac1    : Interpolation factor for Boundary Data 1
!

    type(date_time_type),   intent(in ) :: TIME0, TIME1, TIME2
    real,             intent(out) :: FAC1
    real,             intent(out) :: FAC2

    real        :: TimeDif1
    real        :: TimeDif
!    print *,'Interpolation'
!    print *,'TIME1:', time1
!    print *,'TIME0:', time0
!    print *,'TIME2:', time2

    TimeDif1 = real(datetime2_minus_datetime1(TIME0,TIME2))
    TimeDif  = real(datetime2_minus_datetime1(TIME1,TIME2))
       
    FAC1 = TimeDif1/TimeDif

    FAC2 = 1.-FAC1

!    print *,fac1,fac2

END SUBROUTINE Time_Interp_Fac


! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
SUBROUTINE process_gswp2_veg (nc,nr,regrid,vname,gridnamer,merge)

implicit none
    integer, intent(in) :: nc,nr
    real, dimension (:), allocatable :: catforc,vecforc,catcount	  
    integer, allocatable, target, dimension (:,:) :: gswp2_mask
    integer, allocatable, dimension (:,:) :: tile_id
    integer, parameter :: MAX_NOF_GRID = 15238
    REAL, ALLOCATABLE :: mon_climate(:,:)	
    integer :: ierr, ncid,iret, maxcat
    integer :: i1,k1,n,iv,year,smon,imon,mon,i,j,status
    REAL, PARAMETER :: undef = 1.e+20,UNDEF_GSWP2=-9999.	
    integer :: k,ncatch
    integer :: yr,mn,yr1,mn1
    logical :: regrid
    integer, pointer :: Raster(:,:)
    character(*) :: vname,gridnamer
    character *100 :: fname
    integer, intent(in), optional :: merge
 
     open (10,file=trim(gridnamer)//'.rst',status='old',action='read',  &
          form='unformatted',convert='little_endian')

    allocate (gswp2_mask (1:i_raster,1:j_raster))
    allocate (tile_id    (1:nc,1:nr))         

     do j=1,nr
        read(10)tile_id(:,j)
     end do  
     close (10,status='keep')

    open (10,file=trim(c_data)//'gswp2_mask_2.5.rst',&
             form='unformatted',status='old',action='read',convert='little_endian')

    do j =1,j_raster
         read (10) gswp2_mask(:,j)
    end do
    close (10,status='keep')

    if(regrid) then
       allocate(raster(nc,nr),stat=STATUS); VERIFY_(STATUS)
    else
       raster => gswp2_mask
    end if

    if(regrid) then
       call RegridRaster(gswp2_mask,raster)
    endif 

    fname='clsm/catchment.def'
    open (10,file=fname,status='old',action='read',form='formatted')
	 read (10,*)maxcat
    close (10,status='keep')
  
    allocate(vecforc(1:MAX_NOF_GRID))
    allocate(catforc(maxcat))
    allocate(catcount(maxcat))
    allocate(mon_climate(1:maxcat,1:12))
    mon_climate(:,:)=0.
    catforc=0.

    mon_climate(:,:)=0.

    iret = NF_OPEN(trim(c_data)//trim(vname)//'_uk.nc',NF_NOWRITE, ncid)

    ASSERT_(iret==NF_NOERR)

    if (present (merge)) then
       open (31,file='clsm/lai.dat.gswp2',  &
            form='unformatted',status='unknown',convert='little_endian')
    else
       if(trim(vname) == 'LAI') open (31,file='clsm/lai.dat',  &
            form='unformatted',status='unknown',convert='little_endian')
       if(trim(vname) == 'grnFrac') open (31,file='clsm/green.dat',  &
            form='unformatted',status='unknown',convert='little_endian')
    endif
    
    do year=82,98
        
       smon=(year-82)*12
       imon=0
       do mon=smon+1,smon+12
          imon=imon+1
	  
	  iret = NF_GET_VARA_REAL(ncid, 6,(/1,mon/),(/MAX_NOF_GRID,1/),vecforc)
	  ASSERT_(iret==NF_NOERR)		          
	  catforc =1.e-20
	  catcount=0
          DO j =1,nr
	    DO I = 1,nc
             if((tile_id(i,j).gt.0).and.(tile_id(i,j).le.maxcat)) then
	     if ((Raster(i,j).ge.1).and.(Raster(i,j).le.MAX_NOF_GRID)) then
	        catforc(tile_id(i,j)) = catforc(tile_id(i,j)) +   &
		  vecforc(Raster(i,j))  
		catcount(tile_id(i,j)) = catcount(tile_id(i,j)) + 1. 
	     endif
	     endif
	    END DO
          END DO
	  
	  do i = 1, maxcat
	     if(catcount(i).gt.0.) catforc(i) = catforc (i) /catcount(i)
	  end do
          mon_climate(:,imon)=mon_climate(:,imon)+catforc(:)/17.

       END DO  
    END DO ! Year
    iret = NF_CLOSE(ncid)
    ASSERT_(iret==NF_NOERR)

    fname='clsm/catchment.def'

    open (10,file=fname,status='old',action='read',form='formatted')
    read (10,*) ncatch
    close(10,status='keep')

    do K=0,13
     yr = (k+11)/12
     mn = mod(k+11,12)+1
     yr1= (k+12)/12
     mn1= mod(k+12,12)+1
     write(31) float((/yr,mn,1,0,0,0,yr1,mn1,1,0,0,0,ncatch,1/))
     write(31) mon_climate(:,mod(k+11,12)+1)
    end do 
    close(31,status='keep')

    deallocate(catforc)
    deallocate(mon_climate) 
    deallocate(vecforc)
    deallocate(catcount)  
    deallocate(gswp2_mask)
    deallocate(tile_id)
    if(regrid) then
       deallocate(raster)
    endif 

END SUBROUTINE process_gswp2_veg


! ---------------------------------------------------------------------
!----------------------------------------------------------------------

SUBROUTINE modis_lai (nx,ny,regrid,gfile)

implicit none
type (date_time_type) :: before_time,after_time,end_time
character*300 :: fout,fname
character(*) :: gfile
integer :: i, n, k,j,ncatch
integer :: yr,mn,dy,yr1,mn1,dy1
real, dimension (:,:), target, allocatable :: lai_grid
real, dimension (:), allocatable :: lai,count
character*5 :: mmdd
integer, allocatable, dimension (:,:) :: tile_id
integer :: i_sib,j_sib
integer :: nx,ny,status
logical :: regrid
real, pointer :: Raster(:,:)	

allocate(tile_id(1:nx,1:ny))
i_sib = i_raster
j_sib = j_raster

fname=trim(gfile)//'.rst'
open (10,file=fname,status='old',action='read',  &
      form='unformatted',convert='little_endian')
          
 do j=1,ny
      read(10)tile_id(:,j)
 end do
          
 close (10,status='keep')

before_time%year   =2001
before_time%month  =1
before_time%day    =1
before_time%hour   =0            
before_time%min    =0            
before_time%sec    =0            
before_time%pentad =1            
before_time%dofyr  =1            

end_time%year      =2001
end_time%month     =12
end_time%day       =31
end_time%hour      =23            
end_time%min       =59            
end_time%sec       =59            
end_time%pentad    =73            
end_time%dofyr     =365 
           
after_time = before_time

do j =1,8
   call augment_date_time (86400, after_time)
end do

fname='clsm/catchment.def'

open (10,file=fname,status='old',action='read',form='formatted')
read (10,*) ncatch
close(10,status='keep')

allocate (lai_grid (1:i_sib,1:j_sib))
allocate (lai   (1:ncatch))
allocate (count (1:ncatch))

fname = 'data/CATCH/MODIS_8-DayClim/' &
//'MOD15A2.YYYY.12.27.global_2.5min.data'

!write (*,'(a120)')trim(fname)
!write (*,'(2(f2.0,2f3.0,3f2.0),f9.0,f2.0)') float((/0,12,27,0,0,0,1,1,1,0,0,0,ncatch,1/))
open (20,file=trim(fname),form='unformatted',convert='little_endian', &
      action='read',status='old')

do j =1,j_raster
    read (20) lai_grid (:,j)
end do
close(20,status='keep')
lai   = 0.
count = 0.

if(regrid) then
    allocate(raster(nx,ny),stat=STATUS); VERIFY_(STATUS)
else
    raster => lai_grid
end if

if(regrid) then
    call RegridRasterReal(lai_grid,raster)
endif


do j=1,ny
   do i=1,nx
      if((tile_id(i,j).gt.0).and.(tile_id(i,j).le.ncatch)) then
         if((raster (i,j).ge.0.).and.(raster (i,j).le.10.))   then
            lai  (tile_id(i,j)) = &
                 lai(tile_id(i,j)) + raster(i,j)
            count(tile_id(i,j)) = &
                 count(tile_id(i,j)) + 1. 
         endif
      endif
   end do
end do

DO n =1,ncatch
	if(count(n)/=0.) lai(n)=lai(n)/count(n)	
END DO

fout = 'clsm/lai.dat'
open (30,file=trim(fout),form='unformatted',convert='little_endian', &
      action='write',status='unknown')

write(30) float((/0,12,27,0,0,0,1,1,1,0,0,0,ncatch,1/))
write (30) lai

do while (datetime_le_refdatetime(before_time,end_time))     

   yr = before_time%year -2000 !(k+11)/12
   mn = before_time%month      !mod(k+11,12)+1
   dy = before_time%day
   yr1= after_time%year -2000    !(k+12)/12
   mn1= after_time%month         !mod(k+12,12)+1
   dy1= after_time%day
   write (mmdd,'(i2.2,a1,i2.2)'),mn,'.',dy

   fname ='data/CATCH/MODIS_8-DayClim/' &
	 //'MOD15A2.YYYY.'//mmdd//'.global_2.5min.data'

   open (20,file=trim(fname),form='unformatted',convert='little_endian', &
        action='read',status='old')
	do j =1,j_raster
	   read (20) lai_grid (:,j)
	end do
   close (20,status='keep')

   if(regrid) then
      call RegridRasterReal(lai_grid,raster)      
   else
       raster => lai_grid
   endif

   lai   = 0.
   count = 0.

   do j=1,ny
      do i=1,nx
	 if((tile_id(i,j).gt.0).and.(tile_id(i,j).le.ncatch)) then
		if((raster (i,j).ge.0.).and.(raster (i,j).le.10.))   then
		     lai  (tile_id(i,j)) = &
		          lai(tile_id(i,j)) + raster(i,j)
		     count(tile_id(i,j)) = &
		          count(tile_id(i,j)) + 1. 
		endif
	endif
      end do
   end do

DO n =1,ncatch
	if(count(n)/=0.) lai(n)=lai(n)/count(n)	
END DO
   if(mmdd.eq.'12.27') dy1 = 1 
!   write(*,'(2(f2.0,2f3.0,3f2.0),f9.0,f2.0)') float((/yr,mn,dy,0,0,0,yr1,mn1,dy1,0,0,0,ncatch,1/))
   write(30) float((/yr,mn,dy,0,0,0,yr1,mn1,dy1,0,0,0,ncatch,1/))
   write (30) lai

   do j =1,8
      call augment_date_time (86400, before_time)
      call augment_date_time (86400, after_time )
   end do
end do

fname = 'data/CATCH/MODIS_8-DayClim/' &
//'MOD15A2.YYYY.01.01.global_2.5min.data'

!write (*,'(a120)')trim(fname)
open (20,file=trim(fname),form='unformatted',convert='little_endian', &
      action='read',status='old')

do j =1,j_raster
    read (20) lai_grid (:,j)
end do
close(20,status='keep')

if(regrid) then
     call RegridRasterReal(lai_grid,raster)
 else
     raster => lai_grid
endif

lai   = 0.
count = 0.

do j=1,ny
   do i=1,ny
      if((tile_id(i,j).gt.0).and.(tile_id(i,j).le.ncatch)) then
         if((raster (i,j).ge.0.).and.(raster (i,j).le.10.))   then
            lai  (tile_id(i,j)) = &
                 lai(tile_id(i,j)) + raster(i,j)
            count(tile_id(i,j)) = &
                 count(tile_id(i,j)) + 1. 
         endif
      endif
   end do
end do

DO n =1,ncatch
	if(count(n)/=0.) lai(n)=lai(n)/count(n)	
END DO

! write(*,'(2(f2.0,2f3.0,3f2.0),f9.0,f2.0)') float((/2,1,1,0,0,0,2,1,9,0,0,0,ncatch,1/))
write(30) float((/2,1,1,0,0,0,2,1,9,0,0,0,ncatch,1/))
write(30) lai

    if(regrid) then
       deallocate(raster)
    endif 
close(30,status='keep')

END SUBROUTINE modis_lai

!----------------------------------------------------------------------  

  SUBROUTINE soil_para_high (nx,ny,regrid,gfile,F25Tag)

    implicit none
      real, dimension(12) :: lbee,lpsis,lporo,lcond,lwpwet, &
           atau2,btau2,atau5,btau5
      REAL, ALLOCATABLE :: soildepth (:)
      INTEGER :: soil_class_top,soil_class_com,soil_gswp,swit
      REAL :: BEE, PSIS, POROS,COND,WPWET
      integer :: n,maxcat,count,k1,i1,i,j
      character*100 :: path,fname,fout,metpath
      character(*) :: gfile
      character*10 :: dline
      CHARACTER*20 :: version,resoln,continent
      integer :: iret,ncid
      real, allocatable, target, dimension (:,:) :: SOIL_HIGH
      integer, allocatable, dimension (:,:) :: tile_id
      REAL, ALLOCATABLE :: count_soil(:)
      integer :: tindex, pfafindex,i_sib,j_sib
      integer :: nx,ny,status
      real, allocatable, dimension(:) :: soildepth_gswp2 
      integer, allocatable, dimension (:) :: land_gswp2
      logical :: regrid
      real, pointer :: Raster(:,:)
      logical, intent (in), optional :: F25Tag

      data lbee /3.30, 3.80, 4.34, 5.25, 3.63, 5.96, 7.32,       &
                 8.41, 8.34, 9.70, 10.78, 12.93/
      data lpsis /-0.05, -0.07, -0.16, -0.65, -0.84, -0.24,      &
                 -0.12, -0.63, -0.28, -0.12, -0.58, -0.27/
      data lporo /0.373, 0.386, 0.419, 0.476, 0.471, 0.437,      &
                  0.412, 0.478, 0.447, 0.415, 0.478, 0.450/
      data lcond /2.45e-05, 1.75e-05, 8.35e-06, 2.36e-06,        &
                  1.1e-06, 4.66e-06, 6.31e-06, 1.44e-06,         &
                  2.72e-06, 4.25e-06, 1.02e-06, 1.33e-06/
      data lwpwet /0.033,0.051,0.086,0.169,0.045,0.148,0.156,    &
                  0.249,0.211,0.199,0.286,0.276/

      data atau2/0.0030065,0.0276075,0.0200614,0.0165152,   &
                0.0165152,0.0168748,0.0308809,0.0329365,    &
                0.0437085,0.0466403,0.0956670,0.1257360/

      data btau2/0.0307900,0.0196558,0.0299702,0.0443406,   &
                0.0443406,0.0359961,0.0234851,0.0370919,    &
                0.0312746,0.0249973,0.0222786,0.0193874/

      data atau5/0.0067424,0.0766189,0.0540989,0.0439714,   &
                0.0439714,0.0457011,0.0589881,0.0885157,    &
                0.1175960,0.0692305,0.1348880,0.1535540/

      data btau5/0.0569718,0.0492634,0.0678898,0.0786387,   &
                0.0786387,0.0737872,0.0713841,0.0742609,    &
                0.0693533,0.0745496,0.0732726,0.0718882/

      i_sib = i_raster
      j_sib = j_raster
      fname='clsm/catchment.def'

      open (10,file=fname,status='old',action='read',form='formatted')
      read(10,*) maxcat
      close (10,status='keep')

	  allocate(soildepth(maxcat))
          allocate(soil_high(1:i_raster,1:j_raster))  
          allocate(count_soil(1:maxcat))  
          allocate(tile_id(1:nx,1:ny))

          soil_high =-9999.
          fname=trim(gfile)//'.rst'
          
          open (10,file=fname,status='old',action='read',  &
               form='unformatted',convert='little_endian')
          
          do j=1,ny
             read(10)tile_id(:,j)
          end do
          
          close (10,status='keep')
          
          if (present(F25Tag)) then 

             iret = NF_OPEN('data/CATCH/SoilDepth.nc',NF_NOWRITE, ncid)
             ASSERT_(iret==NF_NOERR)
             allocate (soildepth_gswp2(1: ncat_gswp2))
             allocate (land_gswp2     (1: ncat_gswp2)) 
             iret = NF_GET_VARA_INT (ncid, 3,(/1/),(/ncat_gswp2/),land_gswp2)
	     ASSERT_(iret==NF_NOERR)	
             iret = NF_GET_VARA_REAL(ncid, 4,(/1/),(/ncat_gswp2/),soildepth_gswp2)
	     ASSERT_(iret==NF_NOERR)		          
             iret = NF_CLOSE(ncid)
             ASSERT_(iret==NF_NOERR)

             k1 = i_raster/360

             do n = 1,ncat_gswp2

                j = (land_gswp2(n)-1)/360  + 1
                i = land_gswp2(n) - (j - 1)*360
                j = 181 - j
                soil_high((i-1)*k1+1:i*k1,(j-1)*k1+1:j*k1) = soildepth_gswp2(n)
                
             end do
              deallocate (soildepth_gswp2,land_gswp2)
          else
                          
             open (10,file='data/CATCH/soil_depth_2.5.rst',&
                  form='unformatted',status='old',action='read',convert='little_endian')
             
             do j =1,j_raster
                read (10) soil_high(:,j)
             end do
             close (10,status='keep')

          endif
             
	  if(regrid) then
	      allocate(raster(nx,ny),stat=STATUS); VERIFY_(STATUS)
          else
              raster => soil_high
          end if

          if(regrid) then
               call RegridRasterReal(soil_high,raster)
          endif	  

          soildepth =0.
          count_soil = 0.

          do j=1,ny
             do i=1,nx
                if((tile_id(i,j).gt.0).and.(tile_id(i,j).le.maxcat)) then
                   if(raster(i,j).eq.-9999.)   then
!		   write (*,*)'soil_high UNDEF',i,j,tile_id(i,j),raster(i,j) 
                      !                           stop
                   endif
                   if (raster(i,j).gt.0.) then
 
                      soildepth(tile_id(i,j)) = &
                           soildepth(tile_id(i,j)) + raster(i,j)
                      count_soil(tile_id(i,j)) = &
                           count_soil(tile_id(i,j)) + 1. 
                   endif
                endif
             end do
	  end do

          DO n =1,maxcat
		if(count_soil(n)/=0.) soildepth(n)=soildepth(n)/count_soil(n)	
                if (present(F25Tag)) then
                   soildepth(n) = max(soildepth(n),1.)
                else
                   soildepth(n) = max(soildepth(n),1.334)
                endif
          END DO

          soildepth = soildepth*1000. 

!     Openning files

      fname='clsm/soil_text.top'
          open (10,file=fname,status='old',action='read',form='formatted')
      fname='clsm/soil_text.com'
          open (11,file=fname,status='old',action='read',form='formatted')
      fout='clsm/soil_param.first'          
      open (21,file=fout,status='unknown',action='write',form='formatted')
      fout='clsm/tau_param.dat'          
      open (22,file=fout,status='unknown',action='write',form='formatted')

      swit =0
      DO n=1 , maxcat
         read (10,*) tindex,pfafindex, soil_class_top
         write (22,'(i8,i8,4f10.7)')tindex,pfafindex,atau2(soil_class_top), &
              btau2(soil_class_top),atau5(soil_class_top),btau5(soil_class_top)
              read (11,*) tindex,pfafindex, soil_class_com

        !if (soil_class_com.eq.4) then
        !   soil_gswp = 5
        !elseif (soil_class_com.eq.5) then
        !   soil_gswp = 6
        !elseif (soil_class_com.eq.6) then
        !   soil_gswp = 4
        !elseif (soil_class_com.eq.8) then
        !   soil_gswp = 9
        !elseif (soil_class_com.eq.9) then
        !   soil_gswp = 8
        !else
        !   soil_gswp = soil_class_com
        !endif
        
        soil_gswp = soil_class_com
        
        cond=lcond(soil_gswp)/exp(-1.*zks*gnu)
        wpwet=lwpwet(soil_gswp)/lporo(soil_gswp) 
        write (21,'(i8,i8,i4,i4,3f8.4,f12.8,f7.4,f10.3)')tindex,pfafindex,   &
        soil_class_top,soil_class_com,lBEE(soil_gswp), lPSIS(soil_gswp),          &
        lPORO(soil_gswp),COND,WPWET,soildepth(n)
         
      END DO
      close (10,status='delete')
      close (11,status='delete')
      close (21,status='keep')
      close (22,status='keep')
      deallocate (soildepth, soil_high)
    if(regrid) then
       deallocate(raster)
    endif 
  END SUBROUTINE soil_para_high
!
! ====================================================================
!
  SUBROUTINE remove_tiny_tiles (                         &
       dateline,poles,gout) 

    IMPLICIT NONE
    INTEGER :: ip,ip2,nc_gcm,nr_gcm,nc_ocean,nr_ocean,pick_val,k,nc,nr
    INTEGER :: typ,pfs,ig,jg,indx,indx_old,j_dum,ierr,n,count,count_remain,i_dum
    REAL :: lat,lon,mx_frac,da,tarea
    REAL_ :: fr_gcm,fr_ocean,fr_cat,lats,dx,dy,d2r
    INTEGER :: im,jm,i,j,jk,ik,jx
    INTEGER :: l,imn,imx,jmn,jmx
    CHARACTER*30 :: version
    CHARACTER*128 :: fname,gname,gout,gpath
    character*300 :: string1, string2
    integer(kind=4), allocatable, dimension(:,:) :: grid
    integer(kind=4), allocatable, dimension(:,:) :: grida
    REAL (kind=8), PARAMETER :: threshold=0.01,RADIUS=MAPL_RADIUS,pi= RASTER_PI !3.1415926535898, 6371000.	
    real(kind=8), allocatable, dimension(:) :: tile_frac,total_area,pfaf,tile_area(:),lon_c(:),lat_c(:),int_c(:)
    character*2 :: dateline,poles
    integer, allocatable, dimension(:) :: rev_indx
    real, allocatable, dimension(:,:):: tile_frac_2d
    integer(kind=4),allocatable :: GRIDX(:,:)
    !
    nc=i_raster 
    nr=j_raster
    dx  = 360._8/nc
    dy  = 180._8/nr
    d2r = PI/180._8

    print *,'Revised tile space..:','clsm/'//trim(gout)//'-Pfaf.notiny'
    
    gname='til/'//trim(gout)//'-Pfafstetter'
    fname=  trim(gname)//'.til'

    print *,'Any tile whose geographic area is <',threshold
    print *,'of the AGCM grid box will be dissolved and'
    print *,'the largest geographic neighbor will annex it!'
    print *,'----------------------------------------------'

    open (10,file=trim(fname),status='old',action='read',form='formatted')
    read (10,*)ip
    read (10,*)j_dum
    read (10,'(a)')version
    read (10,*)nc_gcm
    read (10,*)nr_gcm
    read (10,'(a)')version
    read (10,*)nc_ocean
    read (10,*)nr_ocean
 
    count=0

    do n = 1,ip

      read(10,'(I10,3E20.12,9(2I10,E20.12,I10))',IOSTAT=ierr)     &    
            typ,tarea,lon,lat,ig,jg,fr_gcm,indx_old,pfs,j_dum,fr_cat,j_dum
 
       if (typ == 100) ip2 = n
       if ((typ == 100).and.(fr_gcm < threshold)) count =count +1 
       if(ierr /= 0)write (*,*)'Problem reading',fname

    end do
    
    write (*,*)'# of small catchments to be removed  ', count
    if (count < ip2/100) then 

    print *,'Too few tiny tiles, thus exiting .............'
    print *,'CLSM parameters will be generated for ........'
    print *,trim(gname)
    string1 ='til/'//trim(gout)//'-Pfafstetter.til'//' '//&
    'clsm/'//trim(gout)//'-Pfaf.notiny.til'
    call system ('cp '//trim(string1))
    string1 ='rst/'//trim(gout)//'-Pfafstetter.rst'//' '//&
    'clsm/'//trim(gout)//'-Pfaf.notiny.rst'
    call system ('cp '//trim(string1))
    print *,'and, copied those those files to clsm/.'

    stop
    endif

    !
    IM = nc/nc_gcm      ! i-Pixels in GCM box
    JM = nr/(nr_gcm-1)  ! j-Pixels in interior GCM box. Pole boxes have half as many.
    if (index(poles,'PE')/=0) JM = nr/(nr_gcm)  ! pole edge case
    allocate(GRID (nc,jm)) ! Enough space for all pixels in non-pole GCM latitude
    allocate(GRIDA (nc,jm)) ! Enough space for all pixels in non-pole GCM latitude
    allocate(tile_frac(ip))
    grid=0
    grida=0
    
    fname='rst/'//trim(gout)//'-Pfafstetter.rst'
    open (10,file=trim(fname),status='old',action='read',form='unformatted',convert='little_endian')
    fname='clsm/'//trim(gout)//'-catchs_nosmall_rst'
    open (11,file=trim(fname),status='unknown',action='write',form='unformatted',convert='little_endian')
    
    do j=1,nr_gcm !  loop over GCM latitudes
       if(j==1.or.j==nr_gcm) then
          jx=jm/2  !  pole latitudes are half as large 
          if (index(poles,'PE')/=0) jx=jm 
       else
          jx=jm
       endif
       
       allocate(tile_frac_2d(im,jx))
       do jk=1,jx  ! Read raster data for one row of atmos grid points
          if (index(dateline,'DE')/=0) then
             read (10) grid(:,jk)
          else
             read (10) grid(im/2+1:,jk),grid(1:im/2,jk)
          endif
          if(maxval(grid(:,jk)).gt.ip) print *,'MAX EXCEED',maxval(grid(:,jk)),ip,jk
       enddo

       grida=grid
       
       do i=1,nc_gcm   !  loop over GCM longitudes
          tile_frac = 0
          tile_frac_2d = 0
          
          allocate(gridx(im,jx))
          gridx(1:im,1:jx)= grid(1+(I-1)*IM:I*IM,1:JX)
          !
          ! We don't touch ocean, ice and lakes pixels
          do jk=1,jx
             do ik = 1,im
                if(gridx(ik,jk) > ip2)gridx(ik,jk)=0
             end do
          end do
          
          !         We don't have to process 100% ocean, lake or ice pixels
          if (sum(gridx) /= 0) then
             !
             do jk=1,jx
                !             do ik=1+(I-1)*IM,I*IM
                do ik = 1,im 
                   if(gridx(ik,jk) /= 0) then
                      tile_frac(gridx(ik,jk)) = tile_frac(gridx(ik,jk)) + &
                           1./FLOAT(im*jx) 
                   endif
                end do
             end do
             !
             do n= 1,ip2
                if (tile_frac(n) > threshold) then
                   do jk=1,jx
                      do ik = 1,im
                         if(gridx(ik,jk) == n)then   
                            tile_frac_2d(ik,jk) = tile_frac(n) 
                         endif
                      end do
                   end do
                end if
             end do
             !
             if(sum(tile_frac_2d)>0.)then
                do n= 1,ip2
                   if ((tile_frac(n) > 0.).and.(tile_frac(n) <threshold ))then
                      do jk=1,jx
                         do ik = 1,im
                            
                            if (gridx(ik,jk) == n) then
                               mx_frac=0.
                               pick_val=0
                               l=1
                               do 
                                  imx=ik+l
                                  imn=ik-l
                                  jmn=jk-l
                                  jmx=jk+l
                                  imn=MAX(imn,1)
                                  jmn=MAX(jmn,1)
                                  imx=MIN(imx,im)
                                  jmx=MIN(jmx,jx) 
                                  !
                                  do k=imn,imx
                                     if(tile_frac_2d(k,jmn) > mx_frac) then
                                        mx_frac = tile_frac_2d(k,jmn)
                                        pick_val = gridx(k,jmn)
                                     endif
                                  end do
                                  !
                                  do k=imn,imx
                                     if(tile_frac_2d(k,jmx) > mx_frac) then
                                        mx_frac = tile_frac_2d(k,jmx)
                                        pick_val = gridx(k,jmx)
                                     endif
                                  end do
                                  !
                                  do k=jmn,jmx
                                     if(tile_frac_2d(imn,k) > mx_frac) then
                                        mx_frac = tile_frac_2d(imn,k)
                                        pick_val = gridx(imn,k)
                                     endif
                                  end do
                                  !
                                  do k=jmn,jmx
                                     if(tile_frac_2d(imx,k) > mx_frac) then
                                        mx_frac = tile_frac_2d(imx,k)
                                        pick_val = gridx(imx,k)
                                     endif
                                  end do
                                  !
                                  if(pick_val >0) grida ((I-1)*IM+ik ,jk) = &
                                       pick_val
                                  if(pick_val >0) exit
                                  l =l+1
                               end do
                            endif
                         end do
                      end do
                   endif
                end do
             endif
          endif !         We don't have to process 100% ocean, lake or ice pixels  
          
          deallocate(gridx)
       end do           !  loop over GCM longitudes
       deallocate(tile_frac_2d)

 !      print *,maxval(grid),minval(grid)
 !      print *,maxval(grida(:,1:jx)),minval(grida(:,1:jx)),jx       

       do jk=1,jx  ! Read raster data for one row of atmos grid points
          write (11) grida(:,jk)
       enddo
    end do              !  loop over GCM latitudes
    !
    close (10,status='keep')
    close (11,status='keep')

    open (11,file=trim(fname),status='unknown',action='read',form='unformatted',convert='little_endian')
    tile_frac=0.

    do j=1,nr_gcm !  loop over GCM latitudes
       grid=0
       
       if(j==1.or.j==nr_gcm) then
          jx=jm/2  !  pole latitudes are half as large 
          if (index(poles,'PE')/=0) jx=jm 
       else
          jx=jm
       endif
       
       do jk=1,jx  ! Read raster data for one row of atmos grid points
          read (11) grid(:,jk)
       enddo

       do i=1,nc_gcm   !  loop over GCM longitudes
          do jk=1,jx
             do ik=1+(I-1)*IM,I*IM

                tile_frac(grid(ik,jk)) = tile_frac(grid(ik,jk)) + 1./FLOAT(im*jx)

             end do
          enddo
       end do
    enddo
    !
    close (11,status='keep')
    !     
    count=0
    count_remain=0
    allocate(rev_indx(ip))
    rev_indx=0
    do n=1,ip
       if(tile_frac(n) > 0.) then
          count = count + 1
          rev_indx(n) = count
          if((n > ip1).and.(n <= ip2))then
             if ( tile_frac(n) < threshold) count_remain =count_remain + 1
          endif
       end if
    end do
    !
    write(*,*)'# of small catchments after merging',count_remain
    write(*,*)'# of tiles in the before removing tiny tiles :',ip
    write(*,*)'# of tiles in the after removing tiny tiles  :',count
    open (11,file=trim(fname),status='unknown',action='read',form='unformatted',convert='little_endian')
    fname='clsm/'//trim(gout)//'-Pfaf.notiny.rst'
    open (12,file=trim(fname),status='unknown',action='write',form='unformatted',convert='little_endian')
    
    deallocate (grid,grida)
    allocate(GRID (nc,1)) ! Enough space for all pixels in non-pole GCM latitude
    allocate(GRIDA (nc,1)) ! Enough space for all pixels in non-pole GCM latitude
    
    grid=0
    grida=0
    
    allocate(total_area(ip))
    allocate(tile_area(ip))
    allocate(lon_c(ip))
    allocate(lat_c(ip))
    allocate(int_c(ip))
    
    lon_c=0.
    lat_c=0.
    int_c=0.
    tile_area=0.
    total_area=0.
    da = radius*radius*pi*pi/24./24./180./180./1000000.
    
    do jk =1,nr
       lats = -90._8 + (jk - 0.5_8)*dy
       read (11) grid(:,1)
       do ik = 1,nc
          grida(ik,1)=rev_indx(grid(ik,1))
          total_area(rev_indx(grid(ik,1)))=total_area(rev_indx(grid(ik,1))) +1.
          tile_area(rev_indx(grid(ik,1)))=tile_area(rev_indx(grid(ik,1))) + &  
               (sin(d2r*(lats+0.5*dy)) - &
               sin(d2r*(lats-0.5*dy))   )*(dx*d2r)

!               da*cos((-90.+float(jk)/24. -1./48.)*pi/180.)
               
          lat_c(rev_indx(grid(ik,1)))=lat_c(rev_indx(grid(ik,1))) +&
               (-90.+float(jk)/24. -1./48.)
          
          if (index(dateline,'DE')/=0) then
             lon_c(rev_indx(grid(ik,1)))=lon_c(rev_indx(grid(ik,1))) + &
                  (-180.+float(ik)/24. -1./48.)
          else
             if(ik.le.im/2)then

                lon_c(rev_indx(grid(ik,1)))=lon_c(rev_indx(grid(ik,1))) +&
                     (-360.-180.+float(nc-im/2+ik)/24. -1./48.)
             else
                lon_c(rev_indx(grid(ik,1)))=lon_c(rev_indx(grid(ik,1))) + &
                     (-180.+float(ik-im/2)/24. -1./48.)
             endif
          endif
          int_c(rev_indx(grid(ik,1)))=int_c(rev_indx(grid(ik,1))) + 1.
       end do
       if (index(dateline,'DE')/=0) then
          write (12) grida(:,1)   
       else
          write (12) grida(im/2+1:,1),grida(1:im/2,1)    
       endif
    end do
    close (11,status='delete')
    close (12,status='keep')
    
    do n=1,ip
       if(rev_indx(n)>0)then
          lat_c(rev_indx(n))=lat_c(rev_indx(n))/int_c(rev_indx(n))
          lon_c(rev_indx(n))=lon_c(rev_indx(n))/int_c(rev_indx(n))
          if(lon_c(rev_indx(n)).lt.-180.)lon_c(rev_indx(n))=lon_c(rev_indx(n))+360.
       endif
    enddo
    !
    gname='til/'//trim(gout)//'-Pfafstetter'
    fname=  trim(gname)//'.til'
    !
    open (10,file=trim(fname),status='old',action='read',form='formatted')
    read (10,*)ip
    read (10,*)j_dum
    read (10,'(a)')version
    read (10,*)nc_gcm
    read (10,*)nr_gcm
    read (10,'(a)')version
    read (10,*)nc_ocean
    read (10,*)nr_ocean

    allocate(pfaf(36716))
    pfaf=0.

    do n = 1,ip

      read(10,'(I10,3E20.12,9(2I10,E20.12,I10))',IOSTAT=ierr)     &    
            typ,tarea,lon,lat,ig,jg,fr_gcm,j_dum,pfs,indx_old,fr_cat,j_dum

       if(n <= ip2) then
          if (rev_indx(n)>0)pfaf(indx_old)=pfaf(indx_old) +&
               total_area(rev_indx(n))
       endif
    end do
    close (10,status='keep')
    fname=  trim(gname)//'.til'    
    !
    open (10,file=trim(fname),status='old',action='read',form='formatted')
    fname='clsm/'//trim(gout)//'-Pfaf.notiny.til'
    open (20,file=trim(fname),status='unknown',action='write',form='formatted')
    read (10,*)ip
    write (20,*)COUNT, 8640,4320
    read (10,*)j_dum
    write (20,*)j_dum
    read (10,'(a)')version
    write(20,'(a)')version
    read (10,*)nc_gcm
    write (20,*)nc_gcm
    read (10,*)nr_gcm
    write (20,*)nr_gcm
    read (10,'(a)')version
    write (20,'(a)')version
    read (10,*)nc_ocean
    write (20,*)nc_ocean
    read (10,*)nr_ocean
    write (20,*)nr_ocean
    
    do n = 1,ip
       
     read(10,'(I10,3E20.12,9(2I10,E20.12,I10))',IOSTAT=ierr)     &    
            typ,tarea,lon,lat,ig,jg,fr_gcm,indx_old,pfs,i_dum,fr_cat,j_dum

       if(n <= ip2)then
          if (rev_indx(n)>0) then
             
             write(20,'(I10,3E20.12,9(2I10,E20.12,I10))',IOSTAT=ierr) &
                  typ,tile_area(rev_indx(n)),lon_c(rev_indx(n)),lat_c(rev_indx(n)),ig,jg,   &
                  tile_frac(n),indx_old,pfs,i_dum,total_area(rev_indx(n))/pfaf(i_dum),rev_indx(n)

          endif
       else
          if (rev_indx(n)>0)write(20,'(I10,3E20.12,9(2I10,E20.12,I10))',IOSTAT=ierr) &
            typ,tile_area(rev_indx(n)),lon_c(rev_indx(n)),lat_c(rev_indx(n)),ig,jg,   &
            tile_frac(n),indx_old,pfs,i_dum,fr_cat,rev_indx(n)

       endif
    end do

    write(*,*)'Surface Area of the Earth',sum(tile_area)
    write(*,*)'Land area of the Earth',sum(tile_area(rev_indx(ip1+1):rev_indx(ip2)))
    close (10,status='keep')
    close (20,status='keep')
    
  END SUBROUTINE remove_tiny_tiles


! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------

  SUBROUTINE modis_alb_on_tiles (nx,ny,ease_grid,regrid,gfilet,gfiler)
    
    implicit none
    CHARACTER*20 :: version,resoln,continent
    character*100 :: path,fname,fout,metpath
    character (*) :: gfilet,gfiler
    character*10 :: dline
    integer :: n,ip,maxcat,count,k1,i1,i
    integer :: nc_gcm,nr_gcm,nc_ocean,nr_ocean
    REAL :: lat,lon,fr_gcm,fr_cat,tarea
    INTEGER :: typ,pfs,ig,jg,j_dum,ierr,indx_dum,indr1,indr2,indr3 ,ip2
    INTEGER :: laiid,year,mon,smon,imon,iret  
    integer,allocatable :: tile_id(:,:)
    integer :: ialbt,ialbs,yy,j,month
    character*2 :: bw
    character*5 :: cyy
    character*300 :: albtype, albspec
    real, allocatable, target, dimension (:,:) :: alb_in
    real, allocatable, dimension (:) :: alb_count,alb_out
    character*300 :: ifile,ofile
    integer :: nx,ny,status
    logical :: regrid, ease_grid
    real,pointer :: raster (:,:)

    fname=trim(gfilet)//'.til'

    open (10,file=fname,status='old',action='read',form='formatted')
    read (10,*)ip
    read (10,*)j_dum
    read (10,'(a)')version
    read (10,*)nc_gcm
    read (10,*)nr_gcm
    read (10,'(a)')version
    read (10,*)nc_ocean
    read (10,*)nr_ocean
    
    do n = 1,ip
      if (ease_grid) then     
	 read(10,*,IOSTAT=ierr) typ,pfs,lon,lat,ig,jg,fr_gcm
      else
      read(10,'(I10,3E20.12,9(2I10,E20.12,I10))',IOSTAT=ierr)     &    
            typ,tarea,lon,lat,ig,jg,fr_gcm,j_dum,pfs,j_dum,fr_cat,j_dum
       endif
       if (typ == 100) ip2 = n
       if(ierr /= 0)write (*,*)'Problem reading'
    end do

    close (10,status='keep')

    maxcat = ip2

    fname=trim(gfiler)//'.rst'
  
    open (10,file=fname,status='old',action='read',  &
         form='unformatted',convert='little_endian')
    allocate(tile_id(1:nx,1:ny))

    do j=1,ny
       read(10)tile_id(:,j)
    end do

    close (10,status='keep')

    allocate(alb_in(1:i_raster,1:j_raster))  
    allocate(alb_out(1:maxcat)) 
    allocate(alb_count(1:maxcat)) 
    if(regrid) then
	allocate(raster(nx,ny),stat=STATUS); VERIFY_(STATUS)
     else
        raster => alb_in
    end if

    do  ialbt = 2,2
       do ialbs = 1,2
          do yy = 2005,2005
             if (yy.eq.2005)cyy='00-04'
             if(ialbt.eq.1)albtype='BlackSky/'
             if(ialbt.eq.2)albtype='WhiteSky/'
             if(ialbt.eq.1)bw='BS'
             if(ialbt.eq.2)bw='WS'
             
             if(ialbs.eq.1)albspec='0.3_0.7/'
             if(ialbs.eq.2)albspec='0.7_5.0/'           
             ifile=trim(c_data)//'AlbMap.'//bw//'.2x5.'//trim(cyy)//        &
                  '.monthly.'//albspec(1:index(albspec,'/')-1)//'.dat' 
             ofile='clsm/AlbMap.'//bw//'.2x5.'//trim(cyy)//'.monthly-tile.' &
                  //albspec(1:index(albspec,'/')-1)//'.dat'

           open (20,file=trim(ifile),form='unformatted',&
                  convert='big_endian', &
                  action='read',status='old')
           open (30,file=trim(ofile),form='unformatted', &
                convert='big_endian', &
                action='write',status='unknown')

           do month =1,12
              read (20) alb_in
	      if(regrid) then
	         call RegridRasterReal(alb_in,raster)
	      else
                 raster = alb_in
              endif

	      alb_out = 0.
              alb_count = 0.
              do j=1,ny
                 do i=1,nx
                    if((tile_id(i,j).gt.ip1).and.(tile_id(i,j).le.ip2)) then
                       if(raster(i,j).eq.undef)   then
!                           write (*,*)'raster UNDEF',i,j,month,albtype,albspec
!                           stop
                        endif
                        if ((raster(i,j).gt.0.).and.(raster(i,j).ne.undef)) then
                           alb_out(tile_id(i,j)-ip1) = &
                                alb_out(tile_id(i,j)-ip1) + raster(i,j)
                           alb_count(tile_id(i,j)-ip1) = &
                                alb_count(tile_id(i,j)-ip1) + 1. 
                        endif
                     endif
                  end do
               end do

               do n = 1,maxcat
                  if (alb_count(n).gt.0)then
                     alb_out(n) = alb_out(n)/alb_count(n)
                  else
!                     print *,'No albedo for the tile :',n
                     alb_out(n) = alb_out(n-1)
                  endif
               end do
               write (30) alb_out
            end do
            close (20,status='keep')
            close (30,status='keep')           
            
         end do
      end do
   end do

   deallocate (tile_id,alb_in,alb_out,alb_count)
    if(regrid) then
       deallocate(raster)
    endif 
  END SUBROUTINE modis_alb_on_tiles

!----------------------------------------------------------------------

  SUBROUTINE modis_scale_para (ease_grid,gfile)
    
    implicit none
    type (date_time_type) :: gf_green_time,af_green_time,end_time, &
              bf_lai_time,af_lai_time,date_time_new
    logical :: ease_grid
    CHARACTER*20 :: version,resoln,continent
    integer :: nc_gcm,nr_gcm,nc_ocean,nr_ocean
    REAL :: latt,lont,fr_gcm,fr_cat,tsteps,zth, slr,tarea
    INTEGER :: typ,pfs,ig,jg,j_dum,ierr,indx_dum,indr1,indr2,indr3 ,ip2
    character*100 :: path,fname,fout,metpath
    character (*) :: gfile
    integer :: n,maxcat,ip
    integer :: ialbt,ialbs,yy,j,month,unit1,unit2,unit3
    character*2 :: bw
    character*5 :: cyy
    character*300 :: albtype, albspec
    character*30, dimension (2,2) :: sibname
    character*30, dimension (2,2) :: geosname
    integer, allocatable, dimension (:) :: vegcls 
    real, allocatable, dimension (:) :: &
         modisalb,scale_fac,albvr,albnr,albvf,albnf,lat,lon, &
         green,lai,sunang,snw,lai_before,lai_after,grn_before,grn_after
    real, allocatable, dimension (:) :: &
         calbvr,calbnr,calbvf,calbnf
    character*300 :: ifile1,ifile2,ofile
    integer, dimension(12), parameter :: days_in_month_nonleap = &
         (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
    integer :: day, hour, min ,secs_in_day,k
    real :: yr,mn,dy,yr1,mn1,dy1,dum, slice1,slice2

    data sibname /'albvr','albnr',  &
                'albvf','albnf'/
    data geosname /'visdr','nirdr',  &
                'visdf','nirdf'/

    fname='clsm/catchment.def'
    open (10,file=fname,status='old',action='read',form='formatted')
    read (10,*)maxcat

    allocate (albvr    (1:maxcat))
    allocate (albvf    (1:maxcat))
    allocate (albnf    (1:maxcat))
    allocate (calbvf   (1:maxcat))
    allocate (calbnf   (1:maxcat))
    allocate (modisalb (1:maxcat))
    allocate (lai      (1:maxcat))
    allocate (green    (1:maxcat))
    allocate (lai_before (1:maxcat))
    allocate (grn_before (1:maxcat))
    allocate (lai_after  (1:maxcat))
    allocate (grn_after  (1:maxcat))
    allocate (vegcls   (1:maxcat))
    allocate (sunang   (1:maxcat))
    allocate (snw      (1:maxcat))
    close (10,status='keep')

    date_time_new%year   =2002
    date_time_new%month  =1
    date_time_new%day    =1
    date_time_new%hour   =0            
    date_time_new%min    =0            
    date_time_new%sec    =0            
    date_time_new%pentad =1            
    date_time_new%dofyr  =1   

    gf_green_time   = date_time_new
    af_green_time   = date_time_new
    end_time        = date_time_new
    bf_lai_time     = date_time_new
    af_lai_time     = date_time_new   

    fname=trim(gfile)//'.til'

    open (10,file=fname,status='old',action='read',form='formatted')
    fname='clsm/mosaic_veg_typs_fracs'
    open (20,file=fname,status='old',action='read',form='formatted')

    read (10,*)ip
    read (10,*)j_dum
    read (10,'(a)')version
    read (10,*)nc_gcm
    read (10,*)nr_gcm
    read (10,'(a)')version
    read (10,*)nc_ocean
    read (10,*)nr_ocean

    do n = 1,ip
      if (ease_grid) then     
	 read(10,*,IOSTAT=ierr) typ,pfs,lon,lat,ig,jg,fr_gcm
      else
      read(10,'(I10,3E20.12,9(2I10,E20.12,I10))',IOSTAT=ierr)     &    
            typ,tarea,lont,latt,ig,jg,fr_gcm,indx_dum,pfs,j_dum,fr_cat,j_dum
      endif
       if (typ == 100) then
          ip2 = n 
          read (20,'(i8,i8,2(2x,i3),2(2x,f6.4))')     &
            indr1,indr1,vegcls(ip2),indr1,fr_gcm,fr_gcm
       endif
       if(ierr /= 0)write (*,*)'Problem reading'
    end do
    close (10,status='keep')
    close (20,status='keep')

    cyy='00-04'
    albvr    =0.
    albvf    =0.
    albnf    =0.
    calbvf   =0.
    calbnf   =0.
    modisalb =0.
    snw      =0.
    sunang   =0.
    unit1 =10
    unit2 =20
    unit3 =30

    do  ialbt = 2,2
       do ialbs = 1,2
        
        if(ialbt.eq.1)albtype='BlackSky/'
        if(ialbt.eq.2)albtype='WhiteSky/'
        if(ialbt.eq.1)bw='BS'
        if(ialbt.eq.2)bw='WS'
        if(ialbs.eq.1)albspec='0.3_0.7/'
        if(ialbs.eq.2)albspec='0.7_5.0/'         
        ifile1='clsm/AlbMap.'//bw//'.2x5.'//trim(cyy)//'.monthly-tile.'   &
             //albspec(1:index(albspec,'/')-1)//'.dat'
!        write (*,*) 'MODIS file: ',  unit1,trim(ifile1) 
!        write (*,*) '-----------------------------'

        ifile2='clsm/sibalb1.'//trim(sibname(ialbs,ialbt))//'.climatology'
!        write (*,*) 'SiB file: ', unit2, trim(ifile2)

        ofile='clsm/modis_scale_factor.'//trim(sibname(ialbs,ialbt))//'.clim'

!        write (*,*) 'Scale factor: ', unit3, trim(ofile)

        open (unit1,file=trim(ifile1),form='unformatted',convert='big_endian', &
                action='read',status='old')
        open (unit2,file=trim(ifile2),form='unformatted',convert='big_endian', &
                action='write',status='unknown')
        open (unit3,file=trim(ofile),form='unformatted',convert='big_endian', &
             action='write',status='unknown')

        unit1 = unit1 + 1
        unit2 = unit2 + 1
        unit3 = unit3 + 1
     end do
  end do

    fname='clsm/lai.dat'
    open (40,file=fname,status='old',action='read',form='unformatted', &
         convert='little_endian')
    read(40) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
    read(40) lai_before
    call Get_MidTime(yr,mn,dy,yr1,mn1,dy1,bf_lai_time)
    
    read(40) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
    read(40) lai_after
    call Get_MidTime(yr,mn,dy,yr1,mn1,dy1,af_lai_time)

    fname='clsm/green.dat'
    open (41,file=fname,status='old',action='read',form='unformatted', &
         convert='little_endian')
    read(41) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
    read(41) grn_before
    call Get_MidTime(yr,mn,dy,yr1,mn1,dy1,gf_green_time)
 
    read(41) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
    read(41) grn_after
    call Get_MidTime(yr,mn,dy,yr1,mn1,dy1,af_green_time)

  do month=1,12

     write (*,'(a48,i3)') '     Computing MODIS scale parameters for month: ',month

    calbvf   =0.
    calbnf   =0.
    albvr    =0.
    albnr    =0.
    albvf    =0.
    albnf    =0.
    tsteps   =0.
     
    do day = 1,days_in_month_nonleap(month)

    if (datetime_le_refdatetime(date_time_new,af_lai_time)) then

    else
         lai_before = lai_after
	 bf_lai_time = af_lai_time
	 read(40) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
	 read(40) lai_after
         call Get_MidTime(yr,mn,dy,yr1,mn1,dy1,af_lai_time)
    endif
       call Time_Interp_Fac (date_time_new, bf_lai_time, af_lai_time, slice1, slice2)
       lai    = (slice1*lai_before + slice2*lai_after)
!       print *,'LAI'
!       print *,bf_lai_time
!       print *,af_lai_time
!       print *,slice1,slice2
!       print *,minval(lai),maxval(lai)

    if (datetime_le_refdatetime(date_time_new,af_green_time)) then

    else
         grn_before = grn_after
	 gf_green_time = af_green_time
	 read(41) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
	 read(41) grn_after
         call Get_MidTime(yr,mn,dy,yr1,mn1,dy1,af_green_time)

    endif
       call Time_Interp_Fac (date_time_new, gf_green_time, af_green_time, slice1, slice2)
       green  = (slice1*grn_before + slice2*grn_after)
 !      print *,'GREEN'
 !      print *,gf_green_time
 !      print *,af_green_time
 !      print *,slice1,slice2
 !      print *,minval(green),maxval(green)

       call augment_date_time(86400,date_time_new)
        
        tsteps = tsteps + 1.

              call sibalb(                                    &
                   albvr,albvr,albvf,albnf,                   &
                   lai, green, 0.0, snw, vegcls, maxcat)  

              calbvf = calbvf + albvf
              calbnf = calbnf + albnf
              
     end do

     calbvf = calbvf/tsteps
     calbnf = calbnf/tsteps

     unit1 =10
     unit2 =20
     unit3 =30
     
     do  ialbt = 2,2
        do ialbs = 1,2
           
           read (unit1) (modisalb(n),n=1,maxcat)
           if(unit2==20)write (unit2) (calbvf(n),n=1,maxcat)
           if(unit2==21)write (unit2) (calbnf(n),n=1,maxcat)
        
           if(unit2==20) modisalb = modisalb/calbvf
           if(unit2==21) modisalb = modisalb/calbnf

	   do n =1, maxcat
	   if(modisalb(n).le.0)then 
	      print *,'Negative MODIS scale param at cell',n, modisalb(n)
	      print *,'Set to 1'
	      modisalb(n)=1
	   endif

	   if(modisalb(n).gt.100)then 
	      print *,'Too large MODIS scale param',n, modisalb(n)
	      print *,'Set to 1'
	      modisalb(n)=1
	   endif

	   enddo	   

           write (unit3) (modisalb(n),n=1,maxcat) 

           unit1 = unit1 + 1
           unit2 = unit2 + 1
           unit3 = unit3 + 1        
     end do
  end do

  end do
    
  deallocate (modisalb,albvr,albvf,albnf)
  deallocate (green,lai,sunang)
  deallocate (vegcls)
  deallocate (calbvf,calbnf)
  
  unit1 =10
  unit2 =20
  unit3 =30

  do  ialbt = 2,2
     do ialbs = 1,2
        
        close (unit1, status='keep')
        close (unit2, status='keep')
        close (unit3, status='keep')

        unit1 = unit1 + 1
        unit2 = unit2 + 1
        unit3 = unit3 + 1
     end do
  end do

  close (40, status='keep')
  close (41, status='keep')
  
END SUBROUTINE modis_scale_para
 
!----------------------------------------------------------------------

  SUBROUTINE make_75 (nx,ny,regrid,path,gfile)
    implicit none
    integer nc,nr,i,j,i1,i2,j1,j2,cls,ip, ii, jj,xc,xr
    integer, allocatable ::  catid(:,:),catold(:,:),cat75(:,:)
    integer sam(3,3)
    character*100 filename,path
    character (*) :: gfile
    integer :: nx,ny
    logical :: regrid

    nc = i_raster
    nr = j_raster

    filename=trim(path)//'global.cat_id.catch.DL'
    open (9,file=filename,form='formatted',status='old')
    
    filename=trim(gfile)//'.rst'

    open (10,file=filename,convert='little_endian',   &
         form='unformatted',status='old',action='read')
    
    allocate(catid(nc,nr))
    allocate(catold(nc,nr))
    catid=0
    catold=0
    
    do j=1,nr
       read (9,*)(catold(i,j),i=1,nc)
       read (10)(catid(i,j),i=1,nc)
       do i=1,nc
          if((catold(i,j).eq.0).or.(catold(i,j).gt.5999900))catid(i,j)=0
       end do
    end do

    close(9,status='keep')
    close(10,status='keep')

    deallocate(catold)
    allocate(cat75(nc/3,nr/3))
    
    cat75=0

      filename=trim(gfile)//'.7.5.rst'

      open (11,file=filename,convert='little_endian',form='unformatted',status='unknown')
      
      do j=1,1440
         j2=J*3
         j1=j2-2
         do i=1,2880
            i2=i*3
            i1=i2-2
            sam(1:3,1:3)=catid(i1:i2,j1:j2)
            call pick_cat(sam,cat75(i,j))
            !         write(*,*)cat75(i,j)
            !         pause
         end do
         write (11)(cat75(i,j),i=1,2880)
      end do
      deallocate(catid)
      deallocate(cat75)

    end SUBROUTINE make_75

!----------------------------------------------------------------------

    subroutine pick_cat(sam,clr)
      implicit none

      integer sam(9),num_val(9),i,j,cls(1),clr(1)
      num_val(1:9)=0
      clr=0
      do i=1,9
         do j=1,9
            if(sam(i).eq.sam(j))num_val(i)=num_val(i)+1
         end do
      end do
      clr=sam(maxloc(num_val))
      
    end subroutine pick_cat
        

!----------------------------------------------------------------------

  SUBROUTINE catchment_def (nx,ny,regrid,dateline,gfilet,gfiler)

    implicit none

    INTEGER, allocatable, dimension(:) :: CATID  
    integer :: n,ip,maxcat,count,k1,i1,i,j,i_sib,j_sib
    INTEGER, allocatable, dimension (:) :: id,I_INDEX,J_INDEX 
    integer :: nc_gcm,nr_gcm,nc_ocean,nr_ocean
    REAL :: lat,lon,fr_gcm,fr_cat,tarea
    INTEGER :: typ,pfs,ig,jg,j_dum,ierr,indx_dum,indr1,indr2,indr3 ,ip2
    REAL (kind=8), PARAMETER :: RADIUS=MAPL_RADIUS,pi= RASTER_PI !6371000.	
    character*100 :: path,fname,fout,metpath,gtopo30
    character (*) :: gfilet,gfiler
    character*10 :: dline
    CHARACTER*20 :: version,resoln,continent
    REAL, ALLOCATABLE :: limits(:,:)
    REAL :: mnx,mxx,mny,mxy,dx,dy,d2r,lats,sum1,sum2,dx_gcm,dy_gcm
    REAL, dimension (:), allocatable :: tile_ele, tile_area,tile_area_land  
    integer :: nx,ny,status
    logical :: regrid
    real, pointer :: Raster(:,:)

    character*2 :: dateline
    real*4, allocatable , target :: q0 (:,:)

    gtopo30   = 'data/CATCH/srtm30_withKMS_2.5x2.5min.data'
    allocate (q0(1:i_raster,1:j_raster))

    i_sib = nx
    j_sib = ny

    dx  = 360._8/i_sib
    dy  = 180._8/j_sib
    d2r = PI/180._8

    open (10,file=trim(gtopo30),form='unformatted',status='old')
    read (10) q0
    close (10,status='keep')

    if(regrid) then
       allocate(raster(nx,ny),stat=STATUS); VERIFY_(STATUS)
    else
       raster => q0
    end if

    if(regrid) then
       call RegridRasterReal(q0,raster)
    endif

    allocate (catid(1:i_sib))
    catid=0
    fname=trim(gfilet)//'.til'

    open (10,file=fname,status='old',action='read',form='formatted')
    read (10,*)ip
    allocate(id(ip))
    allocate(i_index(ip))
    allocate(j_index(ip))
    allocate(tile_area(ip))  
    id=0
    read (10,*)j_dum
    read (10,'(a)')version
    read (10,*)nc_gcm
    read (10,*)nr_gcm
    read (10,'(a)')version
    read (10,*)nc_ocean
    read (10,*)nr_ocean

    dx_gcm = 360./float(nc_gcm)
    dy_gcm = 180./float(nr_gcm)    

    do n = 1,ip
 
      read(10,'(I10,3E20.12,9(2I10,E20.12,I10))',IOSTAT=ierr)     &    
            typ,tarea,lon,lat,ig,jg,fr_gcm,indx_dum,pfs,j_dum,fr_cat,j_dum
	    tile_area(n) = tarea
            id(n)=pfs
	    i_index(n) = ig 
	    j_index(n) = jg 

       if (typ == 100) ip2 = n
       if(ierr /= 0)write (*,*)'Problem reading'
    end do
    close (10,status='keep')

    maxcat=ip2-ip1

! Tile elevation
    allocate(tile_ele(1:maxcat))
    allocate(tile_area_land(1:maxcat))
    tile_ele = 0.
    tile_area_land = 0.

    fname=trim(gfiler)//'.rst'    
    open (10,file=fname,status='old',action='read',form='unformatted',convert='little_endian')
    
    do j=1,j_sib

       lats = -90._8 + (j - 0.5_8)*dy
       read (10)(catid(i),i=1,i_sib)

       do i=1,i_sib          
             if((catid(i) > ip1).and.(catid(i) <= ip2))then
	        tile_ele(catid(i)-ip1) = tile_ele(catid(i)-ip1) + raster(i,j)*   &	
			(sin(d2r*(lats+0.5*dy)) -sin(d2r*(lats-0.5*dy)))*(dx*d2r)
		tile_area_land(catid(i)-ip1) = tile_area_land(catid(i)-ip1) +  &	
			(sin(d2r*(lats+0.5*dy)) -sin(d2r*(lats-0.5*dy)))*(dx*d2r)		 
	     endif
	enddo	
    enddo
    tile_ele = tile_ele/tile_area_land
    close (10, status='keep')  

    ! adjustment Global Mean Topography to 614.649 (615.662 GTOPO 30) m
    ! --------------------------------------------
    sum1=0.
    sum2=0. 
    do j=1,maxcat	 
         sum1 = sum1 + tile_ele(j)*tile_area(j)
    enddo
    if(sum1/sum(tile_area(1:maxcat)).ne. 614.649D0 ) then
!	print *,sum1/sum(tile_area(1:maxcat))
	tile_ele =tile_ele*(614.649D0 / (sum1/sum(tile_area(1:maxcat))))				  	
	sum1=0.
	sum2=0. 
	do j=1,maxcat	 
	   sum1 = sum1 + tile_ele(j)*tile_area(j)
	enddo
!	print *,sum1/sum(tile_area(1:maxcat))
    endif
    

! catchment def file
! ------------------
    allocate(limits(1:maxcat,1:4))
    limits(:,1)=360.
    limits(:,2)=-360.
    limits(:,3)=90.
    limits(:,4)=-90.

    fname=trim(gfiler)//'.rst'    
    open (10,file=fname,status='old',action='read',form='unformatted',convert='little_endian')
    
    do j=1,j_sib
       mny=-90. + float(j-1)*180./float(j_sib)
       mxy=-90. + float(j)  *180./float(j_sib)
       read (10)(catid(i),i=1,i_sib)       
       if (index(dateline,'DE')/=0) then
          
          do i=1,i_sib          
             if((catid(i) > ip1).and.(catid(i) <= ip2))then
                mnx =-180. + float(i-1)*360./float(i_sib)
                mxx =-180. + float(i)  *360./float(i_sib)
                if(mnx .lt.limits(catid(i)-ip1,1))limits(catid(i)-ip1,1)=mnx 
                if(mxx .gt.limits(catid(i)-ip1,2))limits(catid(i)-ip1,2)=mxx 
                if(mny .lt.limits(catid(i)-ip1,3))limits(catid(i)-ip1,3)=mny 
                if(mxy .gt.limits(catid(i)-ip1,4))limits(catid(i)-ip1,4)=mxy 
             endif
          end do
       else
         do i=1,i_sib- i_sib/nc_gcm/2         
            if((catid(i) > ip1).and.(catid(i) <= ip2))then
               mnx =-180. + float(i-1)*360./float(i_sib)
               mxx =-180. + float(i)  *360./float(i_sib)
               if(mnx .lt.limits(catid(i)-ip1,1))limits(catid(i)-ip1,1)=mnx 
               if(mxx .gt.limits(catid(i)-ip1,2))limits(catid(i)-ip1,2)=mxx 
               if(mny .lt.limits(catid(i)-ip1,3))limits(catid(i)-ip1,3)=mny 
               if(mxy .gt.limits(catid(i)-ip1,4))limits(catid(i)-ip1,4)=mxy 
            endif
         end do
         do i=i_sib- i_sib/nc_gcm/2  +1,i_sib       
            if((catid(i) > ip1).and.(catid(i) <= ip2))then               
               mnx =-360. -180. + float(i-1)*360./float(i_sib)
               mxx =-360. -180. + float(i)  *360./float(i_sib)
               if(mnx < -180. ) mnx = mnx + 360.
               if(mxx <= -180.) mxx = mxx + 360.
               if(mnx .lt.limits(catid(i)-ip1,1))limits(catid(i)-ip1,1)=mnx 
               if(mxx .gt.limits(catid(i)-ip1,2))limits(catid(i)-ip1,2)=mxx 
               if(mny .lt.limits(catid(i)-ip1,3))limits(catid(i)-ip1,3)=mny 
               if(mxy .gt.limits(catid(i)-ip1,4))limits(catid(i)-ip1,4)=mxy 
            endif
         end do
 
      endif

    end do
    close(10,status='keep')
    
    open (10,file='clsm//catchment.def',  &
         form='formatted',status='unknown')
    write (10,*)maxcat

    do j=1,maxcat
       if(limits(j,1).lt.-180.) limits(j,1)= limits(j,1)+360.
       if(limits(j,2).le.-180.) limits(j,2)= limits(j,2)+360.  
 !      if(trim(dateline)=='DC')then
 !         limits(j,1) = max(limits(j,1),(i_index(j)-1)*dx_gcm -180. - dx_gcm/2.)       
 !         limits(j,2) = min(limits(j,2),(i_index(j)-1)*dx_gcm -180. + dx_gcm/2.)  
 !      endif
       write (10,'(i8,i8,5(2x,f9.4))')j+ip1,id(j+ip1),limits(j,1),   &
            limits(j,2),limits(j,3),limits(j,4),tile_ele(j)       
    end do

    deallocate (limits)
    deallocate (catid)
    deallocate (q0)
    if(regrid) then
       deallocate(raster)
    endif 
  END SUBROUTINE catchment_def
 
!----------------------------------------------------------------------
  
  SUBROUTINE create_soil_types_files (nx,ny,ease_grid,gfilet,gfiler)

    implicit none

    !     This program reads global 5'x5' soil texture classification,
    !     then find the dominant Soil Classes for the GCM
    !     http://www.ngdc.noaa.gov/seg/eco/cdroms/reynolds/reynolds/reynolds.htm
    !     http://www.ngdc.noaa.gov/ecosys/cdroms/reynolds/reynolds/reynolds.htm
    !     Soil texture classification follows USDA classification (13 classes)
    !     min. value  : 0
    !     max. value  : 12
    !     legend cats : 13
    !     category 0  : Ocean/No Data
    !     category 1  : Sand
    !     category 2  : Loamy Sand
    !     category 3  : Sandy Loam
    !     category 4  : Silt Loam
    !     category 5  : Silt
    !     category 6  : Loam
    !     category 7  : Sandy Clay Loam
    !     category 8  : Silty Clay Loam
    !     category 9  : Clay Loam
    !     category 10 : Sandy Clay
    !     category 11 : Silty Clay
    !     category 12 : Clay
    !     

    INTEGER col,row,i,j,k,ii,n,ip
    INTEGER colsib,rowsib,isol,jsol,cls,clr1(1),clr2(1)
    PARAMETER(col=4320,row=2160)
    
    INTEGER, allocatable :: SIB_LAY(:,:)
    INTEGER, allocatable, target :: CATID(:,:)
    INTEGER, allocatable :: SOIL1(:,:)
    INTEGER, allocatable :: SOIL2(:,:)
    INTEGER tem1 (13),tem2(13),tem3(13)
    INTEGER, ALLOCATABLE :: TOP(:,:),COM(:,:)
    INTEGER IDVAL,STEX
    INTEGER (kind=1), allocatable :: gtext(:,:)
    INTEGER irrecs, c1,c2,r1,r2
    CHARACTER*100 ifile,ifile2,ofile1,ofile2,fname
    CHARACTER (*) :: gfiler,gfilet
    character*10 :: dline
    CHARACTER*20 :: version,resoln    
    INTEGER, allocatable, dimension (:) :: id !indx,id,indx_old
    integer :: nc_gcm,nr_gcm,nc_ocean,nr_ocean
    REAL :: lat,lon,fr_gcm,fr_cat,tarea
    INTEGER :: typ,pfs,ig,jg,j_dum,ierr,indx_dum,indr1,indr2,indr3 ,ip2
    integer :: nx,ny,status
    logical :: ease_grid 
    colsib = nx
    rowsib = ny
    !
    !	Compute the number of input records per row.
    irrecs = nint (col / 4.0)      
    !
    allocate(catid(1:nx,1:ny))
    catid =0

    ifile=trim(c_data)//'dtex_tp1.bin'
    ifile2=trim(c_data)//'dtex_sb1.bin'
    ofile1='clsm/soil_text.top'
    ofile2='clsm/soil_text.com'     

    fname=trim(gfilet)//'.til'

    open (10,file=fname,status='old',action='read',form='formatted')
    read (10,*)ip
    read (10,*)j_dum
    read (10,'(a)')version
    read (10,*)nc_gcm
    read (10,*)nr_gcm
    read (10,'(a)')version
    read (10,*)nc_ocean
    read (10,*)nr_ocean

    allocate(id(ip))
    id=0
 
    do n = 1,ip
       if (ease_grid) then     
	 read(10,*,IOSTAT=ierr) typ,pfs,lon,lat,ig,jg,fr_gcm
       else 
       read(10,'(I10,3E20.12,9(2I10,E20.12,I10))',IOSTAT=ierr)     &    
            typ,tarea,lon,lat,ig,jg,fr_gcm,indx_dum,pfs,j_dum,fr_cat,j_dum
	endif
         id(n)=pfs
         if (typ == 100) ip2 = n

         if(ierr /= 0)write (*,*)'Problem reading'
      end do

      close (10,status='keep')
!
      fname=trim(gfiler)//'.rst'

      open (1,file=fname,form='unformatted',status='old', convert='little_endian')

      do j=1,ny
         read (1)(catid(i,j),i=1,nx)
      end do

      close(1,status='keep')
      ! write(*,*)'Finished reading CAT_IDs'

      !     Top layer soil classification 0-30cm
      !
!      open (unit=11, file=ifile, form='unformatted', status='old',access='direct',recl=1,  &
!           convert = 'little_endian')
      open (unit=11, file=ifile, form='unformatted', status='old',  &
           convert = 'big_endian')

      !
      allocate(gtext(1:col,1:row))
      allocate(SIB_LAY(1:nx,1:ny))
      gtext(:,:)=0
      SIB_LAY(:,:)=0
      k=0
      do j=row,1,-1
!           do i=1,irrecs
!              k=k+1
!              c1 = (4*i)-3
!              c2 = (4*i)
!              read (unit=11, rec=k) (gtext(ii,j), ii=c1,c2)
              read (unit=11) (gtext(i,j), i=1,col)
!           end do
        end do

        close (11,status='keep')
!
        do j=1,rowsib
           jsol=CEILING(j/(ny/real(row)))
           do i=1,colsib
              isol=CEILING(i/(nx/real(col)))
              SIB_LAY(i,j)=gtext(isol,jsol)
           end do
        end do

        deallocate(gtext)
        !
        !     Top layer on 2x2.5
        allocate(soil1(ip2,1:13))
        soil1(:,:)=0
        do j=1,rowsib
           do i=1,colsib
              if((catid(i,j) > ip1).and.(catid(i,j) <= ip2))then
                 IDVAL=catid(i,j) 
                 STEX=SIB_LAY(i,j)
                 SOIL1(IDVAL,STEX+1)=SOIL1(IDVAL,STEX+1)+1
              end if
           end do
        end do
        !
        !  write(*,*)'Finished reading top layer'
        deallocate(sib_lay)
        !
        !     Bottom layer soil classification 30-100cm
        !
! 	open (unit=11, file=ifile2, form='unformatted', status='old',access='direct',recl=1,  &
!             convert = 'big_endian')
 	open (unit=11, file=ifile2, form='unformatted', status='old',  &
             convert = 'big_endian')

        !
        allocate(gtext(1:col,1:row))
        allocate(SIB_LAY(1:colsib,1:rowsib))
        gtext(:,:)=0
        SIB_LAY(:,:)=0
        k=0
        do j=row,1,-1
!           do i=1,irrecs
!              k=k+1
!              c1 = (4*i)-3
!              c2 = (4*i)
!              read (unit=11, rec=k) (gtext(ii,j), ii=c1,c2)
              read (unit=11) (gtext(i,j), i=1,col)
!           end do
        end do
        !
        close (11,status='keep')

        do j=1,rowsib
           jsol=CEILING(j/(ny/real(row)))
           do i=1,colsib
              isol=CEILING(i/(nx/real(col)))
              SIB_LAY(i,j)=gtext(isol,jsol)
           end do
        end do
        deallocate(gtext)
        ! write(*,*)'Finished reading bottom layer'
        !
        !     Bottom layer on 2x2.5
        allocate(soil2(ip2,1:13))
        soil2(:,:)=0
        do j=1,rowsib
           do i=1,colsib
              if((catid(i,j) > ip1).and.(catid(i,j) <= ip2))then
                 IDVAL=catid(i,j)  

                 STEX=SIB_LAY(i,j)
                 SOIL2(IDVAL,STEX+1)=SOIL2(IDVAL,STEX+1)+1
              endif 
           end do
        end do
        deallocate(sib_lay)
        !
!         write(*,*)'Finished counting pixels for each catchment'
        k=0
        allocate(top(ip2,2))
        allocate(com(ip2,2))
        top=0
        com=0
        do j=1,ip2
           tem1(1:13)=SOIL1(j,1:13)
           tem2(1:13)=SOIL2(j,1:13)

           tem3(:)=3*tem1(:)+7*tem2(:)
           if((sum(tem3).gt.0).and.(sum(tem1).eq.0))then
              tem1(:)=tem3(:)
              write(*,*)'Filled from the bottom layer',j
           end if
            if(sum(tem1).gt.0)then
!              k=k+1
!              ! 
!             clr1=maxloc(tem1)
!             clr2=maxloc(tem3)
!             top(k,1)=j 
!             top(k,2)=clr1(1)-1
!             com(k,1)=j 
!             com(k,2)=clr2(1)-1             
              k=k+1
              ! 
             clr1=maxloc(tem1)
             clr2=maxloc(tem3)
             top(j,1)=j 
             top(j,2)=clr1(1)-1
             com(j,1)=j 
             com(j,2)=clr2(1)-1      
          end if
       end do
       !
       open (unit=11, file=ofile1, form='formatted', status='unknown')
       open (unit=12, file=ofile2, form='formatted', status='unknown')

       !
       if(top(1,2).eq.0)top(1,2)= 3
       if(com(1,2).eq.0)com(1,2)= 9

       do j=1,ip2

          if(top(j,2).eq.0)top(j,2)=top(j-1,2)
          if(com(j,2).eq.0)com(j,2)=com(j-1,2)
 
 !         if(com(j,1).gt.0)then
 !            if(j.gt.1)then
 !               if(top(j,2).eq.0)top(j,2)=top(j-1,2)
 !               if(com(j,2).eq.0)com(j,2)=com(j-1,2)
 !            end if
             !
             write(11,*)j,id(j),top(j,2)
             write(12,*)j,id(j),com(j,2)

      end do
      close(11)
      close(12) 
      deallocate (CATID,soil1,soil2,top,com,id)

  END SUBROUTINE create_soil_types_files

!----------------------------------------------------------------------

  SUBROUTINE compute_mosaic_veg_types (nx,ny,ease_grid,regrid,gfilet,gfiler)

    implicit none
    integer*1, allocatable , dimension (:,:)  :: sib_veg2
    integer,   allocatable , target , dimension (:,:)  :: sib_veg   
    integer, allocatable ::  mos_veg(:,:)
    real, allocatable :: veg_frac(:,:)
    integer :: j,k,mcls
    INTEGER CATID(nx),cls,cls2,bcls
    REAL, allocatable :: veg(:,:),bare_frac(:),zdep2_g(:,:)
    REAL :: fmax0,dummy,tem(6),mfrac,sfrac,bfrac
    
    integer :: n,ip,maxcat,count,k1,i1,i
    INTEGER, allocatable, dimension (:) :: id ! indx,id,indx_old
    integer :: nc_gcm,nr_gcm,nc_ocean,nr_ocean
    REAL :: lat,lon,fr_gcm,fr_cat,tarea
    INTEGER :: typ,pfs,ig,jg,i_dum,j_dum,ierr,indx_dum,indr1,indr2,indr3 ,ip2
    character*100 :: fname,fout
    character (*) :: gfiler,gfilet
    character*10 :: dline
    CHARACTER*20 :: version,resoln,continent
    character*2 :: chyear
    integer :: mon,smon,imon,year
    integer :: nx,ny,status
    logical :: regrid,ease_grid
    integer, pointer :: Raster(:,:)
    real, pointer, dimension (:)  :: z2, z0
    real, dimension (6) :: VGZ2 = (/35.0, 20.0, 17.0, 0.6, 0.5, 0.6/) ! Dorman and Sellers (1989)
    logical :: jpl_height = .true.
 
    fname=trim(gfilet)//'.til'
    open (10,file=fname,status='old',action='read',form='formatted')
    read (10,*)ip
    allocate(id(1:ip))
    
    read (10,*)j_dum
    read (10,'(a)')version
    read (10,*)nc_gcm
    read (10,*)nr_gcm
    read (10,'(a)')version
    read (10,*)nc_ocean
    read (10,*)nr_ocean
    
    do n = 1,ip
      if (ease_grid) then     
	 read(10,*,IOSTAT=ierr) typ,pfs,lon,lat,ig,jg,fr_gcm
      else
      read(10,'(I10,3E20.12,9(2I10,E20.12,I10))',IOSTAT=ierr)     &    
            typ,tarea,lon,lat,ig,jg,fr_gcm,indx_dum,pfs,i_dum,fr_cat,j_dum
      endif
       id(n)=pfs      
       if (typ == 100) ip2 = n
       if(ierr /= 0)write (*,*)'Problem reading'
    end do
    close (10,status='keep')
    maxcat=ip2
    
    allocate(sib_veg2(1:i_raster,1:j_raster))
    allocate(sib_veg (1:i_raster,1:j_raster))

    open (10,file=trim(c_data)//'sib22.5_v2.0.dat',form='unformatted',      &
         status='old',action='read',convert='big_endian')

    READ(10)sib_veg2
    
    close (10,status='keep')
    sib_veg = sib_veg2
    if(regrid) then
       allocate(raster(nx,ny),stat=STATUS); VERIFY_(STATUS)
    else
       raster => sib_veg
    end if
 
    if(regrid) then
       call RegridRaster(sib_veg,raster)
    endif  

    fname=trim(gfiler)//'.rst'  

    open (10,file=fname,status='old',action='read',form='unformatted',convert='little_endian')
    
    allocate(veg(1:maxcat,1:6))
    allocate(zdep2_g(1:maxcat,1:1))
    
    veg=0.
    zdep2_g=0.
    
    n=1

    do j=1,ny
       
       read (10)(catid(i),i=1,nx)
       
       do i=1,nx
          
          if((catid(i) > ip1).and.(catid(i) <= ip2))then
             zdep2_g(catid(i)-ip1,1)=zdep2_g(catid(i)-ip1,1)+1.
             if(raster(i,j).eq.0) then
!                write (*,*)'Warning : SiB2 =0, an ocean pixel found !'
             elseif (raster(i,j).eq.1) then
                veg(catid(i)-ip1,1)=veg(catid(i)-ip1,1) + 1.
             elseif (raster(i,j).eq.2) then
                veg(catid(i)-ip1,2)=veg(catid(i)-ip1,2) + 1.
             elseif (raster(i,j).eq.3) then
                veg(catid(i)-ip1,2)=veg(catid(i)-ip1,2) + 0.5
                veg(catid(i)-ip1,3)=veg(catid(i)-ip1,3) + 0.5                  
             elseif (raster(i,j).eq.4) then
                veg(catid(i)-ip1,3)=veg(catid(i)-ip1,3) + 1.
             elseif (raster(i,j).eq.5) then
                veg(catid(i)-ip1,3)=veg(catid(i)-ip1,3) + 1.
             elseif (raster(i,j).eq.6) then
                veg(catid(i)-ip1,4)=veg(catid(i)-ip1,4) + 1.
             elseif (raster(i,j).eq.7) then
                veg(catid(i)-ip1,5)=veg(catid(i)-ip1,5) + 1.
             elseif (raster(i,j).eq.8) then
!               if (j >= NINT(float(ny)*(140./180.))) then 
!                   veg(catid(i)-ip1,6)=veg(catid(i)-ip1,6) + 1.
!                else
!                   veg(catid(i)-ip1,5)=veg(catid(i)-ip1,5) + 1.   
!                endif
               if ((j > NINT(float(ny)*(40./180.))).and.(j < NINT(float(ny)*(140./180.)))) then 
                   veg(catid(i)-ip1,5)=veg(catid(i)-ip1,5) + 1.
                else
                   veg(catid(i)-ip1,6)=veg(catid(i)-ip1,6) + 1.   
                endif
             elseif (raster(i,j).eq.9) then
                veg(catid(i)-ip1,4)=veg(catid(i)-ip1,4) + 1.
             elseif (raster(i,j).eq.10) then
!                write (*,*)'Warning : SiB2 =10, a water pixel found !'
             elseif (raster(i,j).eq.11) then
!                write (*,*)'Warning : SiB2 =11, an ice pixel found !'
             elseif (raster(i,j).eq.100) then
!                write (*,*)'Warning : SiB2 =100, NODATA pixel found !'
             endif
          endif
       enddo
    enddo
    
    close(10,status='keep')
    
    allocate(mos_veg(1:maxcat,1:2))
    allocate(veg_frac(1:maxcat,1:3))
    mos_veg=0
    veg_frac=0.
    
    k=0
    do j=1,maxcat
       tem(1:6)=veg(j,1:6)
  
       if(sum(tem).le.0.)write(*,*) 'Warning no veg types',j
!       if(sum(tem).le.0.) stop
       if(sum(tem).gt.0)then
          
          k=k+1
          
          mfrac=-10.
          sfrac=-10.
          bfrac=0.
          cls=100
          cls2=100
          do n=1,6
             
             if(mfrac.le.tem(n))then
                sfrac=mfrac
                cls2=cls
                mfrac=tem(n)
                cls=n
             elseif(sfrac.le.tem(n)) then
                if(tem(n).lt.mfrac)then
                   sfrac=tem(n)
                   cls2=n
                endif
             endif
          enddo
          
          mos_veg(k,1)=cls
          mos_veg(k,2)=cls2
          veg_frac(k,1)=mfrac/zdep2_g(k,1)
          veg_frac(k,2)=sfrac/zdep2_g(k,1)
          veg_frac(k,3)=1.-mfrac-sfrac
       else
          k = k + 1
          mos_veg(k,1) = mos_veg(k-1,1)
          mos_veg(k,2) = mos_veg(k-1,2)
          veg_frac(k,1)= veg_frac(k-1,1)
          veg_frac(k,2)= veg_frac(k-1,2)
          veg_frac(k,3)= veg_frac(k-1,3)
       endif
       if(veg_frac(k,1).eq.0.)then
          write(*,*)'Checking for 100% desert soil tiles',k,mos_veg(k,1)&
               ,mos_veg(k,2),veg_frac(k,1),veg_frac(k,2),veg_frac(k,3)
          if(veg_frac(k,3).ne.1.)then
             write(*,*)'it is not 100% desert soil either'
          endif
          
       endif
       if(mos_veg(k,1).eq.100) then
          write(*,*)'Checking for 100% desert soil tiles',k,mos_veg(k,1)&
               ,mos_veg(k,2),veg_frac(k,1),veg_frac(k,2),veg_frac(k,3)
          write(*,*) 'Prob1'
       endif
       if(mos_veg(k,2).eq.100) then
          write(*,*) 'Prob1'
          mos_veg(k,2)=7
          veg_frac(k,2)=0.
          write(*,*)k,tem
          write(*,*)mos_veg(j,1),mos_veg(j,2),veg_frac(j,1),veg_frac(j,2),veg_frac(j,3)
       endif       
    end do
    deallocate(veg)

 ! Canopy height and ASCAT roughness length

    call ascat_r0 (nx,ny,gfiler, z0)

    if(jpl_height) then
       call jpl_canoph (nx,ny,gfiler, z2)
    else
       allocate (z2(1:maxcat))       
    endif   
    
    open (10,file='clsm/mosaic_veg_typs_fracs',  &
         form='formatted',status='unknown')
    do j=1,maxcat
       if (mos_veg(j,1) == 0) then
            if(.not.jpl_height) z2(j) = VGZ2(mos_veg(j,1))
            mos_veg(j,1) = mos_veg(j-1,1)
            mos_veg(j,2) = mos_veg(j-1,2)
            write (10,'(i8,i8,2(2x,i3),2(2x,f6.2),2x,f6.3,2x,f10.7)')     &
                 j+ip1,id(j+ip1),mos_veg(j-1,1),mos_veg(j-1,2),veg_frac(j,1),veg_frac(j,2),z2(j), z0 (j)           
       else
       if(.not.jpl_height) z2(j) = VGZ2(mos_veg(j,1))
       write (10,'(i8,i8,2(2x,i3),2(2x,f6.2),2x,f6.3,2x,f10.7)')     &
            j+ip1,id(j+ip1),mos_veg(j,1),mos_veg(j,2),veg_frac(j,1),veg_frac(j,2),z2(j), z0 (j)
       endif
    end do
    close(10,status='keep')

    open (20,file='clsm/vegdyn.data',status='unknown',action='write',form='unformatted', &
         convert='little_endian')      

    write (20) real(mos_veg(:,1))
    write (20) z2 (:)
    write (20) z0 (:)    
    close (20)   

    deallocate (sib_veg2,sib_veg,mos_veg,veg_frac,zdep2_g,id,  z0, z2)
    if(regrid) then
       deallocate(raster)
    endif     

  END SUBROUTINE compute_mosaic_veg_types
  
!----------------------------------------------------------------------

  SUBROUTINE cti_stat_file (ease_grid,gfile, MaskFile)

    IMPLICIT NONE
    INTEGER, PARAMETER :: nbcat=36716,nofvar=6, SRTM_maxcat = 291284
    INTEGER :: n,i,ip, itext(SRTM_maxcat,2),ix, jx,ip2, maxcat
    INTEGER :: typ,pfs,ig,jg,j_dum,ierr,indx_dum,indr1,indr2,indr3
    INTEGER*8 :: idum8
    INTEGER :: ncat,i_dum
    INTEGER, dimension(:), allocatable :: colin2cat 
    INTEGER, allocatable, dimension (:) :: id,indx_old
    integer :: idum1,idum2
    real :: dum1,dum2,dum3,dum4,dum5,dum6
    integer :: nc_gcm,nr_gcm,nc_ocean,nr_ocean
    REAL :: lat,lon,fr_gcm,fr_cat,tarea
    REAL :: fr
    REAL, allocatable, dimension (:,:) :: var
    REAL, allocatable, dimension (:) :: dummy
    logical :: ease_grid
    CHARACTER*20 :: version
    character*100 :: fname
    character(*) :: gfile
    character(*) :: MaskFile
    !

    fname=trim(gfile)//'.til'
    open (10,file=fname,status='old',action='read',form='formatted')
    read (10,*)ip
    allocate(indx_old(ip))
    allocate(id(ip))
    indx_old=0
    id=0
 
    read (10,*)j_dum
    read (10,'(a)')version
    read (10,*)nc_gcm
    read (10,*)nr_gcm
    read (10,'(a)')version
    read (10,*)nc_ocean
    read (10,*)nr_ocean
    
    do n = 1,ip
      if (ease_grid) then  
         read(10,*,IOSTAT=ierr) typ,pfs,lon,lat,ig,jg,fr_gcm,i_dum
      else
	 read(10,'(I10,3E20.12,9(2I10,E20.12,I10))',IOSTAT=ierr)     &    
            typ,tarea,lon,lat,ig,jg,fr_gcm,indx_dum,pfs,i_dum,fr_cat,j_dum
      endif

       id(n)=pfs
       indx_old(n) = i_dum
       if (index(MaskFile,'GEOS5_10arcsec_mask') /= 0) indx_old(n) = pfs
       if (typ == 100) ip2 = n
       if(ierr /= 0)write (*,*)'Problem reading',fname
       if(ierr /= 0) stop
    end do
    
    close (10,status='keep')
!
    allocate(colin2cat(1:6000000))
    colin2cat=0
    if (index(MaskFile,'GEOS5_10arcsec_mask') /= 0) then
       open (10,file=trim(c_data)//'/SRTM-TopoData/Pfafcatch-routing.dat',   &
            form='formatted', status='old',action='read')       
    else
       open (10,file=trim(c_data)//'/catchment.def',   &
            form='formatted', status='old',action='read')
    endif

    read (10,*) ncat
    do n=1,ncat
        if (index(MaskFile,'GEOS5_10arcsec_mask') /= 0) then
           read (10,*)indx_dum,idum8
        else
           read (10,*)j_dum,indx_dum
        endif
       colin2cat(indx_dum)=n
    end do
    close (10,status='keep')

    if (index(MaskFile,'GEOS5_10arcsec_mask') /= 0) then
       open (10,file=trim(c_data)//'/SRTM-TopoData/SRTM_cti_stats.dat',       &
            form='formatted', status='old',action='read')
    else
       open (10,file=trim(c_data)//'/cti_stats.dat',       &
            form='formatted', status='old',action='read')
    endif

    fname='clsm/cti_stats.dat'
    open (20,file=fname,form='formatted', status='unknown')
    write (20,*)ip2

    read (10,*)ncat
    allocate(var(1:ncat,1:nofvar))
    var=0.
    pfs=1
    do i=1,ncat
       if (index(MaskFile,'GEOS5_10arcsec_mask') /= 0) then
          read(10,'(i8,i15,5(1x,f8.4),i5,e18.3)')n,idum8,dum1,dum2,dum3,dum4,dum5,idum2,dum6
          idum1 = n
       else
          read(10,'(i8,i8,5(1x,f8.4),i5,e18.3)')n,idum1,dum1,dum2,dum3,dum4,dum5,idum2,dum6
       endif
       if(colin2cat(idum1).gt.0)then
          itext(pfs,1)=idum1
          var(pfs,1)=dum1
          var(pfs,2)=dum2
          var(pfs,3)=dum3
          var(pfs,4)=dum4
          var(pfs,5)=dum5
          itext(pfs,2)=idum2
          var(pfs,6)=dum6
          pfs=pfs+1
       endif
    end do

    close (10,status='keep')
!
    do i=1,ip

       if((i > ip1).and.(i <= ip2))then
          if(((id(i).ge.5000142).and.(id(i).le.5025829)))then
             write(20,'(i8,i8,5(1x,f8.4),i5,e18.3)')i,id(i),var(indx_old(i),1)*11.1/9.1,var(indx_old(i),2), &
                  var(indx_old(i),3),var(indx_old(i),4),var(indx_old(i),5),itext(indx_old(i),2),            &
                  var(indx_old(i),6)
          else

             write(20,'(i8,i8,5(1x,f8.4),i5,e18.3)')i,id(i),var(indx_old(i),1),var(indx_old(i),2), & 
                  var(indx_old(i),3),var(indx_old(i),4),var(indx_old(i),5),itext(indx_old(i),2),   &
                  var(indx_old(i),6)
          endif
       endif
       
    end do
    close (20,status='keep')
    deallocate (colin2cat,var,id,indx_old)

  END SUBROUTINE cti_stat_file
 
!---------------------------------------------------------------------

  SUBROUTINE create_model_para (MaskFile)

    implicit none
      integer i,n,k, tindex1,pfaf1,nbcatch
      integer soil_gswp
      real meanlu,stdev,minlu,maxlu,coesk,rzdep
      real minlat,maxlat,minlon,maxlon
      real,allocatable, dimension (:) ::   &
             BEE, PSIS,POROS,COND,WPWET,soildepth, tile_lon, tile_lat
      REAL, allocatable, dimension(:) :: TOPMEAN, TOPVAR, TOPSKEW
      REAL ST(NAR), AC(NAR),COESKEW
      REAL, allocatable, dimension (:) ::   &
	    ARS1,ARS2,ARS3, ARA1,ARA2,ARA3,ARA4, &
	    ARW1,ARW2,ARW3,ARW4,bf1, bf2, bf3,   &
	    tsa1, tsa2,tsb1, tsb2,               &
	    taberr1,taberr2,normerr1,normerr2,   &
	    taberr3,taberr4,normerr3,normerr4
      integer, dimension(12) :: tile_pick
      integer, allocatable, dimension (:) :: soil_class_top,soil_class_com,tindex2,pfaf2
      real watdep(nwt,nrz),wan(nwt,nrz),rzexcn(nwt,nrz),frc(nwt,nrz)
      real, allocatable, dimension  (:,:,:,:) :: &
      gwatdep,gwan,grzexcn,gfrc
      real :: wtdep,wanom,rzaact,fracl,profdep,dist_save,tile_distance 
      character*100 :: pathout,fname,fout,losfile
      character*10 :: dline
      CHARACTER*20 :: version,resoln,continent
      character*6 rdep,ext
      integer :: iwt,irz,group
      character(*) :: MaskFile
      logical :: picked
 
! --------- VARIABLES FOR *OPENMP* PARALLEL ENVIRONMENT ------------
!
! NOTE: "!$" is for conditional compilation
!
logical :: running_omp = .false.
!
!$ integer :: omp_get_thread_num, omp_get_num_threads
!
integer :: n_threads=1, li, ui
!
integer, dimension(:), allocatable :: low_ind, upp_ind
!
! ------------------------------------------------------------------
        
  ! ----------- OpenMP PARALLEL ENVIRONMENT ----------------------------
  !
  ! FIND OUT WHETHER -omp FLAG HAS BEEN SET DURING COMPILATION
  !
  !$ running_omp = .true.         ! conditional compilation
  !
  ! ECHO BASIC OMP VARIABLES
  !
  !$OMP PARALLEL DEFAULT(NONE) SHARED(running_omp,n_threads) 
  !
  !$OMP SINGLE
  !
  !$ n_threads = omp_get_num_threads()
  !
  !$ write (*,*) 'running_omp = ', running_omp
  !$ write (*,*)
  !$ write (*,*) 'parallel OpenMP with ', n_threads, 'threads'
  !$ write (*,*)
  !$OMP ENDSINGLE
  !
  !$OMP CRITICAL
  !$ write (*,*) 'thread ', omp_get_thread_num(), ' alive'
  !$OMP ENDCRITICAL
  !
  !$OMP BARRIER
  !
  !$OMP ENDPARALLEL
      
!c-------------------------------------------------------------------------
      
      losfile =trim(c_data)//'GSWP2_loss_perday/loss_perday'
!c     opening files


      allocate (gwatdep(1:nwt,1:nrz,1:12,1:2))
      allocate (gwan   (1:nwt,1:nrz,1:12,1:2))
      allocate (grzexcn(1:nwt,1:nrz,1:12,1:2))
      allocate (gfrc   (1:nwt,1:nrz,1:12,1:2))

      do n =1,12
        if(n.lt.10)write(ext,'(i1.1)')n
        if(n.ge.10)write(ext,'(i2.2)')n
        do i =1,2
	   if (i==1) rdep='.rz75.'
	   if (i==2) rdep='.rz1.'
	   open (120,file=trim(losfile)//trim(rdep)//trim(ext),  &
           form='formatted',status='old')

	   do iwt=1,nwt
	      do irz=1,nrz
		 read(120,2000) wtdep,wanom,rzaact,fracl
 2000		 format(1x,4e16.8)
		 gwatdep(iwt,irz,n,i)=wtdep
		 gwan(iwt,irz,n,i)=wanom
		 grzexcn(iwt,irz,n,i)=rzaact
		 gfrc(iwt,irz,n,i)=amin1(fracl,1.)
	    enddo
	  enddo
      close (120,status='keep')	   

	end do
      end do    
      fname='clsm/soil_param.first'       
       open (10,file=fname,action='read',       &
          form='formatted',status='old')              
                                                 
     fname='clsm/cti_stats.dat'           
      open (11,file=fname,action='read',        &
          form='formatted',status='old')              

     fname='clsm/catchment.def'           
      open (12,file=fname,action='read',        &
          form='formatted',status='old')                                                             

     fout='clsm/ar.new'               
      open (20,file=fout,action='write',        &
          form='formatted',status='unknown')          
                                                 
     fout='clsm//bf.dat'               
      open (30,file=fout,action='write',        &
          form='formatted',status='unknown')          
                                                 
     fout='clsm//ts.dat'               
      open (40,file=fout,action='write',        &
          form='formatted',status='unknown')        

     if (error_file) then 
	fout='clsm/ar_rmse.dat'           
	open (21,file=fout,action='write',        &
	     form='formatted',status='unknown')

	fout='clsm/bf_rmse.dat'           
	open (31,file=fout,action='write',        &
	     form='formatted',status='unknown')

        fout='clsm/bad_sat_param.tiles'
 	open (41,file=fout,action='write',        &
	     form='formatted',status='unknown')    
     endif 
        fout='clsm/soil_param.dat'
 	open (42,file=fout,action='write',        &
	     form='formatted',status='unknown')    
     read (11,*)nbcatch
     read (12,*)nbcatch

     allocate (tile_lon(1:nbcatch))
     allocate (tile_lat(1:nbcatch))
     allocate (TOPMEAN (1:nbcatch))
     allocate (TOPVAR  (1:nbcatch))
     allocate (TOPSKEW (1:nbcatch))
     allocate (ARS1 (1:nbcatch))
     allocate (ARS2 (1:nbcatch))
     allocate (ARS3 (1:nbcatch))
     allocate (ARA1 (1:nbcatch))
     allocate (ARA2 (1:nbcatch))
     allocate (ARA3 (1:nbcatch))
     allocate (ARA4 (1:nbcatch))
     allocate (ARW1 (1:nbcatch))
     allocate (ARW2 (1:nbcatch))
     allocate (ARW3 (1:nbcatch))
     allocate (ARW4 (1:nbcatch))
     allocate (BF1 (1:nbcatch))
     allocate (BF2 (1:nbcatch))
     allocate (BF3 (1:nbcatch))
     allocate (TSA1 (1:nbcatch))
     allocate (TSA2 (1:nbcatch))
     allocate (TSB1 (1:nbcatch))
     allocate (TSB2 (1:nbcatch))
     allocate (TABERR1 (1:nbcatch))
     allocate (TABERR2 (1:nbcatch))
     allocate (TABERR3 (1:nbcatch))
     allocate (TABERR4 (1:nbcatch))
     allocate (NORMERR1 (1:nbcatch))
     allocate (NORMERR2 (1:nbcatch))
     allocate (NORMERR3 (1:nbcatch))
     allocate (NORMERR4 (1:nbcatch))
     allocate (BEE        (1:nbcatch))
     allocate (PSIS      (1:nbcatch))
     allocate (POROS     (1:nbcatch))
     allocate (COND      (1:nbcatch))
     allocate (WPWET     (1:nbcatch))
     allocate (soildepth (1:nbcatch))
     allocate (soil_class_top (1:nbcatch))
     allocate (soil_class_com (1:nbcatch))
     allocate (tindex2        (1:nbcatch))
     allocate (pfaf2          (1:nbcatch))
 
     do n=1,nbcatch
 
        read(11,'(i8,i8,5(1x,f8.4))') tindex1,pfaf1,meanlu,stdev  &
             ,minlu,maxlu,coesk                                  
        read(10,*) tindex2(n),pfaf2(n),soil_class_top(n),soil_class_com(n),    &
             BEE(n), PSIS(n),POROS(n),COND(n),WPWET(n),soildepth(n)                        

        if(tindex1.ne.tindex2(n))then
           write(*,*)'Warnning 1: tindex mismatched'                        
           stop
        endif

        read (12,*) tindex1,pfaf1,minlon,maxlon,minlat,maxlat
        tile_lon(n) = (minlon + maxlon)/2.
        tile_lat(n) = (minlat + maxlat)/2.

        if(pfaf1.ne.pfaf2(n)) then
           write(*,*)'Warnning 1: pfafstetter mismatched' 
           stop
        endif

        if (index(MaskFile,'GEOS5_10arcsec_mask') /= 0) then
           TOPMEAN(n) = meanlu
        else
           TOPMEAN(n) = 0.961*meanlu-1.957                       
        endif

        TOPVAR(n)  = stdev*stdev                                
        TOPSKEW(n) = coesk*stdev*stdev*stdev                   
        
        if ( TOPVAR(n) .eq. 0. .or. coesk .eq. 0.            &
             .or. topskew(n) .eq. 0.) then                       
           write(*,*) 'Problem: undefined values:'         
           write(*,*) TOPMEAN(n),TOPVAR(n),coesk,            &
                minlu,maxlu
           stop
        endif
     END DO

     rewind(10)

     allocate(low_ind(n_threads))
     allocate(upp_ind(n_threads))
     low_ind(1)         = 1
     upp_ind(n_threads) = nbcatch

     if (running_omp)  then
      do i=1,n_threads-1
        
         upp_ind(i)   = low_ind(i) + (nbcatch/n_threads) - 1 
         low_ind(i+1) = upp_ind(i) + 1
        
      end do 
     end if


!$OMP PARALLELDO DEFAULT(NONE)                          &
!$OMP SHARED( BEE, PSIS,POROS,COND,WPWET,soildepth,     &
!$OMP        TOPMEAN, TOPVAR, TOPSKEW,                  &
!$OMP        ARS1,ARS2,ARS3, ARA1,ARA2,ARA3,ARA4,       &
!$OMP        ARW1,ARW2,ARW3,ARW4,bf1, bf2, bf3,         &
!$OMP        tsa1, tsa2,tsb1, tsb2,                     &
!$OMP        taberr1,taberr2,normerr1,normerr2,         &
!$OMP        taberr3,taberr4,normerr3,normerr4,         &
!$OMP        gwatdep,gwan,grzexcn,gfrc,soil_class_com,  &
!$OMP        n_threads, low_ind, upp_ind )              &
!$OMP PRIVATE(k,li,ui,n,i,watdep,wan,rzexcn,frc,ST,AC,  &
!$OMP COESKEW,profdep)

     do k=1,n_threads

        li = low_ind(k)
        ui = upp_ind(k)

      do n=li,ui
!      if ((n == 877).or.(n == 880).or.(n == 881)) then
!         print *,n
!      endif
!      print *,n
!      pause
!        c Gamma distribution
      CALL TGEN (                              &
          TOPMEAN(n),TOPVAR(n),TOPSKEW(n),     &
          ST,AC,COESKEW)
      
!c      write(*,*) 'tgen4 ok'

!c Areal fractioning parameters
!      print *,'tileid:' ,n
      CALL SAT_PARAM(                                              &
                     BEE(n),PSIS(n),POROS(n),COND(n),              & 
                     WPWET(n), ST, AC, COESKEW,n,                  &
                     soildepth(n),                                 &
                     ars1(n),ars2(n),ars3(n),                      &
                     ara1(n),ara2(n),ara3(n),ara4(n),              &
                     arw1(n),arw2(n),arw3(n),arw4(n),              &
                     taberr1(n),taberr2(n),taberr3(n),taberr4(n),  &
                     normerr1(n),normerr2(n),normerr3(n),normerr4(n))


      CALL BASE_PARAM(                                  &
          BEE(n),PSIS(n),POROS(n),COND(n),              &
          ST, AC,                                       &
          bf1(n),bf2(n),bf3(n),                         &
          taberr1(n),taberr2(n),normerr1(n),normerr2(n) &
          )
            
	    profdep=soildepth(n)/1000.
            profdep=amax1(1.,profdep)
            if (grzdep .gt. .75*profdep) then
              i=1
              else
              i=2 
            end if

	  watdep (:,:) =  gwatdep (:,:,soil_class_com(n),i)   
	  wan    (:,:) =  gwan    (:,:,soil_class_com(n),i)
	  rzexcn (:,:) =  grzexcn (:,:,soil_class_com(n),i)
	  frc    (:,:) =  gfrc    (:,:,soil_class_com(n),i)

      CALL TS_PARAM(                       &
          BEE(n),PSIS(n),POROS(n),         &
          ST, AC,                          &
          watdep,wan,rzexcn,frc,           &
          tsa1(n),tsa2(n),tsb1(n),tsb2(n)  &
          )

     END DO
     END DO
          !$OMP ENDPARALLELDO
     tile_pick = 0
 
     DO n=1,nbcatch
        if((arw1(n).ne.9999.).and.(ars1(n).ne.9999.))then
           if(tile_pick(soil_class_com(n)) == 0)  tile_pick(soil_class_com(n)) = n
        endif
     end do

     DO n=1,nbcatch
        !c Third subroutine for the parameters related to the transfers 
        !c to the water table
        !
        ! Writing the parameters, in the same order as in catchment.def
        !      if((ars1(n).lt.0.).and.(ars2(n).le.0.3).and.(ars3(n).le.0.04).and.(arw1(n).ne.9999.))then
        if((arw1(n).ne.9999.).and.(ars1(n).ne.9999.))then   
           write(20,'(i8,i8,f5.2,11(2x,e14.7))')   &
                tindex2(n),pfaf2(n),gnu,   &
                ars1(n),ars2(n),ars3(n),                   &
                ara1(n),ara2(n),ara3(n),ara4(n),           &
                arw1(n),arw2(n),arw3(n),arw4(n) 
           write(30,'(i8,i8,f5.2,3(2x,e13.7))')tindex2(n),pfaf2(n),gnu,bf1(n),bf2(n),bf3(n)
           write(40,'(i8,i8,f5.2,4(2x,e13.7))')tindex2(n),pfaf2(n),gnu,    &
                tsa1(n),tsa2(n),tsb1(n),tsb2(n)
           
           write(42,'(i8,i8,i4,i4,3f8.4,f12.8,f7.4,f10.3)') tindex2(n),pfaf2(n),soil_class_top(n),soil_class_com(n),    &
                BEE(n), PSIS(n),POROS(n),COND(n),WPWET(n),soildepth(n)  
        else
 
       if(preserve_soiltype) then    
          picked=.false.
          ! Group3
          !     category 1  : Sand
          !     category 2  : Loamy Sand
          !     category 3  : Sandy Loam
          !     category 8  : Silty Clay Loam
          ! Group2
          !     category 4  : Silt Loam
          !     category 5  : Silt
          !     category 6  : Loam
          !     category 7  : Sandy Clay Loam    
          ! Group1
          !     category 9  : Clay Loam
          !     category 10 : Sandy Clay
          !     category 11 : Silty Clay
          !     category 12 : Clay
          
          if ((soil_class_com(n)>=9).and.(soil_class_com(n)<=12)) then	
             group=1
          else if ((soil_class_com(n)>=4).and.(soil_class_com(n)<=7)) then
             group=2
          else
             group=3
          endif
          
          if(tile_pick(soil_class_com(n)) > 0) then
             k = tile_pick(soil_class_com(n))
             picked=.true.
             if (error_file) then
                write (41,*)n,k
             endif
             write(20,'(i8,i8,f5.2,11(2x,e14.7))')   &
                  tindex2(n),pfaf2(n),gnu,   &
                  ars1(k),ars2(k),ars3(k),                   &
                  ara1(k),ara2(k),ara3(k),ara4(k),           &
                  arw1(k),arw2(k),arw3(k),arw4(k) 
             ars1(n)=ars1(k)
             ars2(n)=ars2(k)
             ars3(n)=ars3(k)
             ara1(n)=ara1(k)
             ara2(n)=ara2(k)
             ara3(n)=ara3(k)
             ara4(n)=ara4(k)
             arw1(n)=arw1(k)
             arw2(n)=arw2(k)
             arw3(n)=arw3(k)
             arw4(n)=arw4(k)
             
             write(30,'(i8,i8,f5.2,3(2x,e13.7))')tindex2(n),pfaf2(n),gnu,bf1(k),bf2(k),bf3(k)
             write(40,'(i8,i8,f5.2,4(2x,e13.7))')tindex2(n),pfaf2(n),gnu,    &
                  tsa1(k),tsa2(k),tsb1(k),tsb2(k)             
             write(42,'(i8,i8,i4,i4,3f8.4,f12.8,f7.4,f10.3)') tindex2(n),pfaf2(n),soil_class_top(k),soil_class_com(k),    &
                  BEE(k), PSIS(k),POROS(k),COND(k),WPWET(k),soildepth(k)  
             
          else
             
             do k =n-1,1,-1 
                
                if (group == 1) then
                   if ((soil_class_com(k)>=9).and.(soil_class_com(k)<=12))picked=.true.
                endif
                
                if (group == 2) then
                   if ((soil_class_com(k)>=4).and.(soil_class_com(k)<=7))	picked=.true.
                endif
                
                if (group == 3) then
                   if (((soil_class_com(k)>=1).and.(soil_class_com(k)<=3)).or. &
                        (soil_class_com(k)==8)) picked=.true.
                endif
                
                if (picked) then
                   if (error_file) then
                      write (41,*)n,k
                   endif
                   
                   write(20,'(i8,i8,f5.2,11(2x,e14.7))')   &
                        tindex2(n),pfaf2(n),gnu,   &
                        ars1(k),ars2(k),ars3(k),                   &
                        ara1(k),ara2(k),ara3(k),ara4(k),           &
                        arw1(k),arw2(k),arw3(k),arw4(k) 
                   write(30,'(i8,i8,f5.2,3(2x,e13.7))')tindex2(n),pfaf2(n),gnu,bf1(k),bf2(k),bf3(k)
                   write(40,'(i8,i8,f5.2,4(2x,e13.7))')tindex2(n),pfaf2(n),gnu,    &
                        tsa1(k),tsa2(k),tsb1(k),tsb2(k)             
                   write(42,'(i8,i8,i4,i4,3f8.4,f12.8,f7.4,f10.3)') tindex2(n),pfaf2(n),soil_class_top(k),soil_class_com(k),    &
                        BEE(k), PSIS(k),POROS(k),COND(k),WPWET(k),soildepth(k)  
                   ars1(n)=ars1(k)
                   ars2(n)=ars2(k)
                   ars3(n)=ars3(k)
                   ara1(n)=ara1(k)
                   ara2(n)=ara2(k)
                   ara3(n)=ara3(k)
                   ara4(n)=ara4(k)
                   arw1(n)=arw1(k)
                   arw2(n)=arw2(k)
                   arw3(n)=arw3(k)
                   arw4(n)=arw4(k)
                   exit
                endif
                
                if((k==1).and.not(picked)) then
                   print *,'Warning ar.new is bad at n=',n
                   print *,'Call Sarith ......'
                   stop
                endif
             end do
          endif
 
       
!        write(30,'(i8,i8,f5.2,3(2x,e13.7))')tindex2(n),pfaf2(n),gnu,bf1(n),bf2(n),bf3(n)
!        write(40,'(i8,i8,f5.2,4(2x,e13.7))')tindex2(n),pfaf2(n),gnu,    &
!             tsa1(n),tsa2(n),tsb1(n),tsb2(n)
          else

         dist_save = 1000000.
         k = 0
         do i = 1,nbcatch
            if(i /= n) then
               if((ars1(i).ne.9999.).and.(arw1(i).ne.9999.)) then

                  tile_distance = (tile_lon(i) - tile_lon(n)) * (tile_lon(i) - tile_lon(n)) + &
                                  (tile_lat(i) - tile_lat(n)) * (tile_lat(i) - tile_lat(n))
                  if(tile_distance < dist_save) then
                     k = i
                     dist_save = tile_distance
                  endif
               endif
            endif
         enddo
         write (41,*)n,k
         write(20,'(i8,i8,f5.2,11(2x,e14.7))')   &
              tindex2(n),pfaf2(n),gnu,   &
              ars1(k),ars2(k),ars3(k),                   &
              ara1(k),ara2(k),ara3(k),ara4(k),           &
              arw1(k),arw2(k),arw3(k),arw4(k) 
         write(30,'(i8,i8,f5.2,3(2x,e13.7))')tindex2(n),pfaf2(n),gnu,bf1(k),bf2(k),bf3(k)
         write(40,'(i8,i8,f5.2,4(2x,e13.7))')tindex2(n),pfaf2(n),gnu,    &
              tsa1(k),tsa2(k),tsb1(k),tsb2(k)
         write(42,'(i8,i8,i4,i4,3f8.4,f12.8,f7.4,f10.3)') tindex2(n),pfaf2(n),soil_class_top(k),soil_class_com(k),    &
              BEE(k), PSIS(k),POROS(k),COND(k),WPWET(k),soildepth(k)  

      endif
   endif
   

        if (error_file) then
           write(21,*)tindex2(n),pfaf2(n),taberr1(n),taberr2(n),taberr3(n),taberr4(n), &
		normerr1(n),normerr2(n),normerr3(n),normerr4(n)
           write(31,*)tindex2(n),pfaf2(n),taberr1(n),taberr2(n),normerr1(n),normerr2(n)
        endif
        
     END DO
      
!      Write(*,*) 'END COMPUTING MODEL PARA'

      close(10,status='keep')
      close(20,status='keep')
      close(30,status='keep')
      close(40,status='keep')
!       close(11,status='delete')

      if (error_file) then
	 close(21,status='delete')
	 close(31,status='delete')
         close(41,status='keep')
      endif	 

  END SUBROUTINE create_model_para

!--------------------------------------------------------------------

  SUBROUTINE create_model_para_woesten (Maskfile)

    implicit none
      real, allocatable, dimension (:)  :: a_sand,a_clay,a_silt,a_oc,  &
          atile_sand,atile_clay, tile_lon, tile_lat, grav_vec, soc_vec,&
	  poc_vec,a_sand_surf,a_clay_surf,wpwet_surf,poros_surf

      real, allocatable, dimension (:,:) :: good_clay, good_sand
      integer, allocatable, dimension (:,:) :: tile_add, tile_pick
      type (mineral_perc) :: min_percs
      integer :: CF1, CF2, CF3, CF4
      integer i,j,n,k, tindex1,pfaf1,nbcatch
      integer soil_gswp
      real meanlu,stdev,minlu,maxlu,coesk,rzdep
      real minlat,maxlat,minlon,maxlon
      real,allocatable, dimension (:) ::   &
             BEE, PSIS,POROS,COND,WPWET,soildepth
      REAL, allocatable, dimension(:) :: TOPMEAN, TOPVAR, TOPSKEW
      REAL ST(NAR), AC(NAR),COESKEW
      REAL, allocatable, dimension (:) ::   &
	    ARS1,ARS2,ARS3, ARA1,ARA2,ARA3,ARA4, &
	    ARW1,ARW2,ARW3,ARW4,bf1, bf2, bf3,   &
	    tsa1, tsa2,tsb1, tsb2,               &
	    taberr1,taberr2,normerr1,normerr2,   &
	    taberr3,taberr4,normerr3,normerr4

      integer, allocatable, dimension (:) :: soil_class_com,tindex2,pfaf2, &
                                             soil_class_top
      real watdep(nwt,nrz),wan(nwt,nrz),rzexcn(nwt,nrz),frc(nwt,nrz)
      real, allocatable, dimension  (:,:,:) :: &
      gwatdep,gwan,grzexcn,gfrc
      real :: wtdep,wanom,rzaact,fracl,profdep,dist_save,     &
             ncells_top, ncells_top_pro,ncells_sub_pro,tile_distance
      character*100 :: pathout,fname,fout,losfile
      character*10 :: dline
      CHARACTER*20 :: version,resoln,continent
      character*6 rdep,ext
      character (*) :: MaskFile
      integer :: iwt,irz,group
      logical :: picked

 
! --------- VARIABLES FOR *OPENMP* PARALLEL ENVIRONMENT ------------
!
! NOTE: "!$" is for conditional compilation
!
logical :: running_omp = .false.
!
!$ integer :: omp_get_thread_num, omp_get_num_threads
!
integer :: n_threads=1, li, ui
!
integer, dimension(:), allocatable :: low_ind, upp_ind
!
! ------------------------------------------------------------------
        
  ! ----------- OpenMP PARALLEL ENVIRONMENT ----------------------------
  !
  ! FIND OUT WHETHER -omp FLAG HAS BEEN SET DURING COMPILATION
  !
  !$ running_omp = .true.         ! conditional compilation
  !
  ! ECHO BASIC OMP VARIABLES
  !
  !$OMP PARALLEL DEFAULT(NONE) SHARED(running_omp,n_threads) 
  !
  !$OMP SINGLE
  !
  !$ n_threads = omp_get_num_threads()
  !
  !$ write (*,*) 'running_omp = ', running_omp
  !$ write (*,*)
  !$ write (*,*) 'parallel OpenMP with ', n_threads, 'threads'
  !$ write (*,*)
  !$OMP ENDSINGLE
  !
  !$OMP CRITICAL
  !$ write (*,*) 'thread ', omp_get_thread_num(), ' alive'
  !$OMP ENDCRITICAL
  !
  !$OMP BARRIER
  !
  !$OMP ENDPARALLEL
      
!c-------------------------------------------------------------------------
      fname = trim(c_data)//'SoilClasses-SoilHyd-TauParam.dat' 
      open (11, file=trim(fname), form='formatted',status='old', &
           action = 'read')
      read (11,'(a)')fout           
      losfile =trim(c_data)//'/Woesten_SoilParam/loss_pd_top/loss_perday_rz1m_'

      allocate (a_sand (1:n_SoilClasses))
      allocate (a_silt (1:n_SoilClasses))
      allocate (a_clay (1:n_SoilClasses))
      allocate (a_oc   (1:n_SoilClasses))
      allocate (gwatdep(1:nwt,1:nrz,1:n_SoilClasses))
      allocate (gwan   (1:nwt,1:nrz,1:n_SoilClasses))
      allocate (grzexcn(1:nwt,1:nrz,1:n_SoilClasses))
      allocate (gfrc   (1:nwt,1:nrz,1:n_SoilClasses))
      
      do n =1,n_SoilClasses
         read (11,'(4f7.3)')a_sand(n),a_clay(n),a_silt(n),a_oc(n)
      write (fout,'(i2.2,i2.2,i4.4)')nint(a_sand(n)),nint(a_clay(n)),nint(100*a_oc(n))
	   open (120,file=trim(losfile)//trim(fout),  &
           form='formatted',status='old')

	   do iwt=1,nwt
	      do irz=1,nrz
		 read(120,2000) wtdep,wanom,rzaact,fracl
 2000		 format(1x,4e16.8)
		 gwatdep(iwt,irz,n)= wtdep
		 gwan(iwt,irz,n)   = wanom
		 grzexcn(iwt,irz,n)= rzaact
		 gfrc(iwt,irz,n)   = amin1(fracl,1.)
	      enddo
	  enddo
          close (120,status='keep')	   
      end do  
      close (11,status='keep')  
      deallocate (a_sand,a_silt,a_clay,a_oc)

     fname='clsm/soil_param.first'       
       open (10,file=fname,action='read',       &
          form='formatted',status='old')              
                                                 
     fname='clsm/cti_stats.dat'           
      open (11,file=fname,action='read',        &
          form='formatted',status='old')  
            
     fname='clsm/catchment.def'           
      open (12,file=fname,action='read',        &
          form='formatted',status='old')                                                             

     fout='clsm/ar.new'               
      open (20,file=fout,action='write',        &
          form='formatted',status='unknown')          
                                                 
     fout='clsm//bf.dat'               
      open (30,file=fout,action='write',        &
          form='formatted',status='unknown')          
                                                 
     fout='clsm//ts.dat'               
      open (40,file=fout,action='write',        &
          form='formatted',status='unknown')        

     if (error_file) then 
	fout='clsm/ar_rmse.dat'           
	open (21,file=fout,action='write',        &
	     form='formatted',status='unknown')

	fout='clsm/bf_rmse.dat'           
	open (31,file=fout,action='write',        &
	     form='formatted',status='unknown')

        fout='clsm/bad_sat_param.tiles'
 	open (41,file=fout,action='write',        &
	     form='formatted',status='unknown')  
 
     endif 
       fout='clsm/soil_param.dat'
 	open (42,file=fout,action='write',        &
	     form='formatted',status='unknown')       

     read (11,*)nbcatch
     read (12,*)nbcatch

     allocate (tile_lon(1:nbcatch))
     allocate (tile_lat(1:nbcatch))
     allocate (TOPMEAN (1:nbcatch))
     allocate (TOPVAR  (1:nbcatch))
     allocate (TOPSKEW (1:nbcatch))
     allocate (ARS1 (1:nbcatch))
     allocate (ARS2 (1:nbcatch))
     allocate (ARS3 (1:nbcatch))
     allocate (ARA1 (1:nbcatch))
     allocate (ARA2 (1:nbcatch))
     allocate (ARA3 (1:nbcatch))
     allocate (ARA4 (1:nbcatch))
     allocate (ARW1 (1:nbcatch))
     allocate (ARW2 (1:nbcatch))
     allocate (ARW3 (1:nbcatch))
     allocate (ARW4 (1:nbcatch))
     allocate (BF1 (1:nbcatch))
     allocate (BF2 (1:nbcatch))
     allocate (BF3 (1:nbcatch))
     allocate (TSA1 (1:nbcatch))
     allocate (TSA2 (1:nbcatch))
     allocate (TSB1 (1:nbcatch))
     allocate (TSB2 (1:nbcatch))
     allocate (TABERR1 (1:nbcatch))
     allocate (TABERR2 (1:nbcatch))
     allocate (TABERR3 (1:nbcatch))
     allocate (TABERR4 (1:nbcatch))
     allocate (NORMERR1 (1:nbcatch))
     allocate (NORMERR2 (1:nbcatch))
     allocate (NORMERR3 (1:nbcatch))
     allocate (NORMERR4 (1:nbcatch))
     allocate (BEE        (1:nbcatch))
     allocate (PSIS      (1:nbcatch))
     allocate (POROS     (1:nbcatch))
     allocate (COND      (1:nbcatch))
     allocate (WPWET     (1:nbcatch))
     allocate (soildepth (1:nbcatch))
     allocate (soil_class_top (1:nbcatch))
     allocate (soil_class_com (1:nbcatch))
     allocate (tindex2        (1:nbcatch))
     allocate (pfaf2          (1:nbcatch))
     allocate (atile_clay     (1:nbcatch))
     allocate (atile_sand     (1:nbcatch))
     allocate (grav_vec       (1:nbcatch))
     allocate (soc_vec        (1:nbcatch))
     allocate (poc_vec        (1:nbcatch))
     allocate (a_sand_surf    (1:nbcatch))
     allocate (a_clay_surf    (1:nbcatch))
     allocate (wpwet_surf     (1:nbcatch))
     allocate (poros_surf     (1:nbcatch))
     allocate (good_clay     (1:100,4))
     allocate (good_sand     (1:100,4))
     allocate (tile_add      (1:100,4))   
     allocate (tile_pick     (1:100,4)) 
     tile_add = 0
     tile_pick= 0
     good_clay =0.
     good_sand =0.
 
     do n=1,nbcatch
        read(11,'(i8,i8,5(1x,f8.4))') tindex1,pfaf1,meanlu,stdev  &
             ,minlu,maxlu,coesk                                  
 
        read(10,'(i8,i8,i4,i4,3f8.4,f12.8,f7.4,f10.4,3f7.3,4f7.3,2f10.4)') &
	     tindex2(n),pfaf2(n),soil_class_top(n),soil_class_com(n),      &
             BEE(n), PSIS(n),POROS(n),COND(n),WPWET(n),soildepth(n),       &
	     grav_vec(n),soc_vec(n),poc_vec(n),                            &
	     a_sand_surf(n),a_clay_surf(n),atile_sand(n),atile_clay(n)                     
        if(tindex1.ne.tindex2(n))then
           write(*,*)'Warnning 1: tindex mismatched'                        
           stop
        endif

        read (12,*) tindex1,pfaf1,minlon,maxlon,minlat,maxlat
        tile_lon(n) = (minlon + maxlon)/2.
        tile_lat(n) = (minlat + maxlat)/2.

        if(pfaf1.ne.pfaf2(n)) then
           write(*,*)'Warnning 1: pfafstetter mismatched' 
           stop
        endif

        if (index(MaskFile,'GEOS5_10arcsec_mask') /= 0) then
           TOPMEAN(n) = meanlu
        else
           TOPMEAN(n) = 0.961*meanlu-1.957                       
        endif

        TOPVAR(n)  = stdev*stdev                                
        TOPSKEW(n) = coesk*stdev*stdev*stdev                   
        
        if ( TOPVAR(n) .eq. 0. .or. coesk .eq. 0.            &
             .or. topskew(n) .eq. 0.) then                       
           write(*,*) 'Problem: undefined values:'         
           write(*,*) TOPMEAN(n),TOPVAR(n),coesk,            &
                minlu,maxlu
           stop
        endif
     END DO

     rewind(10)

     allocate(low_ind(n_threads))
     allocate(upp_ind(n_threads))
     low_ind(1)         = 1
     upp_ind(n_threads) = nbcatch

     if (running_omp)  then
      do i=1,n_threads-1
        
         upp_ind(i)   = low_ind(i) + (nbcatch/n_threads) - 1 
         low_ind(i+1) = upp_ind(i) + 1
        
      end do 
     end if


!$OMP PARALLELDO DEFAULT(NONE)                          &
!$OMP SHARED( BEE, PSIS,POROS,COND,WPWET,soildepth,     &
!$OMP        TOPMEAN, TOPVAR, TOPSKEW,                  &
!$OMP        ARS1,ARS2,ARS3, ARA1,ARA2,ARA3,ARA4,       &
!$OMP        ARW1,ARW2,ARW3,ARW4,bf1, bf2, bf3,         &
!$OMP        tsa1, tsa2,tsb1, tsb2,                     &
!$OMP        taberr1,taberr2,normerr1,normerr2,         &
!$OMP        taberr3,taberr4,normerr3,normerr4,         &
!$OMP        gwatdep,gwan,grzexcn,gfrc,soil_class_com,  &
!$OMP        n_threads, low_ind, upp_ind )              &
!$OMP PRIVATE(k,li,ui,n,i,watdep,wan,rzexcn,frc,ST,AC,  &
!$OMP COESKEW,profdep)

     do k=1,n_threads

        li = low_ind(k)
        ui = upp_ind(k)

      do n=li,ui

      CALL TGEN (                              &
          TOPMEAN(n),TOPVAR(n),TOPSKEW(n),     &
          ST,AC,COESKEW)
      
!c Areal fractioning parameters

      CALL SAT_PARAM(                                              &
                     BEE(n),PSIS(n),POROS(n),COND(n),              & 
                     WPWET(n), ST, AC, COESKEW,n,                  &
                     soildepth(n),                                 &
                     ars1(n),ars2(n),ars3(n),                      &
                     ara1(n),ara2(n),ara3(n),ara4(n),              &
                     arw1(n),arw2(n),arw3(n),arw4(n),              &
                     taberr1(n),taberr2(n),taberr3(n),taberr4(n),  &
                     normerr1(n),normerr2(n),normerr3(n),normerr4(n))


      CALL BASE_PARAM(                                  &
          BEE(n),PSIS(n),POROS(n),COND(n),              &
          ST, AC,                                       &
          bf1(n),bf2(n),bf3(n),                         &
          taberr1(n),taberr2(n),normerr1(n),normerr2(n) &
          )
            

	  watdep (:,:) =  gwatdep (:,:,soil_class_com(n))   
	  wan    (:,:) =  gwan    (:,:,soil_class_com(n))
	  rzexcn (:,:) =  grzexcn (:,:,soil_class_com(n))
	  frc    (:,:) =  gfrc    (:,:,soil_class_com(n))

      CALL TS_PARAM(                       &
          BEE(n),PSIS(n),POROS(n),         &
          ST, AC,                          &
          watdep,wan,rzexcn,frc,           &
          tsa1(n),tsa2(n),tsb1(n),tsb2(n)  &
          )

     END DO
     END DO
          !$OMP ENDPARALLELDO

     CF1 =0
     CF2 =0
     CF3 =0
     CF4 =0

     DO n=1,nbcatch

     	if((ars1(n).ne.9999.).and.(arw1(n).ne.9999.))then

	   if ((soil_class_com(n)>=1).and.(soil_class_com(n)<=84)) then	
      	      group=1
           else if ((soil_class_com(n) > 84).and.(soil_class_com(n)<=168)) then
      	      group=2
           else if ((soil_class_com(n) >168).and.(soil_class_com(n)< N_SoilClasses)) then
              group=3
	   else
              group=4
           endif

           min_percs%clay_perc = atile_clay(n)
           min_percs%sand_perc = atile_sand(n) 
           min_percs%silt_perc = 100. - min_percs%clay_perc - min_percs%sand_perc

           if(tile_pick(soil_class (min_percs),group) == 0) then
             tile_pick(soil_class (min_percs),group) = n 

             select case (group)

             case (1)
                 
                CF1 = CF1 + 1
                good_clay (CF1,group) = atile_clay(n)
                good_sand (CF1,group) = atile_sand(n)
                tile_add  (CF1,group) = n

             case (2)
                CF2 = CF2 + 1
                good_clay (CF2,group) = atile_clay(n)
                good_sand (CF2,group) = atile_sand(n)
                tile_add  (CF2,group) = n

             case (3)
                CF3 = CF3 + 1
                good_clay (CF3,group) = atile_clay(n)
                good_sand (CF3,group) = atile_sand(n)
                tile_add  (CF3,group) = n

             case (4)
                CF4 = CF4 + 1
                good_clay (CF4,group) = atile_clay(n)
                good_sand (CF4,group) = atile_sand(n)
                tile_add  (CF4,group) = n

             end select
          endif
	endif	
     END DO	

     DO n=1,nbcatch
        read(10,'(i8,i8,i4,i4,3f8.4,f12.8,f7.4,f10.4,3f7.3,4f7.3,2f10.4)') &
	     tindex2(n),pfaf2(n),soil_class_top(n),soil_class_com(n),         &
             BEE(n), PSIS(n),POROS(n),COND(n),WPWET(n),soildepth(n),       &
	     grav_vec(n),soc_vec(n),poc_vec(n),                            &
	     a_sand_surf(n),a_clay_surf(n),atile_sand(n),atile_clay(n) ,   &
	     wpwet_surf(n),poros_surf(n)
     if((ars1(n).ne.9999.).and.(arw1(n).ne.9999.))then   
      write(20,'(i8,i8,f5.2,11(2x,e14.7))')         &
                     tindex2(n),pfaf2(n),gnu,       &
                     ars1(n),ars2(n),ars3(n),         &
                     ara1(n),ara2(n),ara3(n),ara4(n), &
                     arw1(n),arw2(n),arw3(n),arw4(n) 

      write(30,'(i8,i8,f5.2,3(2x,e13.7))')tindex2(n),pfaf2(n),gnu,bf1(n),bf2(n),bf3(n)
      write(40,'(i8,i8,f5.2,4(2x,e13.7))')tindex2(n),pfaf2(n),gnu,    &
          tsa1(n),tsa2(n),tsb1(n),tsb2(n)

      write(42,'(i8,i8,i4,i4,3f8.4,f12.8,f7.4,f10.4,3f7.3,4f7.3,2f10.4)')  &
	     tindex2(n),pfaf2(n),soil_class_top(n),soil_class_com(n),      &
             BEE(n), PSIS(n),POROS(n),COND(n),WPWET(n),soildepth(n),       &
	     grav_vec(n),soc_vec(n),poc_vec(n),                            &
	     a_sand_surf(n),a_clay_surf(n),atile_sand(n),atile_clay(n),    &
	     wpwet_surf(n),poros_surf(n)

      else
      if(preserve_soiltype) then 
         if ((soil_class_com(n)>=1).and.(soil_class_com(n)<=84)) then	
            group=1
         else if ((soil_class_com(n)>  84).and.(soil_class_com(n)<=168)) then
            group=2
         else if ((soil_class_com(n)> 168).and.(soil_class_com(n)< N_SoilClasses)) then
            group=3
         else
            group=4
         endif
         
         min_percs%clay_perc = atile_clay(n)
         min_percs%sand_perc = atile_sand(n) 
         min_percs%silt_perc = 100. - min_percs%clay_perc - min_percs%sand_perc
         if(tile_pick(soil_class (min_percs),group) > 0) then
            k = tile_pick(soil_class (min_percs),group) 
            
         else
            select case (group)
               
            case (1)
               j = center_pix (good_clay(1:CF1,group),good_sand(1:CF1,group),       &
                    min_percs%clay_perc,min_percs%sand_perc,min_percs%silt_perc,.true.)
               k = tile_add  (j,group) 
            case (2)
               j = center_pix (good_clay(1:CF2,group),good_sand(1:CF2,group),       &
                    min_percs%clay_perc,min_percs%sand_perc,min_percs%silt_perc,.true.)
               k = tile_add  (j,group)   
            case (3)
               j = center_pix (good_clay(1:CF3,group),good_sand(1:CF3,group),       &
                    min_percs%clay_perc,min_percs%sand_perc,min_percs%silt_perc,.true.)
               k = tile_add  (j,group) 
            case (4)
               j = center_pix (good_clay(1:CF4,group),good_sand(1:CF4,group),       &
                    min_percs%clay_perc,min_percs%sand_perc,min_percs%silt_perc,.true.)
               k = tile_add  (j,group)   
            end select
            print *,'NO Similar SoilClass :',soil_class (min_percs),group,n,k
            
         endif
         if (error_file) then
            write (41,*)n,k
            !           write (41,*)tindex2(n),pfaf2(n),soil_class_top(n),soil_class_com(n),    &
            !            BEE(n), PSIS(n),POROS(n),COND(n),WPWET(n),soildepth(n)  
            !           write (41,*)tindex2(k),pfaf2(k),soil_class_top,soil_class_com(k),    &
            !                BEE(k), PSIS(k),POROS(k),COND(k),WPWET(k),soildepth(k) 
         endif
         
         write(20,'(i8,i8,f5.2,11(2x,e14.7))')   &
              tindex2(n),pfaf2(n),gnu,   &
              ars1(k),ars2(k),ars3(k),                   &
              ara1(k),ara2(k),ara3(k),ara4(k),           &
              arw1(k),arw2(k),arw3(k),arw4(k) 
        write(30,'(i8,i8,f5.2,3(2x,e13.7))')tindex2(n),pfaf2(n),gnu,bf1(k),bf2(k),bf3(k)
        write(40,'(i8,i8,f5.2,4(2x,e13.7))')tindex2(n),pfaf2(n),gnu,    &
          tsa1(k),tsa2(k),tsb1(k),tsb2(k)

        write(42,'(i8,i8,i4,i4,3f8.4,f12.8,f7.4,f10.4,3f7.3,4f7.3,2f10.4)') &
              tindex2(n),pfaf2(n),soil_class_top(k),soil_class_com(k),      &
              BEE(k), PSIS(k),POROS(k),COND(k),WPWET(k),soildepth(k),       &
              grav_vec(k),soc_vec(k),poc_vec(k),                            &
              a_sand_surf(k),a_clay_surf(k),atile_sand(k),atile_clay(k) ,   &
	      wpwet_surf(k),poros_surf(k)
     else 

         dist_save = 1000000.
         k = 0
         do i = 1,nbcatch
            if(i /= n) then
               if((ars1(i).ne.9999.).and.(arw1(i).ne.9999.)) then

                  tile_distance = (tile_lon(i) - tile_lon(n)) * (tile_lon(i) - tile_lon(n)) + &
                                  (tile_lat(i) - tile_lat(n)) * (tile_lat(i) - tile_lat(n))
                  if(tile_distance < dist_save) then
                     k = i
                     dist_save = tile_distance
                  endif
               endif
            endif
         enddo
         write (41,*)n,k
         write(20,'(i8,i8,f5.2,11(2x,e14.7))')   &
              tindex2(n),pfaf2(n),gnu,   &
              ars1(k),ars2(k),ars3(k),                   &
              ara1(k),ara2(k),ara3(k),ara4(k),           &
              arw1(k),arw2(k),arw3(k),arw4(k) 
         write(30,'(i8,i8,f5.2,3(2x,e13.7))')tindex2(n),pfaf2(n),gnu,bf1(k),bf2(k),bf3(k)
         write(40,'(i8,i8,f5.2,4(2x,e13.7))')tindex2(n),pfaf2(n),gnu,    &
              tsa1(k),tsa2(k),tsb1(k),tsb2(k)
         write(42,'(i8,i8,i4,i4,3f8.4,f12.8,f7.4,f10.4,3f7.3,4f7.3,2f10.4)')&
              tindex2(n),pfaf2(n),soil_class_top(k),soil_class_com(k),         &
              BEE(k), PSIS(k),POROS(k),COND(k),WPWET(k),soildepth(k),       &
              grav_vec(k),soc_vec(k),poc_vec(k),                            &
              a_sand_surf(k),a_clay_surf(k),atile_sand(k),atile_clay(k) ,   &
	      wpwet_surf(k),poros_surf(k)

      endif
   endif
   
      if (error_file) then
	 write(21,*)tindex2(n),pfaf2(n),taberr1(n),taberr2(n),taberr3(n),taberr4(n), &
		normerr1(n),normerr2(n),normerr3(n),normerr4(n)
	 write(31,*)tindex2(n),pfaf2(n),taberr1(n),taberr2(n),normerr1(n),normerr2(n)
      endif 

      END DO
      
!      Write(*,*) 'END COMPUTING MODEL PARA'

      close(10,status='keep')
      close(11,status='keep')
      close(12,status='keep')
      close(20,status='keep')
      close(30,status='keep')
      close(40,status='keep')
      close(42,status='keep')


      if (error_file) then
	 close(21,status='delete')
	 close(31,status='delete')
         close(41,status='keep')
      endif	 

  END SUBROUTINE create_model_para_woesten


!---------------------------------------------------------------------

      SUBROUTINE TS_PARAM(                                &
                           BEE,PSIS,POROS,                &
                           VALX, PX,                      &
                           watdep,wan,rzexcn,frc,         &
                           tsa1,tsa2,tsb1,tsb2            &
                          )

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                                         c
!c Given pre-computed 1-D relationships between a "local" root zone excess c 
!c and a "local" catchment deficit, the timescale of the bulk vertical     c
!c transfer between the two bulk prognostic variables is computed using    c
!c the distribution of the local deficit established from the distribution c
!c of the topographic index, then an approximated function of catdef and   c
!c rzex is derived.                                                        c
!c                                                                         c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      INTEGER NAR0
      REAL, intent (in) :: BEE, PSIS, POROS
      REAL, intent (in) :: VALX(NAR), PX(NAR)
      real, intent (inout) :: watdep(nwt,nrz),wan(nwt,nrz),  &
                      rzexcn(nwt,nrz),frc(nwt,nrz)
      real, intent (out) ::  tsa1, tsa2 ,tsb1, tsb2

      integer :: tex,iwt,irz,n,idep,k, index1,i0
      REAL VALX0(NAR), PX0(NAR),sumta,sumta2,timean,zbar, rzw
      REAL :: term1, term2, sumdef, suma, frcsat,rzexc, rzact
      real zdep(nar),def(nar),wrz(nar),wbin(500),rze(nar)
      real catd(2,2),tsc(2,2), satfrc,sumfrac,sumz,frac
      real, parameter ::  frcmax = .041
      real  wtdep,wanom,rzaact,fracl,profdep,rzdep

!      logical bug

!c----------------------------------------------------------------
!c Is loss.dat compatible with rzdep = 0.49 ???

      rzdep = grzdep

!c Convert fractions to "per-hour" values
      do iwt=1,nwt
         do irz=1,nrz
            frc(iwt,irz)=1.-((1.-frc(iwt,irz))**(1./24.))
         enddo
      enddo

         nar0=0
         do n=1,nar
            if (px(n) .ne. 0.) then
               nar0=nar0+1
               valx0(nar0)=valx(n)
               px0(nar0)=px(n)
            endif
         enddo

         sumta=0.
         sumta2=0.
         suma=0.
         do n=1,nar0
            sumta=sumta+px0(n)*valx0(n)
            sumta2=sumta2+px0(n)*valx0(n)*valx0(n)
            suma=suma+px0(n)
         enddo

         timean=sumta/suma

!c**** Loop over two water table depths
         do idep=1,2
            if(idep.eq.1) zbar=1.5 ! zbar in meters
            if(idep.eq.2) zbar=2.0

!c**** Compute array of water table depths:
            do k=1,nar0
               term1=(1/gnu)*(valx0(k)-timean)
               zdep(k)=zbar-term1
               if(zdep(k) .lt. 0.) zdep(k)=0.
            enddo
!c            write(*,*)"  End water table depth"
!c**** Compute array of moisture deficits:
            do k=1,nar0
               term1=(psis-zdep(k))/psis
               term1=term1**(1.-1./bee)
               term2=-psis*(bee/(bee-1.))*(term1-1.)
               def(k)=poros*(zdep(k)-term2)
            enddo

!c**** Add deficits to produce catdef:
            sumdef=0.
            do k=1,nar0
               sumdef=sumdef+def(k)*px0(k)*1000.
            enddo
!c            write(*,*)"  End catchment deficit"
!c**** Compute array of root zone moisture (degree of wetness in root zone):
            do k=1,nar0

               if(zdep(k).eq.0.) then
                  wrz(k)=1.
               elseif(zdep(k)-rzdep.lt.0.) then
                  term1=((psis-zdep(k))/psis)**(1.-1./bee)
                  wrz(k)=(-psis/zdep(k))*(bee/(bee-1.))   &
                      *(term1-1.)
                  frcsat=1.-zdep(k)/rzdep
                  wrz(k)=(1.-frcsat)*wrz(k)+frcsat*1.
               else
                  term1=((psis-zdep(k))/psis)**(1.-1./bee)
                  term2=((psis-zdep(k)+rzdep)/psis)    &
                      **(1.-1./bee)
                  wrz(k)=(-psis/rzdep)*(bee/(bee-1.))  &
                      *(term1-term2)
               endif
            enddo

!c       Loop over two root zone excess values:
            do irz=1,2
               if(irz.eq.1) rzexc=-0.1*poros
               if(irz.eq.2) rzexc=0.1*poros

!c       Determine actual root zone excess
               rzact=0.
               do k=1,nar0
                  rze(k)=rzexc
                  rzw=wrz(k)*poros
                  if(rzw+rze(k) .gt. poros) rze(k)=poros-rzw
                  if(rzw+rze(k) .lt. 0.) rze(k)=rzw
                  rzact=rzact+rze(k)*px0(k)
               enddo
!c            write(*,*)"  End root zone excess"
!c       Compute the average timescale

               satfrc=0.
               do k=1,nar0
                  if(zdep(k).lt.0.) satfrc=satfrc+px0(k)
               enddo

               sumfrac=0.
               sumz=0.
               do k=1,nar0
                  sumz=sumz+zdep(k)*px0(k)
                  if(zdep(k) .lt. 1.) frac=frcmax
                  if(zdep(k) .ge. 1.) then
                     index1=1+int(((zdep(k)*100.)-99)/5.)
                     if(index1.gt.nwt) index1 = nwt
                     frac=amin1(frc(index1,1),frcmax)
                     do i0=2,nrz
                        if(rze(k) .ge. rzexcn(index1,i0))  &
                            frac=amin1(frc(index1,i0),frcmax)
                     enddo
                  endif
                  sumfrac=sumfrac+frac*px0(k)
               enddo
!c            write(*,*)"  End average time scale"
               catd(idep,irz)=sumdef
               tsc(idep,irz)=sumfrac

            enddo
         enddo

         tsb1=(alog(tsc(2,2))-alog(tsc(1,2)))/(catd(2,2)-catd(1,2))
         tsb2=(alog(tsc(2,1))-alog(tsc(1,1)))/(catd(2,1)-catd(1,1))
         tsa1=alog(tsc(2,2))-tsb1*catd(2,2)
         tsa2=alog(tsc(2,1))-tsb2*catd(2,1)

       END SUBROUTINE TS_PARAM

!*********************************************************************

      SUBROUTINE BASE_PARAM(                                  &
                           BEE,PSIS,POROS,COND,               &
                           VALX, PX,                          &
                           bf1,bf2,bf3,                       &
                           taberr1,taberr2,normerr1,normerr2  &
                          )

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                                   c
!c New way to get baseflow: we parametrize the relationship between  c 
!c catdef and zbar (two parameters bf1 and bf2).                     c
!c Then, in the LSM/catchment.f/base.f, we use the original relation c 
!c from TOPMODEL to infer baseflow from catdef and the mean of the   c
!c topographic index (topmean=bf3, a third parameter).               c
!c                                                                   c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      INTEGER IDMAX,i1,i2,i,icount

      REAL, intent (in) :: BEE, PSIS,POROS,COND,VALX(NAR),PX(NAR)
      real zbar(nbdep),catdef(nbdep),bflow(nbdep)
      real, intent (out) :: bf1,bf2,bf3,taberr1,taberr2,normerr1,normerr2
      integer :: n,idep
      real suma,sumta,timean

      real catfit(nbdep),bfit(nbdep),dfit(nbdep),catmean,bfmean
      real catref(nbdep),bref(nbdep)
      real err1, err2
!      logical, intent (in) :: bug

         sumta=0.
         suma=0.
         do n=1,nar
            sumta=sumta+px(n)*valx(n)
            suma=suma+px(n)
         enddo
         timean=sumta/suma
         bf3 = timean

!c**** Loop over water table depths

         do idep=1,nbdep

!c           write(*,*) 'idep=',idep

            CALL BASIDEP(                  &
                IDEP,                      &
                BEE,PSIS,POROS,COND,       &
                VALX,PX,TIMEAN,SUMA,       &
                ZBAR,CATDEF,BFLOW)

         enddo


         i1=10   ! zbar= 0 m
         i2=35   ! zbar= 2.5 m

         bf2=zbar(i2)*SQRT(catdef(i1))               &
             /(SQRT(catdef(i2))-SQRT(catdef(i1)))
         bf1=catdef(i1)/(bf2*bf2)

         if (bf1 .le. 0) write(*,*) 'bf1 le 0 for i=',i
         if (bf2 .le. 0) write(*,*) 'bf2 le 0 for i=',i

!c Errors: Root mean square errors: only for points where catdef GT 0.5mm

         do idep=1,nbdep
            catref(idep)=0.
            bref(idep)=0.
         enddo
         catmean=0.
         bfmean=0.
         icount=0
         do idep=1,nbdep
            if (catdef(idep) .gt. lim) then
               icount=icount+1
               catref(icount)=catdef(idep)
               bref(icount)=bflow(idep)
               catfit(icount)=bf1*(zbar(idep)+bf2)             &
                   *(zbar(idep)+bf2)
               dfit(icount)=SQRT(catdef(idep)/bf1)-bf2
               bfit(icount)=cond*exp(-timean-gnu*dfit(icount)) &
                   /gnu
               catmean=catmean+catdef(idep)
               bfmean=bfmean+bflow(idep)
            endif
         enddo
         catmean=catmean/icount
         bfmean=bfmean/icount
	 if (icount.gt.1) then
         call RMSE(catref,catfit,icount,err1)
         call RMSE(bref,bfit,icount,err2)

         taberr1=err1
         taberr2=err2
         normerr1=err1/catmean
         normerr2=err2/bfmean
	 endif
!c---------------------------------------------------------------------
         
       END SUBROUTINE BASE_PARAM

! ************************************************************************

      SUBROUTINE BASIDEP(                          &
                         IDEP,                     &
                         BEE,PSIS,POROS,COND,      &
                         VALX,PX,TIMEAN,SUMA,      &
                         ZBAR,CATDEF,BFLOW)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                                      c
!c This program returns the eight parameters for the areal fractioning  c
!c                                                                      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
     implicit none
      INTEGER, intent (in) :: idep
      integer nref, nind,nmax,indmin,locmax,shift,ord,locmin,ordref,width,k
      REAL, intent (in) :: BEE, PSIS, POROS, COND,VALX(NAR), PX(NAR), &
           suma,timean
      real :: dx,sumdef,dz
      real, intent (out) ::  catdef(nbdep),bflow(nbdep),zbar(idep)
      real term1,term2,sum
      real zdep(nar),locdef(nar)
!      logical bug

!c-------------------------------------------------------------------------
!c integral(f(x)dx)=1. for a pdf
!c here px=f(x)dx

      dx=valx(1)-valx(2)

      if (bug) write(*,*) 'IDEP=',IDEP,' dx=',dx, 'gnu=',gnu

!c the loops over idmax and nbdep are initiated in sta_params4.f

      zbar(idep)=float(idep-10)*slice ! zdep in meters
            
!c**** Compute array of water table depths:
      do k=1,nar
         term1=(1/gnu)*(valx(k)-timean)
         zdep(k)=AMAX1(0.,zbar(idep)-term1)
      enddo
            
!c variable change must be reflected in dx
      dz=dx/gnu

      if (bug) write(*,*) 'basidep: ok1'
      
!c**** Compute array of moisture deficits:
      do k=1,nar
         term1=(psis-zdep(k))/psis
         term1=term1**(1.-1./bee)
         term2=-psis*(bee/(bee-1.))*(term1-1.)
         locdef(k)=zdep(k)-term2
      enddo

!c**** Add deficits to produce catdef:
      sumdef=0.
      do k=1,nar
         sumdef=sumdef+locdef(k)*px(k)
      enddo
      catdef(idep)=poros*1000.*sumdef/suma

      if (bug) write(*,*) 'basidep: ok2'

      bflow(idep)=cond*exp(-timean-gnu*zbar(idep))/gnu

      if (bug) write(*,*) 'basidep: ok3'

    END SUBROUTINE BASIDEP

!*****************************************************************************

      SUBROUTINE SAT_PARAM(                                                   &
                           BEE,PSIS,POROS,COND,                               &
                           WPWET,VALX, PX, COESKEW,PFC,                       &
                           soildepth,                                         &
                           ARS1,ARS2,ARS3,                                    &
                           ARA1,ARA2,ARA3,ARA4,                               &
                           ARW1,ARW2,ARW3,ARW4,                               &
                           taberr1,taberr2,taberr3,taberr4,                   &
                           normerr1,normerr2,normerr3,normerr4,               &
                           DBG_UNIT)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                                      c
!c This program returns the eleven parameters for the areal fractioning c
!c                                                                      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      IMPLICIT NONE

      INTEGER, intent (in) :: pfc
      REAL, intent (in) :: BEE,PSIS,POROS,COND,WPWET, &
           VALX(NAR), PX(NAR)
      REAL, intent (in) :: soildepth, COESKEW
      REAL, intent (inout) :: ARS1,ARS2,ARS3,                                 &
                           ARA1,ARA2,ARA3,ARA4,                               &
                           ARW1,ARW2,ARW3,ARW4,                               &
                           taberr1,taberr2,taberr3,taberr4,                   &
                           normerr1,normerr2,normerr3,normerr4
      INTEGER idep,n,k,i,icount,iref
      integer nar0
      integer nref, nind,nmax,indmin,locmax,shift,ord,locmin
      integer loc1,loc2,loc3,loc0,flag
      REAL VALX0(NAR), PX0(NAR)
      integer :: adjust,loc2save,inc,dec
      real sumta,suma,timean,upval,loval,profdep
      real rjunk,rjunk2
      integer, intent (in), optional :: DBG_UNIT
      real catdef(nbdep),wmin(nbdep),ar1(nbdep),aa(nbdep),aabis(nbdep)
      real ar2(nbdep),ar3(nbdep),swsrf2(nbdep),swsrf3(nbdep),rzeq(nbdep)
      real zbar0,catdef0,wmin0,RZDEP,wminsave(nbdep)

      real x1,x2,x3,x4,w1,w1_0,w2,w3,w4,ref1
      real y0,f1,f2,f3,g1,g2,g3,df,dg,dx,bf,bg,delta,z1,z2

      real nar1(nbdep),nar2(nbdep),nmean2(nbdep),neq(nbdep)
      real shape, nwm, area1,cdi,nar3(nbdep),nmean3
      real err1,err2,err3,err4,sum
      real tabact(nbdep),tabfit(nbdep)

      integer :: mp,isvd,j,first_loop
!      REAL*8, allocatable :: A(:,:),AP(:,:)
!      REAL*8, allocatable :: B(:)
      REAL*8, allocatable, target :: A(:,:)
      REAL*8, allocatable, target :: B(:)
      REAL*8, pointer             :: AP(:,:)
      REAL*8, pointer             :: BP(:)
      REAL*8 V(3,3),W(3),ANS(3),sdmax,sdmin,wbrac

      real :: cdcr1,cdcr2,term1,term2,zmet
      logical :: smooth,ars_svd_loop
      logical, parameter ::  bug=.false.
      logical, parameter :: SingValDecomp = .true.
      integer, parameter :: nl=4, nr=4, m=4, NP=50
      real :: savgol_coeff(NP)  
      integer :: savgol_ind(NP)
      integer :: nbdepl,istart

      ref1 = 100.
!      print *,'PFC', pfc   
      if (bug) write(*,*) 'starting sat_param'

      if(SingValDecomp) then
           savgol_ind(1)=0 
	   j=3
           do i=2, nl+1
              savgol_ind(i)=i-j
	      j=j+2
           end do

           j=2
           do i=nl+2, nl+nr+1
              savgol_ind(i)=i-j
	      j=j+2
           end do   
         call savgol(savgol_coeff,nl+nr+1,nl,nr,0,m)
      endif

      profdep = soildepth
      rzdep =grzdep
      profdep=profdep/1000.
      profdep=amax1(1.,profdep)
      if (rzdep .gt. .75*profdep) then
        rzdep=0.75*profdep
      end if
      
      zmet=profdep
      term1=-1.+((psis-zmet)/psis)**  &
           ((bee-1.)/bee)
      term2=psis*bee/(bee-1)
      cdcr1=1000.*poros*(zmet-(-term2*term1))
      cdcr2=(1-wpwet)*poros*1000.*zmet
!c mean of the topographic index distribution

      nar0=0
      do n=1,nar
         if (px(n) .ne. 0.) then
            nar0=nar0+1
            valx0(nar0)=valx(n)
            px0(nar0)=px(n)
         endif
      enddo

      sumta=0.
      suma=0.
      do n=1,nar0
         sumta=sumta+px0(n)*valx0(n)
         suma=suma+px0(n)
      enddo
      timean=sumta/suma

      if (bug) write(*,*) 'ok 0: sumta,suma,nar0=',sumta,suma,nar0

!c**** Loop over water table depths

         do idep=1,nbdep

            CALL FUNCIDEP(                                       &
                         NAR0,IDEP,                              &
                         BEE,PSIS,POROS,COND,RZDEP,WPWET,        &
                         VALX0,PX0,COESKEW,TIMEAN,SUMA,          &
                         CATDEF,AR1,WMIN,AA,AABIS,               &
                         AR2,AR3,SWSRF2,SWSRF3,RZEQ)             
         enddo

         nbdepl = 100
         if(catdef(50) > cdcr1 + 20.) nbdepl = 50
         if(soildepth > 6500.)  nbdepl = nbdep

         if (bug) write(*,*) 'funcidep loop ok'

!c**** for wmin's adjustment, we need an estimate of its limit toward INF
         adjust =0
         ZBAR0=10.
         CALL FUNCZBAR(                                        &
                         NAR0,ZBAR0,                           &
                         BEE,PSIS,POROS,COND,RZDEP,WPWET,      &
                         VALX0,PX0,COESKEW,TIMEAN,SUMA,        &
                         CATDEF0,WMIN0)

         if (bug) write(*,*) 'funczbar ok'

	 if (wmin0 == 0.9999900) then
	       do idep=1,nbdep-1
	         if(catdef(idep).le.cdcr1+10.) then
		 if((wmin(idep) - wmin(idep +1)) > -0.01) then 
	           wmin0=wmin(idep)
                 endif
		 endif  
	       enddo
              wmin0 = 0.1*(nint(wmin0*100000.)/10000) -0.02		    
	 endif
        
       if(present(dbg_unit)) then
          write (dbg_unit,*) nbdep,nbdepl,wmin0,cdcr1,cdcr2
          write (dbg_unit,*) catdef
          write (dbg_unit,*) ar1
          write (dbg_unit,*) wmin
       endif

!c**** AR1 adjustment: 3 points + limit in INF = 0.

         if (bug) write(*,*) 'STARTING AR1'

         ! Singular value decomposition               
         loc1=1
         loc3=nbdepl
	 loc2=loc3	 

	 do idep = 1,loc2
	    if(ar1(idep) < 1.e-10) then	
	       loc3 = idep - 1
	       exit
	    endif
	 end do	 
  
         first_loop = 0      
         ars_svd_loop = .TRUE.
         DO while (ars_svd_loop)
         
         first_loop = first_loop + 1
         mp = loc3-loc1+1
           
         allocate(A(mp,3))
         allocate(AP(mp,3))
         allocate(B(mp))
              
         a=0.
         ap=0.
         b=0.
         v=0.
         w=0.
         ans=0.
              
         do isvd=loc1,loc3
            A(isvd-loc1+1,1)=catdef(isvd)
            A(isvd-loc1+1,2)=-catdef(isvd)*ar1(isvd)
            A(isvd-loc1+1,3)=-ar1(isvd)*((catdef(isvd))**2.)
            B(isvd-loc1+1)=ar1(isvd)-1.
         end do
              
         ap = a
         call svdcmp(ap,mp,3,w,v)
         sdmax=0.
         do j=1,3
            if(w(j).gt.sdmax)sdmax=w(j)
         end do
         sdmin=sdmax*1.0e-6
         do j=1,3
            if(w(j).lt.sdmin)w(j)=0.
         end do

         call svbksb(ap,w,v,mp,3,b,ans)
              
         ars1 = real(ans(1))
         ars2 = real(ans(2))
         ars3 = real(ans(3))  

         flag=0
         call curve1(ars1,ars2,ars3,cdcr2,flag)
         deallocate (A, AP, B)

         IF(FLAG == 1) THEN
            LOC3 = NBDEP
            LOC1 =1
            IF(first_loop > 1) ars_svd_loop=.FALSE.
            ELSE
            ars_svd_loop=.FALSE.  
         ENDIF
         END DO

         IF (FLAG.EQ.1) then
  
         flag=0
	 loc1=1
         do idep=1,nbdepl
            if (catdef(idep) .le. 20.) loc1=idep
         enddo

         loc3=1
         do idep=1,nbdepl -1
            if ((ar1(idep) >= 0.0001).and.(catdef(idep) <= cdcr1)) loc3=idep + 1
         enddo

         if (loc3.le.loc1+1) then
            loc1=MIN(loc3-4,loc1-4)
            loc1=MAX(1,loc1)
         endif

!c below is what was used for no regression, but it's not equivalent to the 
!c IDL program
         loc2=AINT(loc1-1+(loc3-loc1)*3./5.)+1
         
         w1=ar1(loc1)
         w2=ar1(loc2)
         w3=ar1(loc3)

         if(w3.eq.0.)then
 95         loc3=loc3-1
            if(loc3.eq.loc2)loc2=loc2-1
               w3=ar1(loc3)
               w2=ar1(loc2)
               if(w3.eq.0.)goto 95
            endif
            w4=0.

            if((loc1.ge.loc2).or.(loc2.ge.loc3))then
               loc1=10
               loc2=14
               loc3=18
            endif

 115        x1=catdef(loc1)
            x2=catdef(loc2)
            x3=catdef(loc3)
            w1=ar1(loc1)
            w2=ar1(loc2)
            w3=ar1(loc3)
               
            if (bug) then
               write(*,*) 'loc1,loc2,loc3=',loc1,loc2,loc3
               write(*,*) 'x1,x2,x3=',x1,x2,x3
               write(*,*) 'w1,w2,w3=',w1,w2,w3
            endif

            y0=w4
            f1=(1.-w1)/(w1-y0)/x1
            f2=(1.-w2)/(w2-y0)/x2
            f3=(1.-w3)/(w3-y0)/x3
            g1=(1.-y0)/(w1-y0)
            g2=(1.-y0)/(w2-y0)
            g3=(1.-y0)/(w3-y0)
            df=f2-f1
            dg=g2-g1
            dx=x2-x1
            bf=f1-x1*df/dx
            bg=g1-x1*dg/dx

            ars1 = -(f3-bf-x3*df/dx)/(g3-bg-x3*dg/dx + 1.e-10)
            ars2 = bf+ars1*bg
            ars3 = (df+ars1*dg)/dx

            delta=ars2*ars2-4*ars3
            upval=1.+200.*ars1
            loval=1.+200.*ars2+40000.*ars3
	    z1=0.
            z2=0.

            if (delta .ge. 0.) then !if 8
               z1=(-ars2-SQRT(delta))/2./ars3
               z2=(-ars2+SQRT(delta))/2./ars3 
            endif

            if ((z1 .gt. 0. .and. z1 .lt. cdcr1) .or.   &
               (z2 .gt. 0. .and. z1 .lt. cdcr1) .or.  &
               ((upval/loval).lt.-.01)) then   !if 7
               z1=0.
               z2=0.
               if (loc1 .eq. 10) then 
                  loc1=1
1              else  
                  loc1=1
                  do idep=1,nbdepl
                     if (catdef(idep) .gt. 60.) then
                        loc1=idep
                        if(loc1.ge.loc3-1)then
                        !                           write(*,*)'Loc1 exceeded loc3 in 2nd attempt'
                           loc1=loc3-5
                        endif
                        goto 46
                      endif
                  enddo
            endif
46          loc2=loc1+AINT(float(loc3-loc1)*3./5.)+1
            if(loc2.ge.loc3)loc2=loc3-1
            loc2save=loc2
            INC=1
            DEC=0

47          w1=ar1(loc1)
            w2=ar1(loc2)
            x1=catdef(loc1)
            x2=catdef(loc2)
            
            if (bug) then
               write(*,*) 'z1,z2=',z1,z2,' -> ar1, 2nd try'
               write(*,*) 'loc1,loc2,loc3=',loc1,loc2,loc3
               write(*,*) 'x1,x2,x3=',x1,x2,x3
               write(*,*) 'w1,w2,w3=',w1,w2,w3
            endif
            
            f1=(1.-w1)/(w1-y0)/(x1 + 1.e-20)
            f2=(1.-w2)/(w2-y0)/(x2 + 1.e-20)
            g1=(1.-y0)/(w1-y0 + 1.e-20 )
            g2=(1.-y0)/(w2-y0 + 1.e-20)
            df=f2-f1
            dg=g2-g1
            dx=x2-x1
            bf=f1-x1*df/dx
            bg=g1-x1*dg/dx
            
            ars1 = -(f3-bf-x3*df/dx)/(g3-bg-x3*dg/dx  + 1.e-10)
            ars2 = bf+ars1*bg
            ars3 = (df+ars1*dg)/dx
            delta=ars2*ars2-4*ars3
            upval=1.+200.*ars1
            loval=1.+200.*ars2+40000.*ars3

            if (delta .ge. 0.) then   !if 6
               z1=(-ars2-SQRT(delta))/2./ars3
               z2=(-ars2+SQRT(delta))/2./ars3
            end if

            if ((z1 .gt. 0. .and. z1 .lt. cdcr1) .or.   &
                 (z2 .gt. 0. .and. z1 .lt. cdcr1) .or.   &
                 ((upval/loval).lt.-.01)) then  !if 5
               !c Sarith ---
               z1=0.
               z2=0.
               IF(INC.EQ.1)loc2=loc2+1
               IF(DEC.EQ.1)LOC2=LOC2-1
               if(inc.eq.1)then   !if 4
                  if(loc2.ge.loc3)then   !if 3
                     !                     WRITE(*,*)'INCREASING LOC2 FAILED'
                     INC=0
                     DEC=1
                     loc2=loc2save
                  else
                     adjust=ADJUST+1
                     goto 47
                  end if    !if 3
               endif    !if 4

               if(dec.eq.1)then   !if 2  
                  if(loc2.eq.loc1)then  !if 1
                     !                     WRITE(*,*)'Decreasing too failed'
                     INC=1
                     DEC=0
                     ars1=9999. !ars1old
                     ars2=9999. !ars2old
                     ars3=9999. !ars3old
                     !                     write(*,*) 'AR1: PROBLEM for pfc=',pfc
                  else
                     adjust=ADJUST+1
                     !c                        write(*,*)'ADJUSTING AR1 CYCLE =',ADJUST
                     goto 47
                  end if   !if 1
               endif  !if 2 
            endif     !if 5
            !c               endif    !if 6
         endif            !if 7
         
         !c         endif  !if 8
         flag=0
         call curve1(ars1,ars2,ars3,cdcr2,flag)

         IF (FLAG.EQ.1)then
            !            WRITE(*,*)'Curve problem in the catchment pfc=',pfc
            ars1=9999. 
            ars2=9999. 
            ars3=9999. 
            !                     write(*,*) 'Pick values from icatch-1'
            flag=0
         end if
      endif

         adjust=0

         if (bug) write(*,*) 'ar1 adjustment ok'

!c**** WMIN adjustment: 3 points + limit in INF = wmin0

         if (bug) write(*,*) 'STARTING WMIN'

         w4=wmin0
         y0=w4

!         write(*,*) 'wmin=',(wmin(idep),idep=1,50)

            loc1=1
            do idep=1,nbdepl
               if (catdef(idep) <= 10.) loc1=idep
             enddo
                  
            loc3=1
            do idep=1,nbdepl - 2
               if ((wmin(idep) >= wmin0).and.(catdef(idep) <= cdcr1)) loc3=idep + 2
            enddo
 
            loc2=loc1 + 2
            do idep=1,nbdepl -1 
               if ((wmin(idep) >= wmin0).and.(catdef(idep) <= cdcr1/2.))loc2=idep + 1
            enddo 

!c For global catch         
            INC=1
            DEC=0

          if(loc3.eq.loc2)loc2=loc2-2
          if(loc2 <= loc1) loc1= loc1-2
 44         loc2save=loc2
            if(loc1 < 1) then
               loc1 =1
               loc2 =2 
               loc3 =3
            endif

            w1=wmin(loc1)
            w2=wmin(loc2)
            w3=wmin(loc3)
            x1=catdef(loc1)
            x2=catdef(loc2)
            x3=catdef(loc3)
         
            if (bug) then
               write(*,*) 'loc1,loc2,loc3=',loc1,loc2,loc3
               write(*,*) 'x1,x2,x3=',x1,x2,x3
               write(*,*) 'w1,w2,w3,w4=',w1,w2,w3,w4
            endif

            f1=(1.-w1)/(w1-y0)/x1
            f2=(1.-w2)/(w2-y0)/x2
            f3=(1.-w3)/(w3-y0)/x3
            g1=(1.-y0)/(w1-y0)
            g2=(1.-y0)/(w2-y0)
            g3=(1.-y0)/(w3-y0)
            df=f2-f1
            dg=g2-g1
            dx=x2-x1
            bf=f1-x1*df/dx
            bg=g1-x1*dg/dx
            
            arw1 = -(f3-bf-x3*df/dx)/(g3-bg-x3*dg/dx + 1.e-10)
            arw2 = bf+arw1*bg
            arw3 = (df+arw1*dg)/dx
            arw4 = y0

!c wmin=arw4+(1.-arw4)*(1.+arw1*catdef(idep))
!c     /(1.+arw2*catdef(idep)+arw3*catdef(idep)*catdef(idep))
!c we want to check the roots of the denominator

            delta=arw2*arw2-4*arw3

            if (delta .ge. 0.) then !if 8

               z1=(-arw2-SQRT(delta))/2./arw3
               z2=(-arw2+SQRT(delta))/2./arw3

               if ((z1 .gt. 0. .and. z1 .lt. cdcr1) .or.           &
                   (z2 .gt. 0. .and. z1 .lt. cdcr1)) then !if 7

                  w1_0=w1
                  w1=(1.+w1_0)/2.
                  x1=x1/4.

!                  if (gnu .eq. 3.26/1.5) then 
!                        w1=(1.+w1_0)/3.               ! already difficult
!                        w3=wmin(nint(cdcr1))                ! with gnu=3.26
!                        x3=catdef(nint(cdcr1))
!                        f3=(1.-w3)/(w3-y0)/x3
!                        g3=(1.-y0)/(w3-y0)
!                  endif
                 
                  f1=(1.-w1)/(w1-y0)/x1
                  g1=(1.-y0)/(w1-y0)
                  df=f2-f1
                  dg=g2-g1
                  dx=x2-x1
                  bf=f1-x1*df/dx
                  bg=g1-x1*dg/dx

                  if (bug) then
                     write(*,*) 'z1,z2=',z1,z2,' -> wmin, 2nd try'
                     write(*,*) 'loc1,loc2,loc3=',loc1,loc2,loc3
                     write(*,*) 'x1,x2,x3=',x1,x2,x3
                     write(*,*) 'w1,w2,w3=',w1,w2,w3
                     write(*,*) 'wmin0=',wmin0
                  endif
               
                  arw1 = -(f3-bf-x3*df/dx)/(g3-bg-x3*dg/dx + 1.e-10)
                  arw2 = bf+arw1*bg
                  arw3 = (df+arw1*dg)/dx
                  arw4 = y0
                  
                  delta=arw2*arw2-4*arw3
                  
                  if (delta .ge. 0.) then  !if 6 
                     z1=(-arw2-SQRT(delta))/2./arw3
                     z2=(-arw2+SQRT(delta))/2./arw3

                     if ((z1 .gt. 0. .and. z1 .lt. cdcr1) .or.         &
                         (z2 .gt. 0. .and. z1 .lt. cdcr1)) then    !if 5  
!c Sarith ---
                     IF(INC.EQ.1)loc2=loc2+1
                     IF(DEC.EQ.1)LOC2=LOC2-1
                     if(inc.eq.1)then   !if 4
                     if(loc2.eq.loc3)then   !if 3
!                     WRITE(*,*)'INCREASING LOC2 FAILED: WMIN'
                     INC=0
                     DEC=1
                     loc2=loc2save
                     else
                        adjust=ADJUST+1
!c                        write(*,*)'ADJUSTING AR1 CYCLE =',ADJUST
                        goto 44
                        end if    !if 3
                        endif    !if 4
                     if(dec.eq.1)then   !if 2  
                     if(loc2.eq.loc1)then  !if 1
!                     WRITE(*,*)'Decreasing too failed: WMIN'
                     INC=1
                     DEC=0

                     arw1=9999. 
                     arw2=9999. 
                     arw3=9999. 
                     arw4=9999. 

                     else
                        adjust=ADJUST+1
!c                        write(*,*)'ADJUSTING AR1 CYCLE =',ADJUST
                        goto 44
                        end if   !if 1
                     endif  !if 2                         
                     endif !if 5     
                  endif   !if 6
                  
               endif !if 7
            endif !if 8
         adjust=0
!         endif ! pfc=12821
         flag=0

         call curve2(arw1,arw2,arw3,arw4,cdcr1,flag)

         IF (FLAG.EQ.1) THEN
	     arw1=9999. !arw1old
             arw2=9999. !arw2old
             arw3=9999. !arw3old
             arw4=9999. !arw4old
             flag=0
         endif

         if(arw1==9999.) then 
! Singular Value Decomposition

         w4=wmin0
         y0=w4

            loc1=1
            loc3=nbdepl

         mp = loc3-loc1+1   

         if(mp.lt.3)then

            write(*,*)'WMIN Note: not sufficient points MP = ',mp
	    print *,w4,cdcr1,catdef(loc3),wmin(loc3)
                     arw1 = 9999.
                     arw2 = 9999.
                     arw3 = 9999.
                     arw4 = 9999.
            else
	    
            mp = 1
            istart =1
            w4 = wmin(istart)

            if(w4 <=0) then
               do idep=2,nbdepl
                  if(wmin(idep) > 0.) istart = idep
                  if(wmin(idep) > 0.) exit
               enddo
            endif

            w4 = wmin(istart)

  	    do idep=istart+1,nbdepl
!	       if(wmin(idep).lt.w4) then
               if((wmin(idep) - w4).lt.0.0005) then
	       	w4 = wmin(idep)
		mp = mp +1
	       endif		    		    	          
	    enddo             
	    loc3 = mp   
 	    allocate(A(mp,3))
            allocate(AP(mp,3))
            allocate(B(mp))
            allocate(BP(mp))             
	       smooth = .false.
	       do idep=istart,nbdepl-1
	         if(catdef(idep).le.cdcr1+10.) then
		 if((wmin(idep) - wmin(idep +1)) < -0.01) smooth = .true.   
		 endif  
	       enddo
	       if(smooth) then
	       wminsave = wmin
               ! Apply filter to input data
               do i=istart, nbdepl-nr
    	            wmin(i)=0.
                    do j=1, nl+nr+1
	                if (i+savgol_ind(j).gt.0) then  !skip left points that do not exist
	    	           wmin(i)=wmin(i)+savgol_coeff(j)*wminsave(i+savgol_ind(j))
                        endif
                    end do
               enddo	         
	       wmin (istart:istart+4) = wminsave (istart:istart+4)

	       endif

	       j = 1
	       w4 = wmin(istart)
               do isvd=1,size(wmin)
                  if (j <= mp) then 
                     if(isvd == 1) then
                        wbrac=(wmin(isvd + istart -1)-y0)/(1.-y0 + 1.e-20)
                        A(j,1)=catdef(isvd + istart -1)
                        A(j,2)=-catdef(isvd + istart -1)*wbrac
                        A(j,3)=-wbrac*((catdef(isvd + istart -1))**2.)
                        B(j)=wbrac-1.
                        j = j + 1
                     else
                        if((wmin(isvd + istart -1).lt.w4).and.(wmin(isvd + istart -1).gt.y0)) then
                           wbrac=(wmin(isvd + istart -1)-y0)/(1.-y0 + 1.e-20)
                           A(j,1)=catdef(isvd + istart -1)
                           A(j,2)=-catdef(isvd + istart -1)*wbrac
                           A(j,3)=-wbrac*((catdef(isvd + istart -1))**2.)
                           B(j)=wbrac-1.
                           w4 = wmin(isvd + istart -1)
                           j = j + 1 
                        endif
                     endif
                  endif
               end do

               j = j -1 
               mp = j
               ap => a (1:j,:)
               bp => b (1:j)
               ap(j,1) = catdef(nbdep)
               ap(j,2) = 0.
               ap(j,3) = 0.
               bp (j) = -1.

               call svdcmp(ap,mp,3,w,v)

               sdmax=0.
               do j=1,3
                  if(w(j).gt.sdmax)sdmax=w(j)
               end do

               sdmin=sdmax*1.0e-6
               do j=1,3
                  if(w(j).lt.sdmin)w(j)=0.
               end do

               call svbksb(ap,w,v,mp,3,bp,ans)

               arw1 = real(ans(1))
               arw2 = real(ans(2))
               arw3 = real(ans(3))
               arw4 = y0
               
!c wmin=arw4+(1.-arw4)*(1.+arw1*catdef(idep))
!c     /(1.+arw2*catdef(idep)+arw3*catdef(idep)*catdef(idep))
!c we want to check the roots of the denominator
               
               adjust=0         
               flag=0
               
               call curve2(arw1,arw2,arw3,arw4,cdcr1,flag)
               
               IF (FLAG.EQ.1) THEN
                  !            WRITE(*,*)'Curve2 problem in the catchment:pfc=',pfc
                  
                  arw1 = 9999.
                  arw2 = 9999.
                  arw3 = 9999.
                  arw4 = 9999.
                  
                  flag=0
               end if
               deallocate (A,  B )
               NULLIFY    (AP, BP)
            end if
         endif
         
         if(present(dbg_unit)) then
             write (dbg_unit,*) ars1,ars2,ars3
             write (dbg_unit,*) arw1,arw2,arw3,arw4 
         endif

         if (bug) write(*,*) 'wmin adjustment ok'
         
!c**** SHAPE PARAMETER ADJUSTMENT: with a straight if coeskew > 0.25
!c                                 with 2 segments if not

         if (bug) write(*,*) 'STARTING SHAPE'

         x3=catdef(nbdepl)
         w3=aa(nbdepl)
         x1=0.

         if (coeskew .lt. 0.25) then
            w1=0.1
            loc2=20
            do idep=1,nbdepl
               if (catdef(idep) .gt. ref1) then
                  loc2=idep
                  goto 45
               endif
            enddo
 45         x2=catdef(loc2)
            w2=aabis(loc2)
            ara1 = (w1-w2)/(x1-x2)
            ara2 = w1-ara1*x1
            ara3 = (w2-w3)/(x2-x3)
            ara4 = w2-ara3*x2
         else
            w1=1.
            x2=x1
            w2=w1
            ara3 = (w2-w3)/(x2-x3)
            ara4 = w2-ara3*x2
            ara1 = ara3
            ara2 = ara4       
         endif

         if (bug) write(*,*) 'x1,w1,x2,w2,x3,w3',x1,w1,x2,w2,x3,w3

!**** RMSE checking: on ar1, ar2, swsrf2 and rzeq

         do idep=1,nbdepl
            if(catdef(idep) <= cdcr1) then
            nar1(idep)=AMIN1(1.,AMAX1(0.,(1.+ars1*catdef(idep)) &
                /(1.+ars2*catdef(idep)                          &
                +ars3*catdef(idep)*catdef(idep))))                    
                                                                 
           nwm=AMIN1(1.,AMAX1(0.,arw4+(1.-arw4)*                &
                (1.+arw1*catdef(idep))                          &
                /(1.+arw2*catdef(idep)                          &
                +arw3*catdef(idep)*catdef(idep))))

!c we have to first determine if there is one or two segments
            if (ara1 .ne. ara3) then
               cdi=(ara4-ara2)/(ara1-ara3)
            else
               cdi=0.
            endif

            if (catdef(idep) .ge. cdi) then
               shape=ara3*catdef(idep)+ara4
            else
               shape=ara1*catdef(idep)+ara2
            endif
	    shape =AMIN1(40.,shape)
            area1=exp(-shape*(1.-nwm))*(shape*(1.-nwm)+1.)

!c the threshold for truncation problems is higher than the "usual"
!c E-8 to E-10, because it plays together with the uncertainties coming 
!c from the approximation of the parameters nwm, nar1 and shape.
            if (area1 .ge. 1.-1.E-8) then
		nar1(idep)=1.
		nar2(idep)=0.
		nar3(idep)=0.
		nmean2(idep)=0.  
		nmean3=0.       
                neq(idep)=1.
             else
                
                if (nwm .gt. wpwet) then
                   nar2(idep)=1.-nar1(idep)
                else
                   nar2(idep)=AMAX1(0.,((shape*(wpwet-nwm)+1.)        &
                       *exp(-shape*(wpwet-nwm))                       &
     			- (shape*(1.-nwm)+1.)*exp(-shape*(1.-nwm)))   &
     			* (1.-nar1(idep))/(1.-area1))                  
               endif                                                        
                                                                       
               nar3(idep)=1.-nar1(idep)-nar2(idep)                          
                                                                       
               if (nar3(idep) .lt. 1.E-8) then ! for nwm le wpwet           
                                                                       
                  nmean2(idep)=AMAX1(0.,AMIN1(1.,(nwm + 2./shape +    &
                       shape*exp(-shape*(1.-nwm))*                    &
     			(nwm+nwm/shape-1.-2./shape-2./(shape*shape))) &
     			/(1.-area1)))
                   nmean3=0.

                else

!c WARNING: I think the two values below are false. 
!c But it is never used in this context, because nwm > wpwet !!
                   nmean2(idep)=AMAX1(0.,AMIN1(1.,-shape*(exp(-shape*&
                       (wpwet-nwm))* (nwm*wpwet                      &
                       +nwm/shape-wpwet*wpwet                        &
                       -2.*wpwet/shape-2./(shape*shape))             &
                       - exp(-shape*(1.-nwm))*                       &
     			(nwm+nwm/shape-1.-2./shape-2./(shape*shape)))& 
     			* (1.-nar1(idep))/(1.-area1) / (nar2(idep)+1.e-20)))        
                                                                      
                   nmean3=AMAX1(0.,AMIN1(1.,(nwm+2./shape +          &
                       shape*exp(-shape*(wpwet-nwm))*                &
     			(nwm*wpwet+nwm/shape-wpwet                   &
                       *wpwet-2.*wpwet/shape                         &
                       -2./(shape*shape))) * (1.-nar1(idep))         &
                       /(1.-area1)/(nar3(idep) + 1.e-20)))                       
               endif                                                       
                                                                      
               neq(idep)=nar1(idep)+nar2(idep)*nmean2(idep)          &
                    +nar3(idep)*nmean3
                
                if (area1 .ge. 1.-1.E-5) then
                   nmean2(idep)=1.  
                   nmean3=0.       
                   neq(idep)=1.
                endif

             endif
             endif
          enddo

          if (bug) write(*,*) 'shape adjustment ok'
!c
!c RMSE

!c ERR1
          icount=0
          iref=0
          sum=0.
          do i=1,nbdepl
             if(catdef(i) <= cdcr1) then
             tabact(i)=0.
             tabfit(i)=0.
             endif
          enddo

          do i=1,nbdepl
            if(catdef(i) <= cdcr1) then 
             if (catdef(i) .gt. lim) then
               icount=icount+1
               sum=sum+ar1(i)
               tabfit(icount)=nar1(i)
               tabact(icount)=ar1(i)
             endif
             endif
          enddo

	  if(icount.gt.1) then
          sum=sum/icount
          call RMSE(tabact,tabfit,icount,err1)
          taberr1=err1
          normerr1=err1/sum
	  endif
!c ERR2
          icount=0
          iref=0
          sum=0.
          do i=1,nbdepl
             if(catdef(i) <= cdcr1) then
             tabact(i)=0.
             tabfit(i)=0.
             endif
          enddo

          do i=1,nbdepl
             if(catdef(i) <= cdcr1) then
             if (catdef(i) .gt. lim) then
               icount=icount+1
               sum=sum+ar2(i)
               tabfit(icount)=nar2(i)
               tabact(icount)=ar2(i)
             endif
             endif
          enddo

	  if(icount.gt.1) then
          sum=sum/icount
          call RMSE(tabact,tabfit,icount,err2)
          taberr2=err2
          normerr2=err2/sum
	  endif

!c ERR3
          icount=0
          iref=0
          sum=0.
          do i=1,nbdep
             if(catdef(i) <= cdcr1) then
             tabact(i)=0.
             tabfit(i)=0.
             endif
          enddo

          do i=1,nbdepl
             if(catdef(i) <= cdcr1) then
             if (catdef(i) .gt. lim) then
               icount=icount+1
               sum=sum+swsrf2(i)
               tabfit(icount)=nmean2(i)
               tabact(icount)=swsrf2(i)
             endif
             endif
          enddo

	  if(icount.gt.1) then
          sum=sum/icount
          call RMSE(tabact,tabfit,icount,err3)
          taberr3=err3
          normerr3=err3/sum
	  endif
!c ERR4
          icount=0
          iref=0
          sum=0.
          do i=1,nbdepl
             tabact(i)=0.
             tabfit(i)=0.
          enddo

          do i=1,nbdepl
             if(catdef(i) <= cdcr1) then
             if (catdef(i) .gt. lim) then
               icount=icount+1
               sum=sum+rzeq(i)
               tabfit(icount)=neq(i)
               tabact(icount)=rzeq(i)
             endif
             endif
          enddo

	  if(icount.gt.1) then 
          sum=sum/icount
          call RMSE(tabact,tabfit,icount,err4)
          taberr4=err4
          normerr4=err4/sum
	  endif
     END SUBROUTINE SAT_PARAM
!

! ******************************************************************

!c
      SUBROUTINE CURVE1(ars1,ars2,ars3,cdcr1,flag)
      REAL ars1,ars2,ars3,y,x,yp,cdcr1
      INTEGER i,flag
!c
      yp=1.
      if (abs(ars1+ars2+ars3).le.1.e25) then
      do i=0,nint(cdcr1)
         x=float(i)
         y=(1.+ars1*x)/(1.+ars2*x+ars3*x*x + 1.e-20)
         if((y.gt.0.0).and.(((yp -y) .lt. -1.e-4).or.(y.gt.1.)))then
            flag=1
            goto 99
         endif   
         yp=y
      end do
 99   continue
      else
	flag=1
      endif

    end SUBROUTINE CURVE1


! ******************************************************************

      SUBROUTINE CURVE2(arw1,arw2,arw3,arw4,cdcr1,flag)
      REAL arw1,arw2,arw3,arw4,y,x,yp,cdcr1
      INTEGER i,flag
!c
      yp=1.
      if (abs(arw1+arw2+arw3+arw4).le.1.e25) then
      do i=0,nint(cdcr1)
         x=float(i)
         y=arw4+(1.-arw4)*(1.+arw1*x)/(1.+arw2*x+arw3*x*x + 1.e-20)

         if((y.ge.arw4).and.(((yp -y) .lt. -1.e-4).or.(y.gt.1.)))then
            flag=1
            goto 99
         endif  
         yp=y
      end do
99    continue
      else
      flag=1
      endif
    end SUBROUTINE CURVE2


! ******************************************************************

      subroutine tgen (                         &
          TOPMEAN,TOPVAR,TOPSKEW,               &
          STO,ACO,COESKEW)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                         c
! The difference between tgen4 and tgen3 is that tgen4 deals with arrays  c
! of topmean, topvar and topskew and 2-dim arrays of st and ac.           c
!                                                                         c
! This routine determine the theoretical gamma distribution for the       c
! soil-topographic indexes (Sivapalan et al., 1987), knowing the three    c
! first moments, the min and the max of the observed topographic indexes  c
! in a given catchment.                                                   c
!                                                                         c
! Routine from Dave Wolock.                                               c
! Modified by Agnes (11-06-98): we don't use min and max anymore, and     c
! this strongly improves the behavior for negative skewnesses. It also    c
! improves in general the matching of the moments.                        c
!                                                                         c
! We also add a correction on the skewness to have gamma distributions    c
! that start and end from the x-axis. It is based on the fact that if     c
! TOPETA=1, the gamma is an exponential distribution, and if TOPETA<1,    c
! then the gamma distribution increases towards the infinite when x       c
! decreases towards 0.                                                    c
! To eliminate some numerical pb due to teh discretization of the gamma   c
! distribution, we choose skewness=MAX(MIN(1.9, skewness),-1.6)           c
!                                                                         c
! WE MAY NEED TO COMPUTE IN DOUBLE RESOLUTION !!!! BECAUSE OF THE SMALL   c
! BIN WIDTH
!                                                                         c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  

      IMPLICIT NONE

      real, parameter :: VALMAX=50.
      REAL, intent (in) :: TOPMEAN,TOPVAR,TOPSKEW
      REAL, intent (out) :: COESKEW
      REAL, dimension (NAR), intent (out) :: STO,ACO

      INTEGER I
      REAL ST(NAR),AC(NAR)
      REAL TOPETA,TOPLAM,TOPSCAL,GAMLN,SCALE,ACLN
      real cumac, cum2,cum3

!-------------------------------------------------------------------------

! topmean is the mean of the ln(a/tanB) distribution
! topvar is the variance (2nd moment centerd around the mean) of the ...
! topskew is the skew (3rd moment centerd around the mean) of the ...
! compute the coefficient of skew or skewness (coeskew)

         COESKEW=TOPSKEW/TOPVAR**1.5
         if (coeskew .ge. 0.) then
            COESKEW=AMAX1(0.005, AMIN1(1.9, COESKEW))
         else
            COESKEW=AMAX1(-1.6, AMIN1(-0.005, COESKEW))
         endif

! compute the gamma parameters, eta (topeta) and lambda (toplam), and topscal
! which is the translation parameter

         TOPETA=4./COESKEW**2
         TOPLAM=SQRT(TOPETA)/SQRT(TOPVAR)
         TOPSCAL=TOPMEAN-TOPETA/TOPLAM

! evaluate the gamma function

         CALL GAMMLN(TOPETA,GAMLN)

         CUMAC=0.0

! compute the frequency distribution of ln(a/tanB)
! st(i) are the values of ln(a/tanB)
! ac(i) are the relative frequency values (they should sum to 1)

         DO I=1,NAR
         
            ST(I)=(FLOAT(I)-0.95)*(VALMAX-TOPSCAL)/FLOAT(NAR)+TOPSCAL
            SCALE=ST(I)-TOPSCAL

! below is the logarithmic form of the gamma distribution; this is required 
! because the numerical estimate of the logarithm of the gamma function 
! is more stable than the one of the gamma function.
          
            ACLN=TOPETA*ALOG(TOPLAM)+(TOPETA-1.)*ALOG(SCALE)  &
                -TOPLAM*SCALE-GAMLN
                      
            IF(ACLN.LT.-10.) THEN
               AC(I)=0.
            ELSE
               AC(I)=EXP(ACLN)
            ENDIF

            CUMAC=CUMAC+AC(I)

         ENDDO

! we want the relative frequencies to sum 1.

         IF (CUMAC.eq.0.) THEN
!            write(*,*) 'distrib sum=',CUMAC
            stop
         endif
         CUM2=0.
         DO I=1,NAR
            AC(I) = AC(I) / CUMAC
            CUM2=CUM2+AC(I)
         ENDDO
      
! if the real distribution of the topographic indices is negativeley skewed, 
! we symetrize the gamma distribution (depending on coeskew**2 and always 
! positively skewed), centering on topmean, which preserves topmean and
! topvar, and re-establishes a negative skewness.

         IF (COESKEW.LT.0.) then

            do i=1,nar
               STO(I)=2.*TOPMEAN-ST(I)
               ACO(I)=AC(I)

            enddo
         ELSE
!            if (n .eq. idmax) then
!               write(*,*) 'last catchment'
!            endif
            do i=1,nar
               STO(I)=ST(-I+NAR+1)
               ACO(I)=AC(-I+NAR+1)
            enddo
         ENDIF

!         sum=0.
!         do i=1,nar
!            sum=sum+sto(i)*aco(i)
!         end do

!         sum=0.
!         do i=1,nar
!            sum=sum+aco(i)
!         end do


    END subroutine tgen

  
  ! ********************************************************************

    SUBROUTINE GAMMLN(XX,GAMLN)
      
      implicit none
      DOUBLE PRECISION :: COF(6),STP,HALF,ONE,FPF,X,TMP,SER
      REAL, intent(in) :: XX 
      REAL, intent(out) :: GAMLN
      integer :: j

      DATA COF /76.18009173D0,-86.50532033D0,24.01409822D0,    &
         -1.231739516D0,.120858003D-2,-.536382D-5/
      STP = 2.50662827465D0
      HALF= 0.5D0
      ONE = 1.0D0
      FPF = 5.5D0
      
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE

      DO  J=1,6
         X=X+ONE
         SER=SER+COF(J)/X
      END DO

      GAMLN=TMP+LOG(STP*SER)
      
    END SUBROUTINE GAMMLN
  
  ! ********************************************************************

    SUBROUTINE FUNCIDEP(                                         &
                         NAR0,IDEP,                              &!I
                         BEE,PSIS,POROS,COND,RZDEP,WPWET,        &!I
                         VALX,PX,COESKEW,TIMEAN,SUMA,            &!I
                         CATDEF,AR1,WMIN,AA,AABIS,               &!O
                         AR2,AR3,SWSRF2,SWSRF3,RZEQ)             

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                                      c
!c This program returns the eight parameters for the areal fractioning  c
!c                                                                      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      integer, intent (in) :: NAR0,idep 
      REAL, intent (in) :: BEE, PSIS, POROS, COND, RZDEP, WPWET, COESKEW
      REAL, intent (inout) ::  VALX(NAR), PX(NAR),TIMEAN,SUMA
!      logical, intent(in) :: bug
      real, dimension (nbdep), intent (inout) :: CATDEF,AR1,WMIN,AA,   &
           AABIS,AR2,AR3,SWSRF2,SWSRF3,RZEQ 
      INTEGER :: width, nref, nind,nmax,indmin,locmax,shift,ord,locmin,ordref
      integer :: indimax10,indmin0,k,n,n1,n2
      real dx,zbar

      real test,term1,term2,sum
      real zdep(nar),locdef(nar),wrz(nar),frcunsat
      real valtest(nbdep,nar),ptest(nbdep,nar),denstest(nbdep,nar)
      real dtest(nbdep,nar),cump
      real x1,x2,y1,y2,wa,wb
      real densaux(nar),densaux2(nar),densmax,aux10
      real :: dz, sumdef
!c-------------------------------------------------------------------------

!c integral(f(x)dx)=1. for a pdf
!c here px=f(x)dx
      dx=valx(1)-valx(2)

      if (bug) write(*,*) 'IDEP=',IDEP,' dx=',dx

!c the loops over idmax and nbdep are initiated in sta_params4.f

      zbar=float(idep-10)*slice ! zdep in meters
            
!c**** Compute array of water table depths:
      do k=1,nar0
         term1=(1/gnu)*(valx(k)-timean)
         zdep(k)=AMAX1(0.,zbar-term1)
      enddo
            
!c variable change must be reflected in dx
      dz=dx/gnu

      if (bug) write(*,*) 'funcidep: ok1'
      
!c**** Compute array of moisture deficits:
      do k=1,nar0
         term1=(psis-zdep(k))/psis
         term1=term1**(1.-1./bee)
         term2=-psis*(bee/(bee-1.))*(term1-1.)
         locdef(k)=zdep(k)-term2
      enddo

!c**** Add deficits to produce catdef:
      sumdef=0.
      do k=1,nar0
         sumdef=sumdef+locdef(k)*px(k)
      enddo
      catdef(idep)=poros*1000.*sumdef/suma

      if (bug) write(*,*) 'funcidep: ok2'

!c**** Compute array of root zone moisture (degree of wetness in root zone):
      do k=1,nar0              
         term1=((psis-zdep(k))/psis)              &
             **(1.-1./bee)
        if(zdep(k).le.0.) then
           wrz(k)=1.
        elseif(zdep(k)-rzdep.lt.0.) then
           term2=(-psis/zdep(k))*(bee/(bee-1.))  &
                *(term1-1.)
           frcunsat=zdep(k)/rzdep
           wrz(k)=frcunsat*term2+(1.-frcunsat)*1.
        else
           term2=((psis-zdep(k)+rzdep)           &
                /psis)**(1.-1./bee)
           wrz(k)=(-psis/rzdep)*(bee/            &
                (bee-1.))*(term1-term2)
         endif

      enddo

      if (bug) write(*,*) 'funcidep: ok3'

!c**** compute the densities and dx
!c**** we use a usefull property that is due to the construction of the 
!c**** gamma distribution in tgen3.f : this distribution is continuous, 
!c**** with decreasing values on ln(a/tanb) when n goes from 1 to nar0

!c first we gather in the same bin all the bins with values ge 1 
      nref=1
      nind=1
      ptest(idep,1)=0.
      do k=1,nar0
         if (wrz(k) .eq. 1.) then
            nref=nref+1
            ptest(idep,1) = ptest(idep,1) + px(k)
         endif
      enddo
      if (nref .gt. 1) then
         nind=2
         valtest(idep,1)=1.
      endif
      nmax=nar0-nref+nind
      if (bug) write(*,*) 'nmax,nind,nar0,nref=',nmax,nind,nar0,nref
      
!c definition of the probabilities ptest
      if (nmax .eq. 1) then     ! all the bins have values ge 1 
         dtest(idep,1) = 0.0001
         ptest(idep,1) = 1.
      else                      ! distribution in ar2/ar3
         do n=0,nmax-nind
            valtest(idep,nind+n)=wrz(nref+n)
            ptest(idep,nind+n)=px(nref+n)
         enddo
         
!c we have to define dtest, the size of each bin
         if (nmax .eq. 2) then
            dtest(idep,2) = valtest(idep,1)-valtest(idep,2)
            dtest(idep,1) = dtest(idep,2)/2.
         else                   ! nmax .gt. 2
            do n=2,nmax-1                     
               dtest(idep,n)=(valtest(idep,n-1)-valtest(idep,n+1))/2.   
            enddo
               dtest(idep,1) = dtest(idep,2)/2.
            dtest(idep,nmax) = dtest(idep,nmax-1)
         endif
      endif

      if (bug) write(*,*) 'funcidep: ok4'

!c we can now define the probability density: denstest=ptest/dtest
!c where ptest is the probability and dtest the size of the bin
      do n=1,nmax
         if (ptest(idep,n) .eq. 0.) then
            denstest(idep,n)=0.
         else
            denstest(idep,n)=ptest(idep,n)/dtest(idep,n)
         endif
      enddo

      if (bug) write(*,*) 'funcidep: ok5'

!c NOW we can estimate the parameters for the approximated distrib
!c from the actual distrib

!c 1. AR1=saturated area and AR2 and AR3 + averages of the RZ wetness 
!c    in the different fractions

      ar1(idep)=0.
      ar2(idep)=0.
      ar3(idep)=0.
      swsrf3(idep)=0.
      swsrf2(idep)=0.
      rzeq(idep)=0.
    
      if(valtest(idep,1).eq.1.) ar1(idep)=dtest(idep,1)*denstest(idep,1)
      
      if (nmax .gt. 1) then 
         do n=nind,nmax           
            if (valtest(idep,n) .lt. wpwet) then
               ar3(idep)=ar3(idep)+denstest(idep,n)*dtest(idep,n)
               swsrf3(idep)=swsrf3(idep)+valtest(idep,n)*         &
                    denstest(idep,n)*dtest(idep,n)
            else
               ar2(idep)=ar2(idep)+denstest(idep,n)*dtest(idep,n)
               swsrf2(idep)=swsrf2(idep)+valtest(idep,n)*         &
                   denstest(idep,n)*dtest(idep,n)
            endif
         enddo
      endif
         
      test=ar1(idep)+ar2(idep)+ar3(idep)
      if (test .gt. 1.+1.e-5 .or. test .lt. 1.-1.e-5) then
!         write(*,*) 'PROBLEM at depth ',zbar
!         write(*,*) '  ar1+ar2+ar3=',test
!         write(*,*) '  ar1=',ar1(idep),' ar2=',ar2(idep),' ar3=', &
!             ar3(idep)
      endif
         
      ar1(idep)=ar1(idep)/test
      ar2(idep)=ar2(idep)/test
      ar3(idep)=ar3(idep)/test
      if (ar2(idep) .ne. 0.) swsrf2(idep)=swsrf2(idep)/ar2(idep)
      if (ar3(idep) .ne. 0.) swsrf3(idep)=swsrf3(idep)/ar3(idep)

      rzeq(idep)=ar1(idep)+ar2(idep)*swsrf2(idep)+ar3(idep)*swsrf3(idep)
      
      if (bug) write(*,*) 'funcidep: ok6'

!c 2. Maximum density -> shape parameter 
!c                    -> wmin 

      locmax=3
      shift=15
      ordref=1
      do n=1,nmax
         densaux2(n)=denstest(idep,n)
      enddo
         
      if (nmax .ge. shift*2) then
               
!c we start with sliding mean to facilitate the search for the maximum
        
         ord=MIN(ordref,nmax/shift)
	 
         call smtot(densaux2,nmax,ord,densaux) 
!	 print *,nmax,ord,shift,densaux(shift-14),shift-14,size(densaux)
         do n=nmax,shift,-1
            if (densaux(n) .gt. densaux(n-1) .and.     &
                densaux(n) .gt. densaux(n-2) .and.     &
                densaux(n) .gt. densaux(n-3) .and.     &
                densaux(n) .gt. densaux(n-4) .and.     &
                densaux(n) .gt. densaux(n-5) .and.     &
                densaux(n) .gt. densaux(n-6) .and.     &
                densaux(n) .gt. densaux(n-7) .and.     &
                densaux(n) .gt. densaux(n-8) .and.     &
                densaux(n) .gt. densaux(n-9) .and.     &
                densaux(n) .gt. densaux(n-10) .and.    &
                densaux(n) .gt. densaux(n-11) .and.    &
                densaux(n) .gt. densaux(n-12) .and.    &
                densaux(n) .gt. densaux(n-13) .and.    &
                densaux(n) .gt. densaux(n-14))then ! .and.    &
!                densaux(n) .gt. densaux(n-15)) then
               locmax=n
               goto 30
            endif
         enddo
         
      else

         aux10=-9999.
         indimax10=3
         do n=1,nmax
            if (densaux2(n) .gt. aux10) then
               aux10=densaux2(n)
               indimax10=n
            endif
         enddo
         locmax=MAX(3,indimax10)

      endif      ! if (nmax .ge. shift+1) 
                  
 30   densmax=denstest(idep,locmax)
      aa(idep)=exp(1.)*densmax

      if (bug) write(*,*) 'funcidep: ok7'

!c WMIN=lowest value where the density is strictly gt densmax/100.

      indmin=1
      indmin0=0
      do n=1,nmax
         if (denstest(idep,n) .gt. 0.) indmin0=n
         if (denstest(idep,n) .gt. densmax/100. .and.    &
             valtest(idep,n) .lt. valtest(idep,locmax)) indmin=n
      enddo
      if (indmin .eq.0) indmin=indmin0

      if (indmin .le. 2) then
         wmin(idep) = 0.99999
      else
         x1=valtest(idep,indmin)
         wmin(idep)=x1
      endif

      if (bug) write(*,*) 'funcidep: ok8; first wmin=',wmin(idep)

!c for negative or low coeskew the previous wmin doesn't give good results...
!c wmin is higher !!!

      if (coeskew .lt. 1. ) then

         if (locmax .gt. 3 .and. indmin .ge. locmax+4) then
            n2=MAX(locmax+1,(indmin-locmax)/2+locmax)
            x2=valtest(idep,n2)
            y2=denstest(idep,n2)
            n1=locmax
            x1=valtest(idep,n1)
            y1=denstest(idep,n1)
            wa=(y2-y1)/(x2-x1)
            wb=y1-wa*x1
            wmin(idep)=AMAX1(wmin(idep),-wb/wa)
         endif

!c wmin is even higher in some cases !!!
         if (coeskew .lt. 0.2 ) wmin(idep)=wmin(idep)+0.01
         
      endif  

      if (bug) write(*,*) 'funcidep: ok9; 2nd wmin=',wmin(idep)
      
      if (valtest(idep,locmax) .le. wmin(idep)) then ! doesn't make sense
         wmin(idep)=valtest(idep,locmax)-dx
      endif
      aabis(idep)=1./(valtest(idep,locmax)-wmin(idep)+1.e-20)

      if (bug) write(*,*) 'funcidep: ok10'

    END SUBROUTINE FUNCIDEP
  
  ! ********************************************************************

      SUBROUTINE FUNCZBAR(                                     &   
                         NAR0,ZBAR,                            &
                         BEE,PSIS,POROS,COND,RZDEP,WPWET,      &
                         VALX,PX,COESKEW,TIMEAN,SUMA,          &
                         CATDEF,WMIN)                          

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                                      c
!c This program returns the eight parameters for the areal fractioning  c
!c                                                                      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none                   
      INTEGER , intent (in) :: NAR0
      integer nref,nind,nmax,indmin,locmax,shift,ord,locmin,ordref
      integer indimax10,indmin0
      REAL, intent (in) :: BEE, PSIS, POROS, COND, RZDEP, WPWET, COESKEW
      REAL, intent (inout) ::  VALX(NAR), PX(NAR),TIMEAN,SUMA,zbar
      real, intent (inout) :: catdef,wmin

      REAL  dx,dz,sumdef
      real term1,term2
      real zdep(nar),locdef(nar),wrz(nar),frcunsat
      real valtest(nar),ptest(nar),denstest(nar),dtest(nar)
      real x1,x2,y1,y2,wa,wb
      integer n1,n2,k,n
      real densaux(nar),densaux2(nar),densmax,aux10

!c-------------------------------------------------------------------------
!c integral(f(x)dx)=1. for a pdf
!c here px=f(x)dx
      dx=valx(1)-valx(2)
            
!c**** Compute array of water table depths:
      do k=1,nar0
         term1=(1/gnu)*(valx(k)-timean)
         zdep(k)=AMAX1(0.,zbar-term1)
      enddo

!c variable change must be reflected in dx
      dz=dx/gnu
      
!c**** Compute array of moisture deficits:
      do k=1,nar0
         term1=(psis-zdep(k))/psis
         term1=term1**(1.-1./bee)
         term2=-psis*(bee/(bee-1.))*(term1-1.)
         locdef(k)=zdep(k)-term2
      enddo

!c**** Add deficits to produce catdef:
      sumdef=0.
      do k=1,nar0
         sumdef=sumdef+locdef(k)*px(k)
      enddo
      catdef=poros*1000.*sumdef/suma

!c**** Compute array of root zone moisture (degree of wetness in root zone):
      do k=1,nar0              
         term1=((psis-zdep(k))/psis)  &
             **(1.-1./bee)
         if(zdep(k).le.0.) then
            wrz(k)=1.
         elseif(zdep(k)-rzdep.lt.0.) then
            term2=(-psis/zdep(k))*(bee/(bee-1.))   &
                *(term1-1.)
            frcunsat=zdep(k)/rzdep
            wrz(k)=frcunsat*term2+(1.-frcunsat)*1.
         else
            term2=((psis-zdep(k)+rzdep)     &
                /psis)**(1.-1./bee)
            wrz(k)=(-psis/rzdep)*(bee/      &
                (bee-1.))*(term1-term2)
         endif
      enddo

!c**** compute the densities and dx
!c**** we use a usefull property that is due to the construction of the 
!c**** gamma distribution in tgen3.f : this distribution is continuous, 
!c**** with decreasing values on ln(a/tanb) when n goes from 1 to nar0
!c first we gather in the same bin all the bins with values ge 1 
      nref=1
      nind=1
      ptest(1)=0.
      do k=1,nar0
         if (wrz(k) .eq. 1.) then
            nref=nref+1
            ptest(1) = ptest(1) + px(k)
         endif
      enddo
      if (nref .gt. 1) then
         nind=2
         valtest(1)=1.
      endif
      nmax=nar0-nref+nind
      
!c definition of the probabilities ptest
      if (nmax .eq. 1) then     ! all the bins have values ge 1 
         dtest(1) = 0.0001
         ptest(1) = 1.
      else                      ! distribution in ar2/ar3
         do n=0,nmax-nind
            valtest(nind+n)=wrz(nref+n)
            ptest(nind+n)=px(nref+n)
         enddo
         
!c we have to define dtest, the size of each bin
         if (nmax .eq. 2) then
            dtest(2) = valtest(1)-valtest(2)
            dtest(1) = dtest(2)/2.
         else                   ! nmax .gt. 2
            do n=2,nmax-1
               dtest(n)=(valtest(n-1)-valtest(n+1))/2.            
            enddo             
            dtest(1) = dtest(2)/2.
            dtest(nmax) = dtest(nmax-1)
         endif
      endif

!c we can now define the probability density: denstest=ptest/dtest
!c where ptest is the probability and dtest the size of the bin
      do n=1,nmax
         if (ptest(n) .eq. 0.) then
            denstest(n)=0.
         else
            denstest(n)=ptest(n)/dtest(n)
         endif
      enddo

!c NOW we can estimate the parameters for the approximated distrib
!c from the actual distrib

!c 2. Maximum density -> shape parameter 
!c                    -> wmin 

      locmax=3
      shift=15
      ordref=1
      do n=1,nmax
         densaux2(n)=denstest(n)
      enddo

      if (nmax .ge. shift*2) then
               
!c we start with sliding mean to facilitate the search for the maximum
         
         ord=MIN(ordref,nmax/shift)
         call smtot(densaux2,nmax,ord,densaux)

         do n=nmax,shift,-1
            if (densaux(n) .gt. densaux(n-1) .and.         &
                densaux(n) .gt. densaux(n-2) .and.         &
                densaux(n) .gt. densaux(n-3) .and.         &
                densaux(n) .gt. densaux(n-4) .and.         &
                densaux(n) .gt. densaux(n-5) .and.         &
                densaux(n) .gt. densaux(n-6) .and.         &
                densaux(n) .gt. densaux(n-7) .and.         &
                densaux(n) .gt. densaux(n-8) .and.         &
                densaux(n) .gt. densaux(n-9) .and.         &
                densaux(n) .gt. densaux(n-10) .and.        &
                densaux(n) .gt. densaux(n-11) .and.        &
                densaux(n) .gt. densaux(n-12) .and.        &
                densaux(n) .gt. densaux(n-13) .and.        &
                densaux(n) .gt. densaux(n-14) .and.        &
                densaux(n) .gt. densaux(n-15)) then
               locmax=n
               goto 30
            endif
         enddo

      else

         aux10=-9999.
         indimax10=3
         do n=1,nmax
            if (densaux2(n) .gt. aux10) then
               aux10=densaux2(n)
               indimax10=n
            endif
         enddo
         locmax=MAX(3,indimax10)

      endif      ! if (nmax .ge. shift+1) 
         
 30   densmax=denstest(locmax)

!c WMIN=lowest value where the density is strictly gt densmax/100.

      indmin=1
      indmin0=0
      do n=1,nmax
         if (denstest(n) .gt. 0.) indmin0=n
         if (denstest(n) .gt. densmax/100. .and.         &
             valtest(n) .lt. valtest(locmax)) indmin=n
      enddo
      if (indmin .eq. 0) indmin=indmin0

      if (indmin .le. 2) then
         wmin = 0.99999
      else
         x1=valtest(indmin)
         wmin=x1
      endif

!c for negative or low coeskew the previous wmin doesn't give good results...
!c wmin is higher !!!

      if (coeskew .lt. 1. ) then
         
         if (locmax .gt. 3 .and. indmin .ge. locmax+4) then
            
            n2=MAX(locmax+1,(indmin-locmax)/2+locmax)
            x2=valtest(n2)
            y2=denstest(n2)
            n1=locmax
            x1=valtest(n1)
            y1=denstest(n1)
            wa=(y2-y1)/(x2-x1)
            wb=y1-wa*x1
            wmin=AMAX1(wmin,-wb/wa)
         endif

!c wmin is even higher in some cases !!!
         if (coeskew .lt. 0.2 ) wmin=wmin+0.01
         
      endif  

    END SUBROUTINE FUNCZBAR

! ******************************************************************

       SUBROUTINE RMSE(XX,YY,LEN,ERROR)

!c---------------------------------------------------------------------------
!c Computes the root-mean square error ERROR between two one-dimensional
!c random variables XX and YY of same length LEN
!c---------------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, intent (in) :: LEN
      REAL, intent (in) ::  XX(LEN),YY(LEN)
      REAL, intent (out) :: ERROR
      INTEGER :: I

!c---------------------------------------------------------------------------     
      error=0.
      do i=1,len
         if(abs(xx(i)-yy(i)) >=1.e-10) then
         error=error+(xx(i)-yy(i))*(xx(i)-yy(i))
         endif
      enddo
      error=SQRT(error/float(len))

    END SUBROUTINE RMSE

! ******************************************************************
       SUBROUTINE SMTOT(XX,LEN,ORD,YY)

!c---------------------------------------------------------------------------
!c Runs a sliding average of order ORD through the one-dimensional array XX
!c of length LEN and returns the smoothed YY
!!c---------------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, intent (in) ::  LEN
      INTEGER :: ORD,WIDTH,i,ini,n,end
      REAL, intent (in) :: XX(NAR)
      REAL, intent (out) :: YY(NAR)

!c---------------------------------------------------------------------------     
      do i=1,nar
         yy(i)=0.
      enddo

      width=ord*2+1
      if (width .gt. len/2) then
         write(*,*) 'the order for the sliding average is too large !!!'
         write(*,*) 'regard with the length of the array to be smoothed'
         stop
      endif

      do i=1,len
         ini=MAX(1,i-ord)
         end=MIN(len,i+ord)
         yy(i)=0.
         do n=ini,end
            yy(i)=yy(i)+xx(n)
         enddo
         yy(i)=yy(i)/(end-ini+1)
      enddo

    END SUBROUTINE SMTOT

! -----------------------------------------------------------------------------------

subroutine RegridRaster(Rin,Rout)

  integer, intent(IN)  :: Rin(:,:)
  integer, intent(OUT) :: Rout(:,:)

  REAL_  :: xx, yy
  integer :: i,j,ii,jj

  xx = size(Rin ,1)/float(size(Rout,1))
  yy = size(Rin ,2)/float(size(Rout,2))

  do j=1,size(Rout,2)
     jj = (j-1)*yy + 1
     do i=1,size(Rout,1)
        ii = (i-1)*xx + 1
        Rout(i,j) = Rin(ii,jj)
     end do
  end do

end subroutine RegridRaster

! -----------------------------------------------------------------------------------

subroutine RegridRaster1(Rin,Rout)

  integer*1, intent(IN)  :: Rin(:,:)
  integer*1, intent(OUT) :: Rout(:,:)

  REAL_  :: xx, yy
  integer :: i,j,ii,jj

  xx = size(Rin ,1)/float(size(Rout,1))
  yy = size(Rin ,2)/float(size(Rout,2))

  do j=1,size(Rout,2)
     jj = (j-1)*yy + 1
     do i=1,size(Rout,1)
        ii = (i-1)*xx + 1
        Rout(i,j) = Rin(ii,jj)
     end do
  end do

end subroutine RegridRaster1


! -----------------------------------------------------------------------------------

subroutine RegridRaster2(Rin,Rout)

  integer(kind=2), intent(IN)  :: Rin(:,:)
  integer(kind=2), intent(OUT) :: Rout(:,:)

  REAL_  :: xx, yy
  integer :: i,j,ii,jj

  xx = size(Rin ,1)/float(size(Rout,1))
  yy = size(Rin ,2)/float(size(Rout,2))

  do j=1,size(Rout,2)
     jj = (j-1)*yy + 1
     do i=1,size(Rout,1)
        ii = (i-1)*xx + 1
        Rout(i,j) = Rin(ii,jj)
     end do
  end do

end subroutine RegridRaster2


! -----------------------------------------------------------------------------------

subroutine RegridRasterReal(Rin,Rout)

  real, intent(IN)  :: Rin(:,:)
  real, intent(OUT) :: Rout(:,:)

  REAL_  :: xx, yy
  integer :: i,j,ii,jj

  xx = size(Rin ,1)/float(size(Rout,1))
  yy = size(Rin ,2)/float(size(Rout,2))

  do j=1,size(Rout,2)
     jj = (j-1)*yy + 1
     do i=1,size(Rout,1)
        ii = (i-1)*xx + 1
        Rout(i,j) = Rin(ii,jj)
     end do
  end do

end subroutine RegridRasterReal

!---------------------------------------------------------------------

      SUBROUTINE svbksb(u,w,v,m,n,b,x) 
        implicit none
        INTEGER m,mp,n,np,NMAX 
        REAL*8 b(m),u(m,n),v(n,n),w(n),x(n) 
        PARAMETER (NMAX=500)  !Maximum anticipated value of n
        !------------------------------------------------------------------------------------------- 
        ! Solves A · X = B for a vector X, where A is specified by the arrays u, w, v as returned by 
        ! svdcmp. m and n are the dimensions of a, and will be equal for square matrices. b(1:m) is 
        ! the input right-hand side. x(1:n) is the output solution vector. No input quantities are 
        ! destroyed, so the routine may be called sequentially with different b’s. 
        !-------------------------------------------------------------------------------------------

        INTEGER i,j,jj 
        REAL*8 s,tmp(NMAX) 
        do j=1,n !Calculate UTB. 
           s=0. 
           if(w(j).ne.0.)then !Nonzero result only if wj is nonzero. 
              do i=1,m 
                 s=s+u(i,j)*b(i) 
              end do
              s=s/(w(j) + 1.d-20) !This is the divide by wj . 
           endif
           tmp(j)=s 
        end do
        do j=1,n !Matrix multiply by V to get answer. 
           s=0. 
           do jj=1,n 
              s=s+v(j,jj)*tmp(jj) 
           end do
           x(j)=s 
        end do
        return 
      END SUBROUTINE svbksb

!---------------------------------------------------------------------

      SUBROUTINE svdcmp(a,m,n,w,v) 
        implicit none
        INTEGER m,n,NMAX 
        REAL*8, intent (inout)  :: a(m,n)
        REAL*8, intent (out) :: v(n,n),w(n) 
        PARAMETER (NMAX=500)  !Maximum anticipated value of n. 
        !-------------------------------------------------------------------------------------- 
        ! Given a matrix A(1:m,1:n), this routine computes its singular value decomposition, 
        ! A = U · W · Vt. The matrix U replaces A on output. The diagonal matrix of singular 
        ! values W is output as a vector W(1:n). The matrix V (not the transpose Vt) is output 
        ! as V(1:n,1:n). 
        !--------------------------------------------------------------------------------------

        INTEGER i,its,j,jj,k,l,nm 
        REAL*8 anorm,c,f,g,h,s,scale,x,y,z,rv1(NMAX) 
        real*8, parameter :: EPS=epsilon(1.0d0)
        g=0.d0  !Householder reduction to bidiagonal form. 
        scale=0.d0 
        anorm=0.d0 
        c =0.d0 
        f =0.d0 
        g =0.d0 
        h =0.d0 
        s =0.d0 
        x =0.d0 
        y =0.d0 
        z =0.d0 
        rv1=0.d0 
        w = 0.d0 
        v = 0.d0 
        do i=1,n 
           l=i+1 
           rv1(i)=scale*g 
           g=0.d0 
           s=0.d0 
           scale=0.d0 
           if(i.le.m)then 
              do k=i,m 
                 scale=scale+abs(a(k,i)) 
              end do
              if(scale.ne.0.d0)then 
                 do k=i,m 
                    a(k,i)=a(k,i)/scale 
                    s=s+a(k,i)*a(k,i) 
                 end do
                 f=a(i,i) 
                 g=-dsign(dsqrt(s),f) 
                 h=f*g-s 
                 a(i,i)=f-g 
                 do j=l,n 
                    s=0.d0 
                    do k=i,m 
                       s=s+a(k,i)*a(k,j) 
                    end do
                    f=s/h 
                    do k=i,m 
                       a(k,j)=a(k,j)+f*a(k,i) 
                    end do
                 end do
                 do k=i,m 
                    a(k,i)=scale*a(k,i) 
                 end do
              endif
           endif
           w(i)=scale *g 
           g=0.d0 
           s=0.d0 
           scale=0.d0 
           if((i.le.m).and.(i.ne.n))then 
              do k=l,n 
                 scale=scale+abs(a(i,k)) 
              end do
              if(scale.ne.0.d0)then 
                 do k=l,n 
                    a(i,k)=a(i,k)/scale 
                    s=s+a(i,k)*a(i,k) 
                 end do
                 f=a(i,l) 
                 g=-sign(sqrt(s),f) 
                 h=f*g-s 
                 a(i,l)=f-g 
                 do k=l,n 
                    rv1(k)=a(i,k)/h 
                 end do
                 do j=l,m 
                    s=0.d0 
                    do k=l,n 
                       s=s+a(j,k)*a(i,k) 
                    end do
                    do k=l,n 
                       a(j,k)=a(j,k)+s*rv1(k) 
                    end do
                 end do
                 do k=l,n 
                    a(i,k)=scale*a(i,k) 
                 end do
              endif
           endif
           anorm=max(anorm,(abs(w(i))+abs(rv1(i)))) 
        end do !do i=1,n
        
        do i=n,1,-1 !Accumulation of right-hand transformations. 
           if(i.lt.n)then 
              if(g.ne.0.d0)then 
                 do j=l,n       !Double division to avoid possible underflow. 
                    v(j,i)=(a(i,j)/a(i,l))/g 
                 end do
                 do j=l,n 
                    s=0.d0 
                    do k=l,n 
                       s=s+a(i,k)*v(k,j) 
                    end do
                    do k=l,n 
                       v(k,j)=v(k,j)+s*v(k,i) 
                    end do
                 end do
              endif
              do j=l,n 
                 v(i,j)=0.d0 
                 v(j,i)=0.d0 
              end do
           endif
           v(i,i)=1.d0
           g=rv1(i) 
           l=i 
        end do
        
        do i=min(m,n),1,-1 !Accumulation of left-hand transformations. 
           l=i+1 
           g=w(i) 
           do j=l,n 
              a(i,j)=0.d0 
           end do
           if(g.ne.0.d0)then 
              g=1.d0/g 
              do j=l,n 
                 s=0.d0 
                 do k=l,m 
                    s=s+a(k,i)*a(k,j) 
                 end do
                 f=(s/a(i,i))*g 
                 do k=i,m 
                    a(k,j)=a(k,j)+f*a(k,i) 
                 end do
              end do
              do j=i,m 
                 a(j,i)=a(j,i)*g 
              end do
           else
              do j= i,m 
                 a(j,i)=0.d0 
              end do
           endif
           a(i,i)=a(i,i)+1.d0 
        end do
        
        do k=n,1,-1 !Diagonalization of the bidiagonal form: Loop over 
           !singular values, and over allowed iterations. 
           do its=1,30 
              do l=k,1,-1 !Test for splitting. 
                 nm=l-1 !Note that rv1(1) is always zero.
                 if( abs(rv1(l)) <= EPS*anorm ) goto 2 
                 if( abs(w(nm) ) <= EPS*anorm ) goto 1  
              end do
1             c=0.d0 !Cancellation of rv1(l), if l > 1. 
              s=1.d0 
              do i=l,k 
                 f=s*rv1(i) 
                 rv1(i)=c*rv1(i) 
                 if( abs(f) <= EPS*anorm ) goto 2 
                 g=w(i) 
                 h=pythag(f,g) 
                 w(i)=h 
                 h=1.d0/h 
                 c= (g*h) 
                 s=-(f*h) 
                 do j=1,m 
                    y=a(j,nm) 
                    z=a(j,i) 
                    a(j,nm)=(y*c)+(z*s) 
                    a(j,i)=-(y*s)+(z*c) 
                 end do
              end do
2             z=w(k) 
              if(l.eq.k)then   !Convergence. 
                 if(z.lt.0.d0)then !Singular value is made nonnegative. 
                    w(k)=-z 
                    do j=1,n 
                       v(j,k)=-v(j,k) 
                    end do
                 endif
                 goto 3 
              endif
              if(its.eq.30) print *, 'no convergence in svdcmp' 
 !             if(its.ge.4)  print *, 'its = ',its
              x=w(l) !Shift from bottom 2-by-2 minor. 
              nm=k-1 
              y=w(nm) 
              g=rv1(nm) 
              h=rv1(k) 
              f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y) 
              g=pythag(f,1.d0) 
              f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x 
              c=1.d0 !Next QR transformation: 
              s=1.d0 
              do j=l,nm 
                 i=j+1 
                 g=rv1(i) 
                 y=w(i) 
                 h=s*g 
                 g=c*g 
                 z=pythag(f,h) 
                 rv1(j)=z 
                 c=f/z 
                 s=h/z 
                 f= (x*c)+(g*s) 
                 g=-(x*s)+(g*c) 
                 h=y*s 
                 y=y*c 
                 do jj=1,n 
                    x=v(jj,j) 
                    z=v(jj,i) 
                    v(jj,j)= (x*c)+(z*s) 
                    v(jj,i)=-(x*s)+(z*c) 
                 end do
                 z=pythag(f,h) 
                 w(j)=z !Rotation can be arbitrary if z = 0. 
                 if(z.ne.0.d0)then 
                    z=1.d0/z 
                    c=f*z 
                    s=h*z 
                 endif
                 f= (c*g)+(s*y) 
                 x=-(s*g)+(c*y) 
                 do jj=1,m 
                    y=a(jj,j) 
                    z=a(jj,i) 
                    a(jj,j)= (y*c)+(z*s) 
                    a(jj,i)=-(y*s)+(z*c) 
                 end do
              end do !j=l;nm 
              rv1(l)=0.d0 
              rv1(k)=f 
              w(k)=x 
           end do !its=1,30
3          continue 
        end do !k=n,1,-1 
        return 
      END SUBROUTINE svdcmp
!
! ________________________________________________________________________________
!     
      REAL*8 FUNCTION pythag(a,b) 
        REAL*8 a,b 
        !Computes sqrt(a**2 + b**2) without destructive underflow or overflow.
        REAL*8 absa,absb 
        absa=abs(a) 
        absb=abs(b) 
        if(absa.gt.absb)then 
           pythag=absa*sqrt(1.+(absb/absa)**2) 
        else 
           if(absb.eq.0.)then 
              pythag=0. 
           else
              pythag=absb*sqrt(1.+(absa/absb)**2) 
           endif
        endif
        return 
      END FUNCTION pythag
!
! ________________________________________________________________________________
!     

      SUBROUTINE savgol(c,np,nl,nr,ld,m)
      implicit none
      INTEGER ld,m,nl,np,nr,MMAX 
      real c(np) 
      PARAMETER (MMAX=6)
!-------------------------------------------------------------------------------------------- 
!USES lubksb,ludcmp given below. 
!Returns in c(1:np), in wrap-around order (see reference) consistent with the argument respns 
!in routine convlv, a set of Savitzky-Golay filter coefficients. nl is the number of leftward 
!(past) data points used, while nr is the number of rightward (future) data points, making 
!the total number of data points used nl +nr+1. ld is the order of the derivative desired 
!(e.g., ld = 0 for smoothed function). m is the order of the smoothing polynomial, also 
!equal to the highest conserved moment; usual values are m = 2 or m = 4. 
!--------------------------------------------------------------------------------------------
INTEGER d,icode,imj,ipj,j,k,kk,mm,indx(MMAX+1) 
real fac,sum,a(MMAX+1,MMAX+1),b(MMAX+1)
if(np.lt.nl+nr+1.or.nl.lt.0.or.nr.lt.0.or.ld.gt.m.or.m.gt.MMAX  & 
	  .or.nl+nr.lt.m) pause ' Bad args in savgol.' 
	do ipj=0,2*m        !Set up the normal equations of the desired leastsquares fit. 
    	sum=0. 
    if(ipj.eq.0) sum=1. 
    do k=1,nr 
      sum=sum+dfloat(k)**ipj 
    end do 
    do k=1,nl 
      sum=sum+dfloat(-k)**ipj 
    end do 
    mm=min(ipj,2*m-ipj) 
    do imj=-mm,mm,2 
      a(1+(ipj+imj)/2,1+(ipj-imj)/2)=sum 
    end do 
  end do

  call ludcmp(a,m+1,MMAX+1,indx,d,icode)    !Solve them: LU decomposition. 

  do j=1,m+1 
    b(j)=0. 
  end do 
  b(ld+1)=1.      !Right-hand side vector is unit vector, depending on which derivative we want. 

  call lubksb(a,m+1,MMAX+1,indx,b)   !Backsubstitute, giving one row of the inverse matrix. 

  do kk=1,np                         !Zero the output array (it may be bigger than the number 
    c(kk)=0.                         !of coefficients).  
  end do 
  do k=-nl,nr                        !Each Savitzky-Golay coefficient is the dot product 
    sum=b(1)                         !of powers of an integer with the inverse matrix row. 
    fac=1. 
    do mm=1,m 
      fac=fac*k 
      sum=sum+b(mm+1)*fac 
    end do 
    kk=mod(np-k,np)+1                !Store in wrap-around order. 
    c(kk)=sum 
  end do
  return 
END SUBROUTINE savgol

!***************************************************************
!* Given an N x N matrix A, this routine replaces it by the LU *
!* decomposition of a rowwise permutation of itself. A and N   *
!* are input. INDX is an output vector which records the row   *
!* permutation effected by the partial pivoting; D is output   *
!* as -1 or 1, depending on whether the number of row inter-   *
!* changes was even or odd, respectively. This routine is used *
!* in combination with LUBKSB to solve linear equations or to  *
!* invert a matrix. Return code is 1, if matrix is singular.   *
!***************************************************************
 Subroutine LUDCMP(A,N,NP,INDX,D,CODE)
INTEGER, PARAMETER :: NMAX=100
REAL, PARAMETER :: TINY=1E-12
 real  AMAX,DUM, SUM, A(NP,NP),VV(NMAX)
 INTEGER CODE, D, INDX(N),NP,N,I,J,K,IMAX

 D=1; CODE=0

 DO I=1,N
   AMAX=0.
   DO J=1,N
     IF (ABS(A(I,J)).GT.AMAX) AMAX=ABS(A(I,J))
   END DO ! j loop
   IF(AMAX.LT.TINY) THEN
     CODE = 1
     RETURN
   END IF
   VV(I) = 1. / AMAX
 END DO ! i loop

 DO J=1,N
   DO I=1,J-1
     SUM = A(I,J)
     DO K=1,I-1
       SUM = SUM - A(I,K)*A(K,J) 
     END DO ! k loop
     A(I,J) = SUM
   END DO ! i loop
   AMAX = 0.
   DO I=J,N
     SUM = A(I,J)
     DO K=1,J-1
       SUM = SUM - A(I,K)*A(K,J) 
     END DO ! k loop
     A(I,J) = SUM
     DUM = VV(I)*ABS(SUM)
     IF(DUM.GE.AMAX) THEN
       IMAX = I
       AMAX = DUM
     END IF
   END DO ! i loop  
   
   IF(J.NE.IMAX) THEN
     DO K=1,N
       DUM = A(IMAX,K)
       A(IMAX,K) = A(J,K)
       A(J,K) = DUM
     END DO ! k loop
     D = -D
     VV(IMAX) = VV(J)
   END IF

   INDX(J) = IMAX
   IF(ABS(A(J,J)) < TINY) A(J,J) = TINY

   IF(J.NE.N) THEN
     DUM = 1. / A(J,J)
     DO I=J+1,N
       A(I,J) = A(I,J)*DUM
     END DO ! i loop
   END IF 
 END DO ! j loop

 RETURN
 END Subroutine LUDCMP


!******************************************************************
!* Solves the set of N linear equations A . X = B.  Here A is     *
!* input, not as the matrix A but rather as its LU decomposition, *
!* determined by the routine LUDCMP. INDX is input as the permuta-*
!* tion vector returned by LUDCMP. B is input as the right-hand   *
!* side vector B, and returns with the solution vector X. A, N and*
!* INDX are not modified by this routine and can be used for suc- *
!* cessive calls with different right-hand sides. This routine is *
!* also efficient for plain matrix inversion.                     *
!******************************************************************
 Subroutine LUBKSB(A,N,NP,INDX,B)
 INTEGER :: II,I,J,LL,N,NP
 real  SUM, A(NP,NP),B(N)
 INTEGER INDX(N)

 II = 0

 DO I=1,N
   LL = INDX(I)
   SUM = B(LL)
   B(LL) = B(I)
   IF(II.NE.0) THEN
     DO J=II,I-1
       SUM = SUM - A(I,J)*B(J)
     END DO ! j loop
   ELSE IF(SUM.NE.0.) THEN
     II = I
   END IF
   B(I) = SUM
 END DO ! i loop

 DO I=N,1,-1
   SUM = B(I)
   IF(I < N) THEN
     DO J=I+1,N
       SUM = SUM - A(I,J)*B(J)
     END DO ! j loop
   END IF
   B(I) = SUM / A(I,I)
 END DO ! i loop

 RETURN
 END Subroutine LUBKSB

!
! ====================================================================
!

INTEGER FUNCTION center_pix (x,y,x0,y0,z0,ext_point)

implicit none

real, dimension (:), intent (in) :: x,y
real, allocatable, dimension (:,:) :: length_m
real, allocatable, dimension (:) :: length
real, intent (inout) :: x0,y0,z0
integer :: i,j,npix,ii
logical, intent(in) :: ext_point
real :: zi, zj

npix = size (x)
allocate (length_m (1:npix,1:npix))
allocate (length   (1:npix))
length_m =0.
length   =0.

do i = 1,npix
   zi = 100. - x(i) - y(i)
   if (.not. ext_point) then
      x0 = x(i)
      y0 = y(i)
      z0 = zi
   endif

   do j = i,npix
      zj = 100. - x(j) - y(j)
!      length_m (i,j) = abs (x(j) - x0) + &
!            abs (y(j) - y0) +  abs (zj - z0)
!
      length_m (i,j) = ((x(j) - x0)*(x(j) - x0) &
                     +  (y(j) - y0)*(y(j) - y0) &
                     +  (zj - z0)*(zj - z0))**0.5
      length_m (j,i) = length_m (i,j)
   end do
   length (i) = sum(length_m (i,:))
end do

center_pix = minloc(length,dim=1)

END FUNCTION center_pix

!
!----------------------------------------------------------
!

INTEGER FUNCTION soil_class (min_perc)

! Function returns a unique soil class [1-100], 

IMPLICIT NONE
type(mineral_perc), intent (in)  :: min_perc
!real, intent (in) :: clay_perc,silt_perc,sand_perc
integer :: clay_row, sand_row, silt_row

clay_row = ceiling((100.- min_perc%clay_perc)/10.)
if(clay_row == 0 ) clay_row = 1
if(clay_row == 11) clay_row = 10

sand_row = ceiling((min_perc%sand_perc)/10.)
if(sand_row == 0 ) sand_row = 1
if(sand_row == 11) sand_row = 10

silt_row = ceiling((min_perc%silt_perc)/10.)
if(silt_row == 0 ) silt_row = 1
if(silt_row == 11) silt_row = 10

if(clay_row == 1) soil_class=1

if(clay_row > 1) soil_class=   &
  (clay_row - 1)*(clay_row - 1) + (clay_row - sand_row) + silt_row

end FUNCTION soil_class

! -----------------------------------------------------------------------------------
! -----------------------------------------------------------------------------------
! -----------------------------------------------------------------------------------

SUBROUTINE REFORMAT_VEGFILES
  implicit none
  character*400 :: tmp_string
  integer :: n_tiles
  real, dimension (:), allocatable :: var_array
  character*40 :: header
  real :: a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14
  integer :: month

  tmp_string = 'mkdir -p '//'clsm/g5fmt'
  call system(tmp_string)
  tmp_string = '/bin/mv '//'clsm/lai.dat ' //'clsm/g5fmt/.'
  call system(tmp_string) 
  tmp_string = '/bin/mv '//'clsm/green.dat ' //'clsm/g5fmt/.'
  call system(tmp_string) 

  open (10,file='clsm/g5fmt/lai.dat'  , form = 'unformatted',   &
       convert='little_endian',status='old',action='read' )
  open (11,file='clsm/g5fmt/green.dat', form = 'unformatted',   &
       convert='little_endian',status='old',action='read' )

  open (20,file='clsm/lai.dat', form = 'unformatted',   &
       convert='big_endian',status='unknown',action='write' )
  open (21,file='clsm/green.dat', form = 'unformatted', &
       convert='big_endian',status='unknown',action='write' )

  open (30,file='clsm/catchment.def', form = 'formatted',status='old',action='read' )
  read (30,*) n_tiles
  close(30,status='keep')

  allocate (var_array (1:n_tiles))

  read (10) header
  read (10) var_array
  read (11) header
  read (11) var_array

  do month  =1,12

     read (10) a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14
     read (10) var_array
     print '(12f3.0,f4.00,f2.0,a6,2f6.2)',a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,'LAI ',minval(var_array),maxval(var_array)
     write (20)var_array(:) 
     read (11) a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14
     read (11) var_array
     print '(12f3.0,f4.00,f2.0,a6,2f6.2)',a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,'GREEN ',minval(var_array),maxval(var_array) 
     write (21)var_array(:) 
  end do

END SUBROUTINE REFORMAT_VEGFILES

!
! --------------------------------------------------------
!

SUBROUTINE compute_stats (ndata,cti_val,mu,sig,sk)

implicit none
integer, intent(in)                   :: ndata
real, intent(inout), dimension(ndata) :: cti_val
real, intent(out)                     :: mu,sig,sk
integer                               :: i,j
real                                  :: del

   mu = sum(cti_val(1:ndata))/float(ndata) 
   sig = 0.
   sk  = 0.
   del = 0.

   do i = 1,ndata
      del = CTI_VAL(i) - mu
      sig = sig + del**2
      sk  = sk  + (del*del*del)
   end do

   sig = sig/float(ndata-1)
   sig = sqrt(sig)
   sk  = sk/(sig**3 + 1.e-10)/float(ndata)
   
END SUBROUTINE compute_stats

! -----------------------------------------------------------------------------------
    
    SUBROUTINE ascat_r0 (nc,nr,gfiler, z0)
      
      implicit none

      ! 1) ASCAT roughness 
      ! /discover/nobackup/adarmeno/projects/k14/arlems-roughness.x3600_y1800_t1.nc4
     
      integer, intent (in)               :: nc, nr
      real, pointer, dimension (:), intent (inout) :: z0
      character(*), intent (in)          :: gfiler
      integer  , parameter               :: N_lon_ascat = 3600, N_lat_ascat = 1800
      integer                            :: i,j, status, varid, ncid
      integer                            :: NTILES
      REAL, ALLOCATABLE, dimension (:)   :: count_pix
      REAL, ALLOCATABLE, dimension (:,:) :: z0_grid, data_grid
      INTEGER, ALLOCATABLE, dimension (:,:) :: tile_id
      character*100                      :: fout

      ! Reading number of tiles
      ! -----------------------
      
      open (20, file = 'clsm/catchment.def', form = 'formatted', status = 'old', action =  'read')
      
      read (20, *) NTILES
      
      close (20, status = 'keep')

      ! READ ASCAT source data and regrid
      ! ---------------------------------

      status  = NF_OPEN ('data/CATCH/arlems-roughness.x3600_y1800_t1.nc4', NF_NOWRITE, ncid)
            
      allocate (z0_grid   (1 : NC         , 1 : NR))
      allocate (data_grid (1 : N_lon_ascat, 1 : N_lat_ascat)) 

      status  = NF_INQ_VARID (ncid,'roughness',VarID) ; VERIFY_(STATUS)
      status  = NF_GET_VARA_REAL (ncid,VarID, (/1,1,1/),(/N_lon_ascat, N_lat_ascat,1/), data_grid) ; VERIFY_(STATUS)

      call RegridRasterReal(data_grid, z0_grid)

      status = NF_CLOSE(ncid)

      ! Grid to tile
      ! ------------
                
      ! Reading tile-id raster file

      allocate(tile_id(1:nc,1:nr))
      
      open (10,file=trim(gfiler)//'.rst',status='old',action='read',  &
           form='unformatted',convert='little_endian')
      
      do j=1,nr
         read(10)tile_id(:,j)
      end do       

      close (10,status='keep')     
      
      allocate (z0        (1:NTILES))
      allocate (count_pix (1:NTILES))
      
      z0        = 0.
      count_pix = 0.

      do j = 1,nr
         do i = 1, nc
            if((tile_id(i,j).gt.0).and.(tile_id(i,j).le.NTILES)) then

               ! z0 0. < 0.1
               if((z0_grid(i,j) >= 2.0e-6).and.(z0_grid(i,j) <= 0.1)) then 
                 z0 (tile_id(i,j)) = z0 (tile_id(i,j)) + z0_grid(i,j)
                 count_pix (tile_id(i,j)) = count_pix (tile_id(i,j)) + 1. 
               endif

            endif
         end do
      end do
      
      where (count_pix > 0.) z0 = z0/count_pix
      where (z0 == 0.)       z0 = 2.0e-6
      
      deallocate (count_pix)
      deallocate (z0_grid)
      deallocate (tile_id)
      
  END SUBROUTINE ascat_r0

   ! ----------------------------------------------------------------------------------------------------------------------------

    SUBROUTINE jpl_canoph (nc,nr,gfiler, z2)

      implicit none

      ! 1) JPL Canopy Height 
      ! /discover/nobackup/rreichle/l_data/LandBCs_files_for_mkCatchParam/V001//Simard_Pinto_3DGlobalVeg_JGR.nc4
     
      integer, intent (in)               :: nc, nr
      real, pointer, dimension (:), intent (inout) :: z2
      character(*), intent (in)          :: gfiler
      integer  , parameter               :: N_lon_jpl = 43200, N_lat_jpl = 21600
      integer                            :: i,j, status, varid, ncid
      integer                            :: NTILES
      REAL, ALLOCATABLE, dimension (:)   :: count_pix
      INTEGER, ALLOCATABLE, dimension (:,:) :: data_grid, z2_grid
      INTEGER, ALLOCATABLE, dimension (:,:) :: tile_id
      character*100                      :: fout

      ! Reading number of tiles
      ! -----------------------
      
      open (20, file = 'clsm/catchment.def', form = 'formatted', status = 'old', action =  'read')
      
      read (20, *) NTILES
      
      close (20, status = 'keep')

      ! READ JPL source data files and regrid
      ! -------------------------------------

      status  = NF_OPEN ('data/CATCH/Simard_Pinto_3DGlobalVeg_JGR.nc4', NF_NOWRITE, ncid)
            
      allocate (z2_grid   (1 : NC         , 1 : NR))
      allocate (data_grid (1 : N_lon_jpl, 1 : N_lat_jpl)) 

      status  = NF_INQ_VARID (ncid,'CanopyHeight',VarID) ; VERIFY_(STATUS)
      status  = NF_GET_VARA_INT (ncid,VarID, (/1,1/),(/N_lon_jpl, N_lat_jpl/), data_grid) ; VERIFY_(STATUS)
      
      call RegridRaster(data_grid, z2_grid)

      status = NF_CLOSE(ncid)

      ! Grid to tile
      ! ------------
                
      ! Reading tile-id raster file

      allocate(tile_id(1:nc,1:nr))
      
      open (10,file=trim(gfiler)//'.rst',status='old',action='read',  &
           form='unformatted',convert='little_endian')
      
      do j=1,nr
         read(10)tile_id(:,j)
      end do       

      close (10,status='keep')     
      
      allocate (z2        (1:NTILES))
      allocate (count_pix (1:NTILES))
      
      z2        = 0.
      count_pix = 0.

      do j = 1,nr
         do i = 1, nc
            if((tile_id(i,j).gt.0).and.(tile_id(i,j).le.NTILES)) then

               if(z2_grid(i,j) >= 0.) then 
                 z2 (tile_id(i,j)) = z2 (tile_id(i,j)) + real (z2_grid(i,j))
                 count_pix (tile_id(i,j)) = count_pix (tile_id(i,j)) + 1. 
               endif

            endif
         end do
      end do
      
      where (count_pix >   0.) z2 = z2/count_pix
      where (z2        < 0.01) z2 = 0.01            ! to ensure Z2 >= MIN_VEG_HEIGHT

      deallocate (count_pix)
      deallocate (z2_grid)
      deallocate (tile_id)
      
    END SUBROUTINE jpl_canoph

! -----------------------------------------------------------------------------------
 
END module rmTinyCatchParaMod
