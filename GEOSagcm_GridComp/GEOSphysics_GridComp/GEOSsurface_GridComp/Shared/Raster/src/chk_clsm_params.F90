PROGRAM chk_clsm_params

!
! Usage = chk_clsm_params -s Y -m MaskFile
!  -s single tile (Y/N)
! This program may be handy when debugging problems in sat_param. 
! The program calls sunroutine sat_param and computes ars1,ars2,ars3,arw1,arw2,arw3,arw4
! it writes out arrays of catdef and theoritical values of ar1,wmin1, and computed
! ars1,ars2,ars3,arw1,arw2,arw3,arw4 for debugging
! DEFAULT: the program will process the tile in question, cti and soil data must be provided
! either through the command line or clsm/dbg_file (when you use totalview)
! After running the program, the next step is cd idl_out and ./plot_curves.csh
!

use rmTinyCatchParaMod
implicit none
character*1          :: opt
character*128        :: ARG
real :: TOPMEAN,TOPVAR,TOPSKEW,COESKEW
real :: cti_mean, cti_std, cti_min, cti_max, cti_skew, bee, &         
     psis, poros, cond, wpwet, soildepth
REAL, allocatable, dimension (:) :: ST, AC
real :: &
     ars1,ars2,ars3,                   &
     ara1,ara2,ara3,ara4,              &
     arw1,arw2,arw3,arw4,              &
     taberr1,taberr2,taberr3,taberr4,  &
     normerr1,normerr2,normerr3,normerr4
character*7 :: char_index
logical :: file_exist
integer, parameter :: NAR = 1000
real,dimension(5) :: ctis =  &
(/17.2355,   3.2612,   7.4278,  25.6578,  -1.5686/)
real,dimension(6) :: soils = &
(/3.3000, -0.0500, 0.3730,  0.00187934, 0.0885,  1741.241/)
integer :: ntiles,n,tindex1,pfaf1,tindex2,pfaf2,soil_class_top,soil_class_com 
character*100, parameter :: file2 = 'clsm/soil_param.dat', &
     file1='clsm/cti_stats.dat'
character*1 :: single_tile ='Y'
character*30 :: MaskFile
logical :: read_file  
character*300 :: input_data

read_file = .true.
call system ('mkdir -p idl_out')
call system ('cd bin/ ; /bin/cp ../src/plot_curves.csh . ; cd .. ; /bin/cp bin/plot_curves.csh idl_out/. ; chmod 755 idl_out/plot_curves.csh')


n = iargc()

if(n < 3) then
       print *, "Usage : chk_clsm_params -s Y(N) -m MaskFile"
       print *, "-s Y processing a single tile (DEFAULT); N processing every tile in tile space"
       print *, "-m GEOS5_10arcsec_mask.nc for new mask or just -m for older mask"
       call exit(1)
end if

n = 1
    call getarg(n,arg)
    do while(arg(1:1)=='-')
       opt=arg(2:2)
       if(len(trim(arg))==2) then
          if(scan(opt,'zvh')==0) then
             n = n + 1
             call getarg(n,arg)
          endif
       else
          arg = arg(3:)
       end if
       select case (opt)
       case ('s')
         single_tile  = trim(arg)
       case ('m')
         MaskFile  = trim(arg)
       end select
       n = n + 1
       call getarg(n,arg)
    end do

if((single_tile == 'Y').or.(single_tile == 'y')) read_file = .false.

ntiles=1

if(read_file) then

! Opening input data files
open (10,file=trim(file1), form = 'formatted',status='old', action ='read')
open (11,file=trim(file2), form = 'formatted',status='old', action ='read')
read (10,*) ntiles

else

inquire(file='clsm/dbg_file', exist=file_exist)

if(file_exist) then

open(20,file='clsm/dbg_file',form='formatted',status='old',action='read')
read (20,*)ctis(1),ctis(2),ctis(3),ctis(4),ctis(5)
read (20,*)soils(1),soils(2),soils(3),soils(4),soils(5),soils(6)
close (20,status='keep')

else

print *,'Copy and paste CTI_mean CTI_std CTI_min CTI_max CTI_skew'
print *,'Columns 3-7 in clsm/cti_stats.dat'
print *,'--------------------------------------------------------'
read (*,'(a)') input_data
read (input_data,*)ctis(1),ctis(2),ctis(3),ctis(4),ctis(5)


print *,'    '

print *,'Copy and paste BEE PSIS POROS COND WPWET SOILDEPTH'
print *,'Columns 5-10 in clsm/soil_param.dat'
print *,'--------------------------------------------------'
read (*,'(a)') input_data
read (input_data,*)soils(1),soils(2),soils(3),soils(4),soils(5),soils(6)

endif
print *,'Input MaskFile : ', trim(MaskFile)
print *,'INPUT Topography data:'
write (*,*)ctis(1),ctis(2),ctis(3),ctis(4),ctis(5)
print *,'INPUT soil data:'
write (*,*)soils(1),soils(2),soils(3),soils(4),soils(5),soils(6)

cti_mean  = ctis(1) 
cti_std   = ctis(2) 
cti_min   = ctis(3) 
cti_max   = ctis(4)  
cti_skew  = ctis(5) 
bee       = soils(1)       
psis      = soils(2)     
poros     = soils(3)     
cond      = soils(4)     
wpwet     = soils(5)     
soildepth = soils(6)     

endif

allocate(ST (1:NAR))
allocate(AC (1:NAR))

do n = 1, ntiles

if(read_file) then

   read(10,'(i8,i8,5(1x,f8.4))') tindex1,pfaf1, cti_mean, cti_std, &
        cti_min, cti_max, cti_skew
   read(11,*) tindex2,pfaf2,soil_class_top,soil_class_com,   &
        BEE, PSIS,POROS,COND,WPWET,soildepth   
         
   if(tindex1 /= tindex2) then
      write(*,*)'Warnning 1: pfafstetter mismatched' 
      stop
   endif

endif

write (char_index,'(i7.7)')n
open (20,file='idl_out/file.'//(char_index),form='formatted',status='unknown',action='write')
write (20,*) cti_mean, cti_std,cti_min, cti_max, cti_skew
write (20,*) BEE, PSIS,POROS,COND,WPWET,soildepth

if (index(MaskFile,'GEOS5_10arcsec_mask') /= 0) then
   TOPMEAN = cti_mean
else
   TOPMEAN = 0.961*cti_mean-1.957
endif
                       
TOPVAR  = cti_std*cti_std
TOPSKEW = cti_skew*cti_std*cti_std*cti_std

ST = 0.
AC = 0.

CALL TGEN (                        &  
     TOPMEAN,TOPVAR,TOPSKEW,       &
     ST,AC,COESKEW)

CALL SAT_PARAM(                        &
     BEE,PSIS,POROS,COND,              & 
     WPWET, ST, AC, COESKEW,1,         &
     soildepth,                        &
     ars1,ars2,ars3,                   &
     ara1,ara2,ara3,ara4,              &
     arw1,arw2,arw3,arw4,              &
     taberr1,taberr2,taberr3,taberr4,  &
     normerr1,normerr2,normerr3,normerr4,&
     DBG_UNIT = 20)

close (20,status = 'keep')

enddo

close (10,status = 'keep')
close (11,status = 'keep')

call system ('cd idl_out/ ; ./plot_curves.csh ; cd ..')

END PROGRAM chk_clsm_params
