#define I_AM_MAIN
#include "MAPL_Generic.h"

program ScaleCatch
  
  use MAPL
  use LSM_ROUTINES,      ONLY:          &
       catch_calc_soil_moist,           &
       catch_calc_tp,                   &
       catch_calc_ght
  
  USE CATCH_CONSTANTS,   ONLY:          &
       N_GT              => CATCH_N_GT, &
       DZGT              => CATCH_DZGT, &
       PEATCLSM_POROS_THRESHOLD

  use CatchmentRstMod
  use CatchmentCNRstMod

  implicit none

  character(256)    :: old_fname, new_fname, scale_fname, cnclm
  integer           :: ntiles, n, nargs
  integer           :: iargc
  real              :: SURFLAY        ! (Ganymed-3 and earlier) SURFLAY=20.0 for Old Soil Params
                                      ! (Ganymed-4 and later  ) SURFLAY=50.0 for New Soil Params
  real              :: WEMIN_IN, WEMIN_OUT
  character*256     :: arg(7)


  class(CatchmentRst), allocatable :: old_catch, new_catch, scale_catch

  real,    allocatable, dimension(:)   :: dzsf, ar1, ar2, ar4
  real,    allocatable, dimension(:,:) :: TP_IN, GHT_IN, FICE, GHT_OUT, TP_OUT
  real,    allocatable, dimension(:)   :: swe_in, depth_in, areasc_in, areasc_out, depth_out

  type(Netcdf4_fileformatter) :: formatter
  type(Filemetadata) :: meta
  integer :: i, rc, filetype
  integer :: status
  character(256) :: Iam = "ScaleCatch"
 
! Usage
! -----
  if (iargc() /= 6) then
     write(*,*) "Usage: ScaleCatch.x <Input_Catch> <Regridded_Catch> <Scaled_Catch> <SURFLAY> <WEMIN_IN> <WEMIN_OUT> <cnclm>"
     call exit(2)
  end if

  do n=1,7
    call getarg(n,arg(n))
  enddo

! Open INPUT and Regridded Catch Files
! ------------------------------------
  read(arg(1),'(a)') old_fname

  read(arg(2),'(a)') new_fname

! Open OUTPUT (Scaled) Catch File
! -------------------------------
  read(arg(3),'(a)') scale_fname
! Get SURFLAY Value
! -----------------
  read(arg(4),*) SURFLAY
  read(arg(5),*) WEMIN_IN
  read(arg(6),*) WEMIN_OUT
! catch or catchcn?
  read(arg(7),'(a)') cnclm

  if (index(cnclm,'40') /=0 .or. index(cnclm,'45') /=0 ) then
     allocate(old_catch, source = CatchmentCNRst(old_fname, cnclm))
     allocate(new_catch, source = CatchmentCNRst(new_fname, cnclm))
     allocate(scale_catch, source = new_catch)
  else
     allocate(old_catch, source = CatchmentRst(old_fname))
     allocate(new_catch, source = CatchmentRst(new_fname))
     allocate(scale_catch, source = new_catch)
  endif

  if (SURFLAY.ne.20 .and. SURFLAY.ne.50) then
     print *, "You must supply a valid SURFLAY value:"
     print *, "(Ganymed-3 and earlier) SURFLAY=20.0 for Old Soil Params"
     print *, "(Ganymed-4 and later  ) SURFLAY=50.0 for New Soil Params"
     call exit(2)
  end if
  print *, 'SURFLAY: ',SURFLAY

! 1) soil moisture prognostics
! ----------------------------
  n =count((old_catch%catdef .gt. old_catch%cdcr1))
  
  write(6,200) n,100*n/ntiles

! Scale rxexc regardless of CDCR1, CDCR2 differences
! --------------------------------------------------
  scale_catch%rzexc  = old_catch%rzexc * ( new_catch%vgwmax / &
                                           old_catch%vgwmax )

! Scale catdef regardless of whether CDCR2 is larger or smaller in the new situation
! ----------------------------------------------------------------------------------
  where (old_catch%catdef .gt. old_catch%cdcr1)
 
      scale_catch%catdef = new_catch%cdcr1 +                   &
                        ( old_catch%catdef-old_catch%cdcr1 ) / &
                        ( old_catch%cdcr2 -old_catch%cdcr1 ) * &
                        ( new_catch%cdcr2 -new_catch%cdcr1 )
  end where

! Scale catdef also for the case where catdef le cdcr1.
! -----------------------------------------------------
  where( (old_catch%catdef .le. old_catch%cdcr1))
      scale_catch%catdef = old_catch%catdef * (new_catch%cdcr1 / old_catch%cdcr1)
  end where

! Sanity Check (catch_calc_soil_moist() forces consistency betw. srfexc, rzexc, catdef)
! ------------
  print *, 'Performing Sanity Check ...'
  allocate (   dzsf(ntiles) )
  allocate (   ar1( ntiles) )
  allocate (   ar2( ntiles) )
  allocate (   ar4( ntiles) )

  dzsf = SURFLAY

  call catch_calc_soil_moist( ntiles, dzsf,                                            &
       scale_catch%vgwmax, scale_catch%cdcr1, scale_catch%cdcr2,                          &
       scale_catch%psis,   scale_catch%bee,   scale_catch%poros, scale_catch%wpwet,        &
       scale_catch%ars1,   scale_catch%ars2,  scale_catch%ars3,                           &
       scale_catch%ara1,   scale_catch%ara2,  scale_catch%ara3,  scale_catch%ara4,         &
       scale_catch%arw1,   scale_catch%arw2,  scale_catch%arw3,  scale_catch%arw4,         &
       scale_catch%bf1,    scale_catch%bf2,                                              &
       scale_catch%srfexc, scale_catch%rzexc, scale_catch%catdef,                         &
       ar1,               ar2,              ar4                                 )

  n = count( scale_catch%catdef .ne. new_catch%catdef )
  write(6,300) n,100*n/ntiles
  n = count( scale_catch%srfexc .ne. new_catch%srfexc )
  write(6,400) n,100*n/ntiles
  n = count( scale_catch%rzexc  .ne. new_catch%rzexc  )
  write(6,400) n,100*n/ntiles

! (2) Ground heat
! ---------------

  allocate (TP_IN  (N_GT, Ntiles))
  allocate (GHT_IN (N_GT, Ntiles))
  allocate (GHT_OUT(N_GT, Ntiles))
  allocate (FICE   (N_GT, NTILES))
  allocate (TP_OUT (N_GT, Ntiles))

  GHT_IN (1,:) = old_catch%ghtcnt1
  GHT_IN (2,:) = old_catch%ghtcnt2
  GHT_IN (3,:) = old_catch%ghtcnt3
  GHT_IN (4,:) = old_catch%ghtcnt4
  GHT_IN (5,:) = old_catch%ghtcnt5
  GHT_IN (6,:) = old_catch%ghtcnt6
  
  call catch_calc_tp ( NTILES, old_catch%poros, GHT_IN, tp_in, FICE)
  GHT_OUT = GHT_IN

!  open (99,file='ght.diff', form = 'formatted')

  do n = 1, ntiles
     do i = 1, N_GT
        call catch_calc_ght(dzgt(i), new_catch%poros(n), tp_in(i,n), fice(i,n),  GHT_IN(i,n))
!        if (i == N_GT) then
!           if (GHT_IN(i,n) /= GHT_OUT(i,n)) write (99,*)n,old_catch%poros(n),new_catch%poros(n),ABS(GHT_IN(i,n)-GHT_OUT(i,n))
!        endif
     end do
  end do

  scale_catch%ghtcnt1 = GHT_IN (1,:)
  scale_catch%ghtcnt2 = GHT_IN (2,:)
  scale_catch%ghtcnt3 = GHT_IN (3,:)
  scale_catch%ghtcnt4 = GHT_IN (4,:)
  scale_catch%ghtcnt5 = GHT_IN (5,:)
  scale_catch%ghtcnt6 = GHT_IN (6,:) 

! Deep soil temp sanity check
! ---------------------------

  call catch_calc_tp ( NTILES, new_catch%poros, GHT_IN, tp_out, FICE)

  print *, 'Percent tiles TP Layer 1 differ : ', 100.* count(ABS(tp_out(1,:) - tp_in(1,:)) > 1.e-5) /float (Ntiles)
  print *, 'Percent tiles TP Layer 2 differ : ', 100.* count(ABS(tp_out(2,:) - tp_in(2,:)) > 1.e-5) /float (Ntiles)
  print *, 'Percent tiles TP Layer 3 differ : ', 100.* count(ABS(tp_out(3,:) - tp_in(3,:)) > 1.e-5) /float (Ntiles)
  print *, 'Percent tiles TP Layer 4 differ : ', 100.* count(ABS(tp_out(4,:) - tp_in(4,:)) > 1.e-5) /float (Ntiles)
  print *, 'Percent tiles TP Layer 5 differ : ', 100.* count(ABS(tp_out(5,:) - tp_in(5,:)) > 1.e-5) /float (Ntiles)
  print *, 'Percent tiles TP Layer 6 differ : ', 100.* count(ABS(tp_out(6,:) - tp_in(6,:)) > 1.e-5) /float (Ntiles)


! SNOW scaling
! ------------

  if(wemin_out /= wemin_in) then

     allocate (swe_in     (Ntiles))
     allocate (depth_in   (Ntiles))
     allocate (depth_out  (Ntiles))
     allocate (areasc_in  (Ntiles))
     allocate (areasc_out (Ntiles))

     swe_in    = new_catch%wesnn1 + new_catch%wesnn2 + new_catch%wesnn3
     depth_in  = new_catch%sndzn1 + new_catch%sndzn2 + new_catch%sndzn3
     areasc_in = min(swe_in/wemin_in, 1.)
     areasc_out= min(swe_in/wemin_out,1.)

     where (swe_in .gt. 0.)
        where (areasc_in .lt. 1. .or. areasc_out .lt. 1.)
           !      density_in= swe_in/(areasc_in *  depth_in + 1.e-20)
           !      depth_out = swe_in/(areasc_out*density_in)
           depth_out = areasc_in *  depth_in/(areasc_out + 1.e-20)
           scale_catch%sndzn1 = depth_out/3.
           scale_catch%sndzn2 = depth_out/3.
           scale_catch%sndzn3 = depth_out/3.
        endwhere
     endwhere

     print *, 'Snow scaling summary'
     print *, '....................'
     print *, 'Percent tiles SNDZ scaled : ', 100.* count (scale_catch%sndzn3 .ne. old_catch%sndzn3) /float (count (scale_catch%sndzn3 > 0.)) 
          
  endif

  ! PEATCLSM - ensure low CATDEF on peat tiles where "old" restart is not also peat
  ! -------------------------------------------------------------------------------

  where ( (old_catch%poros < PEATCLSM_POROS_THRESHOLD) .and. (scale_catch%poros >= PEATCLSM_POROS_THRESHOLD) )
     scale_catch%catdef = 25.
     scale_catch%rzexc  =  0.
     scale_catch%srfexc =  0.
  end where

! Write Scaled Catch
! ------------------
  if (filetype ==0) then
     call formatter%open(new_fname, pFIO_READ, __RC__)
     meta = formatter%read(__RC__)
     call formatter%close()
     call scale_catch%write_nc4(scale_fname, meta, cnclm, __RC__)
  else
     call scale_catch%write_bin(scale_fname, __RC__)
  end if

100 format(1x,'Total  Tiles: ',i10)
200 format(1x,'Scaled Tiles: ',i10,2x,'(',i2.2,'%)')
300 format(1x,'CatDef Tiles: ',i10,2x,'(',i2.2,'%)')
400 format(1x,'SrfExc Tiles: ',i10,2x,'(',i2.2,'%)')
500 format(1x,' Rzexc Tiles: ',i10,2x,'(',i2.2,'%)')

  end program
