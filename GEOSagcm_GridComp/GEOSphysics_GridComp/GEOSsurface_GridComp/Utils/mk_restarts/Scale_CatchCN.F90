#define I_AM_MAIN
#include "MAPL_Generic.h"

program Scale_CatchCN

  use MAPL

  use LSM_ROUTINES,      ONLY:          &
       catch_calc_soil_moist,           &
       catch_calc_tp,                   &
       catch_calc_ght
  
  USE CATCH_CONSTANTS,   ONLY:          &
       N_GT              => CATCH_N_GT, &
       DZGT              => CATCH_DZGT, &
       PEATCLSM_POROS_THRESHOLD
  
  implicit none

  character(256)    :: fname1, fname2, fname3
#ifndef __GFORTRAN__
  integer           :: ftell
  external          :: ftell
#endif
  integer           :: bpos, epos, ntiles, n, nargs
  integer           :: old,  new,  sca
  integer           :: iargc
  real              :: SURFLAY        ! (Ganymed-3 and earlier) SURFLAY=20.0 for Old Soil Params
                                      ! (Ganymed-4 and later  ) SURFLAY=50.0 for New Soil Params
  real              :: WEMIN_IN, WEMIN_OUT
  character*256     :: arg(6)

  integer, parameter :: nveg  = 4
  integer, parameter :: nzone = 3
  integer            :: VAR_COL, VAR_PFT
  integer, parameter :: VAR_COL_CLM40 = 40 ! number of CN column restart variables
  integer, parameter :: VAR_PFT_CLM40 = 74 ! number of CN PFT variables per column
  integer, parameter :: npft    = 19  
  integer, parameter :: VAR_COL_CLM45 = 35 ! number of CN column restart variables
  integer, parameter :: VAR_PFT_CLM45 = 75 ! number of CN PFT variables per column
  
  logical            :: clm45  = .false.
  integer :: un_dim3

  type catch_rst
       real, pointer ::        bf1(:)
       real, pointer ::        bf2(:)
       real, pointer ::        bf3(:)
       real, pointer ::     vgwmax(:)
       real, pointer ::      cdcr1(:)
       real, pointer ::      cdcr2(:)
       real, pointer ::       psis(:)
       real, pointer ::        bee(:)
       real, pointer ::      poros(:)
       real, pointer ::      wpwet(:)
       real, pointer ::       cond(:)
       real, pointer ::        gnu(:)
       real, pointer ::       ars1(:)
       real, pointer ::       ars2(:)
       real, pointer ::       ars3(:)
       real, pointer ::       ara1(:)
       real, pointer ::       ara2(:)
       real, pointer ::       ara3(:)
       real, pointer ::       ara4(:)
       real, pointer ::       arw1(:)
       real, pointer ::       arw2(:)
       real, pointer ::       arw3(:)
       real, pointer ::       arw4(:)
       real, pointer ::       tsa1(:)
       real, pointer ::       tsa2(:)
       real, pointer ::       tsb1(:)
       real, pointer ::       tsb2(:)
       real, pointer ::       atau(:)
       real, pointer ::       btau(:)
       real, pointer ::        ity(:,:)
       real, pointer ::        fvg(:,:)
       real, pointer ::         tc(:,:)
       real, pointer ::         qc(:,:)
       real, pointer ::         tg(:,:)
       real, pointer ::      capac(:)
       real, pointer ::     catdef(:)
       real, pointer ::      rzexc(:)
       real, pointer ::     srfexc(:)
       real, pointer ::    ghtcnt1(:)
       real, pointer ::    ghtcnt2(:)
       real, pointer ::    ghtcnt3(:)
       real, pointer ::    ghtcnt4(:)
       real, pointer ::    ghtcnt5(:)
       real, pointer ::    ghtcnt6(:)
       real, pointer ::      tsurf(:)
       real, pointer ::     wesnn1(:)
       real, pointer ::     wesnn2(:)
       real, pointer ::     wesnn3(:)
       real, pointer ::    htsnnn1(:)
       real, pointer ::    htsnnn2(:)
       real, pointer ::    htsnnn3(:)
       real, pointer ::     sndzn1(:)
       real, pointer ::     sndzn2(:)
       real, pointer ::     sndzn3(:)
       real, pointer ::         ch(:,:)
       real, pointer ::         cm(:,:)
       real, pointer ::         cq(:,:)
       real, pointer ::         fr(:,:)
       real, pointer ::         ww(:,:)
       real, pointer ::     TILE_ID(:)
       real, pointer ::       ndep(:)
       real, pointer ::         t2(:)
       real, pointer ::    BGALBVR(:)
       real, pointer ::    BGALBVF(:)
       real, pointer ::    BGALBNR(:)
       real, pointer ::    BGALBNF(:)
       real, pointer ::    CNCOL(:,:)
       real, pointer ::    CNPFT(:,:)
       real, pointer ::    ABM     (:)
       real, pointer ::    FIELDCAP(:)
       real, pointer ::    HDM     (:)
       real, pointer ::    GDP     (:)
       real, pointer ::    PEATF   (:)       
  endtype catch_rst

  type(catch_rst) catch(3)

  real,    allocatable, dimension(:)   :: dzsf, ar1, ar2, ar4
  real,    allocatable, dimension(:,:) :: TP_IN, GHT_IN, FICE, GHT_OUT, TP_OUT
  real,    allocatable, dimension(:)   :: swe_in, depth_in, areasc_in, areasc_out, depth_out

  type(Netcdf4_fileformatter) :: formatter(3)
  type(Filemetadata) :: cfg(3)
  integer :: i, rc, filetype
  integer :: status  
  character(256) :: Iam = "Scale_CatchCN"
  
! Usage
! -----
  if (iargc() /= 6) then
     write(*,*) "Usage: Scale_CatchCN <Input_Catch> <Regridded_Catch> <Scaled_Catch> <SURFLAY> <WEMIN_IN> <WEMIN_OUT>"
     call exit(2)
  end if

  do n=1,6
  call getarg(n,arg(n))
  enddo

! Open INPUT and Regridded Catch Files
! ------------------------------------
  read(arg(1),'(a)') fname1

  read(arg(2),'(a)') fname2

! Open OUTPUT (Scaled) Catch File
! -------------------------------
  read(arg(3),'(a)') fname3

  call MAPL_NCIOGetFileType(fname1, filetype, __RC__)

  if (filetype == 0) then
     call formatter(1)%open(trim(fname1),pFIO_READ, __RC__)
     call formatter(2)%open(trim(fname2),pFIO_READ, __RC__)
     cfg(1)=formatter(1)%read(__RC__) 
     cfg(2)=formatter(2)%read(__RC__)
 ! else
 !    open(unit=10, file=trim(fname1),  form='unformatted')
 !    open(unit=20, file=trim(fname2),  form='unformatted')
 !    open(unit=30, file=trim(fname3),  form='unformatted')
  end if
  
! Get SURFLAY Value
! -----------------
  read(arg(4),*) SURFLAY
  read(arg(5),*) WEMIN_IN
  read(arg(6),*) WEMIN_OUT

  if (SURFLAY.ne.20 .and. SURFLAY.ne.50) then
     print *, "You must supply a valid SURFLAY value:"
     print *, "(Ganymed-3 and earlier) SURFLAY=20.0 for Old Soil Params"
     print *, "(Ganymed-4 and later  ) SURFLAY=50.0 for New Soil Params"
     call exit(2)
  end if
  print *, 'SURFLAY: ',SURFLAY

  VAR_COL = VAR_COL_CLM40 
  VAR_PFT = VAR_PFT_CLM40

  if (filetype ==0) then

     ntiles = cfg(1)%get_dimension('tile', __RC__)
     un_dim3 = cfg(1)%get_dimension('unknown_dim3', __RC__)
     if(un_dim3 == 105) then
        clm45  = .true.
        VAR_COL = VAR_COL_CLM45 
        VAR_PFT = VAR_PFT_CLM45
        print *, 'Processing CLM45 restarts : ', VAR_COL, VAR_PFT, clm45
     else
        print *, 'Processing CLM40 restarts : ', VAR_COL, VAR_PFT, clm45
     endif
!  else
!
!!    Determine NTILES
!!    ----------------
!     bpos=0
!     read(10)
!     epos = ftell(10)            ! ending position of file pointer
!     ntiles = (epos-bpos)/4-2    ! record size (in 4 byte words; 
!     rewind 10

  end if

  write(6,100) ntiles

! Allocate Catches
! ----------------
  do n=1,3
     call allocatch ( ntiles,catch(n) )
  enddo

! Read INPUT Catches
! ------------------
  old = 1
  new = 2
  
  if (filetype ==0) then
     call readcatchcn_nc4 ( catch(old), formatter(old), cfg(old), __RC__ )
     call readcatchcn_nc4 ( catch(new), formatter(new), cfg(new), __RC__ )
!  else
!     call readcatchcn ( 10,catch(old) )
!     call readcatchcn ( 20,catch(new) )
  end if

! Create Scaled Catch
! -------------------
  sca = 3
  
  catch(sca) = catch(new)

! 1) soil moisture prognostics
! ----------------------------
!  n = count( (catch(old)%catdef .gt. catch(old)%cdcr1) .and. &
!             (catch(new)%cdcr2  .gt. catch(old)%cdcr2) )
!
!  write(6,200) n,100*n/ntiles
!
!  where( (catch(old)%catdef .gt. catch(old)%cdcr1) .and. &
!         (catch(new)%cdcr2  .gt. catch(old)%cdcr2) )
! 
!      catch(sca)%rzexc  = catch(old)%rzexc * ( catch(new)%vgwmax / &
!                                               catch(old)%vgwmax )
!
!      catch(sca)%catdef = catch(new)%cdcr1 +                     &
!                        ( catch(old)%catdef-catch(old)%cdcr1 ) / &
!                        ( catch(old)%cdcr2 -catch(old)%cdcr1 ) * &
!                        ( catch(new)%cdcr2 -catch(new)%cdcr1 )
!  end where

  n =count((catch(old)%catdef .gt. catch(old)%cdcr1))
  
  write(6,200) n,100*n/ntiles

! Scale rxexc regardless of CDCR1, CDCR2 differences
! --------------------------------------------------
  catch(sca)%rzexc  = catch(old)%rzexc * ( catch(new)%vgwmax / &
                                               catch(old)%vgwmax )

! Scale catdef regardless of whether CDCR2 is larger or smaller in the new situation
! ----------------------------------------------------------------------------------
  where (catch(old)%catdef .gt. catch(old)%cdcr1)
 
      catch(sca)%catdef = catch(new)%cdcr1 +                     &
                        ( catch(old)%catdef-catch(old)%cdcr1 ) / &
                        ( catch(old)%cdcr2 -catch(old)%cdcr1 ) * &
                        ( catch(new)%cdcr2 -catch(new)%cdcr1 )
  end where

! Scale catdef also for the case where catdef le cdcr1.
! -----------------------------------------------------
  where( (catch(old)%catdef .le. catch(old)%cdcr1))
      catch(sca)%catdef = catch(old)%catdef * (catch(new)%cdcr1 / catch(old)%cdcr1)
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
       catch(sca)%vgwmax, catch(sca)%cdcr1, catch(sca)%cdcr2,                          &
       catch(sca)%psis,   catch(sca)%bee,   catch(sca)%poros, catch(sca)%wpwet,        &
       catch(sca)%ars1,   catch(sca)%ars2,  catch(sca)%ars3,                           &
       catch(sca)%ara1,   catch(sca)%ara2,  catch(sca)%ara3,  catch(sca)%ara4,         &
       catch(sca)%arw1,   catch(sca)%arw2,  catch(sca)%arw3,  catch(sca)%arw4,         &
       catch(sca)%bf1,    catch(sca)%bf2,                                              &
       catch(sca)%srfexc, catch(sca)%rzexc, catch(sca)%catdef,                         &
       ar1,               ar2,              ar4                                 )
  
  n = count( catch(sca)%catdef .ne. catch(new)%catdef )
  write(6,300) n,100*n/ntiles
  n = count( catch(sca)%srfexc .ne. catch(new)%srfexc )
  write(6,400) n,100*n/ntiles
  n = count( catch(sca)%rzexc  .ne. catch(new)%rzexc  )
  write(6,400) n,100*n/ntiles

! (2) Ground heat
! ---------------

  allocate (TP_IN  (N_GT, Ntiles))
  allocate (GHT_IN (N_GT, Ntiles))
  allocate (GHT_OUT(N_GT, Ntiles))
  allocate (FICE   (N_GT, NTILES))
  allocate (TP_OUT (N_GT, Ntiles))

  GHT_IN (1,:) = catch(old)%ghtcnt1
  GHT_IN (2,:) = catch(old)%ghtcnt2
  GHT_IN (3,:) = catch(old)%ghtcnt3
  GHT_IN (4,:) = catch(old)%ghtcnt4
  GHT_IN (5,:) = catch(old)%ghtcnt5
  GHT_IN (6,:) = catch(old)%ghtcnt6
  
  call catch_calc_tp ( NTILES, catch(old)%poros, GHT_IN, tp_in, FICE)
  GHT_OUT = GHT_IN

!  open (99,file='ght.diff', form = 'formatted')

  do n = 1, ntiles
     do i = 1, N_GT
        call catch_calc_ght(dzgt(i), catch(new)%poros(n), tp_in(i,n), fice(i,n),  GHT_IN(i,n))
!        if (i == N_GT) then
!           if (GHT_IN(i,n) /= GHT_OUT(i,n)) write (99,*)n,catch(old)%poros(n),catch(new)%poros(n),ABS(GHT_IN(i,n)-GHT_OUT(i,n))
!        endif
     end do
  end do

  catch(sca)%ghtcnt1 = GHT_IN (1,:)
  catch(sca)%ghtcnt2 = GHT_IN (2,:)
  catch(sca)%ghtcnt3 = GHT_IN (3,:)
  catch(sca)%ghtcnt4 = GHT_IN (4,:)
  catch(sca)%ghtcnt5 = GHT_IN (5,:)
  catch(sca)%ghtcnt6 = GHT_IN (6,:) 

! Deep soil temp sanity check
! ---------------------------

  call catch_calc_tp ( NTILES, catch(new)%poros, GHT_IN, tp_out, FICE)

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

     swe_in    = catch(new)%wesnn1 + catch(new)%wesnn2 + catch(new)%wesnn3
     depth_in  = catch(new)%sndzn1 + catch(new)%sndzn2 + catch(new)%sndzn3
     areasc_in = min(swe_in/wemin_in, 1.)
     areasc_out= min(swe_in/wemin_out,1.)

     !  catch(sca)%sndzn1=catch(old)%sndzn1
     !  catch(sca)%sndzn2=catch(old)%sndzn2
     !  catch(sca)%sndzn3=catch(old)%sndzn3     
     !     do i = 1, ntiles
     !         if((swe_in(i) > 0.).and. ((areasc_in(i) < 1.).OR.(areasc_out(i) < 1.))) then
     !        print *, i, areasc_in(i), depth_in(i)
     !        density_in(i)= swe_in(i)/(areasc_in(i) *  depth_in(i))
     !        depth_out(i) = swe_in(i)/(areasc_out(i)*density_in(i))
     !        depth_out(i) = areasc_in(i) *  depth_in(i)/(areasc_out(i) + 1.e-20)
     !            print *, catch(sca)%sndzn1(i), catch(old)%sndzn1(i),wemin_out/wemin_in
     !            catch(sca)%sndzn1(i) = catch(new)%sndzn1(i)*wemin_out/wemin_in   ! depth_out(i)/3.
     !            catch(sca)%sndzn2(i) = catch(new)%sndzn2(i)*wemin_out/wemin_in   ! depth_out(i)/3.
     !            catch(sca)%sndzn3(i) = catch(new)%sndzn3(i)*wemin_out/wemin_in   ! depth_out(i)/3.
     !         endif 
     !      end do
     
     where (swe_in .gt. 0.)
        where (areasc_in .lt. 1. .or. areasc_out .lt. 1.)
           !      density_in= swe_in/(areasc_in *  depth_in + 1.e-20)
           !      depth_out = swe_in/(areasc_out*density_in)
           depth_out = areasc_in *  depth_in/(areasc_out + 1.e-20)
           catch(sca)%sndzn1 = depth_out/3.
           catch(sca)%sndzn2 = depth_out/3.
           catch(sca)%sndzn3 = depth_out/3.
        endwhere
     endwhere

     print *, 'Snow scaling summary'
     print *, '....................'
     print *, 'Percent tiles SNDZ scaled : ', 100.* count (catch(sca)%sndzn3 .ne. catch(old)%sndzn3) /float (count (catch(sca)%sndzn3 > 0.)) 
          
  endif

  ! PEATCLSM - ensure low CATDEF on peat tiles where "old" restart is not also peat
  ! -------------------------------------------------------------------------------

  where ( (catch(old)%poros < PEATCLSM_POROS_THRESHOLD) .and. (catch(sca)%poros >= PEATCLSM_POROS_THRESHOLD) )
     catch(sca)%catdef = 25.
     catch(sca)%rzexc  =  0.
     catch(sca)%srfexc =  0.
  end where

! Write Scaled Catch
! ------------------
  if (filetype ==0) then
     cfg(3)=cfg(2)
     call formatter(3)%create(fname3, __RC__)
     call formatter(3)%write(cfg(3), __RC__)
     call writecatchcn_nc4 ( catch(sca), formatter(3) ,cfg(3) )
!  else
!     call writecatchcn ( 30,catch(sca) )
  end if

100 format(1x,'Total  Tiles: ',i10)
200 format(1x,'Scaled Tiles: ',i10,2x,'(',i2.2,'%)')
300 format(1x,'CatDef Tiles: ',i10,2x,'(',i2.2,'%)')
400 format(1x,'SrfExc Tiles: ',i10,2x,'(',i2.2,'%)')
500 format(1x,' Rzexc Tiles: ',i10,2x,'(',i2.2,'%)')

  stop

  contains

    subroutine allocatch (ntiles,catch)
      
      integer ntiles
      
      type(catch_rst) catch

       allocate( catch%        bf1(ntiles) )
       allocate( catch%        bf2(ntiles) )
       allocate( catch%        bf3(ntiles) )
       allocate( catch%     vgwmax(ntiles) )
       allocate( catch%      cdcr1(ntiles) )
       allocate( catch%      cdcr2(ntiles) )
       allocate( catch%       psis(ntiles) )
       allocate( catch%        bee(ntiles) )
       allocate( catch%      poros(ntiles) )
       allocate( catch%      wpwet(ntiles) )
       allocate( catch%       cond(ntiles) )
       allocate( catch%        gnu(ntiles) )
       allocate( catch%       ars1(ntiles) )
       allocate( catch%       ars2(ntiles) )
       allocate( catch%       ars3(ntiles) )
       allocate( catch%       ara1(ntiles) )
       allocate( catch%       ara2(ntiles) )
       allocate( catch%       ara3(ntiles) )
       allocate( catch%       ara4(ntiles) )
       allocate( catch%       arw1(ntiles) )
       allocate( catch%       arw2(ntiles) )
       allocate( catch%       arw3(ntiles) )
       allocate( catch%       arw4(ntiles) )
       allocate( catch%       tsa1(ntiles) )
       allocate( catch%       tsa2(ntiles) )
       allocate( catch%       tsb1(ntiles) )
       allocate( catch%       tsb2(ntiles) )
       allocate( catch%       atau(ntiles) )
       allocate( catch%       btau(ntiles) )
       allocate( catch%        ity(ntiles,4) )
       allocate( catch%        fvg(ntiles,4) )
       allocate( catch%         tc(ntiles,4) )
       allocate( catch%         qc(ntiles,4) )
       allocate( catch%         tg(ntiles,4) )
       allocate( catch%      capac(ntiles) )
       allocate( catch%     catdef(ntiles) )
       allocate( catch%      rzexc(ntiles) )
       allocate( catch%     srfexc(ntiles) )
       allocate( catch%    ghtcnt1(ntiles) )
       allocate( catch%    ghtcnt2(ntiles) )
       allocate( catch%    ghtcnt3(ntiles) )
       allocate( catch%    ghtcnt4(ntiles) )
       allocate( catch%    ghtcnt5(ntiles) )
       allocate( catch%    ghtcnt6(ntiles) )
       allocate( catch%      tsurf(ntiles) )
       allocate( catch%     wesnn1(ntiles) )
       allocate( catch%     wesnn2(ntiles) )
       allocate( catch%     wesnn3(ntiles) )
       allocate( catch%    htsnnn1(ntiles) )
       allocate( catch%    htsnnn2(ntiles) )
       allocate( catch%    htsnnn3(ntiles) )
       allocate( catch%     sndzn1(ntiles) )
       allocate( catch%     sndzn2(ntiles) )
       allocate( catch%     sndzn3(ntiles) )
       allocate( catch%         ch(ntiles,4) )
       allocate( catch%         cm(ntiles,4) )
       allocate( catch%         cq(ntiles,4) )
       allocate( catch%         fr(ntiles,4) )
       allocate( catch%         ww(ntiles,4) )
       allocate( catch%     TILE_ID(ntiles) )
       allocate( catch%       ndep(ntiles) )
       allocate( catch%         t2(ntiles) )
       allocate( catch%    BGALBVR(ntiles) )
       allocate( catch%    BGALBVF(ntiles) )
       allocate( catch%    BGALBNR(ntiles) )
       allocate( catch%    BGALBNF(ntiles) )
       allocate( catch%      CNCOL(ntiles,nzone*VAR_COL))
       allocate( catch%      CNPFT(ntiles,nzone*nveg*VAR_PFT))
       allocate( catch%        ABM(ntiles) )
       allocate( catch%   FIELDCAP(ntiles) )
       allocate( catch%        HDM(ntiles) )
       allocate( catch%        GDP(ntiles) )
       allocate( catch%      PEATF(ntiles) )
       
   return
   end subroutine allocatch

   subroutine readcatchcn_nc4 (catch,formatter,cfg, rc)
      type(catch_rst) catch
      type(Filemetadata) :: cfg
      type(Netcdf4_fileformatter) :: formatter
      integer, optional, intent(out) :: rc
      integer :: j, dim1,dim2
      type(Variable), pointer :: myVariable
      character(len=:), pointer :: dname
      integer :: status
      character(256) :: Iam = "readcatchcn_nc4"

      call MAPL_VarRead(formatter,"BF1",catch%bf1, __RC__)
      call MAPL_VarRead(formatter,"BF2",catch%bf2, __RC__)
      call MAPL_VarRead(formatter,"BF3",catch%bf3, __RC__)
      call MAPL_VarRead(formatter,"VGWMAX",catch%vgwmax, __RC__)
      call MAPL_VarRead(formatter,"CDCR1",catch%cdcr1, __RC__)
      call MAPL_VarRead(formatter,"CDCR2",catch%cdcr2, __RC__)
      call MAPL_VarRead(formatter,"PSIS",catch%psis, __RC__)
      call MAPL_VarRead(formatter,"BEE",catch%bee, __RC__)
      call MAPL_VarRead(formatter,"POROS",catch%poros, __RC__)
      call MAPL_VarRead(formatter,"WPWET",catch%wpwet, __RC__)
      call MAPL_VarRead(formatter,"COND",catch%cond, __RC__)
      call MAPL_VarRead(formatter,"GNU",catch%gnu, __RC__)
      call MAPL_VarRead(formatter,"ARS1",catch%ars1, __RC__)
      call MAPL_VarRead(formatter,"ARS2",catch%ars2, __RC__)
      call MAPL_VarRead(formatter,"ARS3",catch%ars3, __RC__)
      call MAPL_VarRead(formatter,"ARA1",catch%ara1, __RC__)
      call MAPL_VarRead(formatter,"ARA2",catch%ara2, __RC__)
      call MAPL_VarRead(formatter,"ARA3",catch%ara3, __RC__)
      call MAPL_VarRead(formatter,"ARA4",catch%ara4, __RC__)
      call MAPL_VarRead(formatter,"ARW1",catch%arw1, __RC__)
      call MAPL_VarRead(formatter,"ARW2",catch%arw2, __RC__)
      call MAPL_VarRead(formatter,"ARW3",catch%arw3, __RC__)
      call MAPL_VarRead(formatter,"ARW4",catch%arw4, __RC__)
      call MAPL_VarRead(formatter,"TSA1",catch%tsa1, __RC__)
      call MAPL_VarRead(formatter,"TSA2",catch%tsa2, __RC__)
      call MAPL_VarRead(formatter,"TSB1",catch%tsb1, __RC__)
      call MAPL_VarRead(formatter,"TSB2",catch%tsb2, __RC__)
      call MAPL_VarRead(formatter,"ATAU",catch%atau, __RC__)
      call MAPL_VarRead(formatter,"BTAU",catch%btau, __RC__)

      myVariable => cfg%get_variable("ITY")
      dname => myVariable%get_ith_dimension(2)
      dim1 = cfg%get_dimension(dname)
      do j=1,dim1
         call MAPL_VarRead(formatter,"ITY",catch%ity(:,j),offset1=j, __RC__)
         call MAPL_VarRead(formatter,"FVG",catch%fvg(:,j),offset1=j, __RC__)
      enddo

      call MAPL_VarRead(formatter,"TC",catch%tc, __RC__)
      call MAPL_VarRead(formatter,"QC",catch%qc, __RC__)
      call MAPL_VarRead(formatter,"TG",catch%tg, __RC__)
      call MAPL_VarRead(formatter,"CAPAC",catch%capac, __RC__)
      call MAPL_VarRead(formatter,"CATDEF",catch%catdef, __RC__)
      call MAPL_VarRead(formatter,"RZEXC",catch%rzexc, __RC__)
      call MAPL_VarRead(formatter,"SRFEXC",catch%srfexc, __RC__)
      call MAPL_VarRead(formatter,"GHTCNT1",catch%ghtcnt1, __RC__)
      call MAPL_VarRead(formatter,"GHTCNT2",catch%ghtcnt2, __RC__)
      call MAPL_VarRead(formatter,"GHTCNT3",catch%ghtcnt3, __RC__)
      call MAPL_VarRead(formatter,"GHTCNT4",catch%ghtcnt4, __RC__)
      call MAPL_VarRead(formatter,"GHTCNT5",catch%ghtcnt5, __RC__)
      call MAPL_VarRead(formatter,"GHTCNT6",catch%ghtcnt6, __RC__)
      call MAPL_VarRead(formatter,"TSURF",catch%tsurf, __RC__)
      call MAPL_VarRead(formatter,"WESNN1",catch%wesnn1, __RC__)
      call MAPL_VarRead(formatter,"WESNN2",catch%wesnn2, __RC__)
      call MAPL_VarRead(formatter,"WESNN3",catch%wesnn3, __RC__)
      call MAPL_VarRead(formatter,"HTSNNN1",catch%htsnnn1, __RC__)
      call MAPL_VarRead(formatter,"HTSNNN2",catch%htsnnn2, __RC__)
      call MAPL_VarRead(formatter,"HTSNNN3",catch%htsnnn3, __RC__)
      call MAPL_VarRead(formatter,"SNDZN1",catch%sndzn1, __RC__)
      call MAPL_VarRead(formatter,"SNDZN2",catch%sndzn2, __RC__)
      call MAPL_VarRead(formatter,"SNDZN3",catch%sndzn3, __RC__)
      call MAPL_VarRead(formatter,"CH",catch%ch, __RC__)
      call MAPL_VarRead(formatter,"CM",catch%cm, __RC__)
      call MAPL_VarRead(formatter,"CQ",catch%cq, __RC__)
      call MAPL_VarRead(formatter,"FR",catch%fr, __RC__)
      call MAPL_VarRead(formatter,"WW",catch%ww, __RC__)
      call MAPL_VarRead(formatter,"TILE_ID",catch%TILE_ID, __RC__)
      call MAPL_VarRead(formatter,"NDEP",catch%ndep, __RC__)
      call MAPL_VarRead(formatter,"CLI_T2M",catch%t2, __RC__)
      call MAPL_VarRead(formatter,"BGALBVR",catch%BGALBVR, __RC__)
      call MAPL_VarRead(formatter,"BGALBVF",catch%BGALBVF, __RC__)
      call MAPL_VarRead(formatter,"BGALBNR",catch%BGALBNR, __RC__)
      call MAPL_VarRead(formatter,"BGALBNF",catch%BGALBNF, __RC__)
      myVariable => cfg%get_variable("CNCOL")
      dname => myVariable%get_ith_dimension(2)
      dim1 = cfg%get_dimension(dname)
      if(clm45) then          
         call MAPL_VarRead(formatter,"ABM",     catch%ABM, __RC__)
         call MAPL_VarRead(formatter,"FIELDCAP",catch%FIELDCAP, __RC__)
         call MAPL_VarRead(formatter,"HDM",     catch%HDM     , __RC__)
         call MAPL_VarRead(formatter,"GDP",     catch%GDP     , __RC__)
         call MAPL_VarRead(formatter,"PEATF",   catch%PEATF   , __RC__)
      endif
      do j=1,dim1
         call MAPL_VarRead(formatter,"CNCOL",catch%CNCOL(:,j),offset1=j, __RC__)
      enddo
      ! The following three lines were added as a bug fix by smahanam on 5 Oct 2020
      ! (to be merged into the "develop" branch in late 2020):
      ! The length of the 2nd dim of CNPFT differs from that of CNCOL.  Prior to this fix,
      ! CNPFT was not read in its entirety and some elements remained uninitialized (or zero),
      ! resulting in bad values in the "regridded" (re-tiled) restart file. 
      ! This impacted re-tiled restarts for both CNCLM40 and CLCLM45.
      ! - reichle, 23 Nov 2020
      myVariable => cfg%get_variable("CNPFT")
      dname => myVariable%get_ith_dimension(2)
      dim1 = cfg%get_dimension(dname)       
      do j=1,dim1
         call MAPL_VarRead(formatter,"CNPFT",catch%CNPFT(:,j),offset1=j, __RC__)
      enddo
      if (present(rc)) rc =0
      !_RETURN(_SUCCESS)
   end subroutine readcatchcn_nc4

   subroutine readcatchcn (unit,catch)
   integer unit, i,j,n
   type(catch_rst) catch

       read(unit) catch%      bf1
       read(unit) catch%      bf2
       read(unit) catch%      bf3
       read(unit) catch%   vgwmax
       read(unit) catch%    cdcr1
       read(unit) catch%    cdcr2
       read(unit) catch%     psis
       read(unit) catch%      bee
       read(unit) catch%    poros
       read(unit) catch%    wpwet
       read(unit) catch%     cond
       read(unit) catch%      gnu
       read(unit) catch%     ars1
       read(unit) catch%     ars2
       read(unit) catch%     ars3
       read(unit) catch%     ara1
       read(unit) catch%     ara2
       read(unit) catch%     ara3
       read(unit) catch%     ara4
       read(unit) catch%     arw1
       read(unit) catch%     arw2
       read(unit) catch%     arw3
       read(unit) catch%     arw4
       read(unit) catch%     tsa1
       read(unit) catch%     tsa2
       read(unit) catch%     tsb1
       read(unit) catch%     tsb2
       read(unit) catch%     atau
       read(unit) catch%     btau
       read(unit) catch%      ity(:,1)
       read(unit) catch%      ity(:,2)
       read(unit) catch%      ity(:,3)
       read(unit) catch%      ity(:,4)
       read(unit) catch%      fvg(:,1)
       read(unit) catch%      fvg(:,2)
       read(unit) catch%      fvg(:,3)
       read(unit) catch%      fvg(:,4)
       read(unit) catch%       tc
       read(unit) catch%       qc
       read(unit) catch%       tg
       read(unit) catch%    capac
       read(unit) catch%   catdef
       read(unit) catch%    rzexc
       read(unit) catch%   srfexc
       read(unit) catch%  ghtcnt1
       read(unit) catch%  ghtcnt2
       read(unit) catch%  ghtcnt3
       read(unit) catch%  ghtcnt4
       read(unit) catch%  ghtcnt5
       read(unit) catch%  ghtcnt6
       read(unit) catch%    tsurf
       read(unit) catch%   wesnn1
       read(unit) catch%   wesnn2
       read(unit) catch%   wesnn3
       read(unit) catch%  htsnnn1
       read(unit) catch%  htsnnn2
       read(unit) catch%  htsnnn3
       read(unit) catch%   sndzn1
       read(unit) catch%   sndzn2
       read(unit) catch%   sndzn3
       read(unit) catch%       ch
       read(unit) catch%       cm
       read(unit) catch%       cq
       read(unit) catch%       fr
       read(unit) catch%       ww
       read(unit) catch%  TILE_ID
       read(unit) catch%     ndep
       read(unit) catch%       t2
       read(unit) catch%  BGALBVR
       read(unit) catch%  BGALBVF
       read(unit) catch%  BGALBNR
       read(unit) catch%  BGALBNF

       do j = 1,nzone * VAR_COL
          read(unit) catch%    CNCOL (:,j)
       end do

       do i = 1,nzone * nveg * VAR_PFT
          read(unit) catch%    CNPFT (:,i)
       end do
   return
   end subroutine readcatchcn

   subroutine writecatchcn_nc4 (catch,formatter,cfg)
   type(catch_rst) catch
   type(Netcdf4_fileformatter) :: formatter
   type(filemetadata) :: cfg
   integer :: i,j, dim1,dim2
   real, dimension (:), allocatable :: var
   type(Variable), pointer :: myVariable
   character(len=:), pointer :: dname

       call MAPL_VarWrite(formatter,"BF1",catch%bf1)
       call MAPL_VarWrite(formatter,"BF2",catch%bf2)
       call MAPL_VarWrite(formatter,"BF3",catch%bf3)
       call MAPL_VarWrite(formatter,"VGWMAX",catch%vgwmax)
       call MAPL_VarWrite(formatter,"CDCR1",catch%cdcr1)
       call MAPL_VarWrite(formatter,"CDCR2",catch%cdcr2)
       call MAPL_VarWrite(formatter,"PSIS",catch%psis)
       call MAPL_VarWrite(formatter,"BEE",catch%bee)
       call MAPL_VarWrite(formatter,"POROS",catch%poros)
       call MAPL_VarWrite(formatter,"WPWET",catch%wpwet)
       call MAPL_VarWrite(formatter,"COND",catch%cond)
       call MAPL_VarWrite(formatter,"GNU",catch%gnu)
       call MAPL_VarWrite(formatter,"ARS1",catch%ars1)
       call MAPL_VarWrite(formatter,"ARS2",catch%ars2)
       call MAPL_VarWrite(formatter,"ARS3",catch%ars3)
       call MAPL_VarWrite(formatter,"ARA1",catch%ara1)
       call MAPL_VarWrite(formatter,"ARA2",catch%ara2)
       call MAPL_VarWrite(formatter,"ARA3",catch%ara3)
       call MAPL_VarWrite(formatter,"ARA4",catch%ara4)
       call MAPL_VarWrite(formatter,"ARW1",catch%arw1)
       call MAPL_VarWrite(formatter,"ARW2",catch%arw2)
       call MAPL_VarWrite(formatter,"ARW3",catch%arw3)
       call MAPL_VarWrite(formatter,"ARW4",catch%arw4)
       call MAPL_VarWrite(formatter,"TSA1",catch%tsa1)
       call MAPL_VarWrite(formatter,"TSA2",catch%tsa2)
       call MAPL_VarWrite(formatter,"TSB1",catch%tsb1)
       call MAPL_VarWrite(formatter,"TSB2",catch%tsb2)
       call MAPL_VarWrite(formatter,"ATAU",catch%atau)
       call MAPL_VarWrite(formatter,"BTAU",catch%btau)

       myVariable => cfg%get_variable("ITY")
       dname => myVariable%get_ith_dimension(2)
       dim1 = cfg%get_dimension(dname)
       do j=1,dim1
          call MAPL_VarWrite(formatter,"ITY",catch%ity(:,j),offset1=j)
          call MAPL_VarWrite(formatter,"FVG",catch%fvg(:,j),offset1=j)
       enddo

       call MAPL_VarWrite(formatter,"TC",catch%tc)
       call MAPL_VarWrite(formatter,"QC",catch%qc)
       call MAPL_VarWrite(formatter,"TG",catch%TG)
       call MAPL_VarWrite(formatter,"CAPAC",catch%capac)
       call MAPL_VarWrite(formatter,"CATDEF",catch%catdef)
       call MAPL_VarWrite(formatter,"RZEXC",catch%rzexc)
       call MAPL_VarWrite(formatter,"SRFEXC",catch%srfexc)
       call MAPL_VarWrite(formatter,"GHTCNT1",catch%ghtcnt1)
       call MAPL_VarWrite(formatter,"GHTCNT2",catch%ghtcnt2)
       call MAPL_VarWrite(formatter,"GHTCNT3",catch%ghtcnt3)
       call MAPL_VarWrite(formatter,"GHTCNT4",catch%ghtcnt4)
       call MAPL_VarWrite(formatter,"GHTCNT5",catch%ghtcnt5)
       call MAPL_VarWrite(formatter,"GHTCNT6",catch%ghtcnt6)
       call MAPL_VarWrite(formatter,"TSURF",catch%tsurf)
       call MAPL_VarWrite(formatter,"WESNN1",catch%wesnn1)
       call MAPL_VarWrite(formatter,"WESNN2",catch%wesnn2)
       call MAPL_VarWrite(formatter,"WESNN3",catch%wesnn3)
       call MAPL_VarWrite(formatter,"HTSNNN1",catch%htsnnn1)
       call MAPL_VarWrite(formatter,"HTSNNN2",catch%htsnnn2)
       call MAPL_VarWrite(formatter,"HTSNNN3",catch%htsnnn3)
       call MAPL_VarWrite(formatter,"SNDZN1",catch%sndzn1)
       call MAPL_VarWrite(formatter,"SNDZN2",catch%sndzn2)
       call MAPL_VarWrite(formatter,"SNDZN3",catch%sndzn3)
       call MAPL_VarWrite(formatter,"CH",catch%ch)
       call MAPL_VarWrite(formatter,"CM",catch%cm)
       call MAPL_VarWrite(formatter,"CQ",catch%cq)
       call MAPL_VarWrite(formatter,"FR",catch%fr)
       call MAPL_VarWrite(formatter,"WW",catch%ww)
       call MAPL_VarWrite(formatter,"TILE_ID",catch%TILE_ID)
       call MAPL_VarWrite(formatter,"NDEP",catch%NDEP)
       call MAPL_VarWrite(formatter,"CLI_T2M",catch%t2)
       call MAPL_VarWrite(formatter,"BGALBVR",catch%BGALBVR)
       call MAPL_VarWrite(formatter,"BGALBVF",catch%BGALBVF)
       call MAPL_VarWrite(formatter,"BGALBNR",catch%BGALBNR)
       call MAPL_VarWrite(formatter,"BGALBNF",catch%BGALBNF)
       myVariable => cfg%get_variable("CNCOL")
       dname => myVariable%get_ith_dimension(2)
       dim1 = cfg%get_dimension(dname)

       do j=1,dim1
          call MAPL_VarWrite(formatter,"CNCOL",catch%CNCOL(:,j),offset1=j)
       enddo
       myVariable => cfg%get_variable("CNPFT")
       dname => myVariable%get_ith_dimension(2)
       dim1 = cfg%get_dimension(dname)
       do j=1,dim1
          call MAPL_VarWrite(formatter,"CNPFT",catch%CNPFT(:,j),offset1=j)
       enddo

       dim1 = cfg%get_dimension('tile')
       allocate (var (dim1))
       var = 0.

       call MAPL_VarWrite(formatter,"BFLOWM", var)
       call MAPL_VarWrite(formatter,"TOTWATM",var)
       call MAPL_VarWrite(formatter,"TAIRM",  var)
       call MAPL_VarWrite(formatter,"TPM",    var)
       call MAPL_VarWrite(formatter,"CNSUM",  var)
       call MAPL_VarWrite(formatter,"SNDZM",  var)
       call MAPL_VarWrite(formatter,"ASNOWM", var)

       myVariable => cfg%get_variable("TGWM")
       dname => myVariable%get_ith_dimension(2)
       dim1 = cfg%get_dimension(dname)
       do j=1,dim1
          call MAPL_VarWrite(formatter,"TGWM",var,offset1=j)
          call MAPL_VarWrite(formatter,"RZMM",var,offset1=j)
       end do

       if (clm45) then
          do j=1,dim1
             call MAPL_VarWrite(formatter,"SFMM",  var,offset1=j)
          enddo

          call MAPL_VarWrite(formatter,"ABM",     catch%ABM, rc =rc     )
          call MAPL_VarWrite(formatter,"FIELDCAP",catch%FIELDCAP)
          call MAPL_VarWrite(formatter,"HDM",     catch%HDM     )
          call MAPL_VarWrite(formatter,"GDP",     catch%GDP     )
          call MAPL_VarWrite(formatter,"PEATF",   catch%PEATF   )
          call MAPL_VarWrite(formatter,"RHM",     var)
          call MAPL_VarWrite(formatter,"WINDM",   var)
          call MAPL_VarWrite(formatter,"RAINFM",  var)
          call MAPL_VarWrite(formatter,"SNOWFM",  var)
          call MAPL_VarWrite(formatter,"RUNSRFM", var)
          call MAPL_VarWrite(formatter,"AR1M",    var)
          call MAPL_VarWrite(formatter,"T2M10D",  var)
          call MAPL_VarWrite(formatter,"TPREC10D",var)
          call MAPL_VarWrite(formatter,"TPREC60D",var)
       else
          call MAPL_VarWrite(formatter,"SFMCM",  var)          
       endif
       
       myVariable => cfg%get_variable("PSNSUNM")
       dname => myVariable%get_ith_dimension(2)
       dim1 = cfg%get_dimension(dname)
       dname => myVariable%get_ith_dimension(3)
       dim2 = cfg%get_dimension(dname)
       do i=1,dim2 
          do j=1,dim1
             call MAPL_VarWrite(formatter,"PSNSUNM",var,offset1=j,offset2=i)
             call MAPL_VarWrite(formatter,"PSNSHAM",var,offset1=j,offset2=i)
          end do
       end do
       call formatter%close()
   return
   end subroutine writecatchcn_nc4

   subroutine writecatchcn (unit,catch)
   integer unit, i,j,n
   type(catch_rst) catch

       write(unit) catch%      bf1
       write(unit) catch%      bf2
       write(unit) catch%      bf3
       write(unit) catch%   vgwmax
       write(unit) catch%    cdcr1
       write(unit) catch%    cdcr2
       write(unit) catch%     psis
       write(unit) catch%      bee
       write(unit) catch%    poros
       write(unit) catch%    wpwet
       write(unit) catch%     cond
       write(unit) catch%      gnu
       write(unit) catch%     ars1
       write(unit) catch%     ars2
       write(unit) catch%     ars3
       write(unit) catch%     ara1
       write(unit) catch%     ara2
       write(unit) catch%     ara3
       write(unit) catch%     ara4
       write(unit) catch%     arw1
       write(unit) catch%     arw2
       write(unit) catch%     arw3
       write(unit) catch%     arw4
       write(unit) catch%     tsa1
       write(unit) catch%     tsa2
       write(unit) catch%     tsb1
       write(unit) catch%     tsb2
       write(unit) catch%     atau
       write(unit) catch%     btau
       write(unit) catch%      ity(:,1)
       write(unit) catch%      ity(:,2)
       write(unit) catch%      ity(:,3)
       write(unit) catch%      ity(:,4)
       write(unit) catch%      fvg(:,1)
       write(unit) catch%      fvg(:,2)
       write(unit) catch%      fvg(:,3)
       write(unit) catch%      fvg(:,4)
       write(unit) catch%       tc
       write(unit) catch%       qc
       write(unit) catch%       tg
       write(unit) catch%    capac
       write(unit) catch%   catdef
       write(unit) catch%    rzexc
       write(unit) catch%   srfexc
       write(unit) catch%  ghtcnt1
       write(unit) catch%  ghtcnt2
       write(unit) catch%  ghtcnt3
       write(unit) catch%  ghtcnt4
       write(unit) catch%  ghtcnt5
       write(unit) catch%  ghtcnt6
       write(unit) catch%    tsurf
       write(unit) catch%   wesnn1
       write(unit) catch%   wesnn2
       write(unit) catch%   wesnn3
       write(unit) catch%  htsnnn1
       write(unit) catch%  htsnnn2
       write(unit) catch%  htsnnn3
       write(unit) catch%   sndzn1
       write(unit) catch%   sndzn2
       write(unit) catch%   sndzn3
       write(unit) catch%       ch
       write(unit) catch%       cm
       write(unit) catch%       cq
       write(unit) catch%       fr
       write(unit) catch%       ww
       write(unit) catch%  TILE_ID
       write(unit) catch%     ndep
       write(unit) catch%       t2
       write(unit) catch%  BGALBVR
       write(unit) catch%  BGALBVF
       write(unit) catch%  BGALBNR
       write(unit) catch%  BGALBNF

       do j = 1,nzone * VAR_COL
          write(unit) catch%    CNCOL (:,j)
       end do

       do i = 1,nzone * nveg * VAR_PFT
          write(unit) catch%    CNPFT (:,i)
       end do

   return
   end subroutine writecatchcn

  end program

