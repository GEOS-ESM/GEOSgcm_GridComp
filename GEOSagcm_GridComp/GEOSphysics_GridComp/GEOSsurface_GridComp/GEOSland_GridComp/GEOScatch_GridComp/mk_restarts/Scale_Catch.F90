  use MAPL_IOMod
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
  character*256     :: arg(4)

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
       real, pointer ::        ity(:)
       real, pointer ::         tc(:,:)
       real, pointer ::         qc(:,:)
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
  endtype catch_rst

  type(catch_rst) catch(3)

  interface
  subroutine calc_soil_moist( &
       ncat,vegcls,dzsf,vgwmax,cdcr1,cdcr2,wpwet,poros, &
       psis,bee,ars1,ars2,ars3,ara1,ara2, &
       ara3,ara4,arw1,arw2,arw3,arw4, &
       srfexc,rzexc,catdef, &
       sfmc, rzmc, prmc,  &
       werror, sfmcun, rzmcun, prmcun )
    
    implicit none
    
    integer, parameter                   :: KSNGL=4
    integer,                  intent(in) :: ncat
    integer, dimension(ncat), intent(in) :: vegcls
    
    real(KIND=KSNGL), dimension(ncat), intent(in) :: dzsf,vgwmax,cdcr1,cdcr2
    real(KIND=KSNGL), dimension(ncat), intent(in) :: wpwet,poros,psis
    real(KIND=KSNGL), dimension(ncat), intent(in) :: bee,ars1
    real(KIND=KSNGL), dimension(ncat), intent(in) :: ars2,ars3,ara1,ara2,ara3
    real(KIND=KSNGL), dimension(ncat), intent(in) :: ara4,arw1,arw2,arw3,arw4
    
    real(KIND=KSNGL), dimension(ncat), intent(inout) :: srfexc, rzexc, catdef
    
    real(KIND=KSNGL), dimension(ncat), intent(out) :: sfmc, rzmc, prmc
    
    real(KIND=KSNGL), dimension(ncat), intent(out), optional :: werror
    
    real(KIND=KSNGL), dimension(ncat), intent(out), optional :: sfmcun
    real(KIND=KSNGL), dimension(ncat), intent(out), optional :: rzmcun
    real(KIND=KSNGL), dimension(ncat), intent(out), optional :: prmcun
  end subroutine calc_soil_moist
  end interface

  real,    allocatable, dimension(:) :: sfmc, rzmc, prmc, werror, sfmcun, rzmcun, prmcun, dzsf
  integer, allocatable, dimension(:) :: vegcls
  real,    allocatable, dimension(:) :: vegdyn

  type(MAPL_NCIO) :: NCIO(3)
  integer :: i, rc, filetype
    
! Usage
! -----
  if (iargc() /= 4) then
     write(*,*) "Usage: Scale_Catch <Input_Catch> <Regridded_Catch> <Scaled_Catch> <SURFLAY>"
     call exit(2)
  end if

  do n=1,4
  call getarg(n,arg(n))
  enddo

! Open INPUT and Regridded Catch Files
! ------------------------------------
  read(arg(1),'(a)') fname1

  read(arg(2),'(a)') fname2

! Open OUTPUT (Scaled) Catch File
! -------------------------------
  read(arg(3),'(a)') fname3

  call MAPL_NCIOGetFileType(fname1, filetype,rc=rc)

  if (filetype == 0) then
     NCIO(1) = MAPL_NCIOOpen(trim(fname1))
     NCIO(2) = MAPL_NCIOOpen(trim(fname2))
  else
     open(unit=10, file=trim(fname1),  form='unformatted')
     open(unit=20, file=trim(fname2),  form='unformatted')
     open(unit=30, file=trim(fname3),  form='unformatted')
  end if

! Get SURFLAY Value
! -----------------
  read(arg(4),*) SURFLAY

  if (SURFLAY.ne.20 .and. SURFLAY.ne.50) then
     print *, "You must supply a valid SURFLAY value:"
     print *, "(Ganymed-3 and earlier) SURFLAY=20.0 for Old Soil Params"
     print *, "(Ganymed-4 and later  ) SURFLAY=50.0 for New Soil Params"
     call exit(2)
  end if
  print *, 'SURFLAY: ',SURFLAY

  if (filetype ==0) then

     call MAPL_NCIOGetDimSizes(NCIO(1),tile=ntiles)

  else

!    Determine NTILES
!    ----------------
     bpos=0
     read(10)
     epos = ftell(10)            ! ending position of file pointer
     ntiles = (epos-bpos)/4-2    ! record size (in 4 byte words; 
     rewind 10

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
     call readcatch_nc4 ( catch(old), NCIO(old) )
     call readcatch_nc4 ( catch(new), NCIO(new) )
  else
     call readcatch ( 10,catch(old) )
     call readcatch ( 20,catch(new) )
  end if

  allocate( vegcls(ntiles) )
            vegcls(:) = catch(new)%ity(:)

! Create Scaled Catch
! -------------------
  sca = 3
  catch(sca) = catch(new)

  n = count( (catch(old)%catdef .gt. catch(old)%cdcr1) .and. &
             (catch(new)%cdcr2  .gt. catch(old)%cdcr2) )

  write(6,200) n,100*n/ntiles

  where( (catch(old)%catdef .gt. catch(old)%cdcr1) .and. &
         (catch(new)%cdcr2  .gt. catch(old)%cdcr2) )
 
      catch(sca)%rzexc  = catch(old)%rzexc * ( catch(new)%vgwmax / &
                                               catch(old)%vgwmax )

      catch(sca)%catdef = catch(new)%cdcr1 +                     &
                        ( catch(old)%catdef-catch(old)%cdcr1 ) / &
                        ( catch(old)%cdcr2 -catch(old)%cdcr1 ) * &
                        ( catch(new)%cdcr2 -catch(new)%cdcr1 )
  end where

! Sanity Check
! ------------
  print *, 'Performing Sanity Check ...'
  allocate (   dzsf(ntiles) )
  allocate (   sfmc(ntiles) )
  allocate (   rzmc(ntiles) )
  allocate (   prmc(ntiles) )
  allocate ( werror(ntiles) )
  allocate ( sfmcun(ntiles) )
  allocate ( rzmcun(ntiles) )
  allocate ( prmcun(ntiles) )

  dzsf = SURFLAY

  call calc_soil_moist( ntiles,vegcls,dzsf,                                                                            &
       catch(sca)%vgwmax,catch(sca)%cdcr1,catch(sca)%cdcr2,catch(sca)%wpwet,catch(sca)%poros,                          &
       catch(sca)%psis,catch(sca)%bee,catch(sca)%ars1,catch(sca)%ars2,catch(sca)%ars3,catch(sca)%ara1,catch(sca)%ara2, &
       catch(sca)%ara3,catch(sca)%ara4,catch(sca)%arw1,catch(sca)%arw2,catch(sca)%arw3,catch(sca)%arw4,                &
       catch(sca)%srfexc,catch(sca)%rzexc,catch(sca)%catdef,                                                           &
       sfmc, rzmc, prmc, werror, sfmcun, rzmcun, prmcun )
    
  n = count( catch(sca)%catdef .ne. catch(new)%catdef )
  write(6,300) n,100*n/ntiles
  n = count( catch(sca)%srfexc .ne. catch(new)%srfexc )
  write(6,400) n,100*n/ntiles
  n = count( catch(sca)%rzexc  .ne. catch(new)%rzexc  )
  write(6,400) n,100*n/ntiles

! Write Scaled Catch
! ------------------
  if (filetype ==0) then
     NCIO(3) = NCIO(2)
     call MAPL_NCIOClose(NCIO(3))
     call MAPL_NCIOSet(NCIO(3),filename=fname3)
     call MAPL_NCIOCreateFile(NCIO(3))
     call writecatch_nc4 ( catch(sca), NCIO(3) )
  else
     call writecatch ( 30,catch(sca) )
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
       allocate( catch%        ity(ntiles) )
       allocate( catch%         tc(ntiles,4) )
       allocate( catch%         qc(ntiles,4) )
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

   return
   end subroutine allocatch

   subroutine readcatch_nc4 (catch,NCIO)
   type(catch_rst) catch
   type(MAPL_NCIO) :: NCIO

       call MAPL_VarRead(NCIO,"BF1",catch%bf1)
       call MAPL_VarRead(NCIO,"BF2",catch%bf2)
       call MAPL_VarRead(NCIO,"BF3",catch%bf3)
       call MAPL_VarRead(NCIO,"VGWMAX",catch%vgwmax)
       call MAPL_VarRead(NCIO,"CDCR1",catch%cdcr1)
       call MAPL_VarRead(NCIO,"CDCR2",catch%cdcr2)
       call MAPL_VarRead(NCIO,"PSIS",catch%psis)
       call MAPL_VarRead(NCIO,"BEE",catch%bee)
       call MAPL_VarRead(NCIO,"POROS",catch%poros)
       call MAPL_VarRead(NCIO,"WPWET",catch%wpwet)
       call MAPL_VarRead(NCIO,"COND",catch%cond)
       call MAPL_VarRead(NCIO,"GNU",catch%gnu)
       call MAPL_VarRead(NCIO,"ARS1",catch%ars1)
       call MAPL_VarRead(NCIO,"ARS2",catch%ars2)
       call MAPL_VarRead(NCIO,"ARS3",catch%ars3)
       call MAPL_VarRead(NCIO,"ARA1",catch%ara1)
       call MAPL_VarRead(NCIO,"ARA2",catch%ara2)
       call MAPL_VarRead(NCIO,"ARA3",catch%ara3)
       call MAPL_VarRead(NCIO,"ARA4",catch%ara4)
       call MAPL_VarRead(NCIO,"ARW1",catch%arw1)
       call MAPL_VarRead(NCIO,"ARW2",catch%arw2)
       call MAPL_VarRead(NCIO,"ARW3",catch%arw3)
       call MAPL_VarRead(NCIO,"ARW4",catch%arw4)
       call MAPL_VarRead(NCIO,"TSA1",catch%tsa1)
       call MAPL_VarRead(NCIO,"TSA2",catch%tsa2)
       call MAPL_VarRead(NCIO,"TSB1",catch%tsb1)
       call MAPL_VarRead(NCIO,"TSB2",catch%tsb2)
       call MAPL_VarRead(NCIO,"ATAU",catch%atau)
       call MAPL_VarRead(NCIO,"BTAU",catch%btau)
       call MAPL_VarRead(NCIO,"OLD_ITY",catch%ity)
       call MAPL_VarRead(NCIO,"TC",catch%tc)
       call MAPL_VarRead(NCIO,"QC",catch%qc)
       call MAPL_VarRead(NCIO,"OLD_ITY",catch%ity)
       call MAPL_VarRead(NCIO,"CAPAC",catch%capac)
       call MAPL_VarRead(NCIO,"CATDEF",catch%catdef)
       call MAPL_VarRead(NCIO,"RZEXC",catch%rzexc)
       call MAPL_VarRead(NCIO,"SRFEXC",catch%srfexc)
       call MAPL_VarRead(NCIO,"GHTCNT1",catch%ghtcnt1)
       call MAPL_VarRead(NCIO,"GHTCNT2",catch%ghtcnt2)
       call MAPL_VarRead(NCIO,"GHTCNT3",catch%ghtcnt3)
       call MAPL_VarRead(NCIO,"GHTCNT4",catch%ghtcnt4)
       call MAPL_VarRead(NCIO,"GHTCNT5",catch%ghtcnt5)
       call MAPL_VarRead(NCIO,"GHTCNT6",catch%ghtcnt6)
       call MAPL_VarRead(NCIO,"TSURF",catch%tsurf)
       call MAPL_VarRead(NCIO,"WESNN1",catch%wesnn1)
       call MAPL_VarRead(NCIO,"WESNN2",catch%wesnn2)
       call MAPL_VarRead(NCIO,"WESNN3",catch%wesnn3)
       call MAPL_VarRead(NCIO,"HTSNNN1",catch%htsnnn1)
       call MAPL_VarRead(NCIO,"HTSNNN2",catch%htsnnn2)
       call MAPL_VarRead(NCIO,"HTSNNN3",catch%htsnnn3)
       call MAPL_VarRead(NCIO,"SNDZN1",catch%sndzn1)
       call MAPL_VarRead(NCIO,"SNDZN2",catch%sndzn2)
       call MAPL_VarRead(NCIO,"SNDZN3",catch%sndzn3)
       call MAPL_VarRead(NCIO,"CH",catch%ch)
       call MAPL_VarRead(NCIO,"CM",catch%cm)
       call MAPL_VarRead(NCIO,"CQ",catch%cq)
       call MAPL_VarRead(NCIO,"FR",catch%fr)
       call MAPL_VarRead(NCIO,"WW",catch%ww)

   return
   end subroutine readcatch_nc4

   subroutine readcatch (unit,catch)
   integer unit
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
       read(unit) catch%      ity
       read(unit) catch%       tc
       read(unit) catch%       qc
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

   return
   end subroutine readcatch

   subroutine writecatch_nc4 (catch,NCIO)
   type(catch_rst) catch
   type(MAPL_NCIO) :: NCIO

       call MAPL_VarWrite(NCIO,"BF1",catch%bf1)
       call MAPL_VarWrite(NCIO,"BF2",catch%bf2)
       call MAPL_VarWrite(NCIO,"BF3",catch%bf3)
       call MAPL_VarWrite(NCIO,"VGWMAX",catch%vgwmax)
       call MAPL_VarWrite(NCIO,"CDCR1",catch%cdcr1)
       call MAPL_VarWrite(NCIO,"CDCR2",catch%cdcr2)
       call MAPL_VarWrite(NCIO,"PSIS",catch%psis)
       call MAPL_VarWrite(NCIO,"BEE",catch%bee)
       call MAPL_VarWrite(NCIO,"POROS",catch%poros)
       call MAPL_VarWrite(NCIO,"WPWET",catch%wpwet)
       call MAPL_VarWrite(NCIO,"COND",catch%cond)
       call MAPL_VarWrite(NCIO,"GNU",catch%gnu)
       call MAPL_VarWrite(NCIO,"ARS1",catch%ars1)
       call MAPL_VarWrite(NCIO,"ARS2",catch%ars2)
       call MAPL_VarWrite(NCIO,"ARS3",catch%ars3)
       call MAPL_VarWrite(NCIO,"ARA1",catch%ara1)
       call MAPL_VarWrite(NCIO,"ARA2",catch%ara2)
       call MAPL_VarWrite(NCIO,"ARA3",catch%ara3)
       call MAPL_VarWrite(NCIO,"ARA4",catch%ara4)
       call MAPL_VarWrite(NCIO,"ARW1",catch%arw1)
       call MAPL_VarWrite(NCIO,"ARW2",catch%arw2)
       call MAPL_VarWrite(NCIO,"ARW3",catch%arw3)
       call MAPL_VarWrite(NCIO,"ARW4",catch%arw4)
       call MAPL_VarWrite(NCIO,"TSA1",catch%tsa1)
       call MAPL_VarWrite(NCIO,"TSA2",catch%tsa2)
       call MAPL_VarWrite(NCIO,"TSB1",catch%tsb1)
       call MAPL_VarWrite(NCIO,"TSB2",catch%tsb2)
       call MAPL_VarWrite(NCIO,"ATAU",catch%atau)
       call MAPL_VarWrite(NCIO,"BTAU",catch%btau)
       call MAPL_VarWrite(NCIO,"OLD_ITY",catch%ity)
       call MAPL_VarWrite(NCIO,"TC",catch%tc)
       call MAPL_VarWrite(NCIO,"QC",catch%qc)
       call MAPL_VarWrite(NCIO,"OLD_ITY",catch%ity)
       call MAPL_VarWrite(NCIO,"CAPAC",catch%capac)
       call MAPL_VarWrite(NCIO,"CATDEF",catch%catdef)
       call MAPL_VarWrite(NCIO,"RZEXC",catch%rzexc)
       call MAPL_VarWrite(NCIO,"SRFEXC",catch%srfexc)
       call MAPL_VarWrite(NCIO,"GHTCNT1",catch%ghtcnt1)
       call MAPL_VarWrite(NCIO,"GHTCNT2",catch%ghtcnt2)
       call MAPL_VarWrite(NCIO,"GHTCNT3",catch%ghtcnt3)
       call MAPL_VarWrite(NCIO,"GHTCNT4",catch%ghtcnt4)
       call MAPL_VarWrite(NCIO,"GHTCNT5",catch%ghtcnt5)
       call MAPL_VarWrite(NCIO,"GHTCNT6",catch%ghtcnt6)
       call MAPL_VarWrite(NCIO,"TSURF",catch%tsurf)
       call MAPL_VarWrite(NCIO,"WESNN1",catch%wesnn1)
       call MAPL_VarWrite(NCIO,"WESNN2",catch%wesnn2)
       call MAPL_VarWrite(NCIO,"WESNN3",catch%wesnn3)
       call MAPL_VarWrite(NCIO,"HTSNNN1",catch%htsnnn1)
       call MAPL_VarWrite(NCIO,"HTSNNN2",catch%htsnnn2)
       call MAPL_VarWrite(NCIO,"HTSNNN3",catch%htsnnn3)
       call MAPL_VarWrite(NCIO,"SNDZN1",catch%sndzn1)
       call MAPL_VarWrite(NCIO,"SNDZN2",catch%sndzn2)
       call MAPL_VarWrite(NCIO,"SNDZN3",catch%sndzn3)
       call MAPL_VarWrite(NCIO,"CH",catch%ch)
       call MAPL_VarWrite(NCIO,"CM",catch%cm)
       call MAPL_VarWrite(NCIO,"CQ",catch%cq)
       call MAPL_VarWrite(NCIO,"FR",catch%fr)
       call MAPL_VarWrite(NCIO,"WW",catch%ww)

   return
   end subroutine writecatch_nc4

   subroutine writecatch (unit,catch)
   integer unit
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
       write(unit) catch%      ity
       write(unit) catch%       tc
       write(unit) catch%       qc
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

   return
   end subroutine writecatch

  end

  subroutine calc_soil_moist( &
       ncat,vegcls,dzsf,vgwmax,cdcr1,cdcr2,wpwet,poros, &
       psis,bee,ars1,ars2,ars3,ara1,ara2, &
       ara3,ara4,arw1,arw2,arw3,arw4, &
       srfexc,rzexc,catdef, &
       sfmc, rzmc, prmc,  &
       werror, sfmcun, rzmcun, prmcun )
    
    ! Calculate diagnostic soil moisture content from prognostic
    ! excess/deficit variables.
    !
    ! On input, also check validity of prognostic excess/deficit variables
    ! and modify if necessary.  Perturbed or updated excess/deficit variables 
    ! in data assimilation integrations may be unphysical.  
    ! Optional output "werror" contains excess or missing water related
    ! to inconsistency.
    !
    ! Optional outputs "smfcun", "rzmcun", "prmcun" are surface,
    ! root zone, and profile moisture content for unsaturated areas only,
    ! ie. excluding the saturated area of the catchment.    
    !
    ! NOTE: When calling with optional output arguments, use keywords
    !       unless arguments are in proper order!
    !       
    !       Example: 
    !       (don't want "werror" as output, but want "*mcun" output)
    !       
    !       call calc_soil_moist(         & 
    !            ncat, ...                &
    !            sfmc, rzmc, prmc,        &
    !            sfmcun=sfmc_unsat,   &
    !            rzmcun=rzmc_unsat,   & 
    !            prmcun=prmc_unsat )
    !
    ! replaces moisture_sep_22_2003.f (and older moisture.f)
    !
    ! koster+reichle, Feb 5, 2004
    !
    ! revised - koster+reichle, Mar 19, 2004
    !
    ! added optional *un output - koster+reichle, Apr 6, 2004
    !
    ! ----------------------------------------------------------------

    
    implicit none
    
    integer, parameter                   :: KSNGL=4
    integer,                  intent(in) :: ncat
    integer, dimension(ncat), intent(in) :: vegcls
    
    real(KIND=KSNGL), dimension(ncat), intent(in) :: dzsf,vgwmax,cdcr1,cdcr2
    real(KIND=KSNGL), dimension(ncat), intent(in) :: wpwet,poros,psis
    real(KIND=KSNGL), dimension(ncat), intent(in) :: bee,ars1
    real(KIND=KSNGL), dimension(ncat), intent(in) :: ars2,ars3,ara1,ara2,ara3
    real(KIND=KSNGL), dimension(ncat), intent(in) :: ara4,arw1,arw2,arw3,arw4
    
    real(KIND=KSNGL), dimension(ncat), intent(inout) :: srfexc, rzexc, catdef
    
    real(KIND=KSNGL), dimension(ncat), intent(out) :: sfmc, rzmc, prmc
    
    real(KIND=KSNGL), dimension(ncat), intent(out), optional :: werror
    
    real(KIND=KSNGL), dimension(ncat), intent(out), optional :: sfmcun
    real(KIND=KSNGL), dimension(ncat), intent(out), optional :: rzmcun
    real(KIND=KSNGL), dimension(ncat), intent(out), optional :: prmcun
    
    ! ----------------------------
    !    
    ! local variables
    
    integer :: n
    
    real(KIND=KSNGL), parameter :: dtstep_dummy = -9999.
    
    real(KIND=KSNGL), dimension(ncat) :: rzeq, runsrf_dummy, catdef_dummy
    real(KIND=KSNGL), dimension(ncat) :: ar1, ar2, ar4, prmc_orig
    real(KIND=KSNGL), dimension(ncat) :: srfmn, srfmx, swsrf1, swsrf2, swsrf4, rzi
    

    ! --------------------------------------------------------------------
    !
    ! compute soil water storage upon input [mm]
    
    do n=1,ncat
       prmc_orig(n) =                                                 &
            (cdcr2(n)/(1.-wpwet(n))-catdef(n)+rzexc(n)+srfexc(n))
    enddo
       
    ! -----------------------------------
    !
    ! check limits of catchment deficit
    !
    ! increased minimum catchment deficit from 0.01 to 1. to make the 
    ! check work with perturbed parameters and initial condition
    ! reichle, 16 May 01
    !
    ! IT REALLY SHOULD WORK WITH catdef > 0 (rather than >1.) ????
    ! reichle, 5 Feb 2004 
    
    do n=1,ncat     
       catdef(n)=max(1.,min(cdcr2(n),catdef(n)))
    end do
    
    ! ------------------------------------------------------------------
    !     
    ! check limits of root zone excess
    !     
    ! calculate root zone equilibrium moisture for given catchment deficit
    
    call rzequil( &
         ncat, vegcls, catdef, vgwmax, &
         cdcr1, cdcr2, wpwet, &
         ars1, ars2, ars3, ara1, ara2, ara3, ara4, &
         arw1, arw2, arw3, arw4, &
         rzeq)
    
    ! assume srfexc=0 and constrain rzexc appropriately
    ! (iteration would be needed to contrain srfexc and rzexc simultaneously)
    
    do n=1,ncat
       rzexc(n)=max(wpwet(n)*vgwmax(n)-rzeq(n),min(vgwmax(n)-rzeq(n),rzexc(n)))
    end do
    
    ! this translates into:
    !
    ! wilting level < rzmc < porosity
    ! 
    ! or more precisely:  wpwet*vgwmax < rzeq+rzexc < vgwmax
    ! 
    ! NOTE: root zone moisture is not allowed to drop below wilting level 
    
    ! -----------------------------------------------------------------
    !
    ! Call partition() for computation of surface moisture content.
    !
    ! Call to partition() also checks limits of surface excess.
    !
    ! Call partition with dtstep_dummy:
    !  In partition, dtstep is only used for a correction that
    !  puts water into runsrf (for which runsrf_dummy is used here).
    !  Also use catdef_dummy because partition() updates catdef
    !  whenever srfexc exceeds physical bounds, but this is not desired here.
    
    runsrf_dummy = 0.
    catdef_dummy = catdef          
    
    call partition( &
         ncat,dtstep_dummy,vegcls,dzsf,rzexc, &
         rzeq,vgwmax,cdcr1,cdcr2, &
         psis,bee,poros,wpwet, &
         ars1,ars2,ars3, &
         ara1,ara2,ara3,ara4, &
         arw1,arw2,arw3,arw4,.false., &
         srfexc,catdef_dummy,runsrf_dummy, &
         ar1, ar2, ar4,srfmx,srfmn, & 
         swsrf1,swsrf2,swsrf4,rzi &
         )
    
    ! compute surface, root zone, and profile soil moisture
    
    do n=1,ncat

       sfmc(n) = poros(n) *                                           &
            (swsrf1(n)*ar1(n) + swsrf2(n)*ar2(n) + swsrf4(n)*ar4(n))
       
       rzmc(n) = (rzeq(n)+rzexc(n)+srfexc(n))*poros(n)/vgwmax(n)
       
       ! compute revised soil water storage [mm]
       
       prmc(n) =                                                               &
            (cdcr2(n)/(1.-wpwet(n))-catdef(n)+rzexc(n)+srfexc(n))

       ! compute error in soil water storage [mm] (if argument is present)
       
       if (present(werror))  werror(n)=(prmc(n)-prmc_orig(n))
       
       ! convert to volumetric soil moisture
       ! note: dzpr = (cdcr2/(1-wpwet)) / poros 
       
       prmc(n) = prmc(n)*poros(n) / (cdcr2(n)/(1.-wpwet(n)))

       
       ! check for negative soil moisture 
       
       if ( (sfmc(n)<.0) .or. (rzmc(n)<.0) .or. (prmc(n)<.0) ) then
          
          write (*,*) 'FOUND NEGATIVE SOIL MOISTURE CONTENT.... stopping'
          write (*,*) n, sfmc(n), rzmc(n), prmc(n)
          stop
       end if

       ! compute moisture content in unsaturated areas [m3/m3] (if arg present)

       if (ar1(n)<1.) then       

          if (present(prmcun))  prmcun(n)=(prmc(n)-poros(n)*ar1(n))/(1.-ar1(n))
          if (present(rzmcun))  rzmcun(n)=(rzmc(n)-poros(n)*ar1(n))/(1.-ar1(n))
          if (present(sfmcun))  sfmcun(n)=(sfmc(n)-poros(n)*ar1(n))/(1.-ar1(n))

       else          
          
          if (present(prmcun))  prmcun(n)=poros(n)
          if (present(rzmcun))  rzmcun(n)=poros(n)
          if (present(sfmcun))  sfmcun(n)=poros(n)
          
       end if
       
    enddo

  return

  end subroutine calc_soil_moist

      SUBROUTINE PARTITION (                                                   &
                            NCH,DTSTEP,ITYP,DZSF,RZEXC,RZEQ,VGWMAX,CDCR1,CDCR2,&
                            PSIS,BEE,poros,WPWET,                              &
                            ars1,ars2,ars3,ara1,ara2,ara3,ara4,                &
                            arw1,arw2,arw3,arw4,BUG,                           &
                            srfexc,catdef,runsrf,                              &
                            AR1, AR2, AR4, srfmx, srfmn,                       &
                            SWSRF1,SWSRF2,SWSRF4,RZI                           &
                           )

      IMPLICIT NONE

! -------------------------------------------------------------------
      INTEGER, INTENT(IN) :: NCH
      INTEGER, INTENT(IN), DIMENSION(NCH) :: ITYP

      REAL, INTENT(IN) :: DTSTEP
      REAL, INTENT(IN), DIMENSION(NCH) :: DZSF,RZEXC,RZEQ,VGWMAX,CDCR1,CDCR2,  &
                                          PSIS,BEE,poros,WPWET,                &
                                          ars1,ars2,ars3,ara1,ara2,ara3,ara4,  &
                                          arw1,arw2,arw3,arw4

      LOGICAL, INTENT(IN) :: BUG
! -------------------------------------------------------------------
      REAL, INTENT(INOUT), DIMENSION(NCH) :: srfexc,catdef,runsrf
! -------------------------------------------------------------------
      REAL, INTENT(OUT), DIMENSION(NCH) :: AR1, AR2, AR4, srfmx, srfmn,        &
                                           SWSRF1, SWSRF2, SWSRF4, RZI
! -------------------------------------------------------------------
      INTEGER :: N

      REAL :: cor, A150, W150, WMIN, AX, WMNEW, WRZ, TERM1, TERM2, TERM3,      &
              AREA0, AREA1, AREA2, AREA3, AREA4, ASCALE, WILT, D1, D2, CDI,    &
              DELTA1, DELTA2, DELTA4, MULTAR, CATDEFX, RZEQX, RZEQW, FACTOR,   &
              X0, RZEQY, CATDEFW, AR1W, ASUM, RZEQYI, RZEQWI, RZEQXI, AR20,    &
              ARG1, EXPARG1, ARG2, EXPARG2, ARG3, EXPARG3  !, surflay

      LOGICAL :: LSTRESS


      DATA LSTRESS/.FALSE./    !,surflay/20./

!****
!**** --------------------------------------------------       

!rr   next line for debugging, sep 23, 2003, reichle
!rr
!rr   write (*,*) 'entering partition()'

      DO N=1,NCH

        WILT=WPWET(N)
        WRZ=RZEXC(N)/VGWMAX(N)
        CATDEFX=AMIN1( CATDEF(N) , CDCR1(N) )

! CDI DEFINES IF THE SHAPE PARAMETER IS ADJUSTED IN ONE OR TWO SEGMENTS
        if (ara1(n) .ne. ara3(n)) then
            cdi=(ara4(n)-ara2(n))/(ara1(n)-ara3(n))
          else
            cdi=0.
          endif

        AR1(N)= AMIN1(1.,AMAX1(0.,(1.+ars1(n)*CATDEFX)                         &
                 /(1.+ars2(n)*CATDEFX+ars3(n)*CATDEFX*CATDEFX))) 

        if (CATDEFX .ge. cdi) then
            ax=ara3(n)*CATDEFX+ara4(n)
          else
            ax=ara1(n)*CATDEFX+ara2(n)
          endif

        WMIN=AMIN1(1.,AMAX1(0.,arw4(n)+(1.-arw4(n))*                           &
                 (1.+arw1(n)*CATDEFX)                                          &
                 /(1.+arw2(n)*CATDEFX+arw3(n)*CATDEFX*CATDEFX)))

!**** CRITICAL VALUE 1: AVERAGE MOISTURE IN ROOT ZONE AT WMIN
!**** ASSOCIATED WITH CATDEF.
        ARG1=AMAX1(-40., AMIN1(40., -AX*(1.-WMIN)))
        EXPARG1=EXP(ARG1)
        RZEQX=(WMIN-1.-(2./AX))*EXPARG1 + WMIN + (2./AX)
        RZEQXI=AX*EXPARG1 *                                                    &
          ( -1. -(2./AX) - (2./(AX*AX)) + WMIN + (WMIN/AX) )                   &
          + WMIN + 2./AX
        AR20=1.+(-AX-1.+AX*WMIN)*EXPARG1
        RZEQXI=RZEQXI/(AR20+1.E-20)

!**** CRITICAL VALUE 2: AVERAGE MOISTURE IN ROOT ZONE WHEN WMIN
!**** IS EXACTLY AT WILTING POINT.
        ARG2=AMAX1(-40., AMIN1(40., -AX*(1.-WILT)))
        EXPARG2=EXP(ARG2)
        RZEQW=(WILT-1.-(2./AX))*EXPARG2 + WILT + (2./AX)
        RZEQWI=AX*EXPARG2 *                                                    &
         ( -1. -(2./AX) - (2./(AX*AX)) + WILT + (WILT/AX) )                    &
         + WILT + 2./AX
        AR20=1.+(-AX-1.+AX*WILT)*EXPARG2
        RZEQWI=RZEQWI/(AR20+1.E-20)

!**** SITUATION 1: CATDEF LE CDCR1
        IF(CATDEF(N) .LE. CDCR1(N)) THEN
          RZEQY=RZEQX+WRZ
          RZEQYI=RZEQXI+WRZ
          WMNEW=WMIN+WRZ
          ARG3=AMAX1(-40., AMIN1(40., -AX*(1.-WMNEW)))
          EXPARG3=EXP(ARG3)
          AREA1=(1.+AX-AX*WMIN)*EXPARG1
          AREA2=(1.+AX-AX*WMNEW)*EXPARG3
          IF(WMNEW .GE. WILT) THEN
            AR1(N)=AR1(N)+AREA2-AREA1
            AR2(N)=1.-AR1(N)
            AR4(N)=0.
            ENDIF
          IF(WMNEW .LT. WILT) THEN
            AREA3=(1.+AX-AX*WILT)*EXPARG2
            AR1(N)=AR1(N)+AREA3-AREA1
            AR2(N)=1.-AR1(N)
            FACTOR=(RZEQX+WRZ-WILT)/(RZEQW-WILT)
            AR1(N)=AR1(N)*FACTOR
            AR2(N)=AR2(N)*FACTOR
            AR4(N)=1.-FACTOR
            ENDIF
          ENDIF

!**** SITUATION 2: CATDEF GT CDCR1
        IF(CATDEF(N) .GT. CDCR1(N)) THEN
          FACTOR=(CDCR2(N)-CATDEF(N))/(CDCR2(N)-CDCR1(N))
          RZEQY=WILT+(RZEQX-WILT)*FACTOR+WRZ
          RZEQYI=WILT+(RZEQXI-WILT)*FACTOR+WRZ

          IF(RZEQY .LT. WILT) THEN
            IF(RZEQY .LT. WILT-.001) THEN
!rr                WRITE(*,*) 'RZEXC WAY TOO LOW!  N=',N,' RZEQY=',RZEQY
!rr                WRITE(*,*) 'SRFEXC=',SRFEXC(N),'RZEXC=',RZEXC(N),
!rr     &                     'CATDEF=',CATDEF(N)
!             ELSE
!               WRITE(*,*) 'RZEXC TOO LOW  N=',N
              ENDIF
            RZEQY=WILT
            RZEQYI=WILT
            ENDIF

          IF(RZEQY .GE. RZEQX) THEN  ! RZEXC BRINGS MOISTURE ABOVE CDCR1 POINT
            WMNEW=WMIN+(RZEQY-RZEQX)
            ARG3=AMAX1(-40., AMIN1(40., -AX*(1.-WMNEW)))
            EXPARG3=EXP(ARG3)
            AREA1=(1.+AX-AX*WMIN)*EXPARG1
            AREA2=(1.+AX-AX*WMNEW)*EXPARG3
            AR1(N)=AR1(N)+(AREA2-AREA1)
            AR2(N)=1.-AR1(N)
            AR4(N)=0.
            ENDIF

          IF(RZEQY .LT. RZEQX .AND. RZEQY .GE. RZEQW) THEN
            CATDEFW=CDCR2(N)+((RZEQW-WILT)/(RZEQX-WILT))*(CDCR1(N)-CDCR2(N))
            AR1W= AMIN1(1.,AMAX1(0.,(1.+ars1(n)*CATDEFW)                       &
                 /(1.+ars2(n)*CATDEFW+ars3(n)*CATDEFW*CATDEFW)))
            FACTOR=(RZEQY-RZEQW)/(RZEQX-RZEQW)
            AR1(N)=AR1W+FACTOR*(AR1(N)-AR1W)
            AR2(N)=1.-AR1(N)
            AR4(N)=0.
            ENDIF

          IF(RZEQY .LT. RZEQW) THEN
            CATDEFW=CDCR2(N)+((RZEQW-WILT)/(RZEQX-WILT))*(CDCR1(N)-CDCR2(N))
            AR1W= AMIN1(1.,AMAX1(0.,(1.+ars1(n)*CATDEFW)                       &
                 /(1.+ars2(n)*CATDEFW+ars3(n)*CATDEFW*CATDEFW)))
            AR1(N)=AR1W
            AR2(N)=1.-AR1(N)
            FACTOR=(RZEQY-WILT)/(RZEQW-WILT)
            AR1(N)=AR1(N)*FACTOR
            AR2(N)=AR2(N)*FACTOR
            AR4(N)=1.-FACTOR
            ENDIF

          ENDIF

        RZI(N)=RZEQYI

        SWSRF1(N)=1.
!mjs: changed .001 temporarily because of large bee.
        SWSRF2(N)=AMIN1(1., AMAX1(0.01, RZEQYI))
        SWSRF4(N)=AMIN1(1., AMAX1(0.01, WILT))

!**** EXTRAPOLATION OF THE SURFACE WETNESSES

! 1st step: surface wetness in the unstressed fraction without considering
!           the surface excess; we just assume an equilibrium profile from 
!           the middle of the root zone to the surface.

        SWSRF2(N)=((SWSRF2(N)**(-BEE(N))) - (.5/PSIS(N)))**(-1./BEE(N))
        SWSRF4(N)=((SWSRF4(N)**(-BEE(N))) - (.5/PSIS(N)))**(-1./BEE(N))

! srfmx is the maximum amount of water that can be added to the surface layer
! The choice of defining SWSRF4 like SWSRF2 needs to be better examined.
        srfmx(n)=ar2(n)*(1.-swsrf2(n))*(dzsf(n)*poros(n))
        srfmx(n)=srfmx(n)+ar4(n)*(1.-swsrf4(n))*(dzsf(n)*poros(n))
!**** For calculation of srfmn, assume surface moisture associated with
!**** AR1 is constantly replenished by water table.
        srfmn(n)=-(ar2(n)*swsrf2(n)+ar4(n)*swsrf4(n))*(dzsf(n)*poros(n))

        if(srfexc(n).gt.srfmx(n)) then
            cor=srfexc(n)-srfmx(n)     !  The correction is here
            srfexc(n)=srfmx(n)
            catdef(n)=catdef(n)-cor
            if(catdef(n).lt.0.) then
              runsrf(n)=runsrf(n)-catdef(n)/dtstep
              catdef(n)=0.
              endif
          else if(srfexc(n).lt.srfmn(n)) then
            cor=srfexc(n)-srfmn(n)
            catdef(n)=catdef(n)-cor
            srfexc(n)=srfmn(n)
          else
            cor=0.
          endif
          
        SWSRF2(N)=SWSRF2(N)+SRFEXC(N)/(dzsf(n)*poros(n)*(1.-ar1(n))+1.e-20)
        SWSRF2(N)=AMIN1(1., AMAX1(1.E-5, SWSRF2(N)))
        swsrf4(n)=swsrf4(n)+srfexc(n)/(dzsf(n)*poros(n)*(1.-ar1(n))+1.e-20)
        SWSRF4(N)=AMIN1(1., AMAX1(1.E-5, SWSRF4(N)))

        IF (AR1(N) .ge. 1.-1.E-5) then
          AR1(N)=1.
          AR2(N)=0.
          AR4(N)=0.
          SWSRF2(N)=1.
          SWSRF4(N)=wilt
          ENDIF

        IF (AR1(N) .LT. 0.) then
!rr          IF(AR1(N) .LT. -1.E-3) WRITE(*,*) 'AR1 TOO LOW: AR1=',AR1(N)
          AR1(N)=0.
          ENDIF
        ar1(n)=amax1(0., amin1(1., ar1(n)))
        ar2(n)=amax1(0., amin1(1., ar2(n)))
        ar4(n)=amax1(0., amin1(1., ar4(n)))
        asum=ar1(n)+ar2(n)+ar4(n)
        if(asum .lt. .9999 .or. asum .gt. 1.0001) then
          write(*,*) 'Areas do not add to 1: sum=',asum,'N=',n
       endif


        ENDDO

         
      RETURN
      END SUBROUTINE PARTITION

      SUBROUTINE RZEQUIL (                                                     &
                          NCH,ITYP,CATDEF,VGWMAX,CDCR1,CDCR2,WPWET,            &
                          ars1,ars2,ars3,ara1,ara2,ara3,ara4,                  &
                          arw1,arw2,arw3,arw4,                                 &
                          RZEQ                                                 &
                         )

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NCH
      INTEGER, INTENT(IN), DIMENSION(NCH) :: ITYP
      REAL, INTENT(IN), DIMENSION(NCH) :: CATDEF, VGWMAX, CDCR1, CDCR2,        &
                   WPWET, ars1, ars2, ars3, ara1, ara2, ara3, ara4, arw1,      &
                   arw2, arw3, arw4

      REAL, INTENT(OUT), DIMENSION(NCH) :: RZEQ

      INTEGER N
      REAL AX,WMIN,ASCALE,cdi,wilt,catdefx,factor,ARG1,EXPARG1

! ----------------------------------------------------------------------

      DO N=1,NCH

        WILT=WPWET(N)
        CATDEFX=AMIN1( CATDEF(N) , CDCR1(N) )

! CDI DEFINES IF THE SHAPE PARAMETER IS ADJUSTED IN ONE OR TWO SEGMENTS
        if (ara1(n) .ne. ara3(n)) then
            cdi=(ara4(n)-ara2(n))/(ara1(n)-ara3(n))
          else
            cdi=0.
          endif

        if (CATDEFX .ge. cdi) then
            ax=ara3(n)*CATDEFX+ara4(n)
          else
            ax=ara1(n)*CATDEFX+ara2(n)
          endif

        WMIN=AMIN1(1.,AMAX1(0.,arw4(n)+(1.-arw4(n))*(1.+arw1(n)*CATDEFX)       &
                 /(1.+arw2(n)*CATDEFX+arw3(n)*CATDEFX*CATDEFX)))

        ARG1=AMAX1(-40., AMIN1(40., -AX*(1.-WMIN)))
        EXPARG1=EXP(ARG1)
        RZEQ(N)=(WMIN-1.-(2./AX))*EXPARG1 + WMIN + (2./AX)

        IF(CATDEF(N) .GT. CDCR1(N)) THEN
          FACTOR=(CDCR2(N)-CATDEF(N))/(CDCR2(N)-CDCR1(N))
          RZEQ(N)=WILT+(RZEQ(N)-WILT)*FACTOR
          ENDIF

! scaling:    
        RZEQ(N)=AMIN1(1.,AMAX1(0.,RZEQ(N)))
        RZEQ(N)=RZEQ(N)*VGWMAX(N)

      ENDDO

      RETURN
      END SUBROUTINE RZEQUIL

