#include "MAPL_Generic.h"

module CatchmentRstMod
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
#ifndef __GFORTRAN__
  integer           :: ftell
  external          :: ftell
#endif
  type :: CatchmentRst
     integer :: ntiles
     real, allocatable ::        bf1(:)
     real, allocatable ::        bf2(:)
     real, allocatable ::        bf3(:)
     real, allocatable ::     vgwmax(:)
     real, allocatable ::      cdcr1(:)
     real, allocatable ::      cdcr2(:)
     real, allocatable ::       psis(:)
     real, allocatable ::        bee(:)
     real, allocatable ::      poros(:)
     real, allocatable ::      wpwet(:)
     real, allocatable ::       cond(:)
     real, allocatable ::        gnu(:)
     real, allocatable ::       ars1(:)
     real, allocatable ::       ars2(:)
     real, allocatable ::       ars3(:)
     real, allocatable ::       ara1(:)
     real, allocatable ::       ara2(:)
     real, allocatable ::       ara3(:)
     real, allocatable ::       ara4(:)
     real, allocatable ::       arw1(:)
     real, allocatable ::       arw2(:)
     real, allocatable ::       arw3(:)
     real, allocatable ::       arw4(:)
     real, allocatable ::       tsa1(:)
     real, allocatable ::       tsa2(:)
     real, allocatable ::       tsb1(:)
     real, allocatable ::       tsb2(:)
     real, allocatable ::       atau(:)
     real, allocatable ::       btau(:)
     real, allocatable ::        ity(:)
     real, allocatable ::         tc(:,:)
     real, allocatable ::         qc(:,:)
     real, allocatable ::      capac(:)
     real, allocatable ::     catdef(:)
     real, allocatable ::      rzexc(:)
     real, allocatable ::     srfexc(:)
     real, allocatable ::    ghtcnt1(:)
     real, allocatable ::    ghtcnt2(:)
     real, allocatable ::    ghtcnt3(:)
     real, allocatable ::    ghtcnt4(:)
     real, allocatable ::    ghtcnt5(:)
     real, allocatable ::    ghtcnt6(:)
     real, allocatable ::      tsurf(:)
     real, allocatable ::     wesnn1(:)
     real, allocatable ::     wesnn2(:)
     real, allocatable ::     wesnn3(:)
     real, allocatable ::    htsnnn1(:)
     real, allocatable ::    htsnnn2(:)
     real, allocatable ::    htsnnn3(:)
     real, allocatable ::     sndzn1(:)
     real, allocatable ::     sndzn2(:)
     real, allocatable ::     sndzn3(:)
     real, allocatable ::         ch(:,:)
     real, allocatable ::         cm(:,:)
     real, allocatable ::         cq(:,:)
     real, allocatable ::         fr(:,:)
     real, allocatable ::         ww(:,:)
  contains
     procedure :: read_GEOSldas_rst_bin      
     procedure :: write_nc4      
     procedure :: read_shared_nc4      
     procedure :: write_shared_nc4
     procedure :: add_bcs_to_rst      
     procedure :: allocate_catch
  endtype CatchmentRst

  interface CatchmentRst
     module procedure CatchmentRst_Create
  end interface

contains

  function CatchmentRst_create(filename, rc) result (catch)
    type(CatchmentRst) :: catch
    character(*), intent(in)           :: filename
    integer, optional, intent(out) :: rc
    integer :: status
    character(len=256) :: Iam = "CatchmentRst_create"
    type(Netcdf4_fileformatter) :: formatter
    integer :: filetype, ntiles, unit
    type(FileMetadata) :: meta
    integer :: bpos, epos, n

    call MAPL_NCIOGetFileType(filename, filetype, __RC__)

     if(filetype == 0) then
       ! nc4 format
       call formatter%open(filename, pFIO_READ, __RC__)
       meta   = formatter%read(__RC__)
       ntiles = meta%get_dimension('tile', __RC__)
       this%ntiles = ntiles
       call catch%allocate_catch()
       call catch%read_shared_nc4(formatter, __RC__)
       call MAPL_VarRead(formatter,"OLD_ITY",catch%ity, __RC__)
       call formatter%close()
     else
       !GEOSldas binary
       open(newunit=unit, file=filename,  form='unformatted', action = 'read')
       bpos=0
       read(unit)
       epos = ftell(unit)          ! ending position of file pointer
       close(unit)
       ntiles = (epos-bpos)/4-2    ! record size (in 4 byte words; 
       this%ntiles = ntiles
       call catch%allocate_catch()
       call catch%read_GEOSldas_rst_bin(filename, __RC__)
     endif
     _RETURN(_SUCCESS)
   end function CatchmentRst_Create

   subroutine read_GEOSldas_rst_bin(this, filename, rc)
      class(CatchmentRst), intent(inout) :: this
      character(*), intent(in) :: filename
      integer, optional, intent(out) :: rc
      integer :: unit
      open(newunit=unit, file=filename, form='unformatted', action='read')
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
      close(unit)
      _RETURN(_SUCCESS)
   end subroutine read_GEOSldas_rst_bin

   subroutine read_shared_nc4(this, formatter, rc)
     class(CatchmentRst), intent(inout):: this
     type(Netcdf4_fileformatter),intent(inout) :: formatter
     integer, optional, intent(out):: rc
     integer :: status
     
     call MAPL_VarRead(formatter,"BF1",this%bf1, __RC__)
     call MAPL_VarRead(formatter,"BF2",this%bf2, __RC__)
     call MAPL_VarRead(formatter,"BF3",this%bf3, __RC__)
     call MAPL_VarRead(formatter,"VGWMAX",this%vgwmax, __RC__)
     call MAPL_VarRead(formatter,"CDCR1",this%cdcr1, __RC__)
     call MAPL_VarRead(formatter,"CDCR2",this%cdcr2, __RC__)
     call MAPL_VarRead(formatter,"PSIS",this%psis, __RC__)
     call MAPL_VarRead(formatter,"BEE",this%bee, __RC__)
     call MAPL_VarRead(formatter,"POROS",this%poros, __RC__)
     call MAPL_VarRead(formatter,"WPWET",this%wpwet, __RC__)
     call MAPL_VarRead(formatter,"COND",this%cond, __RC__)
     call MAPL_VarRead(formatter,"GNU",this%gnu, __RC__)
     call MAPL_VarRead(formatter,"ARS1",this%ars1, __RC__)
     call MAPL_VarRead(formatter,"ARS2",this%ars2, __RC__)
     call MAPL_VarRead(formatter,"ARS3",this%ars3, __RC__)
     call MAPL_VarRead(formatter,"ARA1",this%ara1, __RC__)
     call MAPL_VarRead(formatter,"ARA2",this%ara2, __RC__)
     call MAPL_VarRead(formatter,"ARA3",this%ara3, __RC__)
     call MAPL_VarRead(formatter,"ARA4",this%ara4, __RC__)
     call MAPL_VarRead(formatter,"ARW1",this%arw1, __RC__)
     call MAPL_VarRead(formatter,"ARW2",this%arw2, __RC__)
     call MAPL_VarRead(formatter,"ARW3",this%arw3, __RC__)
     call MAPL_VarRead(formatter,"ARW4",this%arw4, __RC__)
     call MAPL_VarRead(formatter,"TSA1",this%tsa1, __RC__)
     call MAPL_VarRead(formatter,"TSA2",this%tsa2, __RC__)
     call MAPL_VarRead(formatter,"TSB1",this%tsb1, __RC__)
     call MAPL_VarRead(formatter,"TSB2",this%tsb2, __RC__)
     call MAPL_VarRead(formatter,"ATAU",this%atau, __RC__)
     call MAPL_VarRead(formatter,"BTAU",this%btau, __RC__)
     call MAPL_VarRead(formatter,"TC",this%tc, __RC__)
     call MAPL_VarRead(formatter,"QC",this%qc, __RC__)
!
!     call MAPL_VarRead(formatter,"OLD_ITY",this%ity, __RC__)
!
     call MAPL_VarRead(formatter,"CAPAC",this%capac, __RC__)
     call MAPL_VarRead(formatter,"CATDEF",this%catdef, __RC__)
     call MAPL_VarRead(formatter,"RZEXC",this%rzexc, __RC__)
     call MAPL_VarRead(formatter,"SRFEXC",this%srfexc, __RC__)
     call MAPL_VarRead(formatter,"GHTCNT1",this%ghtcnt1, __RC__)
     call MAPL_VarRead(formatter,"GHTCNT2",this%ghtcnt2, __RC__)
     call MAPL_VarRead(formatter,"GHTCNT3",this%ghtcnt3, __RC__)
     call MAPL_VarRead(formatter,"GHTCNT4",this%ghtcnt4, __RC__)
     call MAPL_VarRead(formatter,"GHTCNT5",this%ghtcnt5, __RC__)
     call MAPL_VarRead(formatter,"GHTCNT6",this%ghtcnt6, __RC__)
     call MAPL_VarRead(formatter,"TSURF",this%tsurf, __RC__)
     call MAPL_VarRead(formatter,"WESNN1",this%wesnn1, __RC__)
     call MAPL_VarRead(formatter,"WESNN2",this%wesnn2, __RC__)
     call MAPL_VarRead(formatter,"WESNN3",this%wesnn3, __RC__)
     call MAPL_VarRead(formatter,"HTSNNN1",this%htsnnn1, __RC__)
     call MAPL_VarRead(formatter,"HTSNNN2",this%htsnnn2, __RC__)
     call MAPL_VarRead(formatter,"HTSNNN3",this%htsnnn3, __RC__)
     call MAPL_VarRead(formatter,"SNDZN1",this%sndzn1, __RC__)
     call MAPL_VarRead(formatter,"SNDZN2",this%sndzn2, __RC__)
     call MAPL_VarRead(formatter,"SNDZN3",this%sndzn3, __RC__)
     call MAPL_VarRead(formatter,"CH",this%ch, __RC__)
     call MAPL_VarRead(formatter,"CM",this%cm, __RC__)
     call MAPL_VarRead(formatter,"CQ",this%cq, __RC__)
     call MAPL_VarRead(formatter,"FR",this%fr, __RC__)
     call MAPL_VarRead(formatter,"WW",this%ww, __RC__)
 
      _RETURN(_SUCCESS)
   end subroutine

   subroutine write_nc4 (this, filename, meta, cnclm, rc)
     class(CatchmentRst), intent(inout):: this
     character(*), intent(in) :: filename
     type(FileMetaData), intent(inout) :: meta
     character(*), intent(in) :: cnclm
     integer, optional, intent(out):: rc

     type(Netcdf4_fileformatter) :: formatter
     integer :: status
     character(256) :: Iam = "write_nc4"

     call formatter%create(filename, __RC__)
     call formatter%write(meta, __RC__)
     call this%write_shared_nc4(formatter, __RC__)
     call MAPL_VarWrite(formatter,"OLD_ITY",this%ity)
     call formatter%close()
     _RETURN(_SUCCESS)

   end subroutine write_nc4

   subroutine write_shared_nc4(this, formatter, rc)
     class(CatchmentRst), intent(inout):: this
     type(Netcdf4_fileformatter),intent(inout) :: formatter
     integer, optional, intent(out):: rc
     integer :: status
     call MAPL_VarWrite(formatter,"BF1",this%bf1)
     call MAPL_VarWrite(formatter,"BF2",this%bf2)
     call MAPL_VarWrite(formatter,"BF3",this%bf3)
     call MAPL_VarWrite(formatter,"VGWMAX",this%vgwmax)
     call MAPL_VarWrite(formatter,"CDCR1",this%cdcr1)
     call MAPL_VarWrite(formatter,"CDCR2",this%cdcr2)
     call MAPL_VarWrite(formatter,"PSIS",this%psis)
     call MAPL_VarWrite(formatter,"BEE",this%bee)
     call MAPL_VarWrite(formatter,"POROS",this%poros)
     call MAPL_VarWrite(formatter,"WPWET",this%wpwet)
     call MAPL_VarWrite(formatter,"COND",this%cond)
     call MAPL_VarWrite(formatter,"GNU",this%gnu)
     call MAPL_VarWrite(formatter,"ARS1",this%ars1)
     call MAPL_VarWrite(formatter,"ARS2",this%ars2)
     call MAPL_VarWrite(formatter,"ARS3",this%ars3)
     call MAPL_VarWrite(formatter,"ARA1",this%ara1)
     call MAPL_VarWrite(formatter,"ARA2",this%ara2)
     call MAPL_VarWrite(formatter,"ARA3",this%ara3)
     call MAPL_VarWrite(formatter,"ARA4",this%ara4)
     call MAPL_VarWrite(formatter,"ARW1",this%arw1)
     call MAPL_VarWrite(formatter,"ARW2",this%arw2)
     call MAPL_VarWrite(formatter,"ARW3",this%arw3)
     call MAPL_VarWrite(formatter,"ARW4",this%arw4)
     call MAPL_VarWrite(formatter,"TSA1",this%tsa1)
     call MAPL_VarWrite(formatter,"TSA2",this%tsa2)
     call MAPL_VarWrite(formatter,"TSB1",this%tsb1)
     call MAPL_VarWrite(formatter,"TSB2",this%tsb2)
     call MAPL_VarWrite(formatter,"ATAU",this%atau)
     call MAPL_VarWrite(formatter,"BTAU",this%btau)
     call MAPL_VarWrite(formatter,"TC",this%tc)
     call MAPL_VarWrite(formatter,"QC",this%qc)
     call MAPL_VarWrite(formatter,"CAPAC",this%capac)
     call MAPL_VarWrite(formatter,"CATDEF",this%catdef)
     call MAPL_VarWrite(formatter,"RZEXC",this%rzexc)
     call MAPL_VarWrite(formatter,"SRFEXC",this%srfexc)
     call MAPL_VarWrite(formatter,"GHTCNT1",this%ghtcnt1)
     call MAPL_VarWrite(formatter,"GHTCNT2",this%ghtcnt2)
     call MAPL_VarWrite(formatter,"GHTCNT3",this%ghtcnt3)
     call MAPL_VarWrite(formatter,"GHTCNT4",this%ghtcnt4)
     call MAPL_VarWrite(formatter,"GHTCNT5",this%ghtcnt5)
     call MAPL_VarWrite(formatter,"GHTCNT6",this%ghtcnt6)
     call MAPL_VarWrite(formatter,"TSURF",this%tsurf)
     call MAPL_VarWrite(formatter,"WESNN1",this%wesnn1)
     call MAPL_VarWrite(formatter,"WESNN2",this%wesnn2)
     call MAPL_VarWrite(formatter,"WESNN3",this%wesnn3)
     call MAPL_VarWrite(formatter,"HTSNNN1",this%htsnnn1)
     call MAPL_VarWrite(formatter,"HTSNNN2",this%htsnnn2)
     call MAPL_VarWrite(formatter,"HTSNNN3",this%htsnnn3)
     call MAPL_VarWrite(formatter,"SNDZN1",this%sndzn1)
     call MAPL_VarWrite(formatter,"SNDZN2",this%sndzn2)
     call MAPL_VarWrite(formatter,"SNDZN3",this%sndzn3)
     call MAPL_VarWrite(formatter,"CH",this%ch)
     call MAPL_VarWrite(formatter,"CM",this%cm)
     call MAPL_VarWrite(formatter,"CQ",this%cq)
     call MAPL_VarWrite(formatter,"FR",this%fr)
     call MAPL_VarWrite(formatter,"WW",this%ww)
     _RETURN(_SUCCESS)

   end subroutine write_shared_nc4

   subroutine allocate_catch(this,rc)
     class(CatchmentRst), intent(inout) :: this
     integer, optional, intent(out):: rc
     integer :: ntiles
     ntiles = this%ntiles
     allocate( this%        bf1(ntiles) )
     allocate( this%        bf2(ntiles) )
     allocate( this%        bf3(ntiles) )
     allocate( this%     vgwmax(ntiles) )
     allocate( this%      cdcr1(ntiles) )
     allocate( this%      cdcr2(ntiles) )
     allocate( this%       psis(ntiles) )
     allocate( this%        bee(ntiles) )
     allocate( this%      poros(ntiles) )
     allocate( this%      wpwet(ntiles) )
     allocate( this%       cond(ntiles) )
     allocate( this%        gnu(ntiles) )
     allocate( this%       ars1(ntiles) )
     allocate( this%       ars2(ntiles) )
     allocate( this%       ars3(ntiles) )
     allocate( this%       ara1(ntiles) )
     allocate( this%       ara2(ntiles) )
     allocate( this%       ara3(ntiles) )
     allocate( this%       ara4(ntiles) )
     allocate( this%       arw1(ntiles) )
     allocate( this%       arw2(ntiles) )
     allocate( this%       arw3(ntiles) )
     allocate( this%       arw4(ntiles) )
     allocate( this%       tsa1(ntiles) )
     allocate( this%       tsa2(ntiles) )
     allocate( this%       tsb1(ntiles) )
     allocate( this%       tsb2(ntiles) )
     allocate( this%       atau(ntiles) )
     allocate( this%       btau(ntiles) )
     allocate( this%        ity(ntiles) )
     allocate( this%         tc(ntiles,4) )
     allocate( this%         qc(ntiles,4) )
     allocate( this%      capac(ntiles) )
     allocate( this%     catdef(ntiles) )
     allocate( this%      rzexc(ntiles) )
     allocate( this%     srfexc(ntiles) )
     allocate( this%    ghtcnt1(ntiles) )
     allocate( this%    ghtcnt2(ntiles) )
     allocate( this%    ghtcnt3(ntiles) )
     allocate( this%    ghtcnt4(ntiles) )
     allocate( this%    ghtcnt5(ntiles) )
     allocate( this%    ghtcnt6(ntiles) )
     allocate( this%      tsurf(ntiles) )
     allocate( this%     wesnn1(ntiles) )
     allocate( this%     wesnn2(ntiles) )
     allocate( this%     wesnn3(ntiles) )
     allocate( this%    htsnnn1(ntiles) )
     allocate( this%    htsnnn2(ntiles) )
     allocate( this%    htsnnn3(ntiles) )
     allocate( this%     sndzn1(ntiles) )
     allocate( this%     sndzn2(ntiles) )
     allocate( this%     sndzn3(ntiles) )
     allocate( this%         ch(ntiles,4) )
     allocate( this%         cm(ntiles,4) )
     allocate( this%         cq(ntiles,4) )
     allocate( this%         fr(ntiles,4) )
     allocate( this%         ww(ntiles,4) )
     _RETURN(_SUCCESS)
   end subroutine allocate

   ! This subroutine reads BCs from BCSDIR and hydrological varable
   subroutine add_bcs_to_rst(this, surflay, DataDir, rc)
      class(CatchmentRst), intent(inout) :: this
      real, intent(in) :: surflay
      integer, optional, intent(out) :: rc

      real, allocatable ::  DP2BR(:)
      real :: CanopH
      integer, allocatable :: ity(:)
      integer       :: ntiles, STATUS
      integer       :: idum, i,j,n, ib, nv
      real          :: rdum, zdep1, zdep2, zdep3, zmet, term1, term2, bare,fvg(4)
      logical       :: NEWLAND
      logical       :: file_exists

      type(NetCDF4_Fileformatter) :: CatchFmt

      character*256        :: Iam = "add_bcs"

      ntiles = this%ntiles

      allocate (DP2BR(ntiles), ity(ntiles) )

      inquire(file = trim(DataDir)//'/catch_params.nc4', exist=file_exists)
      inquire(file = trim(DataDir)//"CLM_veg_typs_fracs",exist=NewLand )

      if(file_exists) then
        call CatchFmt%Open(trim(DataDir)//'/catch_params.nc4', pFIO_READ, __RC__)
        call MAPL_VarRead ( CatchFmt ,'OLD_ITY', this%ITY, __RC__)
        call MAPL_VarRead ( CatchFmt ,'ARA1', this%ARA1, __RC__)
        call MAPL_VarRead ( CatchFmt ,'ARA2', this%ARA2, __RC__)
        call MAPL_VarRead ( CatchFmt ,'ARA3', this%ARA3, __RC__)
        call MAPL_VarRead ( CatchFmt ,'ARA4', this%ARA4, __RC__)
        call MAPL_VarRead ( CatchFmt ,'ARS1', this%ARS1, __RC__)
        call MAPL_VarRead ( CatchFmt ,'ARS2', this%ARS2, __RC__)
        call MAPL_VarRead ( CatchFmt ,'ARS3', this%ARS3, __RC__)
        call MAPL_VarRead ( CatchFmt ,'ARW1', this%ARW1, __RC__)
        call MAPL_VarRead ( CatchFmt ,'ARW2', this%ARW2, __RC__)
        call MAPL_VarRead ( CatchFmt ,'ARW3', this%ARW3, __RC__)
        call MAPL_VarRead ( CatchFmt ,'ARW4', this%ARW4, __RC__)

        if( SURFLAY.eq.20.0 ) then
          call MAPL_VarRead ( CatchFmt ,'ATAU2', this%ATAU, __RC__)
          call MAPL_VarRead ( CatchFmt ,'BTAU2', this%BTAU, __RC__)
        endif

        if( SURFLAY.eq.50.0 ) then
          call MAPL_VarRead ( CatchFmt ,'ATAU5', this%ATAU, __RC__)
          call MAPL_VarRead ( CatchFmt ,'BTAU5', this%BTAU, __RC__)
        endif

        call MAPL_VarRead ( CatchFmt ,'PSIS', this%PSIS, __RC__)
        call MAPL_VarRead ( CatchFmt ,'BEE', this%BEE, __RC__)
        call MAPL_VarRead ( CatchFmt ,'BF1', this%BF1, __RC__)
        call MAPL_VarRead ( CatchFmt ,'BF2', this%BF2, __RC__)
        call MAPL_VarRead ( CatchFmt ,'BF3', this%BF3, __RC__)
        call MAPL_VarRead ( CatchFmt ,'TSA1', this%TSA1, __RC__)
        call MAPL_VarRead ( CatchFmt ,'TSA2', this%TSA2, __RC__)
        call MAPL_VarRead ( CatchFmt ,'TSB1', this%TSB1, __RC__)
        call MAPL_VarRead ( CatchFmt ,'TSB2', this%TSB2, __RC__)
        call MAPL_VarRead ( CatchFmt ,'COND', this%COND, __RC__)
        call MAPL_VarRead ( CatchFmt ,'GNU', this%GNU, __RC__)
        call MAPL_VarRead ( CatchFmt ,'WPWET', this%WPWET, __RC__)
        call MAPL_VarRead ( CatchFmt ,'DP2BR', DP2BR, __RC__)
        call MAPL_VarRead ( CatchFmt ,'POROS', this%POROS, __RC__)
        call CatchFmt%close()
      else
        open(unit=21, file=trim(DataDir)//'mosaic_veg_typs_fracs',form='formatted')
        open(unit=22, file=trim(DataDir)//'bf.dat'               ,form='formatted')
        open(unit=23, file=trim(DataDir)//'soil_param.dat'       ,form='formatted')
        open(unit=24, file=trim(DataDir)//'ar.new'               ,form='formatted')
        open(unit=25, file=trim(DataDir)//'ts.dat'               ,form='formatted')
        open(unit=26, file=trim(DataDir)//'tau_param.dat'        ,form='formatted')

        do n=1,ntiles
           ! W.J notes: CanopH is not used. If CLM_veg_typs_fracs exists, the read some dummy ???? Ask Sarith  
           if (NewLand) then
              read(21,*) I, j, ITY(N),idum, rdum, rdum, CanopH
           else
              read(21,*) I, j, ITY(N),idum, rdum, rdum
           endif

           read (22, *) i,j, this%GNU(n), this%BF1(n), this%BF2(n), this%BF3(n)

           read (23, *) i,j, idum, idum, this%BEE(n), this%PSIS(n),&
                this%POROS(n), this%COND(n), this%WPWET(n), DP2BR(n)

           read (24, *) i,j, rdum, this%ARS1(n), this%ARS2(n), this%ARS3(n), &
                this%ARA1(n), this%ARA2(n), this%ARA3(n), this%ARA4(n),      &
                this%ARW1(n), this%ARW2(n), this%ARW3(n), this%ARW4(n)

           read (25, *) i,j, rdum, this%TSA1(n), this%TSA2(n), this%TSB1(n), this%TSB2(n)

           if( SURFLAY.eq.20.0 ) read (26, *) i,j, this%ATAU(n), this%BTAU(n), rdum, rdum   ! for old soil params
           if( SURFLAY.eq.50.0 ) read (26, *) i,j, rdum , rdum, this%ATAU(n), this%BTAU(n)  ! for new soil params

        end do
        this%ity = real(ity)
        CLOSE (21, STATUS = 'KEEP')
        CLOSE (22, STATUS = 'KEEP')
        CLOSE (23, STATUS = 'KEEP')
        CLOSE (24, STATUS = 'KEEP')
        CLOSE (25, STATUS = 'KEEP')
        CLOSE (26, STATUS = 'KEEP')
      endif

      do n=1,ntiles
        zdep2=1000.
        zdep3=amax1(1000.,DP2BR(n))

        if (zdep2 .gt.0.75*zdep3) then
           zdep2  =  0.75*zdep3
        end if

        zdep1=20.
        zmet=zdep3/1000.

        term1=-1.+((this%PSIS(n)-zmet)/this%PSIS(n))**((this%BEE(n)-1.)/this%BEE(n))
        term2=this%PSIS(n)*this%BEE(n)/(this%BEE(n)-1)

        this%VGWMAX(n) = this%POROS(n)*zdep2
        this%CDCR1(n)  = 1000.*this%POROS(n)*(zmet-(-term2*term1))
        this%CDCR2(n)  = (1.-this%WPWET(n))*this%POROS(n)*zdep3
      enddo

      _RETURN(_SUCCESS)
   end subroutine add_bcs_to_rst

end module CatchmentRstMod
