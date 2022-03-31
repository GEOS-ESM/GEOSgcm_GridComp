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
     procedure :: write_bin      
     procedure :: write_nc4      
     procedure :: read_shared_nc4      
     procedure :: write_shared_nc4      
     procedure :: allocate
  endtype CatchmentRst

  interface CatchmentRst
     module procedure CatchmentRst_Create
  end interface

contains

  function CatchmentRst_create(filename, rc) result (catch)
    type(CatchmentRst) :: catch
    character(*), intent(in) :: filename
    integer, optional, intent(out) :: rc
    integer :: status
    character(len=256) :: Iam = "CatchmentRst_create"
    type(Netcdf4_fileformatter) :: formatter
    integer :: filetype, ntiles, unit
    type(FileMetadata) :: meta
    integer :: bpos, epos


     call MAPL_NCIOGetFileType(filename, filetype, __RC__)
     if (filetype == 0) then
        ! nc4 format
        call formatter%open(filename, pFIO_READ, __RC__)
        meta   = formatter%read(__RC__)
        ntiles = meta%get_dimension('tile', __RC__)
        call catch%allocate(ntiles)
        call catch%read_shared_nc4(formatter, __RC__)
        call MAPL_VarRead(formatter,"OLD_ITY",catch%ity, __RC__)
        call formatter%close()
  
     else ! binary
        open(newunit=unit, file=filename,  form='unformatted')
        bpos=0
        read(unit)
        epos = ftell(unit)            ! ending position of file pointer
        ntiles = (epos-bpos)/4-2    ! record size (in 4 byte words; 
        rewind(unit)
        call catch%allocate(ntiles)
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
     endif
     _RETURN(_SUCCESS)
   end function CatchmentRst_Create

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

   subroutine write_bin (this, filename, rc)
     class(CatchmentRst), intent(in):: this
     character(*), intent(in) :: filename
     integer, optional, intent(out) :: rc
     integer :: status, unit
     character(256) :: Iam = "write_bin"
     open(newunit=unit, file=filename,  form='unformatted') 
     write(unit) this%      bf1
     write(unit) this%      bf2
     write(unit) this%      bf3
     write(unit) this%   vgwmax
     write(unit) this%    cdcr1
     write(unit) this%    cdcr2
     write(unit) this%     psis
     write(unit) this%      bee
     write(unit) this%    poros
     write(unit) this%    wpwet
     write(unit) this%     cond
     write(unit) this%      gnu
     write(unit) this%     ars1
     write(unit) this%     ars2
     write(unit) this%     ars3
     write(unit) this%     ara1
     write(unit) this%     ara2
     write(unit) this%     ara3
     write(unit) this%     ara4
     write(unit) this%     arw1
     write(unit) this%     arw2
     write(unit) this%     arw3
     write(unit) this%     arw4
     write(unit) this%     tsa1
     write(unit) this%     tsa2
     write(unit) this%     tsb1
     write(unit) this%     tsb2
     write(unit) this%     atau
     write(unit) this%     btau
     write(unit) this%      ity
     write(unit) this%       tc
     write(unit) this%       qc
     write(unit) this%    capac
     write(unit) this%   catdef
     write(unit) this%    rzexc
     write(unit) this%   srfexc
     write(unit) this%  ghtcnt1
     write(unit) this%  ghtcnt2
     write(unit) this%  ghtcnt3
     write(unit) this%  ghtcnt4
     write(unit) this%  ghtcnt5
     write(unit) this%  ghtcnt6
     write(unit) this%    tsurf
     write(unit) this%   wesnn1
     write(unit) this%   wesnn2
     write(unit) this%   wesnn3
     write(unit) this%  htsnnn1
     write(unit) this%  htsnnn2
     write(unit) this%  htsnnn3
     write(unit) this%   sndzn1
     write(unit) this%   sndzn2
     write(unit) this%   sndzn3
     write(unit) this%       ch
     write(unit) this%       cm
     write(unit) this%       cq
     write(unit) this%       fr
     write(unit) this%       ww
     close(unit)
     _RETURN(_SUCCESS)
   end subroutine write_bin

   subroutine allocate(this, ntiles,rc)
     class(CatchmentRst), intent(inout) :: this
     integer, intent(in) :: ntiles
     integer, optional, intent(out):: rc
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

end module CatchmentRstMod
