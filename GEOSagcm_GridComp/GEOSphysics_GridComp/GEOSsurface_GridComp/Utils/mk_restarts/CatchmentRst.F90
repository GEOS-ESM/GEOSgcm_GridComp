#include "MAPL_Generic.h"

module CatchmentRstMod
  use mk_restarts_getidsMod, ONLY:      &
       GetIds,                          &                              
       ReadTileFile_RealLatLon
  use MAPL
  use MAPL_Base,         ONLY: MAPL_UNDEF
  use mpi
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
     type(FileMetadata) :: meta
     character(len=:), allocatable  :: time !yyyymmddhhmm
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
     real, allocatable ::    snowalb(:)
     real, allocatable ::         ch(:,:)
     real, allocatable ::         cm(:,:)
     real, allocatable ::         cq(:,:)
     real, allocatable ::         fr(:,:)
     real, allocatable ::         ww(:,:)
     ! save old values for scale
     real, allocatable :: old_catdef(:)
     real, allocatable :: old_cdcr1 (:)
     real, allocatable :: old_cdcr2 (:)
     real, allocatable :: old_rzexc (:)
     real, allocatable :: old_vgwmax(:)
     real, allocatable :: old_ghtcnt1(:)
     real, allocatable :: old_ghtcnt2(:)
     real, allocatable :: old_ghtcnt3(:)
     real, allocatable :: old_ghtcnt4(:)
     real, allocatable :: old_ghtcnt5(:)
     real, allocatable :: old_ghtcnt6(:)
     real, allocatable :: old_poros(:)
     real, allocatable :: old_sndzn3(:)
     ! intermediate
     real, allocatable, dimension(:) :: lonc,latc,lonn,latt, latg
     integer, allocatable, dimension(:) :: id_glb
  contains
     procedure :: read_GEOSldas_rst_bin      
     procedure :: write_nc4      
     procedure :: read_shared_nc4      
     procedure :: write_shared_nc4
     procedure :: add_bcs_to_rst      
     procedure :: allocate_catch
     procedure :: re_tile
     procedure :: re_scale
     procedure :: set_scale_var
  endtype CatchmentRst

  interface CatchmentRst
     module procedure CatchmentRst_Create
     module procedure CatchmentRst_empty
  end interface

contains

  function CatchmentRst_create(filename, time, rc) result (catch)
     type(CatchmentRst) :: catch
     character(*), intent(in)           :: filename
     character(*), intent(in) :: time ! yyyymmddhhmm format
     integer, optional, intent(out) :: rc
     integer :: status
     character(len=256) :: Iam = "CatchmentRst_create"
     type(Netcdf4_fileformatter) :: formatter
     integer :: filetype, ntiles, unit
     type(FileMetadata) :: meta
     integer :: bpos, epos, n, myid, mpierr

     call MPI_COMM_RANK( MPI_COMM_WORLD, myid, mpierr )
     catch%time = time
     call MAPL_NCIOGetFileType(filename, filetype, __RC__)

     if(filetype == 0) then
       ! nc4 format
       call formatter%open(filename, pFIO_READ, __RC__)
       meta   = formatter%read(__RC__)
       ntiles = meta%get_dimension('tile', __RC__)
       catch%meta   = meta
       catch%ntiles = ntiles
       if (myid ==0) then
          call catch%allocate_catch()
          call catch%read_shared_nc4(formatter, __RC__)
       endif
       call formatter%close()
     else
       !if ( .not. present(time)) then
       !  _ASSERT(.false., 'Please provide time for binary catch, format yyyymmddhhmm')
       !endif
       !GEOSldas binary
       open(newunit=unit, file=filename,  form='unformatted', action = 'read')
       bpos=0
       read(unit)
       epos = ftell(unit)          ! ending position of file pointer
       close(unit)
       ntiles = (epos-bpos)/4-2    ! record size (in 4 byte words; 
       catch%ntiles = ntiles
       catch%meta = create_meta(ntiles, time)
       if (myid ==0) then
          call catch%allocate_catch()
          call catch%read_GEOSldas_rst_bin(filename, __RC__)
       endif
     endif
     _RETURN(_SUCCESS)
   end function CatchmentRst_Create

  function CatchmentRst_empty(meta, time, rc) result (catch)
     type(CatchmentRst) :: catch
     character(*), intent(in) :: time
     type(FileMetadata), intent(in) :: meta
     
     integer, optional, intent(out) :: rc
     integer :: status, myid, mpierr
     character(len=256) :: Iam = "CatchmentRst_create"
     type(Netcdf4_fileformatter) :: formatter

     call MPI_COMM_RANK( MPI_COMM_WORLD, myid, mpierr )
       ! nc4 format
     catch%ntiles = meta%get_dimension('tile', __RC__)
     catch%meta = meta 
     catch%time = time
     if (myid ==0) then
        call catch%allocate_catch()
     endif
     _RETURN(_SUCCESS)
   end function CatchmentRst_empty

   subroutine read_GEOSldas_rst_bin(this, filename, rc)
      class(CatchmentRst), intent(inout) :: this
      character(*), intent(in) :: filename
      integer, optional, intent(out) :: rc
      integer :: unit
      open(newunit=unit, file=filename, form='unformatted', action='read')
      read(unit)                  ! skip     bf1
      read(unit)                  ! skip     bf2
      read(unit)                  ! skip     bf3
      read(unit) this%   vgwmax
      read(unit) this%    cdcr1
      read(unit) this%    cdcr2
      read(unit)                  ! skip    psis
      read(unit)                  ! skip     bee
      read(unit) this%    poros
      read(unit)                  ! skip   wpwet
      read(unit)                  ! skip    cond
      read(unit)                  ! skip     gnu
      read(unit)                  ! skip    ars1
      read(unit)                  ! skip    ars2
      read(unit)                  ! skip    ars3
      read(unit)                  ! skip    ara1
      read(unit)                  ! skip    ara2
      read(unit)                  ! skip    ara3
      read(unit)                  ! skip    ara4
      read(unit)                  ! skip    arw1
      read(unit)                  ! skip    arw2
      read(unit)                  ! skip    arw3
      read(unit)                  ! skip    arw4
      read(unit)                  ! skip    tsa1
      read(unit)                  ! skip    tsa2
      read(unit)                  ! skip    tsb1
      read(unit)                  ! skip    tsb2
      read(unit)                  ! skip    atau
      read(unit)                  ! skip    btau
      read(unit)                  ! skip     ity
      read(unit) this%       tc
      read(unit) this%       qc
      read(unit) this%    capac
      read(unit) this%   catdef
      read(unit) this%    rzexc
      read(unit) this%   srfexc
      read(unit) this%  ghtcnt1
      read(unit) this%  ghtcnt2
      read(unit) this%  ghtcnt3
      read(unit) this%  ghtcnt4
      read(unit) this%  ghtcnt5
      read(unit) this%  ghtcnt6
      read(unit) this%    tsurf
      read(unit) this%   wesnn1
      read(unit) this%   wesnn2
      read(unit) this%   wesnn3
      read(unit) this%  htsnnn1
      read(unit) this%  htsnnn2
      read(unit) this%  htsnnn3
      read(unit) this%   sndzn1
      read(unit) this%   sndzn2
      read(unit) this%   sndzn3
      read(unit) this%       ch
      read(unit) this%       cm
      read(unit) this%       cq
      read(unit) this%       fr
      read(unit) this%       ww
      close(unit)
      _RETURN(_SUCCESS)
   end subroutine read_GEOSldas_rst_bin

   subroutine read_shared_nc4(this, formatter, rc)
     class(CatchmentRst), intent(inout):: this
     type(Netcdf4_fileformatter),intent(inout) :: formatter
     integer, optional, intent(out):: rc
     integer :: status
     
     ! these four (time-invariant) variables are used for rescaling of prognostic variables
     call MAPL_VarRead(formatter,"VGWMAX",this%vgwmax, __RC__)
     call MAPL_VarRead(formatter,"CDCR1",this%cdcr1, __RC__)
     call MAPL_VarRead(formatter,"CDCR2",this%cdcr2, __RC__)
     call MAPL_VarRead(formatter,"POROS",this%poros, __RC__)

     ! Catchment model prognostic variables (and some diagnostics needed in Catch restart for GCM) 
     call MAPL_VarRead(formatter,"TC",this%tc, __RC__)
     call MAPL_VarRead(formatter,"QC",this%qc, __RC__)
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
     if (this%meta%has_variable('TSURF')) then
       call MAPL_VarRead(formatter,"TSURF",this%tsurf, __RC__)
     endif
     call MAPL_VarRead(formatter,"WESNN1",this%wesnn1, __RC__)
     call MAPL_VarRead(formatter,"WESNN2",this%wesnn2, __RC__)
     call MAPL_VarRead(formatter,"WESNN3",this%wesnn3, __RC__)
     call MAPL_VarRead(formatter,"HTSNNN1",this%htsnnn1, __RC__)
     call MAPL_VarRead(formatter,"HTSNNN2",this%htsnnn2, __RC__)
     call MAPL_VarRead(formatter,"HTSNNN3",this%htsnnn3, __RC__)
     call MAPL_VarRead(formatter,"SNDZN1",this%sndzn1, __RC__)
     call MAPL_VarRead(formatter,"SNDZN2",this%sndzn2, __RC__)
     call MAPL_VarRead(formatter,"SNDZN3",this%sndzn3, __RC__)
     if (this%meta%has_variable('CH')) then
       call MAPL_VarRead(formatter,"CH",this%ch, __RC__)
     endif
     if (this%meta%has_variable('CM')) then
       call MAPL_VarRead(formatter,"CM",this%cm, __RC__)
     endif
     if (this%meta%has_variable('CQ')) then
       call MAPL_VarRead(formatter,"CQ",this%cq, __RC__)
     endif
     if (this%meta%has_variable('FR')) then
       call MAPL_VarRead(formatter,"FR",this%fr, __RC__)
     endif
     if (this%meta%has_variable('WW')) then
       call MAPL_VarRead(formatter,"WW",this%ww, __RC__)
     endif
     _RETURN(_SUCCESS)
   end subroutine read_shared_nc4

   subroutine write_nc4 (this, filename, rc)
     class(CatchmentRst), intent(inout):: this
     character(*), intent(in) :: filename
     integer, optional, intent(out):: rc

     type(Netcdf4_fileformatter) :: formatter
     integer :: status
     character(256) :: Iam = "write_nc4"

     call formatter%create(filename, __RC__)
     call formatter%write(this%meta, __RC__)
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
     if (this%meta%has_variable('TSURF')) then
        call MAPL_VarWrite(formatter,"TSURF",this%tsurf)
     endif
     call MAPL_VarWrite(formatter,"WESNN1",this%wesnn1)
     call MAPL_VarWrite(formatter,"WESNN2",this%wesnn2)
     call MAPL_VarWrite(formatter,"WESNN3",this%wesnn3)
     call MAPL_VarWrite(formatter,"HTSNNN1",this%htsnnn1)
     call MAPL_VarWrite(formatter,"HTSNNN2",this%htsnnn2)
     call MAPL_VarWrite(formatter,"HTSNNN3",this%htsnnn3)
     call MAPL_VarWrite(formatter,"SNDZN1",this%sndzn1)
     call MAPL_VarWrite(formatter,"SNDZN2",this%sndzn2)
     call MAPL_VarWrite(formatter,"SNDZN3",this%sndzn3)

     if (this%meta%has_variable('CH')) then
       call MAPL_VarWrite(formatter,"CH",this%ch)
     endif
     if (this%meta%has_variable('CM')) then
       call MAPL_VarWrite(formatter,"CM",this%cm)
     endif
     if (this%meta%has_variable('CQ')) then
       call MAPL_VarWrite(formatter,"CQ",this%cq)
     endif
     if (this%meta%has_variable('FR')) then
       call MAPL_VarWrite(formatter,"FR",this%fr)
     endif
     if (this%meta%has_variable('WW')) then
       call MAPL_VarWrite(formatter,"WW",this%ww)
     endif
     if (this%meta%has_variable('SNOWALB')) then
       call MAPL_VarWrite(formatter,"SNOWALB",this%snowalb)
     endif

     _RETURN(_SUCCESS)

   end subroutine write_shared_nc4

   subroutine allocate_catch(this,rc)
     class(CatchmentRst), intent(inout) :: this
     integer, optional, intent(out):: rc
     integer :: ntiles
     ntiles = this%ntiles
     allocate( this%     vgwmax(ntiles) )
     allocate( this%      cdcr1(ntiles) )
     allocate( this%      cdcr2(ntiles) )
     allocate( this%      poros(ntiles) )

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
     
     if (this%meta%has_variable('TSURF')) then
       allocate( this%    tsurf(ntiles) )
     endif

     allocate( this%     wesnn1(ntiles) )
     allocate( this%     wesnn2(ntiles) )
     allocate( this%     wesnn3(ntiles) )
     allocate( this%    htsnnn1(ntiles) )
     allocate( this%    htsnnn2(ntiles) )
     allocate( this%    htsnnn3(ntiles) )
     allocate( this%     sndzn1(ntiles) )
     allocate( this%     sndzn2(ntiles) )
     allocate( this%     sndzn3(ntiles) )

     if (this%meta%has_variable('CH')) then
       allocate( this%         ch(ntiles,4) )
     endif
     if (this%meta%has_variable('CM')) then
       allocate( this%         cm(ntiles,4) )
     endif
     if (this%meta%has_variable('CQ')) then
       allocate( this%         cq(ntiles,4) )
     endif
     if (this%meta%has_variable('FR')) then
       allocate( this%         fr(ntiles,4) )
     endif
     if (this%meta%has_variable('WW')) then
       allocate( this%         ww(ntiles,4) )
     endif

     _RETURN(_SUCCESS)
   end subroutine allocate_catch

   ! This subroutine reads BCs from BCSDIR and hydrological variable (??)
   subroutine add_bcs_to_rst(this, surflay, DataDir, rc)
      class(CatchmentRst), intent(inout) :: this
      real, intent(in) :: surflay
      character(*), intent(in) :: DataDir
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
      type(Variable)              :: var
      type(FileMetadata)          :: meta_

      character*256               :: Iam = "add_bcs_to_rst"

      ntiles = this%ntiles

      allocate (DP2BR(ntiles), ity(ntiles) )

      inquire(file = trim(DataDir)//'/clsm/catch_params.nc4', exist=file_exists)
      inquire(file = trim(DataDir)//"/clsm/CLM_veg_typs_fracs",exist=NewLand )

      this%ity   = DP2BR
      this%ARA1  = DP2BR      
      this%ARA2  = DP2BR      
      this%ARA3  = DP2BR      
      this%ARA4  = DP2BR      
      this%ARS1  = DP2BR      
      this%ARS2  = DP2BR      
      this%ARS3  = DP2BR      
      this%ARW1  = DP2BR      
      this%ARW2  = DP2BR      
      this%ARW3  = DP2BR      
      this%ARW4  = DP2BR      
      this%ATAU  = DP2BR      
      this%BTAU  = DP2BR      
      this%PSIS  = DP2BR      
      this%BEE   = DP2BR      
      this%BF1   = DP2BR      
      this%BF2   = DP2BR      
      this%BF3   = DP2BR      
      this%TSA1  = DP2BR      
      this%TSA2  = DP2BR      
      this%TSB1  = DP2BR      
      this%TSB2  = DP2BR      
      this%GNU   = DP2BR      
      this%COND  = DP2BR      
      this%WPWET = DP2BR      
      this%POROS = DP2BR
      this%VGWMAX = DP2BR
      this%cdcr1 = DP2BR
      this%cdcr2 = DP2BR

      if(file_exists) then
        call CatchFmt%Open(trim(DataDir)//'/clsm/catch_params.nc4', pFIO_READ, __RC__)
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

        meta_ = CatchFmt%read(__RC__)

        if (meta_%has_variable('SNOWALB')) then
           if ( .not. allocated(this%snowalb)) allocate(this%snowalb(ntiles))
           call MAPL_VarRead ( CatchFmt ,'SNOWALB', this%snowalb, __RC__)
           if ( .not. this%meta%has_variable('SNOWALB')) then
              var = Variable(type=pFIO_REAL32, dimensions='tile')
              call var%add_attribute('long_name', 'snow_reflectivity')
              call var%add_attribute('units', '1')
              call this%meta%add_variable('SNOWALB', var)
           endif
        elseif (this%meta%has_variable('SNOWALB')) then
           call this%meta%remove_variable('SNOWALB')
        endif
        call CatchFmt%close()
      else
        open(unit=21, file=trim(DataDir)//'/clsm/mosaic_veg_typs_fracs',form='formatted')
        open(unit=22, file=trim(DataDir)//'/clsm//bf.dat'               ,form='formatted')
        open(unit=23, file=trim(DataDir)//'/clsm/soil_param.dat'       ,form='formatted')
        open(unit=24, file=trim(DataDir)//'/clsm/ar.new'               ,form='formatted')
        open(unit=25, file=trim(DataDir)//'/clsm/ts.dat'               ,form='formatted')
        open(unit=26, file=trim(DataDir)//'/clsm/tau_param.dat'        ,form='formatted')

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

        if (this%meta%has_variable('SNOWALB'))  call this%meta%remove_variable('SNOWALB')

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

   type(FileMetadata) function create_meta(ntiles, t, rc) result(meta)
     integer, intent(in) :: ntiles
     character(*), intent(in) :: t ! yyyymmddmmhh format
     integer, optional, intent(out) :: rc
     character(64), dimension(58, 3) :: fields
     integer :: n, status
     character(:), allocatable :: s
     type(Variable) :: var
 
     fields(1,:)  = [character(len=64)::"ARA1"     ,  "shape_param_1"                             ,     "m+2 kg-1"]
     fields(2,:)  = [character(len=64)::"ARA2"     ,  "shape_param_2"                             ,     "1"]
     fields(3,:)  = [character(len=64)::"ARA3"     ,  "shape_param_3"                             ,     "m+2 kg-1"]
     fields(4,:)  = [character(len=64)::"ARA4"     ,  "shape_param_4"                             ,     "1"]
     fields(5,:)  = [character(len=64)::"ARS1"     ,  "wetness_param_1"                           ,     "m+2 kg-1"]
     fields(6,:)  = [character(len=64)::"ARS2"     ,  "wetness_param_2"                           ,     "m+2 kg-1"]
     fields(7,:)  = [character(len=64)::"ARS3"     ,  "wetness_param_3"                           ,     "m+4 kg-2"]
     fields(8,:)  = [character(len=64)::"ARW1"     ,  "min_theta_param_1"                         ,     "m+2 kg-1"]
     fields(9,:)  = [character(len=64)::"ARW2"     ,  "min_theta_param_2"                         ,     "m+2 kg-1"]
     fields(10,:) = [character(len=64)::"ARW3"     ,  "min_theta_param_3"                         ,     "m+4 kg-2"]
     fields(11,:) = [character(len=64)::"ARW4"     ,  "min_theta_param_4"                         ,     "1"]
     fields(12,:) = [character(len=64)::"ATAU"     ,  "water_transfer_param_5"                    ,     "1"]
     fields(13,:) = [character(len=64)::"BEE"      ,  "clapp_hornberger_b"                        ,     "1"]
     fields(14,:) = [character(len=64)::"BF1"      ,  "topo_baseflow_param_1"                     ,     "kg m-4"]
     fields(15,:) = [character(len=64)::"BF2"      ,  "topo_baseflow_param_2"                     ,     "m"]
     fields(16,:) = [character(len=64)::"BF3"      ,  "topo_baseflow_param_3"                     ,     "log(m)"]
     fields(17,:) = [character(len=64)::"BTAU"     ,  "water_transfer_param_6"                    ,     "1"]
     fields(18,:) = [character(len=64)::"CAPAC"    ,  "interception_reservoir_capac"              ,     "kg m-2"]
     fields(19,:) = [character(len=64)::"CATDEF"   ,  "catchment_deficit"                         ,     "kg m-2"]
     fields(20,:) = [character(len=64)::"CDCR1"    ,  "moisture_threshold"                        ,     "kg m-2"]
     fields(21,:) = [character(len=64)::"CDCR2"    ,  "max_soil_water_content_above_wilting_point",     "kg m-2"]
     fields(22,:) = [character(len=64)::"COND"     ,  "sfc_sat_hydraulic_conduct"                 ,     "m s-1"]
     fields(23,:) = [character(len=64)::"GHTCNT1"  ,  "soil_heat_content_layer_1"                 ,     "J m-2"]
     fields(24,:) = [character(len=64)::"GHTCNT2"  ,  "soil_heat_content_layer_2"                 ,     "J_m-2"]
     fields(25,:) = [character(len=64)::"GHTCNT3"  ,  "soil_heat_content_layer_3"                 ,     "J m-2"]
     fields(26,:) = [character(len=64)::"GHTCNT4"  ,  "soil_heat_content_layer_4"                 ,     "J m-2"]
     fields(27,:) = [character(len=64)::"GHTCNT5"  ,  "soil_heat_content_layer_5"                 ,     "J m-2"]
     fields(28,:) = [character(len=64)::"GHTCNT6"  ,  "soil_heat_content_layer_6"                 ,     "J m-2"]
     fields(29,:) = [character(len=64)::"GNU"      ,  "vertical_transmissivity"                   ,     "m-1"]
     fields(30,:) = [character(len=64)::"HTSNNN1"  ,  "heat_content_snow_layer_1"                 ,     "J m-2"]
     fields(31,:) = [character(len=64)::"HTSNNN2"  ,  "heat_content_snow_layer_2"                 ,     "J m-2"]
     fields(32,:) = [character(len=64)::"HTSNNN3"  ,  "heat_content_snow_layer_3"                 ,     "J m-2"]
     fields(33,:) = [character(len=64)::"OLD_ITY"  ,  "Placeholder. Used to be vegetation_type."  ,     "1"]
     fields(34,:) = [character(len=64)::"POROS"    ,  "soil_porosity"                             ,     "1"]
     fields(35,:) = [character(len=64)::"PSIS"     ,  "saturated_matric_potential"                ,     "m"]
     fields(36,:) = [character(len=64)::"RZEXC"    ,  "root_zone_excess"                          ,     "kg m-2"]
     fields(37,:) = [character(len=64)::"SNDZN1"   ,  "snow_depth_layer_1"                        ,     "m"]
     fields(38,:) = [character(len=64)::"SNDZN2"   ,  "snow_depth_layer_2"                        ,     "m"]
     fields(39,:) = [character(len=64)::"SNDZN3"   ,  "snow_depth_layer_3"                        ,     "m"]
     fields(40,:) = [character(len=64)::"SRFEXC"   ,  "surface_excess"                            ,     "kg m-2"]
     fields(41,:) = [character(len=64)::"TILE_ID"  ,  "catchment_tile_id"                         ,     "1"]
     fields(42,:) = [character(len=64)::"TSA1"     ,  "water_transfer_param_1"                    ,     "1"]
     fields(43,:) = [character(len=64)::"TSA2"     ,  "water_transfer_param_2"                    ,     "1"]
     fields(44,:) = [character(len=64)::"TSB1"     ,  "water_transfer_param_3"                    ,     "1"]
     fields(45,:) = [character(len=64)::"TSB2"     ,  "water_transfer_param_4"                    ,     "1"]
     fields(46,:) = [character(len=64)::"TSURF"    ,  "mean_catchment_temp_incl_snw"              ,     "K"]
     fields(47,:) = [character(len=64)::"VGWMAX"   ,  "max_rootzone_water_content"                ,     "kg m-2"]
     fields(48,:) = [character(len=64)::"WESNN1"   ,  "snow_mass_layer_1"                         ,     "kg m-2"]
     fields(49,:) = [character(len=64)::"WESNN2"   ,  "snow_mass_layer_2"                         ,     "kg m-2"]
     fields(50,:) = [character(len=64)::"WESNN3"   ,  "snow_mass_layer_3"                         ,     "kg m-2"]
     fields(51,:) = [character(len=64)::"WPWET"    ,  "wetness_at_wilting_point"                  ,     "1"]
     fields(52,:) = [character(len=64)::"CH"       ,  "surface_heat_exchange_coefficient"         ,     "kg m-2 s-1"]
     fields(53,:) = [character(len=64)::"CM"       ,  "surface_momentum_exchange_coefficient"     ,     "kg m-2 s-1"]
     fields(54,:) = [character(len=64)::"CQ"       ,  "surface_moisture_exchange_coffiecient"     ,     "kg m-2 s-1"]
     fields(55,:) = [character(len=64)::"FR"       ,  "subtile_fractions"                         ,     "1"]
     fields(56,:) = [character(len=64)::"QC"       ,  "canopy_specific_humidity"                  ,     "kg kg-1"]
     fields(57,:) = [character(len=64)::"TC"       ,  "canopy_temperature"                        ,     "K"]
     fields(58,:) = [character(len=64)::"WW"       ,  "vertical_velocity_scale_squared"           ,     "m+2 s-2"]

     call meta%add_dimension('tile', ntiles)
     call meta%add_dimension('subtile', 4)
     call meta%add_dimension('time',1)

     do n = 1, 58
        if (n >=52) then
           var = Variable(type=pFIO_REAL32, dimensions='tile,subtile')
        else
           var = Variable(type=pFIO_REAL32, dimensions='tile')
        endif
        call var%add_attribute('long_name', trim(fields(n,2)))
        call var%add_attribute('units', trim(fields(n,3)))
        call meta%add_variable(trim(fields(n,1)), var)
     enddo
     var = Variable(type=pFIO_REAL32, dimensions='time')
     s   = "minutes since "//t(1:4)//"-"//t(5:6)//"-"//t(7:8)//" "//t(9:10)//":00:00"
     call var%add_attribute('units', s)
     call meta%add_variable('time', var)
     _RETURN(_SUCCESS)
   end function create_meta

   subroutine re_tile(this, InTileFile, OutBcsDir, OutTileFile, surflay, rc)
     class(CatchmentRst), intent(inout) :: this
     character(*), intent(in) :: InTileFile
     character(*), intent(in) :: OutBcsDir
     character(*), intent(in) :: OutTileFile
     real, intent(in) :: surflay
     integer, optional, intent(out) :: rc
     integer :: status, in_ntiles, out_ntiles, myid, numprocs
     real, allocatable :: var_out(:), tmp2d(:,:)
     real   , allocatable , dimension (:) :: lonn,latt
     real   , pointer, dimension      (:) :: long, latg, lonc, latc
     integer, allocatable , dimension (:) :: low_ind, upp_ind, nt_local
     integer, allocatable , dimension (:) :: Id_glb, id_loc, tid_offl
     logical :: root_proc
     integer :: mpierr, n, i, k, req
     type(CatchmentRst) :: xgrid
     type(FileMetadata) :: meta
     character(*), parameter :: Iam = "Catchment::Re_tile"

     open (10,file =trim(OutBcsDir)//"/clsm/catchment.def",status='old',form='formatted')
     read (10,*) out_ntiles
     close (10, status = 'keep')

     in_ntiles = this%ntiles

     call this%meta%modify_dimension('tile', out_ntiles, __RC__)

     call MPI_COMM_RANK( MPI_COMM_WORLD, myid, mpierr )
     call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, mpierr )
     root_proc = .false.
     if (myid == 0)  root_proc = .true.

     if(root_proc) then
       print *,'ntiles in target BCs      : ',out_ntiles
       print *,'ntiles in restarts        : ',in_ntiles
     endif

     
     ! Domain decomposition
     ! --------------------

     
     allocate(low_ind (   numprocs))
     allocate(upp_ind (   numprocs))
     allocate(nt_local(   numprocs))
     low_ind (:)    = 1
     upp_ind (:)    = out_ntiles
     nt_local(:)    = out_ntiles

     if (numprocs > 1) then
       do i = 1, numprocs - 1
          upp_ind(i)   = low_ind(i) + (out_ntiles/numprocs) - 1
          low_ind(i+1) = upp_ind(i) + 1
          nt_local(i)  = upp_ind(i) - low_ind(i) + 1
       end do
       nt_local(numprocs) = upp_ind(numprocs) - low_ind(numprocs) + 1
     endif

     allocate (id_loc (nt_local (myid + 1)))
     allocate (lonn   (nt_local (myid + 1)))
     allocate (latt   (nt_local (myid + 1)))
     allocate (lonc   (1:in_ntiles))
     allocate (latc   (1:in_ntiles))
     allocate (tid_offl (in_ntiles))

     if (root_proc) then
        allocate (long   (out_ntiles))
        allocate (latg   (out_ntiles))
        call ReadTileFile_RealLatLon ( OutTileFile, n, xlon=long, xlat=latg)
        _ASSERT( n == out_ntiles, "Out tile number should match")
        this%latg = latg
        call ReadTileFile_RealLatLon ( InTileFile, n, xlon=lonc, xlat=latc)
        _ASSERT( n == in_ntiles, "In tile number should match")
     endif

     do i = 1, in_ntiles
        tid_offl(i)    = i
     end do

     ! create mapping, nearest
     call MPI_Barrier(MPI_COMM_WORLD, STATUS)

     do i = 1, numprocs
       if((I == 1).and.(myid == 0)) then
          lonn(:) = long(low_ind(i) : upp_ind(i))
          latt(:) = latg(low_ind(i) : upp_ind(i))
       else if (I > 1) then
          if(I-1 == myid) then
             ! receiving from root
             call MPI_RECV(lonn,nt_local(i) , MPI_REAL, 0,995,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
             call MPI_RECV(latt,nt_local(i) , MPI_REAL, 0,994,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
          else if (myid == 0) then
             ! root sends
             call MPI_ISend(long(low_ind(i) : upp_ind(i)),nt_local(i),MPI_REAL,i-1,995,MPI_COMM_WORLD,req,mpierr)
             call MPI_WAIT (req,MPI_STATUS_IGNORE,mpierr)
             call MPI_ISend(latg(low_ind(i) : upp_ind(i)),nt_local(i),MPI_REAL,i-1,994,MPI_COMM_WORLD,req,mpierr)
             call MPI_WAIT (req,MPI_STATUS_IGNORE,mpierr)
          endif
       endif
     enddo
     if(root_proc) deallocate (long)

     call MPI_BCAST(lonc,in_ntiles,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
     call MPI_BCAST(latc,in_ntiles,MPI_REAL,0,MPI_COMM_WORLD,mpierr)

    ! --------------------------------------------------------------------------------
    ! Here we create transfer index array to map offline restarts to output tile space
    ! --------------------------------------------------------------------------------   

    ! id_glb for hydrologic variable
     this%lonc = lonc
     this%latc = latc
     this%lonn = lonn
     this%latt = latt

     call GetIds(lonc, latc, lonn, latt, id_loc, tid_offl)
     if(root_proc)  allocate (id_glb  (out_ntiles))

     call MPI_Barrier(MPI_COMM_WORLD, STATUS)

     do i = 1, numprocs
        if((I == 1).and.(myid == 0)) then
           id_glb(low_ind(i) : upp_ind(i)) = Id_loc(:)
        else if (I > 1) then
           if(I-1 == myid) then
              ! send to root
              call MPI_ISend(id_loc,nt_local(i),MPI_INTEGER,0,993,MPI_COMM_WORLD,req,mpierr)
              call MPI_WAIT (req,MPI_STATUS_IGNORE,mpierr)
           else if (myid == 0) then
              ! root receives
              call MPI_RECV(id_glb(low_ind(i) : upp_ind(i)),nt_local(i) , MPI_INTEGER, i-1,993,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
           endif
        endif
     end do

     deallocate (id_loc)

     this%ntiles = out_ntiles
     if (root_proc)  then
       ! regrid
       this%id_glb = id_glb
       var_out    = this%poros(id_glb(:))
       this%poros = var_out

       var_out      = this%vgwmax(id_glb(:))
       this%vgwmax  = var_out

       var_out     = this%cdcr1(id_glb(:))
       this%cdcr1  = var_out

       var_out     = this%cdcr2(id_glb(:))
       this%cdcr2  = var_out

       tmp2d = this%tc
       deallocate(this%tc)
       allocate(this%tc(out_ntiles, 4))
       do k = 1, 4
          this%tc(:,k) = tmp2d(id_glb(:),k)
       enddo

       tmp2d = this%qc
       deallocate(this%qc)
       allocate(this%qc(out_ntiles, 4))
       do k = 1, 4
          this%qc(:,k) = tmp2d(id_glb(:),k)
       enddo

       var_out      = this%capac(id_glb(:))
       this%capac   = var_out

       var_out       = this%catdef(id_glb(:))
       this%catdef   = var_out

       var_out       = this%rzexc(id_glb(:))
       this%rzexc    = var_out

       var_out       = this%SRFEXC(id_glb(:))
       this%SRFEXC   = var_out

       var_out       = this%GHTCNT1(id_glb(:))
       this%GHTCNT1  = var_out

       var_out       = this%GHTCNT2(id_glb(:))
       this%GHTCNT2  = var_out

       var_out       = this%GHTCNT3(id_glb(:))
       this%GHTCNT3  = var_out

       var_out       = this%GHTCNT4(id_glb(:))
       this%GHTCNT4  = var_out

       var_out       = this%GHTCNT5(id_glb(:))
       this%GHTCNT5  = var_out

       var_out       = this%GHTCNT6(id_glb(:))
       this%GHTCNT6  = var_out

       var_out       = this%WESNN1(id_glb(:))
       this%WESNN1   = var_out

       var_out       = this%WESNN2(id_glb(:))
       this%WESNN2   = var_out

       var_out       = this%WESNN3(id_glb(:))
       this%WESNN3   = var_out

       var_out       = this%HTSNNN1(id_glb(:))
       this%HTSNNN1  = var_out

       var_out       = this%HTSNNN2(id_glb(:))
       this%HTSNNN2  = var_out

       var_out       = this%HTSNNN3(id_glb(:))
       this%HTSNNN3  = var_out

       var_out       = this%SNDZN1(id_glb(:))
       this%SNDZN1   = var_out

       var_out       = this%SNDZN2(id_glb(:))
       this%SNDZN2   = var_out

       var_out       = this%SNDZN3(id_glb(:))
       this%SNDZN3   = var_out
    
       !set tsurf to zero
       if (this%meta%has_variable('TSURF')) then
          var_out = this%tsurf(id_glb(:)) 
          this%tsurf = var_out
       endif

       ! CH CM CQ FR WW
       ! WW
       if(allocated(tmp2d)) deallocate(tmp2d)
       allocate(tmp2d(out_ntiles,4))
       !tmp2d = 0.001
       if (this%meta%has_variable('CH')) then
          do k = 1,4
            tmp2d(:,k) = this%ch(id_glb(:),k)
          enddo
          this%ch  = tmp2d
       endif
       if (this%meta%has_variable('CM')) then
          do k = 1,4
            tmp2d(:,k) = this%cm(id_glb(:),k)
          enddo
          this%cm  = tmp2d
       endif
       if (this%meta%has_variable('CQ')) then
          do k = 1,4
            tmp2d(:,k) = this%cq(id_glb(:),k)
          enddo
          this%cq  = tmp2d
       endif
       !tmp2d = 0.25
       if (this%meta%has_variable('FR')) then
          do k = 1,4
            tmp2d(:,k) = this%fr(id_glb(:),k)
          enddo
          this%fr  = tmp2d
       endif
       !tmp2d = 0.1
       if (this%meta%has_variable('WW')) then
          do k = 1,4
            tmp2d(:,k) = this%ww(id_glb(:),k)
          enddo
          this%ww  = tmp2d
       endif

       call this%set_scale_var()
     endif

     if(associated(long)) deallocate(long)
     if(associated(latg)) deallocate(latg)
     if(associated(lonc)) deallocate(lonc)
     if(associated(latc)) deallocate(latc)

     _RETURN(_SUCCESS) 
   end subroutine re_tile

   subroutine re_scale(this, surflay, wemin_in, wemin_out,  rc)
     class(CatchmentRst), intent(inout) :: this
     real, intent(in) :: surflay
     real, intent(in) :: wemin_in
     real, intent(in) :: wemin_out
     integer, optional, intent(out) :: rc
     integer :: ntiles

     real,    allocatable, dimension(:)   :: dzsf, ar1, ar2, ar4
     real,    allocatable, dimension(:,:) :: TP_IN, GHT_IN, FICE, TP_OUT
     real,    allocatable, dimension(:)   :: swe_in, depth_in, areasc_in, areasc_out, depth_out

     type(Netcdf4_fileformatter) :: formatter
     integer :: i, filetype, n
     integer :: status

     ntiles =  this%ntiles

     n =count((this%old_catdef .gt. this%old_cdcr1))

     print*, "Scale tile and pesentile :", n,100*n/ntiles

! Scale rxexc regardless of CDCR1, CDCR2 differences
! --------------------------------------------------
     this%rzexc  = this%old_rzexc * ( this%vgwmax / this%old_vgwmax )

! Scale catdef regardless of whether CDCR2 is larger or smaller in the new situation
! ----------------------------------------------------------------------------------
     where (this%old_catdef .gt. this%old_cdcr1)

        this%catdef = this%cdcr1 +              &
                   ( this%old_catdef - this%old_cdcr1 ) / &
                   ( this%old_cdcr2  - this%old_cdcr1 ) * &
                   ( this%cdcr2 - this%cdcr1)
     end where

! Scale catdef also for the case where catdef le cdcr1.
! -----------------------------------------------------
     where( this%old_catdef .le. this%old_cdcr1)
       this%catdef = this%old_catdef * (this%cdcr1 / this%old_cdcr1)
     end where

! Sanity Check (catch_calc_soil_moist() forces consistency betw. srfexc, rzexc, catdef)
! ------------
     print *, 'Performing Sanity Check ...'

     allocate (   dzsf(ntiles), source = SURFLAY )
     allocate (   ar1( ntiles) )
     allocate (   ar2( ntiles) )
     allocate (   ar4( ntiles) )

     call catch_calc_soil_moist( ntiles, dzsf,                  &
        this%vgwmax, this%cdcr1, this%cdcr2,                    &
        this%psis,   this%bee,   this%poros, this%wpwet,        &
        this%ars1,   this%ars2,  this%ars3,                     &
        this%ara1,   this%ara2,  this%ara3,  this%ara4,         &
        this%arw1,   this%arw2,  this%arw3,  this%arw4,         &
        this%bf1,    this%bf2,                                  &
        this%srfexc, this%rzexc, this%catdef,                   &
        ar1,               ar2,              ar4 )


     allocate (TP_IN  (N_GT, Ntiles))
     allocate (GHT_IN (N_GT, Ntiles))
     allocate (FICE   (N_GT, NTILES))
     allocate (TP_OUT (N_GT, Ntiles))

     GHT_IN (1,:) = this%old_ghtcnt1
     GHT_IN (2,:) = this%old_ghtcnt2
     GHT_IN (3,:) = this%old_ghtcnt3
     GHT_IN (4,:) = this%old_ghtcnt4
     GHT_IN (5,:) = this%old_ghtcnt5
     GHT_IN (6,:) = this%old_ghtcnt6

     call catch_calc_tp ( NTILES, this%old_poros, GHT_IN, tp_in, FICE)

     do n = 1, ntiles
        do i = 1, N_GT
           call catch_calc_ght(dzgt(i), this%poros(n), tp_in(i,n), fice(i,n),  GHT_IN(i,n))
        end do
     end do

     this%ghtcnt1 = GHT_IN (1,:)
     this%ghtcnt2 = GHT_IN (2,:)
     this%ghtcnt3 = GHT_IN (3,:)
     this%ghtcnt4 = GHT_IN (4,:)
     this%ghtcnt5 = GHT_IN (5,:)
     this%ghtcnt6 = GHT_IN (6,:)

! Deep soil temp sanity check
! ---------------------------

     call catch_calc_tp ( NTILES, this%poros, GHT_IN, tp_out, FICE)

     print *, 'Percent tiles TP Layer 1 differ : ', 100.* count(ABS(tp_out(1,:) - tp_in(1,:)) > 1.e-5) /float (Ntiles)
     print *, 'Percent tiles TP Layer 2 differ : ', 100.* count(ABS(tp_out(2,:) - tp_in(2,:)) > 1.e-5) /float (Ntiles)
     print *, 'Percent tiles TP Layer 3 differ : ', 100.* count(ABS(tp_out(3,:) - tp_in(3,:)) > 1.e-5) /float (Ntiles)
     print *, 'Percent tiles TP Layer 4 differ : ', 100.* count(ABS(tp_out(4,:) - tp_in(4,:)) > 1.e-5) /float (Ntiles)
     print *, 'Percent tiles TP Layer 5 differ : ', 100.* count(ABS(tp_out(5,:) - tp_in(5,:)) > 1.e-5) /float (Ntiles)
     print *, 'Percent tiles TP Layer 6 differ : ', 100.* count(ABS(tp_out(6,:) - tp_in(6,:)) > 1.e-5) /float (Ntiles)

! SNOW scaling
! ------------

     if(abs(wemin_out-wemin_in) >1.e-6 ) then

        allocate (swe_in     (Ntiles))
        allocate (depth_in   (Ntiles))
        allocate (depth_out  (Ntiles))
        allocate (areasc_in  (Ntiles))
        allocate (areasc_out (Ntiles))
        swe_in    = this%wesnn1 + this%wesnn2 + this%wesnn3
        depth_in  = this%sndzn1 + this%sndzn2 + this%sndzn3
        areasc_in = min(swe_in/wemin_in, 1.)
        areasc_out= min(swe_in/wemin_out,1.)

        where (swe_in .gt. 0.)
           where (areasc_in .lt. 1. .or. areasc_out .lt. 1.)
           !      density_in= swe_in/(areasc_in *  depth_in + 1.e-20)
           !      depth_out = swe_in/(areasc_out*density_in)
             depth_out = areasc_in *  depth_in/(areasc_out + 1.e-20)
             this%sndzn1 = depth_out/3.
             this%sndzn2 = depth_out/3.
             this%sndzn3 = depth_out/3.
          endwhere
        endwhere

        print *, 'Snow scaling summary'
        print *, '....................'
        print *, 'Percent tiles SNDZ scaled : ', 100.* count (this%sndzn3 .ne. this%old_sndzn3) /float (count (this%sndzn3 > 0.))

     endif

  ! PEATCLSM - ensure low CATDEF on peat tiles where "old" restart is not also peat
  ! -------------------------------------------------------------------------------

     where ( (this%old_poros < PEATCLSM_POROS_THRESHOLD) .and. (this%poros >= PEATCLSM_POROS_THRESHOLD) )
        this%catdef = 25.
        this%rzexc  =  0.
        this%srfexc =  0.
     end where

   end subroutine re_scale

   subroutine set_scale_var(this )
     class(CatchmentRst), intent(inout) :: this
     this%old_poros  = this%poros
     this%old_cdcr1  = this%cdcr1
     this%old_cdcr2  = this%cdcr2
     this%old_vgwmax = this%vgwmax

     this%old_catdef = this%catdef
     this%old_rzexc  = this%rzexc
     this%old_sndzn3 = this%sndzn3
     this%old_ghtcnt1 = this%ghtcnt1
     this%old_ghtcnt2 = this%ghtcnt2
     this%old_ghtcnt3 = this%ghtcnt3
     this%old_ghtcnt4 = this%ghtcnt4
     this%old_ghtcnt5 = this%ghtcnt5
     this%old_ghtcnt6 = this%ghtcnt6
   end subroutine set_scale_var 

end module CatchmentRstMod
