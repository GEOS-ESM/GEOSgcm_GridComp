#include "MAPL_Generic.h"

module CatchmentCNRstMod
  use MAPL
  use LSM_ROUTINES,      ONLY:          &
       catch_calc_soil_moist,           &
       catch_calc_tp,                   &
       catch_calc_ght
  USE CATCH_CONSTANTS,   ONLY:          &
       N_GT              => CATCH_N_GT, &
       DZGT              => CATCH_DZGT, &
       PEATCLSM_POROS_THRESHOLD
  use CatchmentRstMod, only : CatchmentRst
  implicit none

  integer, parameter :: nveg  =  integer, parameter :: nveg    = 4
  integer, parameter :: nzone = 3
  integer, parameter :: VAR_COL_CLM40 = 40 ! number of CN column restart variables
  integer, parameter :: VAR_PFT_CLM40 = 74 ! number of CN PFT variables per column
  integer, parameter :: npft    = 19
  integer, parameter :: npft_clm45    = 19
  integer, parameter :: VAR_COL_CLM45 = 35 ! number of CN column restart variables
  integer, parameter :: VAR_PFT_CLM45 = 75 ! number of CN PFT variables per column 

 type, extends(CatchmentRst) :: CatchmentCNRst
     logical :: isCLM45
     integer :: VAR_COL
     integer :: VAR_PFT
     real, allocatable ::    cnity(:,:)
     real, allocatable ::    fvg(:,:)
     real, allocatable ::    tg(:,:)
     real, allocatable ::    TILE_ID(:)
     real, allocatable ::    ndep(:)
     real, allocatable ::    t2(:)
     real, allocatable ::    BGALBVR(:)
     real, allocatable ::    BGALBVF(:)
     real, allocatable ::    BGALBNR(:)
     real, allocatable ::    BGALBNF(:)
     real, allocatable ::    CNCOL(:,:)
     real, allocatable ::    CNPFT(:,:)
     real, allocatable ::    ABM     (:)
     real, allocatable ::    FIELDCAP(:)
     real, allocatable ::    HDM     (:)
     real, allocatable ::    GDP     (:)
     real, allocatable ::    PEATF   (:) 
  contains
     procedure :: write_nc4
     procedure :: allocate_cn   
     procedure :: add_bcs_to_rst   
  endtype CatchmentCNRst

  interface CatchmentCNRst
     module procedure CatchmentCNRst_Create
  end interface

contains

  function CatchmentCNRst_create(filename, cnclm, rc) result (catch)
    type(CatchmentCNRst) :: catch
    character(*), intent(in) :: filename
    character(*), intent(in) :: cnclm
    integer, optional, intent(out) :: rc
    integer :: status
    type(Netcdf4_fileformatter) :: formatter
    integer :: filetype, ntiles, unit
    integer :: j, dim1,dim2
    type(Variable), pointer :: myVariable
    character(len=:), pointer :: dname
    type(FileMetadata) :: meta
    character(len=256) :: Iam = "CatchmentCNRst_create"

     call MAPL_NCIOGetFileType(filename, filetype, __RC__)
     if (filetype /= 0) then
        _ASSERT( .false., "CatchmentCN only support nc4 file restart")
     endif

     catch%isCLM45 = .false.
     call formatter%open(filename, pFIO_READ, __RC__)
     meta   = formatter%read(__RC__)
     ntiles = meta%get_dimension('tile', __RC__)
     catch%ntiles = ntiles
     if (index(cnclm, '40' /=0) then
        catch%VAR_COL = VAR_COL_CLM40
        catch%VAR_PFT = VAR_PFT_CLM40
     endif
     if (index(cnclm, '45' /=0) then
        catch%VAR_COL = VAR_COL_CLM45
        catch%VAR_PFT = VAR_PFT_CLM45
        catch%isCLM45 = .true.
     endif

     call catch%allocate_cn(__RC__)
     call catch%read_shared_nc4(formatter, __RC__)

     myVariable => meta%get_variable("ITY")
     dname => myVariable%get_ith_dimension(2)
     dim1 = meta%get_dimension(dname)
     do j=1,dim1
        call MAPL_VarRead(formatter,"ITY",catch%cnity(:,j),offset1=j, __RC__)
        call MAPL_VarRead(formatter,"FVG",catch%fvg(:,j),offset1=j, __RC__)
     enddo

     call MAPL_VarRead(formatter,"TG",catch%tg, __RC__)
     call MAPL_VarRead(formatter,"TILE_ID",catch%TILE_ID, __RC__)
     call MAPL_VarRead(formatter,"NDEP",catch%ndep, __RC__)
     call MAPL_VarRead(formatter,"CLI_T2M",catch%t2, __RC__)
     call MAPL_VarRead(formatter,"BGALBVR",catch%BGALBVR, __RC__)
     call MAPL_VarRead(formatter,"BGALBVF",catch%BGALBVF, __RC__)
     call MAPL_VarRead(formatter,"BGALBNR",catch%BGALBNR, __RC__)
     call MAPL_VarRead(formatter,"BGALBNF",catch%BGALBNF, __RC__)

     myVariable => meta%get_variable("CNCOL")
     dname => myVariable%get_ith_dimension(2)
     dim1 = meta%get_dimension(dname)
     if( catch%isCLM45) then
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
     myVariable => meta%get_variable("CNPFT")
     dname => myVariable%get_ith_dimension(2)
     dim1 = meta%get_dimension(dname)
     do j=1,dim1
        call MAPL_VarRead(formatter,"CNPFT",catch%CNPFT(:,j),offset1=j, __RC__)
     enddo
     call formatter%close()
     if (present(rc)) rc =0
   end function CatchmentCNRst_Create

   subroutine write_nc4(this, filename, meta, rc)
     class(CatchmentCNRst), intent(inout):: this
     character(*), intent(in) :: filename
     type(FileMetadata), intent(inout) :: meta
     integer, optional, intent(out):: rc

     type(Netcdf4_fileformatter) :: formatter
     integer :: status
     character(256) :: Iam = "write_nc4"
     integer :: i,j, dim1,dim2
     real, dimension (:), allocatable :: var
     type(Variable), pointer :: myVariable
     character(len=:), pointer :: dname

     call formatter%create(filename, __RC__)
     call formatter%write(meta, __RC__)

     call this%write_shared_nc4(formatter, __RC__)

     myVariable => meta%get_variable("ITY")
     dname => myVariable%get_ith_dimension(2)
     dim1 = meta%get_dimension(dname)
     do j=1,dim1
        call MAPL_VarWrite(formatter,"ITY",this%cnity(:,j),offset1=j)
        call MAPL_VarWrite(formatter,"FVG",this%fvg(:,j),offset1=j)
     enddo

     call MAPL_VarWrite(formatter,"TILE_ID",this%TILE_ID)
     call MAPL_VarWrite(formatter,"NDEP",this%NDEP)
     call MAPL_VarWrite(formatter,"CLI_T2M",this%t2)
     call MAPL_VarWrite(formatter,"BGALBVR",this%BGALBVR)
     call MAPL_VarWrite(formatter,"BGALBVF",this%BGALBVF)
     call MAPL_VarWrite(formatter,"BGALBNR",this%BGALBNR)
     call MAPL_VarWrite(formatter,"BGALBNF",this%BGALBNF)
     myVariable => meta%get_variable("CNCOL")
     dname => myVariable%get_ith_dimension(2)
     dim1 = meta%get_dimension(dname)

     do j=1,dim1
        call MAPL_VarWrite(formatter,"CNCOL",this%CNCOL(:,j),offset1=j)
     enddo
     myVariable => meta%get_variable("CNPFT")
     dname => myVariable%get_ith_dimension(2)
     dim1 = meta%get_dimension(dname)
     do j=1,dim1
        call MAPL_VarWrite(formatter,"CNPFT",this%CNPFT(:,j),offset1=j)
     enddo

     dim1 = meta%get_dimension('tile')
     allocate (var (dim1))
     var = 0.

     call MAPL_VarWrite(formatter,"BFLOWM", var)
     call MAPL_VarWrite(formatter,"TOTWATM",var)
     call MAPL_VarWrite(formatter,"TAIRM",  var)
     call MAPL_VarWrite(formatter,"TPM",    var)
     call MAPL_VarWrite(formatter,"CNSUM",  var)
     call MAPL_VarWrite(formatter,"SNDZM",  var)
     call MAPL_VarWrite(formatter,"ASNOWM", var)

     myVariable => meta%get_variable("TGWM")
     dname => myVariable%get_ith_dimension(2)
     dim1 = meta%get_dimension(dname)
     do j=1,dim1
        call MAPL_VarWrite(formatter,"TGWM",var,offset1=j)
        call MAPL_VarWrite(formatter,"RZMM",var,offset1=j)
     end do

     if (this%isCLM45) then
        do j=1,dim1
           call MAPL_VarWrite(formatter,"SFMM",  var,offset1=j)
        enddo

        call MAPL_VarWrite(formatter,"ABM",     this%ABM, rc =rc     )
        call MAPL_VarWrite(formatter,"FIELDCAP",this%FIELDCAP)
        call MAPL_VarWrite(formatter,"HDM",     this%HDM     )
        call MAPL_VarWrite(formatter,"GDP",     this%GDP     )
        call MAPL_VarWrite(formatter,"PEATF",   this%PEATF   )
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
     myVariable => meta%get_variable("PSNSUNM")
     dname => myVariable%get_ith_dimension(2)
     dim1 = meta%get_dimension(dname)
     dname => myVariable%get_ith_dimension(3)
     dim2 = meta%get_dimension(dname)
     do i=1,dim2
        do j=1,dim1
           call MAPL_VarWrite(formatter,"PSNSUNM",var,offset1=j,offset2=i)
           call MAPL_VarWrite(formatter,"PSNSHAM",var,offset1=j,offset2=i)
        end do
     end do
     call formatter%close()
     
     _RETURN(_SUCCESS)
   end subroutine write_nc4

   subroutine allocate_cn(this, ncol, npft, rc)
     class(CatchmentCNRst), intent(inout) :: this
     integer, intent(in) :: ncol,npft
     integer, optional, intent(out):: rc
     integer :: status

     ntiles = this%ntiles
     ncol = nzone* this%VAR_COL 
     npft = nzone*nveg*this%VAR_PFT

     call this%CatchmentRst%allocate_catch(__RC__)

     allocate(this%cnity(ntiles,4))
     allocate(this%fvg(ntiles,4))
     allocate(this%tg(ntiles,4))
     allocate(this%TILE_ID(ntiles))
     allocate(this%ndep(ntiles))
     allocate(this%t2(ntiles))
     allocate(this%BGALBVR(ntiles))
     allocate(this%BGALBVF(ntiles))
     allocate(this%BGALBNR(ntiles))
     allocate(this%BGALBNF(ntiles))
     allocate(this%CNCOL(ntiles,ncol))
     allocate(this%CNPFT(ntiles,npft))
     allocate(this%ABM(ntiles))
     allocate(this%FIELDCAP(ntiles))
     allocate(this%HDM(ntiles))
     allocate(this%GDP(ntiles))
     allocate(this%PEATF(ntiles))
     _RETURN(_SUCCESS)
   end subroutine allocatecn

   SUBROUTINE add_bcs_to_rst (this, SURFLAY, DataDir,rc)
    class(CatchmentCNRst), intent(inout) :: this
    real, intent (in)                    :: SURFLAY
    character(*), intent (in)            :: DataDir
    integer, optional, intent(out) :: rc
    real, allocatable :: CLMC_pf1(:), CLMC_pf2(:), CLMC_sf1(:), CLMC_sf2(:)
    real, allocatable :: CLMC_pt1(:), CLMC_pt2(:), CLMC_st1(:), CLMC_st2(:)    
    real, allocatable :: CLMC45_pf1(:), CLMC45_pf2(:), CLMC45_sf1(:), CLMC45_sf2(:)
    real, allocatable :: CLMC45_pt1(:), CLMC45_pt2(:), CLMC45_st1(:), CLMC45_st2(:)    
    real, allocatable :: NDEP(:), BVISDR(:), BVISDF(:), BNIRDR(:), BNIRDF(:) 
    real, allocatable :: T2(:), var1(:), hdm(:), fc(:), gdp(:), peatf(:)
    integer, allocatable :: ity(:), abm (:)
    integer       :: STATUS, ntiles, uit27, unit28, unit29, unit30
    integer       :: idum, i,j,n, ib, nv
    real          :: rdum, zdep1, zdep2, zdep3, zmet, term1, term2, bare,fvg(4)
    logical       :: NEWLAND
    logical       :: file_exists

    type(NetCDF4_Fileformatter) :: CatchCNFmt
    character*256        :: Iam = "add_bcs"
 
    call this%CatchmentRst%add_bcs(surflay, DataDir, __RC__)
    ntiles = this%ntiles
    allocate (BVISDR(ntiles),  BVISDF(ntiles),  BNIRDR(ntiles)  )
    allocate (BNIRDF(ntiles),      T2(ntiles),    NDEP(ntiles)  )    
    allocate (CLMC_pf1(ntiles), CLMC_pf2(ntiles), CLMC_sf1(ntiles))
    allocate (CLMC_sf2(ntiles), CLMC_pt1(ntiles), CLMC_pt2(ntiles))
    allocate (CLMC45_pf1(ntiles), CLMC45_pf2(ntiles), CLMC45_sf1(ntiles))
    allocate (CLMC45_sf2(ntiles), CLMC45_pt1(ntiles), CLMC45_pt2(ntiles))
    allocate (CLMC_st1(ntiles), CLMC_st2(ntiles))
    allocate (CLMC45_st1(ntiles), CLMC45_st2(ntiles))
    allocate (hdm(ntiles), fc(ntiles), gdp(ntiles))
    allocate (peatf(ntiles), abm(ntiles), var1(ntiles)

    inquire(file = trim(DataDir)//'/catchcn_params.nc4', exist=file_exists)
    inquire(file = trim(DataDir)//"CLM_veg_typs_fracs"   ,exist=NewLand )
    _ASSERT(Newland, "catchcn should get bc from newland")

    if(file_exists) then
       call CatchCNFmt%Open(trim(DataDir)//'/catchcn_params.nc4', pFIO_READ, __RC__)    
       call MAPL_VarRead ( CatchCNFmt ,'BGALBNF', BNIRDF, __RC__)
       call MAPL_VarRead ( CatchCNFmt ,'BGALBNR', BNIRDR, __RC__)
       call MAPL_VarRead ( CatchCNFmt ,'BGALBVF', BVISDF, __RC__)
       call MAPL_VarRead ( CatchCNFmt ,'BGALBVR', BVISDR, __RC__)
       call MAPL_VarRead ( CatchCNFmt ,'NDEP', NDEP, __RC__)
       call MAPL_VarRead ( CatchCNFmt ,'T2_M', T2, __RC__)
       call MAPL_VarRead(CatchCNFmt,'ITY',CLMC_pt1,offset1=1, __RC__)     !  30
       call MAPL_VarRead(CatchCNFmt,'ITY',CLMC_pt2,offset1=2, __RC__)     !  31
       call MAPL_VarRead(CatchCNFmt,'ITY',CLMC_st1,offset1=3, __RC__)     !  32
       call MAPL_VarRead(CatchCNFmt,'ITY',CLMC_st2,offset1=4, __RC__)     !  33
       call MAPL_VarRead(CatchCNFmt,'FVG',CLMC_pf1,offset1=1, __RC__)     !  34
       call MAPL_VarRead(CatchCNFmt,'FVG',CLMC_pf2,offset1=2, __RC__)     !  35
       call MAPL_VarRead(CatchCNFmt,'FVG',CLMC_sf1,offset1=3, __RC__)     !  36
       call MAPL_VarRead(CatchCNFmt,'FVG',CLMC_sf2,offset1=4, __RC__)     !  37
       call CatchCNFmt%close()
    else

       open(newunit=unit27, file=trim(DataDir)//'CLM_veg_typs_fracs'   ,form='formatted')
       open(newunit=unit28, file=trim(DataDir)//'CLM_NDep_SoilAlb_T2m' ,form='formatted')

       do n=1,ntiles
          read (unit27, *) i,j, CLMC_pt1(n), CLMC_pt2(n), CLMC_st1(n), CLMC_st2(n), &
                CLMC_pf1(n), CLMC_pf2(n), CLMC_sf1(n), CLMC_sf2(n)
             
          read (unit28, *) NDEP(n), BVISDR(n), BVISDF(n), BNIRDR(n), BNIRDF(n), T2(n) ! MERRA-2 Annual Mean Temp is default.
          if(this%isCLM45) then
          endif
       end do
       
       CLOSE (unit27, STATUS = 'KEEP')
       CLOSE (unit28, STATUS = 'KEEP')

    endif

    if (this%isCLM45 ) then

      open(newunit=unit29, file=trim(DataDir)//'CLM4.5_veg_typs_fracs',form='formatted')
      open(newunit=unit30, file=trim(DataDir)//'CLM4.5_abm_peatf_gdp_hdm_fc' ,form='formatted')
      do n=1,ntiles
         read (unit29, *) i,j, CLMC45_pt1(n), CLMC45_pt2(n), CLMC45_st1(n), CLMC45_st2(n), &
                   CLMC45_pf1(n), CLMC45_pf2(n), CLMC45_sf1(n), CLMC45_sf2(n)
         read (unit30, *) i, j, abm(n), peatf(n), &
               gdp(n), hdm(n), fc(n)
      end do
      CLOSE (unit29, STATUS = 'KEEP')
      CLOSE (unit30, STATUS = 'KEEP')
    endif
    
    do n=1,ntiles
      BVISDR(n) = amax1(1.e-6, BVISDR(n))
      BVISDF(n) = amax1(1.e-6, BVISDF(n))
      BNIRDR(n) = amax1(1.e-6, BNIRDR(n))
      BNIRDF(n) = amax1(1.e-6, BNIRDF(n))

      ! convert % to fractions
      
      CLMC_pf1(n) = CLMC_pf1(n) / 100.
      CLMC_pf2(n) = CLMC_pf2(n) / 100.
      CLMC_sf1(n) = CLMC_sf1(n) / 100.
      CLMC_sf2(n) = CLMC_sf2(n) / 100.
      
      fvg(1) = CLMC_pf1(n)
      fvg(2) = CLMC_pf2(n)
      fvg(3) = CLMC_sf1(n)
      fvg(4) = CLMC_sf2(n)
      
      BARE = 1.      
      
      DO NV = 1, NVEG
         BARE = BARE - FVG(NV)! subtract vegetated fractions 
      END DO
      
      if (BARE /= 0.) THEN
         IB = MAXLOC(FVG(:),1)
         FVG (IB) = FVG(IB) + BARE ! This also corrects all cases sum ne 0.
      ENDIF
      
      CLMC_pf1(n) = fvg(1)
      CLMC_pf2(n) = fvg(2)
      CLMC_sf1(n) = fvg(3)
      CLMC_sf2(n) = fvg(4)
    enddo

    if(this%isCLM45) then
       do n =1, ntiles
         CLMC45_pf1(n) = CLMC45_pf1(n) / 100.
         CLMC45_pf2(n) = CLMC45_pf2(n) / 100.
         CLMC45_sf1(n) = CLMC45_sf1(n) / 100.
         CLMC45_sf2(n) = CLMC45_sf2(n) / 100.
         
         fvg(1) = CLMC45_pf1(n)
         fvg(2) = CLMC45_pf2(n)
         fvg(3) = CLMC45_sf1(n)
         fvg(4) = CLMC45_sf2(n)
         
         BARE = 1.      
         
         DO NV = 1, NVEG
            BARE = BARE - fvg(NV)! subtract vegetated fractions 
         END DO
         
         if (BARE /= 0.) THEN
            IB = MAXLOC(fvg(:),1)
            fvg (IB) = fvg(IB) + BARE ! This also corrects all cases sum ne 0.
         ENDIF
         
         CLMC45_pf1(n) = fvg(1)
         CLMC45_pf2(n) = fvg(2)
         CLMC45_sf1(n) = fvg(3)
         CLMC45_sf2(n) = fvg(4)
      enddo
    endif
       
    NDEP = NDEP * 1.e-9
    
    ! prevent trivial fractions
    ! -------------------------
    do n = 1,ntiles
       if(CLMC_pf1(n) <= 1.e-4) then
          CLMC_pf2(n) = CLMC_pf2(n) + CLMC_pf1(n)
          CLMC_pf1(n) = 0.
       endif
       
       if(CLMC_pf2(n) <= 1.e-4) then
          CLMC_pf1(n) = CLMC_pf1(n) + CLMC_pf2(n)
          CLMC_pf2(n) = 0.
       endif
       
       if(CLMC_sf1(n) <= 1.e-4) then
          if(CLMC_sf2(n) > 1.e-4) then
             CLMC_sf2(n) = CLMC_sf2(n) + CLMC_sf1(n)
          else if(CLMC_pf2(n) > 1.e-4) then
             CLMC_pf2(n) = CLMC_pf2(n) + CLMC_sf1(n)
          else if(CLMC_pf1(n) > 1.e-4) then
             CLMC_pf1(n) = CLMC_pf1(n) + CLMC_sf1(n)
          else
             stop 'fveg3'
          endif
          CLMC_sf1(n) = 0.
       endif
       
       if(CLMC_sf2(n) <= 1.e-4) then
          if(CLMC_sf1(n) > 1.e-4) then
             CLMC_sf1(n) = CLMC_sf1(n) + CLMC_sf2(n)
          else if(CLMC_pf2(n) > 1.e-4) then
             CLMC_pf2(n) = CLMC_pf2(n) + CLMC_sf2(n)
          else if(CLMC_pf1(n) > 1.e-4) then
             CLMC_pf1(n) = CLMC_pf1(n) + CLMC_sf2(n)
          else
             stop 'fveg4'
          endif
          CLMC_sf2(n) = 0.
       endif
     enddo 
     if (this%isCLM45) then
        do n = 1, ntiles
          if(CLMC45_pf1(n) <= 1.e-4) then
             CLMC45_pf2(n) = CLMC45_pf2(n) + CLMC45_pf1(n)
             CLMC45_pf1(n) = 0.
          endif
          
          if(CLMC45_pf2(n) <= 1.e-4) then
             CLMC45_pf1(n) = CLMC45_pf1(n) + CLMC45_pf2(n)
             CLMC45_pf2(n) = 0.
          endif
          
          if(CLMC45_sf1(n) <= 1.e-4) then
             if(CLMC45_sf2(n) > 1.e-4) then
                CLMC45_sf2(n) = CLMC45_sf2(n) + CLMC45_sf1(n)
             else if(CLMC45_pf2(n) > 1.e-4) then
                CLMC45_pf2(n) = CLMC45_pf2(n) + CLMC45_sf1(n)
             else if(CLMC45_pf1(n) > 1.e-4) then
                CLMC45_pf1(n) = CLMC45_pf1(n) + CLMC45_sf1(n)
             else
                stop 'fveg3'
             endif
             CLMC45_sf1(n) = 0.
          endif
          
          if(CLMC45_sf2(n) <= 1.e-4) then
             if(CLMC45_sf1(n) > 1.e-4) then
                CLMC45_sf1(n) = CLMC45_sf1(n) + CLMC45_sf2(n)
             else if(CLMC45_pf2(n) > 1.e-4) then
                CLMC45_pf2(n) = CLMC45_pf2(n) + CLMC45_sf2(n)
             else if(CLMC45_pf1(n) > 1.e-4) then
                CLMC45_pf1(n) = CLMC45_pf1(n) + CLMC45_sf2(n)
             else
                stop 'fveg4'
             endif
             CLMC45_sf2(n) = 0.
          endif
       enddo
    endif
     
    this%cnity(:,1) = CLMC_pt1
    this%cnity(:,2) = CLMC_pt2
    this%cnity(:,3) = CLMC_st1
    this%cnity(:,4) = CLMC_st2
    this%fvg(:,1) = CLMC_pf1
    this%fvg(:,2) = CLMC_pf2
    this%fvg(:,3) = CLMC_sf1
    this%fvg(:,4) = CLMC_sf2
    
    this%ndep = ndep
    this%t2   = t2
    this%BGALBVR = BVISDR
    this%BGALBVF = BVISDF
    this%BGALBNR = BNIRDR
    this%BGALBNF = BNIRDF
 
    if(this%isCLM45) then
       this%abm       = real(abm)
       this%fieldcap  = fc
       this%hdm       = hdm
       this%gdp       = gdp
       this%peatf     = peatf
     endif

     deallocate (BVISDR,  BVISDF,  BNIRDR  )
     deallocate (BNIRDF,      T2,    NDEP  )    
     deallocate (CLMC_pf1, CLMC_pf2, CLMC_sf1)
     deallocate (CLMC_sf2, CLMC_pt1, CLMC_pt2)
     deallocate (CLMC_st1,CLMC_st2)

     _RETURN(_SUCCESS)
  END SUBROUTINE add_bcs_to_rst

end module CatchmentCNRstMod
