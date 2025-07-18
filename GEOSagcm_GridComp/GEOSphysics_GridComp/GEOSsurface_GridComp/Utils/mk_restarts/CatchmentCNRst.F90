#include "MAPL_Generic.h"

module CatchmentCNRstMod
  use mk_restarts_getidsMod, ONLY:      &
       GetIds  
  use mpi
  use ESMF
  use MAPL
  use CatchmentRstMod, only : CatchmentRst
  use clm_varpar_shared , only : nzone => NUM_ZON_CN, nveg_40 => NUM_VEG_CN, nveg_51 => NUM_VEG_CN51, &
                                 VAR_COL_40, VAR_PFT_40, VAR_COL_45, VAR_PFT_45, &
                                 VAR_COL_51, VAR_PFT_51, &
                                 npft => numpft_CN, npft_51 => numpft_CN51
  use nanMod         , only : nan
  
  implicit none

  real,    parameter :: fmin= 1.e-4 ! ignore vegetation fractions at or below this value
  integer :: iclass_40(npft) = (/1,1,2,3,3,4,5,5,6,7,8,9,10,11,12,11,12,11,12/)
  integer :: iclass_45(npft) = (/1,1,2,3,3,4,5,5,6,7,8,9,10,11,12,11,12,11,12/)
  integer :: iclass_51(npft_51) = (/1,1,2,3,3,4,5,5,6,7,9,10,11,11,11/)
  integer, dimension(:), allocatable :: iclass

  type, extends(CatchmentRst) :: CatchmentCNRst
     logical :: isCLM45
     logical :: isCLM51
     logical :: isCLM40

     integer :: VAR_COL
     integer :: VAR_PFT
     integer :: NVEG
     real, allocatable ::    cnity(:,:)
     real, allocatable ::    fvg(:,:)
     real, allocatable ::    tg(:,:)
     real, allocatable ::    tgwm(:,:)
     real, allocatable ::    rzmm(:,:)
     real, allocatable ::    sfmm(:,:)
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

     real, allocatable :: bflowm(:)
     real, allocatable :: totwatm(:)
     real, allocatable :: tairm(:)
     real, allocatable :: tpm(:)
     real, allocatable :: cnsum(:)
     real, allocatable :: sndzm(:)
     real, allocatable :: asnowm(:)
     real, allocatable :: ar1m(:)
     real, allocatable :: rainfm(:)
     real, allocatable :: rhm(:)
     real, allocatable :: runsrfm(:)
     real, allocatable :: snowfm(:)
     real, allocatable :: windm(:)
     real, allocatable :: tprec10d(:)
     real, allocatable :: tprec60d(:)
     real, allocatable :: t2m10d(:)
     real, allocatable :: sfmcm(:)
     real, allocatable :: psnsunm(:,:,:)
     real, allocatable :: psnsham(:,:,:)
     real, allocatable :: lmrsunm(:,:,:)
     real, allocatable :: lmrsham(:,:,:)
     real, allocatable :: laisunm(:,:,:)
     real, allocatable :: laisham(:,:,:)
     real, allocatable :: rh30d(:) 
     real, allocatable :: tg10d(:)
     real, allocatable :: t2mmin5d(:)
     real, allocatable :: sndzm5d(:)    
     
  contains
     procedure :: write_nc4
     procedure :: allocate_cn   
     procedure :: add_bcs_to_cnrst   
     procedure :: re_tile
  endtype CatchmentCNRst

  interface CatchmentCNRst
     module procedure CatchmentCNRst_Create
  end interface

contains

  function CatchmentCNRst_create(filename, cnclm, time, rc) result (catch)
    type(CatchmentCNRst) :: catch
    character(*), intent(in) :: filename
    character(*), intent(in) :: cnclm
    character(*), intent(in) :: time
    integer, optional, intent(out) :: rc
    integer :: status
    type(Netcdf4_fileformatter) :: formatter
    integer :: filetype, ntiles, unit
    integer :: j, dim1,dim2, myid, mpierr
    type(Variable), pointer :: myVariable
    character(len=:), pointer :: dname
    type(FileMetadata) :: meta
    character(len=256) :: Iam = "CatchmentCNRst_create"

     call MAPL_NCIOGetFileType(filename, filetype, __RC__)
     if (filetype /= 0) then
        _ASSERT( .false., "CatchmentCN only support nc4 file restart")
     endif
  
     call MPI_COMM_RANK( MPI_COMM_WORLD, myid, mpierr )
 
     catch%isCLM45 = .false.
     catch%isCLM51 = .false.
     catch%isCLM40 = .false.

     call formatter%open(filename, pFIO_READ, __RC__)
     meta  = formatter%read(__RC__)
     ntiles = meta%get_dimension('tile', __RC__)
     catch%ntiles = ntiles
     catch%meta  = meta
     catch%time = time
     if (index(cnclm, '40') /=0) then
        catch%isCLM40 = .true.
        catch%VAR_COL = VAR_COL_40
        catch%VAR_PFT = VAR_PFT_40
        catch%NVEG    = nveg_40
     endif
     if (index(cnclm, '45') /=0) then
        catch%isCLM45 = .true.
        catch%VAR_COL = VAR_COL_45
        catch%VAR_PFT = VAR_PFT_45
        catch%nveg    = nveg_40
     endif
     if (index(cnclm, '51') /=0) then
        catch%VAR_COL = VAR_COL_51
        catch%VAR_PFT = VAR_PFT_51
        catch%isCLM51 = .true.
        catch%nveg    = nveg_51
     endif

     if (myid == 0) then
        call catch%allocate_cn(__RC__)
        call catch%read_shared_nc4(formatter, __RC__)

        call MAPL_VarRead(formatter,"ITY",catch%cnity, __RC__)
        call MAPL_VarRead(formatter,"FVG",catch%fvg, __RC__)

        call MAPL_VarRead(formatter,"TG",catch%tg, __RC__)
        call MAPL_VarRead(formatter,"TILE_ID",catch%TILE_ID, __RC__)
        call MAPL_VarRead(formatter,"NDEP",catch%ndep, __RC__)
        call MAPL_VarRead(formatter,"CLI_T2M",catch%t2, __RC__)
        call MAPL_VarRead(formatter,"BGALBVR",catch%BGALBVR, __RC__)
        call MAPL_VarRead(formatter,"BGALBVF",catch%BGALBVF, __RC__)
        call MAPL_VarRead(formatter,"BGALBNR",catch%BGALBNR, __RC__)
        call MAPL_VarRead(formatter,"BGALBNF",catch%BGALBNF, __RC__)

        if( catch%isCLM40 ) then
           call MAPL_VarRead(formatter,"SFMCM",   catch%sfmcm , __RC__)
        endif

        if( catch%isCLM45 ) then
           call MAPL_VarRead(formatter,"ABM",     catch%ABM, __RC__)
           call MAPL_VarRead(formatter,"FIELDCAP",catch%FIELDCAP, __RC__)
           call MAPL_VarRead(formatter,"HDM",     catch%HDM     , __RC__)
           call MAPL_VarRead(formatter,"GDP",     catch%GDP     , __RC__)
           call MAPL_VarRead(formatter,"PEATF",   catch%PEATF   , __RC__)


           call MAPL_VarRead(formatter,"RHM",       catch%rhm        , __RC__)
           call MAPL_VarRead(formatter,"WINDM",     catch%windm      , __RC__)
           call MAPL_VarRead(formatter,"RAINFM",    catch%rainfm     , __RC__)
           call MAPL_VarRead(formatter,"SNOWFM",    catch%snowfm     , __RC__)
           call MAPL_VarRead(formatter,"RUNSRFM",   catch%runsrfm    , __RC__)
           call MAPL_VarRead(formatter,"AR1M",      catch%ar1m       , __RC__)
           call MAPL_VarRead(formatter,"T2M10D",    catch%T2M10D     , __RC__)
           call MAPL_VarRead(formatter,"TPREC10D",  catch%TPREC10D   , __RC__)
           call MAPL_VarRead(formatter,"TPREC60D",  catch%TPREC60D   , __RC__)
           call MAPL_VarRead(formatter,"SFMM",      catch%sfmm       , __RC__)
        endif

        if( catch%isCLM51) then
           call MAPL_VarRead(formatter,"ABM",     catch%ABM, __RC__)
           call MAPL_VarRead(formatter,"FIELDCAP",catch%FIELDCAP, __RC__)
           call MAPL_VarRead(formatter,"HDM",     catch%HDM     , __RC__)
           call MAPL_VarRead(formatter,"GDP",     catch%GDP     , __RC__)
           call MAPL_VarRead(formatter,"PEATF",   catch%PEATF   , __RC__)
           call MAPL_VarRead(formatter,"RHM",     catch%RHM     , __RC__)
           call MAPL_VarRead(formatter,"WINDM",   catch%WINDM   , __RC__)
           call MAPL_VarRead(formatter,"RAINFM",  catch%RAINFM  , __RC__)
           call MAPL_VarRead(formatter,"SNOWFM",  catch%SNOWFM  , __RC__)
           call MAPL_VarRead(formatter,"RUNSRFM", catch%RUNSRFM, __RC__)
           call MAPL_VarRead(formatter,"AR1M",    catch%AR1M    , __RC__)
           call MAPL_VarRead(formatter,"SNDZM5D", catch%SNDZM5D , __RC__)
           call MAPL_VarRead(formatter,"T2M10D",  catch%T2M10D  , __RC__)
           call MAPL_VarRead(formatter,"T2MMIN5D",catch%T2MMIN5D, __RC__)
           call MAPL_VarRead(formatter,"TG10D",   catch%TG10D   , __RC__)
           call MAPL_VarRead(formatter,"RH30D",   catch%RH30D   , __RC__)
           call MAPL_VarRead(formatter,"TPREC10D",catch%TPREC10D, __RC__)
           call MAPL_VarRead(formatter,"TPREC60D",catch%TPREC60D, __RC__)
        endif

        call MAPL_VarRead(formatter,"CNCOL",catch%CNCOL, __RC__)

        ! The following three lines were added as a bug fix by smahanam on 5 Oct 2020
        ! (to be merged into the "develop" branch in late 2020):
        ! The length of the 2nd dim of CNPFT differs from that of CNCOL.  Prior to this fix,
        ! CNPFT was not read in its entirety and some elements remained uninitialized (or zero),
        ! resulting in bad values in the "regridded" (re-tiled) restart file. 
        ! This impacted re-tiled restarts for both CNCLM40 and CLCLM45.
        ! - reichle, 23 Nov 2020
        call MAPL_VarRead(formatter,"CNPFT",catch%CNPFT, __RC__)

        ! more reading
        call MAPL_VarRead(formatter,  "BFLOWM",  catch%bflowm ,_RC) 
        call MAPL_VarRead(formatter,  "TOTWATM", catch%totwatm,_RC) 
        call MAPL_VarRead(formatter,  "TAIRM",   catch%tairm  ,_RC) 
        call MAPL_VarRead(formatter,  "TPM",     catch%tpm    ,_RC) 
        call MAPL_VarRead(formatter,  "CNSUM",   catch%cnsum  ,_RC) 
        call MAPL_VarRead(formatter,  "SNDZM",   catch%sndzm  ,_RC) 
        call MAPL_VarRead(formatter,  "ASNOWM",  catch%asnowm ,_RC) 
        call MAPL_VarRead(formatter,  "PSNSUNM", catch%psnsunm,_RC) 
        call MAPL_VarRead(formatter,  "PSNSHAM", catch%psnsham,_RC) 
        call MAPL_VarRead(formatter,  "LMRSUNM", catch%lmrsunm,_RC)
        call MAPL_VarRead(formatter,  "LMRSHAM", catch%lmrsham,_RC)
        call MAPL_VarRead(formatter,  "LAISUNM", catch%laisunm,_RC)
        call MAPL_VarRead(formatter,  "LAISHAM", catch%laisham,_RC)
        call MAPL_VarRead(formatter,  "RZMM",    catch%rzmm   ,_RC)
        call MAPL_VarRead(formatter,  "TGWM",    catch%tgwm   ,_RC)
     endif

     call formatter%close()

     if (present(rc)) rc =0
   end function CatchmentCNRst_Create

  function CatchmentCNRst_empty(meta, cnclm, time, rc) result (catch)
    type(CatchmentCNRst) :: catch
    type(FileMetadata), intent(in) :: meta
    character(*), intent(in) :: cnclm
    character(*), intent(in) :: time
    integer, optional, intent(out) :: rc
    integer :: status, myid, mpierr
    character(len=256) :: Iam = "CatchmentCNRst_empty"

     catch%isCLM45 = .false.
     catch%isCLM51 = .false.
     catch%isCLM40 = .false.

     catch%ntiles = meta%get_dimension('tile', __RC__)
     catch%time = time
     catch%meta = meta
     if (index(cnclm, '40') /=0) then
        catch%isCLM40 = .true.
        catch%VAR_COL = VAR_COL_40
        catch%VAR_PFT = VAR_PFT_40
     endif
     if (index(cnclm, '45') /=0) then
        catch%isCLM45 = .true.
        catch%VAR_COL = VAR_COL_45
        catch%VAR_PFT = VAR_PFT_45
     endif
     if (index(cnclm, '51') /=0) then
        catch%VAR_COL = VAR_COL_51
        catch%VAR_PFT = VAR_PFT_51
        catch%isCLM51 = .true.
     endif

     call MPI_COMM_RANK( MPI_COMM_WORLD, myid, mpierr )
     if (myid ==0) call catch%allocate_cn(__RC__)
     if(present(rc)) rc = 0
   end function CatchmentCNRst_empty

   subroutine write_nc4(this, filename, rc)
     class(CatchmentCNRst), intent(inout):: this
     character(*), intent(in) :: filename
     integer, optional, intent(out):: rc

     type(Netcdf4_fileformatter) :: formatter
     integer :: status
     character(256) :: Iam = "write_nc4"
     integer :: i,j, dim1,dim2
     real, dimension (:), allocatable :: var
     type(Variable), pointer :: myVariable
     character(len=:), pointer :: dname
     type(FileMetadata) :: meta

     meta = this%meta
     call formatter%create(filename, __RC__)
     call formatter%write(meta, __RC__)

     call this%write_shared_nc4(formatter, __RC__)

     call MAPL_VarWrite(formatter,"ITY",this%cnity)
     call MAPL_VarWrite(formatter,"FVG",this%fvg)
     call MAPL_VarWrite(formatter,"TG",this%tg)

     call MAPL_VarWrite(formatter,"TILE_ID",this%TILE_ID)
     call MAPL_VarWrite(formatter,"NDEP",this%NDEP)
     call MAPL_VarWrite(formatter,"CLI_T2M",this%t2)
     call MAPL_VarWrite(formatter,"BGALBVR",this%BGALBVR)
     call MAPL_VarWrite(formatter,"BGALBVF",this%BGALBVF)
     call MAPL_VarWrite(formatter,"BGALBNR",this%BGALBNR)
     call MAPL_VarWrite(formatter,"BGALBNF",this%BGALBNF)
     call MAPL_VarWrite(formatter,"CNCOL",this%CNCOL)
     call MAPL_VarWrite(formatter,"CNPFT",this%CNPFT)

     call MAPL_VarWrite(formatter,"BFLOWM", this%bflowm )
     call MAPL_VarWrite(formatter,"TOTWATM",this%totwatm)
     call MAPL_VarWrite(formatter,"TAIRM",  this%tairm  )
     call MAPL_VarWrite(formatter,"TPM",    this%tpm    )
     call MAPL_VarWrite(formatter,"CNSUM",  this%cnsum  )
     call MAPL_VarWrite(formatter,"SNDZM",  this%sndzm  )
     call MAPL_VarWrite(formatter,"ASNOWM", this%asnowm )
     call MAPL_VarWrite(formatter,"TGWM",   this%tgwm)
     call MAPL_VarWrite(formatter,"RZMM",   this%rzmm)

     if (this%isCLM45) then
        call MAPL_VarWrite(formatter,"SFMM",  this%sfmm)

        call MAPL_VarWrite(formatter,"ABM",     this%ABM, rc =rc     )
        call MAPL_VarWrite(formatter,"FIELDCAP",this%FIELDCAP)
        call MAPL_VarWrite(formatter,"HDM",     this%HDM     )
        call MAPL_VarWrite(formatter,"GDP",     this%GDP     )
        call MAPL_VarWrite(formatter,"PEATF",   this%PEATF   )

        call MAPL_VarWrite(formatter,"RHM",     this%rhm      )
        call MAPL_VarWrite(formatter,"WINDM",   this%windm    )
        call MAPL_VarWrite(formatter,"RAINFM",  this%rainfm   )
        call MAPL_VarWrite(formatter,"SNOWFM",  this%snowfm   )
        call MAPL_VarWrite(formatter,"RUNSRFM", this%runsrfm  )
        call MAPL_VarWrite(formatter,"AR1M",    this%ar1m     )
        call MAPL_VarWrite(formatter,"T2M10D",  this%t2m10d   )
        call MAPL_VarWrite(formatter,"TPREC10D",this%tprec10d )
        call MAPL_VarWrite(formatter,"TPREC60D",this%tprec60d )
        call MAPL_VarWrite(formatter,"LMRSUNM", this%LMRSUNM )
        call MAPL_VarWrite(formatter,"LMRSHAM", this%LMRSHAM )

     elseif (this%isCLM51) then

         call MAPL_VarWrite(formatter,"SFMM",  this%sfmm)

          call MAPL_VarWrite(formatter,"ABM",     this%ABM, rc =rc     )
          call MAPL_VarWrite(formatter,"FIELDCAP",this%FIELDCAP)
          call MAPL_VarWrite(formatter,"HDM",     this%HDM     )
          call MAPL_VarWrite(formatter,"GDP",     this%GDP     )
          call MAPL_VarWrite(formatter,"PEATF",   this%PEATF   )
          call MAPL_VarWrite(formatter,"RHM",     this%RHM)
          call MAPL_VarWrite(formatter,"WINDM",   this%WINDM)
          call MAPL_VarWrite(formatter,"RAINFM",  this%RAINFM)
          call MAPL_VarWrite(formatter,"SNOWFM",  this%SNOWFM)
          call MAPL_VarWrite(formatter,"RUNSRFM", this%RUNSRFM)
          call MAPL_VarWrite(formatter,"AR1M",    this%AR1M)
          call MAPL_VarWrite(formatter,"SNDZM5D", this%SNDZM5D)
          call MAPL_VarWrite(formatter,"T2M10D",  this%T2M10D)
          call MAPL_VarWrite(formatter,"T2MMIN5D",this%T2MMIN5D)
          call MAPL_VarWrite(formatter,"TG10D",   this%TG10D)
          call MAPL_VarWrite(formatter,"RH30D",   this%RH30D)
          call MAPL_VarWrite(formatter,"TPREC10D",this%TPREC10D)
          call MAPL_VarWrite(formatter,"TPREC60D",this%TPREC60D)
          call MAPL_VarWrite(formatter,"LMRSUNM", this%LMRSUNM )
          call MAPL_VarWrite(formatter,"LMRSHAM", this%LMRSHAM )
          call MAPL_VarWrite(formatter,"LAISUNM", this%LAISUNM )
          call MAPL_VarWrite(formatter,"LAISHAM", this%LAISHAM )

     endif

     if (this%isCLM40)   call MAPL_VarWrite(formatter,"SFMCM",  this%sfmcm)

     call MAPL_VarWrite(formatter,"PSNSUNM", this%PSNSUNM )
     call MAPL_VarWrite(formatter,"PSNSHAM", this%PSNSHAM )

     call formatter%close()
     
     _RETURN(_SUCCESS)
   end subroutine write_nc4

   subroutine allocate_cn(this,rc)
     class(CatchmentCNRst), intent(inout) :: this
     integer, optional, intent(out):: rc
     integer :: status
     integer  :: ncol,npft, ntiles, nveg
    
     
     nveg = this%NVEG
     ntiles = this%ntiles
     ncol = nzone* this%VAR_COL 
     npft = nzone*nveg*this%VAR_PFT

     call this%CatchmentRst%allocate_catch(__RC__)

     ! W.Jiang notes : some varaiables are not allocated because they are set to zero directly during write
     allocate(this%cnity(ntiles,nveg))
     allocate(this%fvg(ntiles,nveg))
     allocate(this%tg(ntiles,nveg))
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

     allocate(this%bflowm  (ntiles))
     allocate(this%totwatm (ntiles))
     allocate(this%tairm   (ntiles))
     allocate(this%tpm     (ntiles))
     allocate(this%cnsum   (ntiles))
     allocate(this%sndzm   (ntiles))
     allocate(this%asnowm  (ntiles))
     allocate(this%psnsunm(ntiles,nveg,nzone))
     allocate(this%psnsham(ntiles,nveg,nzone))
     allocate(this%lmrsunm(ntiles,nveg,nzone))
     allocate(this%lmrsham(ntiles,nveg,nzone))
     allocate(this%laisunm(ntiles,nveg,nzone))
     allocate(this%laisham(ntiles,nveg,nzone))
     allocate(this%rzmm   (ntiles,nzone))
     allocate(this%tgwm   (ntiles,nzone))

     if (this%isCLM40) then
        allocate(this%sfmcm   (ntiles))
     endif
     if (this%isCLM45) then
        allocate(this%ar1m    (ntiles))
        allocate(this%rainfm  (ntiles))
        allocate(this%rhm     (ntiles))
        allocate(this%runsrfm (ntiles))
        allocate(this%snowfm  (ntiles))
        allocate(this%windm   (ntiles))
        allocate(this%tprec10d(ntiles))
        allocate(this%tprec60d(ntiles))
        allocate(this%t2m10d  (ntiles))
        allocate(this%sfmm    (ntiles,nzone))
     endif
     if (this%isCLM51) then
        allocate(this%ar1m    (ntiles))
        allocate(this%rainfm  (ntiles))
        allocate(this%rhm     (ntiles))
        allocate(this%runsrfm (ntiles))
        allocate(this%snowfm  (ntiles))
        allocate(this%windm   (ntiles))
        allocate(this%tprec10d(ntiles))
        allocate(this%tprec60d(ntiles))
        allocate(this%t2m10d  (ntiles))
        allocate(this%sfmm    (ntiles,nzone))
        allocate(this%rh30d   (ntiles))
        allocate(this%tg10d   (ntiles))
        allocate(this%t2mmin5d(ntiles))
        allocate(this%sndzm5d (ntiles))
     endif

     _RETURN(_SUCCESS)
   end subroutine allocate_cn

   SUBROUTINE add_bcs_to_cnrst (this, SURFLAY, OutBcsDir,rc)
    class(CatchmentCNRst), intent(inout) :: this
    real, intent (in)                    :: SURFLAY
    character(*), intent (in)            :: OutBcsDir
    integer, optional, intent(out) :: rc
    real, allocatable :: CLMC_pf1(:), CLMC_pf2(:), CLMC_sf1(:), CLMC_sf2(:)
    real, allocatable :: CLMC_pt1(:), CLMC_pt2(:), CLMC_st1(:), CLMC_st2(:)    
    real, allocatable :: NDEP(:), BVISDR(:), BVISDF(:), BNIRDR(:), BNIRDF(:) 
    real, allocatable :: T2(:), hdm(:), fc(:), gdp(:), peatf(:)
    integer, allocatable :: ity(:), abm (:)
    integer       :: STATUS, ntiles, unit27, unit28, unit29, unit30
    integer       :: idum, i,j,n, ib, nv, nveg
    real          :: rdum, zdep1, zdep2, zdep3, zmet, term1, term2, bare,fvg(4)
    integer, dimension(npft) :: map_pft
    logical       :: NEWLAND
    logical       :: file_exists

    type(NetCDF4_Fileformatter) :: CatchCNFmt
    character*256        :: Iam = "add_bcs"

    nveg = this%nveg 

    open (10,file =trim(OutBcsDir)//"/clsm/catchment.def",status='old',form='formatted')
    read (10,*) ntiles
    close (10, status = 'keep')
 
    !ntiles = this%ntiles
    !call this%CatchmentRst%add_bcs_to_rst(surflay, OutBcsDir, __RC__)

    allocate (BVISDR(ntiles),  BVISDF(ntiles),  BNIRDR(ntiles)  )
    allocate (BNIRDF(ntiles),      T2(ntiles),    NDEP(ntiles)  )    
    allocate (CLMC_pf1(ntiles), CLMC_pf2(ntiles), CLMC_sf1(ntiles))
    allocate (CLMC_sf2(ntiles), CLMC_pt1(ntiles), CLMC_pt2(ntiles))
    allocate (CLMC_st1(ntiles), CLMC_st2(ntiles))
    allocate (hdm(ntiles), fc(ntiles), gdp(ntiles))
    allocate (peatf(ntiles), abm(ntiles))

    inquire(file = trim(OutBcsDir)//'/clsm/catchcn_params.nc4', exist=file_exists)
    inquire(file = trim(OutBcsDir)//'/clsm/CLM_veg_typs_fracs', exist=NewLand )
    _ASSERT(Newland, "catchcn should get bc from newland")

    if(file_exists) then
       call CatchCNFmt%Open(trim(OutBcsDir)//'/clsm/catchcn_params.nc4', pFIO_READ, __RC__)    
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

       open(newunit=unit27, file=trim(OutBcsDir)//'/clsm/CLM_veg_typs_fracs'   ,form='formatted')
       open(newunit=unit28, file=trim(OutBcsDir)//'/clsm/CLM_NDep_SoilAlb_T2m' ,form='formatted')

       do n=1,ntiles
          read (unit27, *) i,j, CLMC_pt1(n), CLMC_pt2(n), CLMC_st1(n), CLMC_st2(n), &
                CLMC_pf1(n), CLMC_pf2(n), CLMC_sf1(n), CLMC_sf2(n)
             
          read (unit28, *) NDEP(n), BVISDR(n), BVISDF(n), BNIRDR(n), BNIRDF(n), T2(n) ! MERRA-2 Annual Mean Temp is default.
       end do
       
       CLOSE (unit27, STATUS = 'KEEP')
       CLOSE (unit28, STATUS = 'KEEP')

    endif

    if ((this%isCLM45) .or. (this%isCLM51)) then

      open(newunit=unit30, file=trim(OutBcsDir)//'/clsm/CLM4.5_abm_peatf_gdp_hdm_fc' ,form='formatted')
      do n=1,ntiles
         read (unit30, *) i, j, abm(n), peatf(n), &
               gdp(n), hdm(n), fc(n)
      end do
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

     ! if using Catchment-CN5.1, reduce down to 2 PFTs
     ! step 1: map split PFTs to their parent type
     ! step 2: add up area fractions

     if (this%isCLM51) then

        map_pft = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10, 11, 12, 13, 13, 14, 14, 15, 15 /)

        do n = 1,ntiles
   
             ! map split PFTs to parent PFTs
             CLMC_pt1(n) = map_pft(CLMC_pt1(n))
             CLMC_pt2(n) = map_pft(CLMC_pt2(n))
             CLMC_st1(n) = map_pft(CLMC_st1(n))
             CLMC_st2(n) = map_pft(CLMC_st2(n))
    
             ! combine area fractions of same PFTs, 
             ! otherwise retain area fraction of single present PFT

             if (CLMC_pt1(n).eq.CLMC_pt2(n)) then
                CLMC_pf1(n) = CLMC_pf1(n) + CLMC_pf2(n)
                CLMC_pf2(n) = 0.
             else if (CLMC_pt1(n).ne.CLMC_pt2(n)) then
                CLMC_pf1(n) = maxval((/ CLMC_pf1(n), CLMC_pf2(n) /))
             endif

             if (CLMC_st1(n).eq.CLMC_st2(n)) then
                CLMC_sf1(n) = CLMC_sf1(n) + CLMC_sf2(n)
                CLMC_sf2(n) = 0. 
             else if (CLMC_st1(n).ne.CLMC_st2(n)) then
                CLMC_sf1(n) = maxval((/ CLMC_sf1(n), CLMC_sf2(n) /))
             endif
        end do
 
     endif

     if ((this%isCLM40).or.(this%isCLM45)) then
        this%cnity = reshape([CLMC_pt1,CLMC_pt2,CLMC_st1,CLMC_st2],[ntiles,4])
        this%fvg   = reshape([CLMC_pf1,CLMC_pf2,CLMC_sf1,CLMC_sf2],[ntiles,4])
     elseif (this%isCLM51) then
        this%cnity = reshape([CLMC_pt1,CLMC_st1],[ntiles,2])
        this%fvg   = reshape([CLMC_pf1,CLMC_sf1],[ntiles,2])
     endif 

     this%ndep = ndep
     this%t2   = t2
     this%BGALBVR = BVISDR
     this%BGALBVF = BVISDF
     this%BGALBNR = BNIRDR
     this%BGALBNF = BNIRDF
 
     if ((this%isCLM45) .or. (this%isCLM51))then
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
   end subroutine add_bcs_to_cnrst

   subroutine re_tile(this, InTileFile, OutBcsDir, OutTileFile, surflay, rc)
     class(CatchmentCNRst), intent(inout) :: this
     character(*), intent(in) :: InTileFile
     character(*), intent(in) :: OutBcsDir
     character(*), intent(in) :: OutTileFile
     real, intent(in)         :: surflay
     integer, optional, intent(out) :: rc
     
     real   , allocatable, dimension (:)   :: DAYX
     integer, allocatable, dimension (:)   :: low_ind, upp_ind, nt_local
     integer, allocatable, dimension (:,:) :: Id_glb_cn, id_loc_cn
     integer, allocatable, dimension (:)   :: tid_offl, id_loc
     real, allocatable, dimension (:)      :: CLMC_pf1, CLMC_pf2, CLMC_sf1, CLMC_sf2, &
         CLMC_pt1, CLMC_pt2,CLMC_st1,CLMC_st2
     integer                :: AGCM_YY,AGCM_MM,AGCM_DD,AGCM_HR=0,AGCM_DATE, &
                               AGCM_MI, AGCM_S,  dofyr

     real,    allocatable, dimension(:,:) :: fveg_offl,  ityp_offl, tg_tmp, dummy_tmp
     real, allocatable :: var_off_col (:,:,:), var_off_pft (:,:,:,:), var_out(:), var_tmp3d(:,:,:), &
                          var_out_zone(:,:)
     integer :: status, in_ntiles, out_ntiles, numprocs, npft_int
     logical :: root_proc
     integer :: mpierr, n, i, k, tag, req, st, ed, myid, L, iv, nv,nz, var_col, var_pft, nveg
     real, allocatable, dimension(:) :: lat_tmp
     type(MAPL_SunOrbit)         :: ORBIT
     type(ESMF_Time)             :: CURRENT_TIME
     type(ESMF_TimeInterval)     :: timeStep
     type(ESMF_Clock)            :: CLOCK  
     type(ESMF_Config)           :: CF

     character(*), parameter :: Iam = "CatchmentCN::Re_tile"

     nveg = this%NVEG
     in_ntiles = this%ntiles
     var_pft = this%var_pft
     var_col = this%var_col
     call this%CatchmentRst%re_tile(InTileFile, OutBcsDir, OutTileFile, surflay, _RC)

     out_ntiles = this%ntiles 
     call MPI_COMM_RANK( MPI_COMM_WORLD, myid, mpierr )
     call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, mpierr )

     root_proc = .false.
     if (myid == 0)  root_proc = .true.

     allocate(low_ind (numprocs))
     allocate(upp_ind (numprocs))
     allocate(nt_local(numprocs))
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

     allocate (CLMC_pf1(nt_local (myid + 1)))
     allocate (CLMC_pf2(nt_local (myid + 1)))
     allocate (CLMC_sf1(nt_local (myid + 1)))
     allocate (CLMC_sf2(nt_local (myid + 1)))
     allocate (CLMC_pt1(nt_local (myid + 1)))
     allocate (CLMC_pt2(nt_local (myid + 1)))
     allocate (CLMC_st1(nt_local (myid + 1)))
     allocate (CLMC_st2(nt_local (myid + 1)))
     allocate (ityp_offl (in_ntiles,nveg))
     allocate (fveg_offl (in_ntiles,nveg))
     allocate (tid_offl(in_ntiles))
     allocate (id_loc_cn (nt_local (myid + 1),nveg))

     allocate (lat_tmp(in_ntiles))

     do n = 1, in_ntiles
        tid_offl(n) = n
     enddo

     if (root_proc) then

        allocate (DAYX   (out_ntiles))

        READ(this%time(1:8),'(I8)') AGCM_DATE
        AGCM_YY = AGCM_DATE / 10000
        AGCM_MM = (AGCM_DATE - AGCM_YY*10000) / 100
        AGCM_DD = (AGCM_DATE - AGCM_YY*10000 - AGCM_MM*100)
        AGCM_MI = 0
        AGCM_S = 0


        !1) Set current date & time
        ! -----------------------

        call ESMF_CalendarSetDefault ( ESMF_CALKIND_GREGORIAN, rc=status )

        call ESMF_TimeSet  ( CURRENT_TIME, YY = AGCM_YY,       &
                                            MM = AGCM_MM,       &
                                            DD = AGCM_DD,       &
                                            H  = AGCM_HR,       &
                                            M  = AGCM_MI,       &
                                            S  = AGCM_S ,       &
                                            rc=status )
         VERIFY_(STATUS)

        !2) create a clock
        ! time interval value is not critical here, just for a clock

        call ESMF_TimeIntervalSet(TimeStep,  S=450, RC=status)
        clock = ESMF_ClockCreate(TimeStep, startTime = CURRENT_TIME, RC=status)
        VERIFY_(STATUS)
        call ESMF_ClockSet ( clock, CurrTime=CURRENT_TIME, rc=status )

        !3) create an orbit
        CF = ESMF_ConfigCreate(RC=STATUS)
        VERIFY_(status)

        ORBIT = MAPL_SunOrbitCreateFromConfig(CF, CLOCK, .false., RC=status)
        VERIFY_(status) 

        !4) current daylight duration
        lat_tmp = this%latg*MAPL_PI/180.
        call MAPL_SunGetDaylightDuration(ORBIT, lat_tmp, dayx, currTime=CURRENT_TIME,RC=STATUS)
        VERIFY_(STATUS)

        ! save the old vaues dimension (in_ntiles, nv)
        ityp_offl = this%cnity
        fveg_offl = this%fvg

        if ((this%isCLM40) .or. (this%isCLM45)) then
            npft_int = npft
        else if (this%isCLM51) then
            npft_int = npft_51
        endif

        do n = 1, in_ntiles
           do nv = 1,nveg
              if(ityp_offl(n,nv)<0 .or. ityp_offl(n,nv)>npft_int)    stop 'ityp'
              if(fveg_offl(n,nv)<0..or. fveg_offl(n,nv)>1.00001) stop 'fveg'
           end do

           if (nint(this%tile_id(n)) /= n) stop ("cannot assign ity_offl to cnity and fvg_offl to fvg")

           if ((this%isCLM40) .or. (this%isCLM45)) then

              if((ityp_offl(N,3) == 0).and.(ityp_offl(N,4) == 0)) then
                 if(ityp_offl(N,1) /= 0) then
                    ityp_offl(N,3) = ityp_offl(N,1)
                 else
                    ityp_offl(N,3) = ityp_offl(N,2)
                 endif
              endif

              if((ityp_offl(N,1) == 0).and.(ityp_offl(N,2) /= 0)) ityp_offl(N,1) = ityp_offl(N,2)
              if((ityp_offl(N,2) == 0).and.(ityp_offl(N,1) /= 0)) ityp_offl(N,2) = ityp_offl(N,1)
              if((ityp_offl(N,3) == 0).and.(ityp_offl(N,4) /= 0)) ityp_offl(N,3) = ityp_offl(N,4)
              if((ityp_offl(N,4) == 0).and.(ityp_offl(N,3) /= 0)) ityp_offl(N,4) = ityp_offl(N,3)

           elseif (this%isCLM51) then

              if((ityp_offl(N,1) == 0).and.(ityp_offl(N,2) /= 0)) then
                ityp_offl(N,1) = ityp_offl(N,2)
              else if((ityp_offl(N,2) == 0).and.(ityp_offl(N,1) /= 0)) then
                ityp_offl(N,2) = ityp_offl(N,1)
              endif   

           endif
       end do
    endif

    call MPI_BCAST(ityp_offl,size(ityp_offl),MPI_REAL   ,0,MPI_COMM_WORLD,mpierr)
    call MPI_BCAST(fveg_offl,size(fveg_offl),MPI_REAL   ,0,MPI_COMM_WORLD,mpierr)

     if (root_proc ) then

        ! after this call, the cnity and fvg is the dimension of (out_ntiles, nveg)
        call this%add_bcs_to_cnrst(surflay, OutBcsDir, __RC__)

        do i = 1, numprocs -1
           st  = low_ind(i+1)
           l   = nt_local(i+1)
           tag = i*numprocs
           if ((this%isCLM40) .or. (this%isCLM45)) then
              call MPI_send(this%cnity(st,1),l, MPI_REAL, i, tag, MPI_COMM_WORLD, mpierr)
              call MPI_send(this%cnity(st,2),l, MPI_REAL, i, tag+1, MPI_COMM_WORLD, mpierr)
              call MPI_send(this%cnity(st,3),l, MPI_REAL, i, tag+2, MPI_COMM_WORLD, mpierr)
              call MPI_send(this%cnity(st,4),l, MPI_REAL, i, tag+3, MPI_COMM_WORLD, mpierr)
              call MPI_send(this%fvg(st,1),l,   MPI_REAL, i, tag+4, MPI_COMM_WORLD, mpierr)
              call MPI_send(this%fvg(st,2),l,   MPI_REAL, i, tag+5, MPI_COMM_WORLD, mpierr)
              call MPI_send(this%fvg(st,3),l,   MPI_REAL, i, tag+6, MPI_COMM_WORLD, mpierr)
              call MPI_send(this%fvg(st,4),l,   MPI_REAL, i, tag+7, MPI_COMM_WORLD, mpierr)
            else if (this%isCLM51) then
              call MPI_send(this%cnity(st,1),l, MPI_REAL, i, tag, MPI_COMM_WORLD, mpierr)
              call MPI_send(this%cnity(st,2),l, MPI_REAL, i, tag+1, MPI_COMM_WORLD, mpierr)
              call MPI_send(this%fvg(st,1),l,   MPI_REAL, i, tag+2, MPI_COMM_WORLD, mpierr)
              call MPI_send(this%fvg(st,2),l,   MPI_REAL, i, tag+3, MPI_COMM_WORLD, mpierr)
            endif 
        enddo
        st  = low_ind(1)
        l   = nt_local(1)
        ed  = st + l -1
        if ((this%isCLM40) .or. (this%isCLM45)) then
           CLMC_pt1 = this%cnity(st:ed,1)
           CLMC_pt2 = this%cnity(st:ed,2)
           CLMC_st1 = this%cnity(st:ed,3)
           CLMC_st2 = this%cnity(st:ed,4)
           CLMC_pf1 = this%fvg(st:ed,1)
           CLMC_pf2 = this%fvg(st:ed,2)
           CLMC_sf1 = this%fvg(st:ed,3)
           CLMC_sf2 = this%fvg(st:ed,4)
        elseif (this%isCLM51) then
           CLMC_pt1 = this%cnity(st:ed,1)
           CLMC_st1 = this%cnity(st:ed,2)
           CLMC_pf1 = this%fvg(st:ed,1)
           CLMC_sf1 = this%fvg(st:ed,2)
        endif 
     else
        tag = myid*numprocs
        if ((this%isCLM40) .or. (this%isCLM45)) then
           call MPI_RECV(CLMC_pt1,nt_local(myid+1) , MPI_REAL, 0, tag,  MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
           call MPI_RECV(CLMC_pt2,nt_local(myid+1) , MPI_REAL, 0, tag+1, MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
           call MPI_RECV(CLMC_st1,nt_local(myid+1) , MPI_REAL, 0, tag+2, MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
           call MPI_RECV(CLMC_st2,nt_local(myid+1) , MPI_REAL, 0, tag+3, MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
           call MPI_RECV(CLMC_pf1,nt_local(myid+1) , MPI_REAL, 0, tag+4, MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
           call MPI_RECV(CLMC_pf2,nt_local(myid+1) , MPI_REAL, 0, tag+5, MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
           call MPI_RECV(CLMC_sf1,nt_local(myid+1) , MPI_REAL, 0, tag+6, MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
           call MPI_RECV(CLMC_sf2,nt_local(myid+1) , MPI_REAL, 0, tag+7, MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
        elseif (this%isCLM51) then
           call MPI_RECV(CLMC_pt1,nt_local(myid+1) , MPI_REAL, 0, tag,  MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
           call MPI_RECV(CLMC_st1,nt_local(myid+1) , MPI_REAL, 0, tag+1, MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
           call MPI_RECV(CLMC_pf1,nt_local(myid+1) , MPI_REAL, 0, tag+2, MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
           call MPI_RECV(CLMC_sf1,nt_local(myid+1) , MPI_REAL, 0, tag+3, MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
        end if 
    endif

    call MPI_Barrier(MPI_COMM_WORLD, STATUS)
 
    if(root_proc) print*, "GetIDs...."

    call GetIds(this%lonc,this%latc,this%lonn,this%latt,id_loc_cn, tid_offl, &
             CLMC_pf1, CLMC_pf2, CLMC_sf1, CLMC_sf2, CLMC_pt1, CLMC_pt2,CLMC_st1,CLMC_st2, &
             fveg_offl, ityp_offl,this%isCLM51)

    call MPI_Barrier(MPI_COMM_WORLD, STATUS)

    if(root_proc) allocate (id_glb_cn  (out_ntiles,nveg))

    allocate (id_loc (out_ntiles))
    deallocate (CLMC_pf1, CLMC_pf2, CLMC_sf1, CLMC_sf2)
    deallocate (CLMC_pt1, CLMC_pt2, CLMC_st1, CLMC_st2)
    deallocate (lat_tmp)

    do nv = 1, nveg
       call MPI_Barrier(MPI_COMM_WORLD, STATUS)
       do i = 1, numprocs
           if((I == 1).and.(myid == 0)) then
              id_loc(low_ind(i) : upp_ind(i)) = Id_loc_cn(:,nv)
           else if (I > 1) then
              if(I-1 == myid) then
                    ! send to root
                  call MPI_ISend(id_loc_cn(:,nv),nt_local(i),MPI_INTEGER,0,994,MPI_COMM_WORLD,req,mpierr)
                  call MPI_WAIT (req,MPI_STATUS_IGNORE,mpierr)
              else if (myid == 0) then
                    ! root receives
                  call MPI_RECV(id_loc(low_ind(i) : upp_ind(i)),nt_local(i) , MPI_INTEGER, i-1,994,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
              endif
           endif
       end do

       if(root_proc) id_glb_cn (:,nv) = id_loc

    end do

    if(root_proc) then

        allocate (var_off_col (1: in_ntiles, 1 : nzone,1 : var_col))
        allocate (var_off_pft (1: in_ntiles, 1 : nzone,1 : nveg, 1 : var_pft))
        allocate (var_out     (out_ntiles))
        allocate (var_tmp3d   (out_ntiles, nveg, nzone))
        allocate (var_out_zone(out_ntiles, nzone))

        this%tile_id = [(i*1.0, i=1, out_ntiles)]

        allocate (tg_tmp(out_ntiles, 4),source = 0.)
        do i = 1, 3
          tg_tmp(:,i) = this%tg(this%id_glb(:),i)
        enddo        
        this%tg = tg_tmp
        deallocate(tg_tmp)


        var_out = this%bflowm (this%id_glb(:))
        this%bflowm = var_out
        var_out = this%totwatm(this%id_glb(:))
        this%totwatm= var_out
        var_out = this%tairm  (this%id_glb(:))
        this%tairm  = var_out
        var_out = this%tpm    (this%id_glb(:))
        this%tpm    = var_out
        var_out = this%cnsum  (this%id_glb(:))
        this%cnsum  = var_out
        var_out = this%sndzm  (this%id_glb(:))
        this%sndzm  = var_out
        var_out = this%asnowm (this%id_glb(:))
        this%asnowm = var_out
        do nz = 1, nzone
           do nv = 1, nveg
               var_tmp3d(:,nv,nz) = this%psnsunm(this%id_glb(:), nv,nz)
           enddo
        enddo
        this%psnsunm= var_tmp3d

        do nz = 1, nzone
           do nv = 1, nveg
               var_tmp3d(:,nv,nz) = this%psnsham(this%id_glb(:), nv,nz)
           enddo
        enddo
        this%psnsham = var_tmp3d

        do nz = 1, nzone   
           do nv = 1, nveg 
               var_tmp3d(:,nv,nz) = this%lmrsunm(this%id_glb(:), nv,nz)
           enddo
        enddo
        this%lmrsunm= var_tmp3d

        do nz = 1, nzone
           do nv = 1, nveg
               var_tmp3d(:,nv,nz) = this%lmrsham(this%id_glb(:), nv,nz)
           enddo
        enddo
        this%lmrsham = var_tmp3d

        do nz = 1, nzone   
           do nv = 1, nveg 
               var_tmp3d(:,nv,nz) = this%laisunm(this%id_glb(:), nv,nz)
           enddo
        enddo
        this%laisunm= var_tmp3d

        do nz = 1, nzone
           do nv = 1, nveg
               var_tmp3d(:,nv,nz) = this%laisham(this%id_glb(:), nv,nz)
           enddo
        enddo
        this%laisham = var_tmp3d


        do nz = 1, nzone
           var_out_zone(:,nz) = this%rzmm(this%id_glb(:), nz)
        enddo
        this%rzmm = var_out_zone

        do nz = 1, nzone
           var_out_zone(:,nz) = this%tgwm(this%id_glb(:), nz)
        enddo
        this%tgwm = var_out_zone

        if (this%isCLM40) then
           var_out = this%sfmcm (this%id_glb(:))
           this%sfmcm = var_out
        endif
        if (this%isCLM45) then
           var_out = this%ar1m    (this%id_glb(:))
           this%ar1m    = var_out
           var_out = this%rainfm  (this%id_glb(:))
           this%rainfm  = var_out
           var_out = this%rhm     (this%id_glb(:))
           this%rhm     = var_out
           var_out = this%runsrfm (this%id_glb(:))
           this%runsrfm = var_out
           var_out = this%snowfm  (this%id_glb(:))
           this%snowfm  = var_out
           var_out = this%windm   (this%id_glb(:))
           this%windm   = var_out
           var_out = this%tprec10d(this%id_glb(:))
           this%tprec10d= var_out
           var_out = this%tprec60d(this%id_glb(:))
           this%tprec60d= var_out
           var_out = this%t2m10d  (this%id_glb(:))
           this%t2m10d  = var_out
           do nz = 1, nzone
              var_out_zone(:,nz) = this%sfmm(this%id_glb(:), nz)
           enddo
           this%sfmm = var_out_zone
        endif

        i = 1
        do nv = 1,VAR_COL
           do nz = 1,nzone
              var_off_col(:,nz,nv) = this%cncol(:,i)
              i = i + 1
           end do
        end do

        i = 1
        do iv = 1,VAR_PFT
           do nv = 1,nveg
              do nz = 1,nzone
                 var_off_pft(:, nz,nv,iv) = this%cnpft(:,i)
                 i = i + 1
              end do
           end do
        end do

        where(isnan(var_off_pft))  var_off_pft = 0.
        where(var_off_pft /= var_off_pft)  var_off_pft = 0.

        print *, 'calculating regridded carbn'

        if (this%isCLM40) then
           allocate(iclass(1:npft))
           iclass = iclass_40
        elseif (this%isCLM45) then 
           allocate(iclass(1:npft))
           iclass = iclass_45
        elseif (this%isCLM51) then 
           allocate(iclass(1:npft_51))
           iclass = iclass_51
        end if
        
       

        call regrid_carbon (out_NTILES, in_ntiles,id_glb_cn, &
                DAYX, var_off_col,var_off_pft, ityp_offl, fveg_offl, iclass)
        deallocate (var_off_col,var_off_pft)
     endif
     call MPI_Barrier(MPI_COMM_WORLD, STATUS)

    _RETURN(_SUCCESS)

  contains
     SUBROUTINE regrid_carbon (NTILES, in_ntiles, id_glb, &
        DAYX, var_off_col, var_off_pft, ityp_offl, fveg_offl,iclass_in)

     ! write out regridded carbon variables
     implicit none
     integer, intent (in) :: NTILES, in_ntiles,id_glb (ntiles,nveg)
     real, intent (in)    :: DAYX (NTILES), var_off_col(in_ntiles,NZONE,var_col), var_off_pft(in_ntiles,NZONE, NVEG, var_pft)
     real, intent (in),  dimension(in_ntiles,nveg) :: fveg_offl,  ityp_offl
     integer, intent(in), dimension(:) :: iclass_in
     real, allocatable, dimension (:)    :: CLMC_pf1, CLMC_pf2, CLMC_sf1, CLMC_sf2, &
          CLMC_pt1, CLMC_pt2,CLMC_st1,CLMC_st2, var_dum
     real, allocatable :: var_col_out (:,:,:), var_pft_out (:,:,:,:)
     integer           :: N, STATUS, nv, nx, offl_cell, ityp_new, i, j, nz, iv
     real              :: fveg_new
     character(256) :: Iam = "write_regridded_carbon"


     allocate (CLMC_pf1(NTILES))
     allocate (CLMC_pf2(NTILES))
     allocate (CLMC_sf1(NTILES))
     allocate (CLMC_sf2(NTILES))
     allocate (CLMC_pt1(NTILES))
     allocate (CLMC_pt2(NTILES))
     allocate (CLMC_st1(NTILES))
     allocate (CLMC_st2(NTILES))
     allocate (VAR_DUM (NTILES))

     if ((this%isCLM40).or.(this%isCLM45)) then
        CLMC_pt1 = this%cnity(:,1)
        CLMC_pt2 = this%cnity(:,2)
        CLMC_st1 = this%cnity(:,3)
        CLMC_st2 = this%cnity(:,4)
        CLMC_pf1 = this%fvg(:,1)
        CLMC_pf2 = this%fvg(:,2)
        CLMC_sf1 = this%fvg(:,3)
        CLMC_sf2 = this%fvg(:,4)

     elseif (this%isCLM51) then

        CLMC_pt1 = this%cnity(:,1)
        CLMC_st1 = this%cnity(:,3)
        CLMC_st2 = this%cnity(:,4)
        CLMC_pf1 = this%fvg(:,1)
        CLMC_pf2 = this%fvg(:,2)
        CLMC_sf1 = this%fvg(:,3)
        CLMC_sf2 = this%fvg(:,4)

     end if 

     allocate (var_col_out (1: NTILES, 1 : nzone,1 : var_col))
     allocate (var_pft_out (1: NTILES, 1 : nzone,1 : nveg, 1 : var_pft))

     var_col_out = 0.
     var_pft_out = NaN

     OUT_TILE : DO N = 1, NTILES

        ! if(mod (n,1000) == 0) print *, myid +1, n, Id_glb(n,:)

        NVLOOP2 : do nv = 1, nveg

           if ((this%isCLM40).or.(this%isCLM45)) then
             
              if(nv <= 2) then ! index for secondary PFT index if primary or primary if secondary
                 nx = nv + 2
              else
                 nx = nv - 2
              endif

              if (nv == 1) ityp_new = CLMC_pt1(n)
              if (nv == 1) fveg_new = CLMC_pf1(n)
              if (nv == 2) ityp_new = CLMC_pt2(n)
              if (nv == 2) fveg_new = CLMC_pf2(n)
              if (nv == 3) ityp_new = CLMC_st1(n)
              if (nv == 3) fveg_new = CLMC_sf1(n)
              if (nv == 4) ityp_new = CLMC_st2(n)
              if (nv == 4) fveg_new = CLMC_sf2(n)
 
           elseif (this%isCLM51) then
             
              if(nv <= 1) then ! index for secondary PFT index if primary or primary if secondary
                 nx = nv + 1
              else
                 nx = nv - 1
              endif

              if (nv == 1) ityp_new = CLMC_pt1(n)
              if (nv == 1) fveg_new = CLMC_pf1(n)
              if (nv == 2) ityp_new = CLMC_st1(n)
              if (nv == 2) fveg_new = CLMC_sf1(n) 

           end if

           if (fveg_new > fmin) then

              offl_cell    = Id_glb(n,nv)

              if(ityp_new      == ityp_offl (offl_cell,nv) .and. fveg_offl (offl_cell,nv)> fmin) then
                 iv = nv                                     ! same type fraction (primary of secondary)                          
              else if(ityp_new == ityp_offl (offl_cell,nx) .and. fveg_offl (offl_cell,nx)> fmin) then
                 iv = nx                                     ! not same fraction
              else if(iclass_in(ityp_new)==iclass_in(ityp_offl(offl_cell,nv)) .and. fveg_offl (offl_cell,nv)> fmin) then
                 iv = nv                                     ! primary, other type (same class)
              else if(fveg_offl (offl_cell,nx)> fmin) then
                 iv = nx                                     ! secondary, other type (same class)
              endif

              ! Get col and pft variables for the Id_glb(nv) grid cell from offline catchcn_internal_rst
              ! ----------------------------------------------------------------------------------------

              ! call NCDF_reshape_getOput (NCFID,Id_glb(n,nv),var_off_col,var_off_pft,.true.)  

              var_pft_out (n,:,nv,:) = var_off_pft(Id_glb(n,nv), :,iv,:)
              var_col_out (n,:,:)    = var_col_out(n,:,:) + fveg_new * var_off_col(Id_glb(n,nv), :,:) ! gkw: column state simple weighted mean; ! could use "woody" fraction?

              ! Check whether var_pft_out is realistic
              do nz = 1, nzone
                 do j = 1, VAR_PFT
                    if (isnan(var_pft_out (n, nz,nv,j))) print *,j,nv,nz,n,var_pft_out (n, nz,nv,j),fveg_new
                    !if(isnan(var_pft_out (n, nz,nv,69))) var_pft_out (n, nz,nv,69) = 1.e-6
                    !if(isnan(var_pft_out (n, nz,nv,70))) var_pft_out (n, nz,nv,70) = 1.e-6   
                    !if(isnan(var_pft_out (n, nz,nv,73))) var_pft_out (n, nz,nv,73) = 1.e-6
                    !if(isnan(var_pft_out (n, nz,nv,74))) var_pft_out (n, nz,nv,74) = 1.e-6                 
                 end do
              end do
           endif

        end do NVLOOP2

        ! reset carbon if negative < 10g
        ! ------------------------

        NZLOOP : do nz = 1, nzone

           if(var_col_out (n, nz,14) < 10.) then

              var_col_out(n, nz, 1) = max(var_col_out(n, nz, 1), 0.)
              var_col_out(n, nz, 2) = max(var_col_out(n, nz, 2), 0.)
              var_col_out(n, nz, 3) = max(var_col_out(n, nz, 3), 0.)
              var_col_out(n, nz, 4) = max(var_col_out(n, nz, 4), 0.)
              var_col_out(n, nz, 5) = max(var_col_out(n, nz, 5), 0.)
              var_col_out(n, nz,10) = max(var_col_out(n, nz,10), 0.)
              var_col_out(n, nz,11) = max(var_col_out(n, nz,11), 0.)
              var_col_out(n, nz,12) = max(var_col_out(n, nz,12), 0.)
              var_col_out(n, nz,13) = max(var_col_out(n, nz,13),10.)   ! soil4c       
              var_col_out(n, nz,14) = max(var_col_out(n, nz,14), 0.)
              var_col_out(n, nz,15) = max(var_col_out(n, nz,15), 0.)
              var_col_out(n, nz,16) = max(var_col_out(n, nz,16), 0.)
              var_col_out(n, nz,17) = max(var_col_out(n, nz,17), 0.)
              var_col_out(n, nz,18) = max(var_col_out(n, nz,18), 0.)
              var_col_out(n, nz,19) = max(var_col_out(n, nz,19), 0.)
              var_col_out(n, nz,20) = max(var_col_out(n, nz,20), 0.)
              var_col_out(n, nz,24) = max(var_col_out(n, nz,24), 0.)
              var_col_out(n, nz,25) = max(var_col_out(n, nz,25), 0.)
              var_col_out(n, nz,26) = max(var_col_out(n, nz,26), 0.)
              var_col_out(n, nz,27) = max(var_col_out(n, nz,27), 0.)
              var_col_out(n, nz,28) = max(var_col_out(n, nz,28), 1.)
              var_col_out(n, nz,29) = max(var_col_out(n, nz,29), 0.)

              NVLOOP3 : do nv = 1,nveg
                 if ((this%isCLM40).or.(this%isCLM45)) then
                    if (nv == 1) ityp_new = CLMC_pt1(n)
                    if (nv == 1) fveg_new = CLMC_pf1(n)
                    if (nv == 2) ityp_new = CLMC_pt2(n)
                    if (nv == 2) fveg_new = CLMC_pf2(n)
                    if (nv == 3) ityp_new = CLMC_st1(n)
                    if (nv == 3) fveg_new = CLMC_sf1(n)
                    if (nv == 4) ityp_new = CLMC_st2(n)
                    if (nv == 4) fveg_new = CLMC_sf2(n)
                 elseif (this%isCLM51) then
                    if (nv == 1) ityp_new = CLMC_pt1(n)
                    if (nv == 1) fveg_new = CLMC_pf1(n)
                    if (nv == 2) ityp_new = CLMC_st1(n)
                    if (nv == 2) fveg_new = CLMC_sf1(n)
                 end if

                 if(fveg_new > fmin) then
                    var_pft_out(n, nz,nv, 1) = max(var_pft_out(n, nz,nv, 1),0.)
                    var_pft_out(n, nz,nv, 2) = max(var_pft_out(n, nz,nv, 2),0.)
                    var_pft_out(n, nz,nv, 3) = max(var_pft_out(n, nz,nv, 3),0.)
                    var_pft_out(n, nz,nv, 4) = max(var_pft_out(n, nz,nv, 4),0.)

                    if(ityp_new <= 12) then ! tree or shrub deadstemc
                       var_pft_out(n, nz,nv, 5) = max(var_pft_out(n, nz,nv, 5),0.1)
                    else
                       var_pft_out(n, nz,nv, 5) = max(var_pft_out(n, nz,nv, 5),0.0)
                    endif

                    var_pft_out(n, nz,nv, 6) = max(var_pft_out(n, nz,nv, 6),0.)
                    var_pft_out(n, nz,nv, 7) = max(var_pft_out(n, nz,nv, 7),0.)
                    var_pft_out(n, nz,nv, 8) = max(var_pft_out(n, nz,nv, 8),0.)
                    var_pft_out(n, nz,nv, 9) = max(var_pft_out(n, nz,nv, 9),0.)
                    var_pft_out(n, nz,nv,10) = max(var_pft_out(n, nz,nv,10),0.)
                    var_pft_out(n, nz,nv,11) = max(var_pft_out(n, nz,nv,11),0.)
                    var_pft_out(n, nz,nv,12) = max(var_pft_out(n, nz,nv,12),0.)

                    if(ityp_new <=2 .or. ityp_new ==4 .or. ityp_new ==5 .or. ityp_new == 9) then
                       var_pft_out(n, nz,nv,13) = max(var_pft_out(n, nz,nv,13),1.)  ! leaf carbon display for evergreen
                       var_pft_out(n, nz,nv,14) = max(var_pft_out(n, nz,nv,14),0.)
                    else
                       var_pft_out(n, nz,nv,13) = max(var_pft_out(n, nz,nv,13),0.)
                       var_pft_out(n, nz,nv,14) = max(var_pft_out(n, nz,nv,14),1.)  ! leaf carbon storage for deciduous
                    endif

                    var_pft_out(n, nz,nv,15) = max(var_pft_out(n, nz,nv,15),0.)
                    var_pft_out(n, nz,nv,16) = max(var_pft_out(n, nz,nv,16),0.)
                    var_pft_out(n, nz,nv,17) = max(var_pft_out(n, nz,nv,17),0.)
                    var_pft_out(n, nz,nv,18) = max(var_pft_out(n, nz,nv,18),0.)
                    var_pft_out(n, nz,nv,19) = max(var_pft_out(n, nz,nv,19),0.)
                    var_pft_out(n, nz,nv,20) = max(var_pft_out(n, nz,nv,20),0.)
                    var_pft_out(n, nz,nv,21) = max(var_pft_out(n, nz,nv,21),0.)
                    var_pft_out(n, nz,nv,22) = max(var_pft_out(n, nz,nv,22),0.)
                    var_pft_out(n, nz,nv,23) = max(var_pft_out(n, nz,nv,23),0.)
                    var_pft_out(n, nz,nv,25) = max(var_pft_out(n, nz,nv,25),0.)
                    var_pft_out(n, nz,nv,26) = max(var_pft_out(n, nz,nv,26),0.)
                    var_pft_out(n, nz,nv,27) = max(var_pft_out(n, nz,nv,27),0.)
                    var_pft_out(n, nz,nv,41) = max(var_pft_out(n, nz,nv,41),0.)
                    var_pft_out(n, nz,nv,42) = max(var_pft_out(n, nz,nv,42),0.)
                    var_pft_out(n, nz,nv,44) = max(var_pft_out(n, nz,nv,44),0.)
                    var_pft_out(n, nz,nv,45) = max(var_pft_out(n, nz,nv,45),0.)
                    var_pft_out(n, nz,nv,46) = max(var_pft_out(n, nz,nv,46),0.)
                    var_pft_out(n, nz,nv,47) = max(var_pft_out(n, nz,nv,47),0.)
                    var_pft_out(n, nz,nv,48) = max(var_pft_out(n, nz,nv,48),0.)
                    var_pft_out(n, nz,nv,49) = max(var_pft_out(n, nz,nv,49),0.)
                    var_pft_out(n, nz,nv,50) = max(var_pft_out(n, nz,nv,50),0.)
                    var_pft_out(n, nz,nv,51) = max(var_pft_out(n, nz,nv, 5)/500.,0.)
                    var_pft_out(n, nz,nv,52) = max(var_pft_out(n, nz,nv,52),0.)
                    var_pft_out(n, nz,nv,53) = max(var_pft_out(n, nz,nv,53),0.)
                    var_pft_out(n, nz,nv,54) = max(var_pft_out(n, nz,nv,54),0.)
                    var_pft_out(n, nz,nv,55) = max(var_pft_out(n, nz,nv,55),0.)
                    var_pft_out(n, nz,nv,56) = max(var_pft_out(n, nz,nv,56),0.)
                    var_pft_out(n, nz,nv,57) = max(var_pft_out(n, nz,nv,13)/25.,0.)
                    var_pft_out(n, nz,nv,58) = max(var_pft_out(n, nz,nv,14)/25.,0.)
                    var_pft_out(n, nz,nv,59) = max(var_pft_out(n, nz,nv,59),0.)
                    var_pft_out(n, nz,nv,60) = max(var_pft_out(n, nz,nv,60),0.)
                    var_pft_out(n, nz,nv,61) = max(var_pft_out(n, nz,nv,61),0.)
                    var_pft_out(n, nz,nv,62) = max(var_pft_out(n, nz,nv,62),0.)
                    var_pft_out(n, nz,nv,63) = max(var_pft_out(n, nz,nv,63),0.)
                    var_pft_out(n, nz,nv,64) = max(var_pft_out(n, nz,nv,64),0.)
                    var_pft_out(n, nz,nv,65) = max(var_pft_out(n, nz,nv,65),0.)
                    var_pft_out(n, nz,nv,66) = max(var_pft_out(n, nz,nv,66),0.)
                    var_pft_out(n, nz,nv,67) = max(var_pft_out(n, nz,nv,67),0.)
                    var_pft_out(n, nz,nv,68) = max(var_pft_out(n, nz,nv,68),0.)
                    var_pft_out(n, nz,nv,69) = max(var_pft_out(n, nz,nv,69),0.)
                    var_pft_out(n, nz,nv,70) = max(var_pft_out(n, nz,nv,70),0.)
                    var_pft_out(n, nz,nv,73) = max(var_pft_out(n, nz,nv,73),0.)
                    var_pft_out(n, nz,nv,74) = max(var_pft_out(n, nz,nv,74),0.)
                    if(this%isCLM45) var_pft_out(n, nz,nv,75) = max(var_pft_out(n, nz,nv,75),0.)
                    if(this%isCLM51) then
                       var_pft_out(n, nz,nv,76) = max(var_pft_out(n, nz,nv,76),0.)
                       var_pft_out(n, nz,nv,77) = max(var_pft_out(n, nz,nv,77),0.)
                       var_pft_out(n, nz,nv,78) = max(var_pft_out(n, nz,nv,78),0.)
                       var_pft_out(n, nz,nv,79) = max(var_pft_out(n, nz,nv,79),0.)
                       var_pft_out(n, nz,nv,80) = max(var_pft_out(n, nz,nv,80),0.)
                       var_pft_out(n, nz,nv,81) = max(var_pft_out(n, nz,nv,81),0.)
                       var_pft_out(n, nz,nv,82) = max(var_pft_out(n, nz,nv,82),0.)
                       var_pft_out(n, nz,nv,83) = max(var_pft_out(n, nz,nv,83),0.)
                    end if
                 endif
              end do NVLOOP3  ! end veg loop                 
           endif    ! end carbon check         
        end do NZLOOP ! end zone loop

        ! Update dayx variable var_pft_out (:,:,28)

        do j = 28, 28  !  1,VAR_PFT var_pft_out (:,:,:,28)
           do nv = 1,nveg
              do nz = 1,nzone
                 var_pft_out (n, nz,nv,j) = dayx(n)
              end do
           end do
        end do

        ! call NCDF_reshape_getOput (OutID,N,var_col_out,var_pft_out,.false.)  

        ! column vars clm40                         clm45
        ! -----------------                         ---------------------
        !  1 clm3%g%l%c%ccs%col_ctrunc            !  1 ccs%col_ctrunc_vr   (:,1)              
        !  2 clm3%g%l%c%ccs%cwdc            !  2 ccs%decomp_cpools_vr(:,1,4)   ! cwdc        
        !  3 clm3%g%l%c%ccs%litr1c                 !  3 ccs%decomp_cpools_vr(:,1,1)   ! litr1c 
        !  4 clm3%g%l%c%ccs%litr2c             !  4 ccs%decomp_cpools_vr(:,1,2)   ! litr2c 
        !  5 clm3%g%l%c%ccs%litr3c             !  5 ccs%decomp_cpools_vr(:,1,3)   ! litr3c 
        !  6 clm3%g%l%c%ccs%pcs_a%totvegc          !  6 ccs%totvegc_col                        
        !  7 clm3%g%l%c%ccs%prod100c              !  7 ccs%prod100c                           
        !  8 clm3%g%l%c%ccs%prod10c            !  8 ccs%prod10c                            
        !  9 clm3%g%l%c%ccs%seedc                  !  9 ccs%seedc                              
        ! 10 clm3%g%l%c%ccs%soil1c                 ! 10 ccs%decomp_cpools_vr(:,1,5)   ! soil1c 
        ! 11 clm3%g%l%c%ccs%soil2c                 ! 11 ccs%decomp_cpools_vr(:,1,6)   ! soil2c 
        ! 12 clm3%g%l%c%ccs%soil3c                 ! 12 ccs%decomp_cpools_vr(:,1,7)   ! soil3c 
        ! 13 clm3%g%l%c%ccs%soil4c                 ! 13 ccs%decomp_cpools_vr(:,1,8)   ! soil4c 
        ! 14 clm3%g%l%c%ccs%totcolc                ! 14 ccs%totcolc                            
        ! 15 clm3%g%l%c%ccs%totlitc                ! 15 ccs%totlitc                            
        ! 16 clm3%g%l%c%cns%col_ntrunc            ! 16 cns%col_ntrunc_vr   (:,1)              
        ! 17 clm3%g%l%c%cns%cwdn            ! 17 cns%decomp_npools_vr(:,1,4)   ! cwdn        
        ! 18 clm3%g%l%c%cns%litr1n                 ! 18 cns%decomp_npools_vr(:,1,1)   ! litr1n 
        ! 19 clm3%g%l%c%cns%litr2n             ! 19 cns%decomp_npools_vr(:,1,2)   ! litr2n 
        ! 20 clm3%g%l%c%cns%litr3n             ! 20 cns%decomp_npools_vr(:,1,3)   ! litr3n 
        ! 21 clm3%g%l%c%cns%prod100n              ! 21 cns%prod100n                           
        ! 22 clm3%g%l%c%cns%prod10n                ! 22 cns%prod10n                            
        ! 23 clm3%g%l%c%cns%seedn                  ! 23 cns%seedn                              
        ! 24 clm3%g%l%c%cns%sminn                  ! 24 cns%sminn_vr        (:,1)              
        ! 25 clm3%g%l%c%cns%soil1n                 ! 25 cns%decomp_npools_vr(:,1,5)   ! soil1n 
        ! 26 clm3%g%l%c%cns%soil2n                 ! 26 cns%decomp_npools_vr(:,1,6)   ! soil2n 
        ! 27 clm3%g%l%c%cns%soil3n                 ! 27 cns%decomp_npools_vr(:,1,7)   ! soil3n 
        ! 28 clm3%g%l%c%cns%soil4n                 ! 28 cns%decomp_npools_vr(:,1,8)   ! soil4n 
        ! 29 clm3%g%l%c%cns%totcoln                ! 29 cns%totcoln                            
        ! 30 clm3%g%l%c%cps%ann_farea_burned       ! 30 cps%fpg                                
        ! 31 clm3%g%l%c%cps%annsum_counter         ! 31 cps%annsum_counter                     
        ! 32 clm3%g%l%c%cps%cannavg_t2m              ! 32 cps%cannavg_t2m                             
        ! 33 clm3%g%l%c%cps%cannsum_npp              ! 33 cps%cannsum_npp                             
        ! 34 clm3%g%l%c%cps%farea_burned           ! 34 cps%farea_burned                            
        ! 35 clm3%g%l%c%cps%fire_prob             ! 35 cps%fpi_vr          (:,1)              
        ! 36 clm3%g%l%c%cps%fireseasonl              ! OLD           ! 30 cps%altmax                   
        ! 37 clm3%g%l%c%cps%fpg                     ! OLD           ! 31 cps%annsum_counter            
        ! 38 clm3%g%l%c%cps%fpi                     ! OLD           ! 32 cps%cannavg_t2m               
        ! 39 clm3%g%l%c%cps%me                      ! OLD           ! 33 cps%cannsum_npp          
        ! 40 clm3%g%l%c%cps%mean_fire_prob         ! OLD           ! 34 cps%farea_burned         
                                                   ! OLD           ! 35 cps%altmax_lastyear      
                                                   ! OLD           ! 36 cps%altmax_indx          
                                                   ! OLD           ! 37 cps%fpg                  
                                                   ! OLD           ! 38 cps%fpi_vr          (:,1)
                                                   ! OLD           ! 39 cps%altmax_lastyear_indx 

        ! PFT vars CLM40                                CLM45    
        ! --------------                                -----
        !  1 clm3%g%l%c%p%pcs%cpool                    !  1 pcs%cpool                           
        !  2 clm3%g%l%c%p%pcs%deadcrootc             !  2 pcs%deadcrootc                      
        !  3 clm3%g%l%c%p%pcs%deadcrootc_storage      !  3 pcs%deadcrootc_storage              
        !  4 clm3%g%l%c%p%pcs%deadcrootc_xfer              !  4 pcs%deadcrootc_xfer                 
        !  5 clm3%g%l%c%p%pcs%deadstemc              !  5 pcs%deadstemc                       
        !  6 clm3%g%l%c%p%pcs%deadstemc_storage       !  6 pcs%deadstemc_storage               
        !  7 clm3%g%l%c%p%pcs%deadstemc_xfer               !  7 pcs%deadstemc_xfer                  
        !  8 clm3%g%l%c%p%pcs%frootc                      !  8 pcs%frootc                          
        !  9 clm3%g%l%c%p%pcs%frootc_storage               !  9 pcs%frootc_storage                  
        ! 10 clm3%g%l%c%p%pcs%frootc_xfer                  ! 10 pcs%frootc_xfer                     
        ! 11 clm3%g%l%c%p%pcs%gresp_storage                ! 11 pcs%gresp_storage                   
        ! 12 clm3%g%l%c%p%pcs%gresp_xfer             ! 12 pcs%gresp_xfer                      
        ! 13 clm3%g%l%c%p%pcs%leafc                    ! 13 pcs%leafc                           
        ! 14 clm3%g%l%c%p%pcs%leafc_storage             ! 14 pcs%leafc_storage                   
        ! 15 clm3%g%l%c%p%pcs%leafc_xfer             ! 15 pcs%leafc_xfer                      
        ! 16 clm3%g%l%c%p%pcs%livecrootc        ! 16 pcs%livecrootc                      
        ! 17 clm3%g%l%c%p%pcs%livecrootc_storage      ! 17 pcs%livecrootc_storage              
        ! 18 clm3%g%l%c%p%pcs%livecrootc_xfer              ! 18 pcs%livecrootc_xfer                 
        ! 19 clm3%g%l%c%p%pcs%livestemc              ! 19 pcs%livestemc                       
        ! 20 clm3%g%l%c%p%pcs%livestemc_storage       ! 20 pcs%livestemc_storage               
        ! 21 clm3%g%l%c%p%pcs%livestemc_xfer               ! 21 pcs%livestemc_xfer                  
        ! 22 clm3%g%l%c%p%pcs%pft_ctrunc             ! 22 pcs%pft_ctrunc                      
        ! 23 clm3%g%l%c%p%pcs%xsmrpool                     ! 23 pcs%xsmrpool                        
        ! 24 clm3%g%l%c%p%pepv%annavg_t2m                  ! 24 pepv%annavg_t2m                     
        ! 25 clm3%g%l%c%p%pepv%annmax_retransn             ! 25 pepv%annmax_retransn                
        ! 26 clm3%g%l%c%p%pepv%annsum_npp                  ! 26 pepv%annsum_npp                     
        ! 27 clm3%g%l%c%p%pepv%annsum_potential_gpp        ! 27 pepv%annsum_potential_gpp           
        ! 28 clm3%g%l%c%p%pepv%dayl                    ! 28 pepv%dayl                           
        ! 29 clm3%g%l%c%p%pepv%days_active                 ! 29 pepv%days_active                    
        ! 30 clm3%g%l%c%p%pepv%dormant_flag                ! 30 pepv%dormant_flag                   
        ! 31 clm3%g%l%c%p%pepv%offset_counter              ! 31 pepv%offset_counter                 
        ! 32 clm3%g%l%c%p%pepv%offset_fdd                  ! 32 pepv%offset_fdd                     
        ! 33 clm3%g%l%c%p%pepv%offset_flag                 ! 33 pepv%offset_flag                    
        ! 34 clm3%g%l%c%p%pepv%offset_swi                  ! 34 pepv%offset_swi                     
        ! 35 clm3%g%l%c%p%pepv%onset_counter               ! 35 pepv%onset_counter                  
        ! 36 clm3%g%l%c%p%pepv%onset_fdd             ! 36 pepv%onset_fdd      
        ! 37 clm3%g%l%c%p%pepv%onset_flag                  ! 37 pepv%onset_flag                     
        ! 38 clm3%g%l%c%p%pepv%onset_gdd             ! 38 pepv%onset_gdd                      
        ! 39 clm3%g%l%c%p%pepv%onset_gddflag               ! 39 pepv%onset_gddflag                  
        ! 40 clm3%g%l%c%p%pepv%onset_swi             ! 40 pepv%onset_swi                      
        ! 41 clm3%g%l%c%p%pepv%prev_frootc_to_litter       ! 41 pepv%prev_frootc_to_litter          
        ! 42 clm3%g%l%c%p%pepv%prev_leafc_to_litter        ! 42 pepv%prev_leafc_to_litter           
        ! 43 clm3%g%l%c%p%pepv%tempavg_t2m                 ! 43 pepv%tempavg_t2m                    
        ! 44 clm3%g%l%c%p%pepv%tempmax_retransn       ! 44 pepv%tempmax_retransn               
        ! 45 clm3%g%l%c%p%pepv%tempsum_npp                 ! 45 pepv%tempsum_npp                    
        ! 46 clm3%g%l%c%p%pepv%tempsum_potential_gpp       ! 46 pepv%tempsum_potential_gpp          
        ! 47 clm3%g%l%c%p%pepv%xsmrpool_recover       ! 47 pepv%xsmrpool_recover               
        ! 48 clm3%g%l%c%p%pns%deadcrootn             ! 48 pns%deadcrootn                      
        ! 49 clm3%g%l%c%p%pns%deadcrootn_storage      ! 49 pns%deadcrootn_storage              
        ! 50 clm3%g%l%c%p%pns%deadcrootn_xfer              ! 50 pns%deadcrootn_xfer                 
        ! 51 clm3%g%l%c%p%pns%deadstemn              ! 51 pns%deadstemn                       
        ! 52 clm3%g%l%c%p%pns%deadstemn_storage       ! 52 pns%deadstemn_storage               
        ! 53 clm3%g%l%c%p%pns%deadstemn_xfer               ! 53 pns%deadstemn_xfer                  
        ! 54 clm3%g%l%c%p%pns%frootn                      ! 54 pns%frootn                          
        ! 55 clm3%g%l%c%p%pns%frootn_storage               ! 55 pns%frootn_storage                  
        ! 56 clm3%g%l%c%p%pns%frootn_xfer                  ! 56 pns%frootn_xfer                     
        ! 57 clm3%g%l%c%p%pns%leafn                    ! 57 pns%leafn                           
        ! 58 clm3%g%l%c%p%pns%leafn_storage                ! 58 pns%leafn_storage                   
        ! 59 clm3%g%l%c%p%pns%leafn_xfer             ! 59 pns%leafn_xfer                      
        ! 60 clm3%g%l%c%p%pns%livecrootn             ! 60 pns%livecrootn                      
        ! 61 clm3%g%l%c%p%pns%livecrootn_storage      ! 61 pns%livecrootn_storage              
        ! 62 clm3%g%l%c%p%pns%livecrootn_xfer              ! 62 pns%livecrootn_xfer                 
        ! 63 clm3%g%l%c%p%pns%livestemn              ! 63 pns%livestemn                       
        ! 64 clm3%g%l%c%p%pns%livestemn_storage       ! 64 pns%livestemn_storage               
        ! 65 clm3%g%l%c%p%pns%livestemn_xfer               ! 65 pns%livestemn_xfer                  
        ! 66 clm3%g%l%c%p%pns%npool                    ! 66 pns%npool                           
        ! 67 clm3%g%l%c%p%pns%pft_ntrunc        ! 67 pns%pft_ntrunc                      
        ! 68 clm3%g%l%c%p%pns%retransn                    ! 68 pns%retransn                        
        ! 69 clm3%g%l%c%p%pps%elai                     ! 69 pps%elai                            
        ! 70 clm3%g%l%c%p%pps%esai                     ! 70 pps%esai                            
        ! 71 clm3%g%l%c%p%pps%hbot                     ! 71 pps%hbot                            
        ! 72 clm3%g%l%c%p%pps%htop                     ! 72 pps%htop                            
        ! 73 clm3%g%l%c%p%pps%tlai                     ! 73 pps%tlai                            
        ! 74 clm3%g%l%c%p%pps%tsai                     ! 74 pps%tsai                            
                                                           ! 75 pepv%plant_ndemand                  
                                                           ! OLD           ! 75 pps%gddplant        
                                                           ! OLD           ! 76 pps%gddtsoi         
                                                           ! OLD           ! 77 pps%peaklai         
                                                           ! OLD           ! 78 pps%idop            
                                                           ! OLD           ! 79 pps%aleaf           
                                                           ! OLD           ! 80 pps%aleafi          
                                                           ! OLD           ! 81 pps%astem           
                                                           ! OLD           ! 82 pps%astemi          
                                                           ! OLD           ! 83 pps%htmx            
                                                           ! OLD           ! 84 pps%hdidx           
                                                           ! OLD           ! 85 pps%vf              
                                                           ! OLD           ! 86 pps%cumvd           
                                                           ! OLD           ! 87 pps%croplive        
                                                           ! OLD           ! 88 pps%cropplant       
                                                           ! OLD           ! 89 pps%harvdate        
                                                           ! OLD           ! 90 pps%gdd1020         
                                                           ! OLD           ! 91 pps%gdd820          
                                                           ! OLD           ! 92 pps%gdd020          
                                                           ! OLD           ! 93 pps%gddmaturity     
                                                           ! OLD           ! 94 pps%huileaf         
                                                           ! OLD           ! 95 pps%huigrain        
                                                           ! OLD           ! 96 pcs%grainc          
                                                           ! OLD           ! 97 pcs%grainc_storage  
                                                           ! OLD           ! 98 pcs%grainc_xfer     
                                                           ! OLD           ! 99 pns%grainn          
                                                           ! OLD           !100 pns%grainn_storage  
                                                           ! OLD           !101 pns%grainn_xfer     
                                                           ! OLD           !102 pepv%fert_counter   
                                                           ! OLD           !103 pnf%fert            
                                                           ! OLD           !104 pepv%grain_flag     

     end do OUT_TILE

     i = 1
     deallocate(this%cncol)
     allocate(this%cncol(NTILES, nzone*VAR_COL))
     do nv = 1,VAR_COL
        do nz = 1,nzone
           this%cncol(:,i) = var_col_out(:, nz,nv)
           !STATUS = NF_PUT_VARA_REAL(OutID,VarID(OutID,'CNCOL'), (/1,i/), (/NTILES,1 /),var_col_out(:, nz,nv))  ; VERIFY_(STATUS)
           i = i + 1
        end do
     end do

     i = 1
     deallocate(this%cnpft)
     allocate(this%cnpft(NTILES,VAR_PFT*nveg*nzone))
     if(this%isclm45) then
        do iv = 1,VAR_PFT
           do nv = 1,nveg
              do nz = 1,nzone
                 if(iv <= 74) then
                    this%cnpft(:,i) = var_pft_out(:, nz,nv,iv)
                    !STATUS = NF_PUT_VARA_REAL(OutID,VarID(OutID,'CNPFT'), (/1,i/), (/NTILES,1 /),var_pft_out(:, nz,nv,iv))  ; VERIFY_(STATUS)
                 else
                    if((iv == 78) .OR. (iv == 89)) then    ! idop and harvdate
                       var_dum = 999
                       this%cnpft(:,i) = var_dum
                       !STATUS = NF_PUT_VARA_REAL(OutID,VarID(OutID,'CNPFT'), (/1,i/), (/NTILES,1 /),var_dum)  ; VERIFY_(STATUS)
                    else
                       var_dum = 0.
                       this%cnpft(:,i) = var_dum
                       !STATUS = NF_PUT_VARA_REAL(OutID,VarID(OutID,'CNPFT'), (/1,i/), (/NTILES,1 /),var_dum)  ; VERIFY_(STATUS)
                    endif
                 endif
                 i = i + 1
              end do
           end do
        end do
     elseif(this%isclm51) then
        do iv = 1,VAR_PFT
           do nv = 1,nveg
              do nz = 1,nzone
                 this%cnpft(:,i) = var_pft_out(:, nz,nv,iv)
                    !STATUS = NF_PUT_VARA_REAL(OutID,VarID(OutID,'CNPFT'), (/1,i/), (/NTILES,1 /),var_pft_out(:, nz,nv,iv))  ; VERIFY_(STATUS)
                 i = i + 1
              end do
           end do
        end do
     else
        do iv = 1,VAR_PFT
           do nv = 1,nveg
              do nz = 1,nzone
                 this%cnpft(:,i) = var_pft_out(:, nz,nv,iv)
                 !STATUS = NF_PUT_VARA_REAL(OutID,VarID(OutID,'CNPFT'), (/1,i/), (/NTILES,1 /),var_pft_out(:, nz,nv,iv))  ; VERIFY_(STATUS)
                 i = i + 1
              end do
           end do
        end do
     endif

     deallocate (var_col_out,var_pft_out)
     deallocate (CLMC_pf1, CLMC_pf2, CLMC_sf1, CLMC_sf2)
     deallocate (CLMC_pt1, CLMC_pt2, CLMC_st1, CLMC_st2)

    end subroutine regrid_carbon


  end subroutine re_tile

end module CatchmentCNRstMod
