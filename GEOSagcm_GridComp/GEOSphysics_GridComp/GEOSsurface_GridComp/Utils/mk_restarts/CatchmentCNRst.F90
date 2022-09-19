#include "MAPL_Generic.h"

module CatchmentCNRstMod
  use mk_restarts_getidsMod
  use mpi
  use MAPL
  use CatchmentRstMod, only : CatchmentRst
  implicit none

  real, parameter :: ECCENTRICITY  = 0.0167
  real, parameter :: PERIHELION    = 102.0
  real, parameter :: OBLIQUITY     = 23.45
  integer, parameter :: EQUINOX    = 80

  integer, parameter :: nveg    = 4
  integer, parameter :: nzone = 3
  integer, parameter :: VAR_COL_CLM40 = 40 ! number of CN column restart variables
  integer, parameter :: VAR_PFT_CLM40 = 74 ! number of CN PFT variables per column
  integer, parameter :: npft    = 19
  integer, parameter :: npft_clm45    = 19
  integer, parameter :: VAR_COL_CLM45 = 35 ! number of CN column restart variables
  integer, parameter :: VAR_PFT_CLM45 = 75 ! number of CN PFT variables per column 
  real,    parameter :: nan = O'17760000000'
  real,    parameter :: fmin= 1.e-4 ! ignore vegetation fractions at or below this value
  integer :: iclass(npft) = (/1,1,2,3,3,4,5,5,6,7,8,9,10,11,12,11,12,11,12/)

  type, extends(CatchmentRst) :: CatchmentCNRst
     logical :: isCLM45
     integer :: VAR_COL
     integer :: VAR_PFT
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
     ! below is not necessary. It is not read. It is set to 0 during writing
     !real, allocatable :: bflowm(:)
     !real, allocatable :: totwatm(:)
     !real, allocatable :: tairm(:)
     !real, allocatable :: tpm(:)
     !real, allocatable :: cnsum(:)
     !real, allocatable :: sndzm(:)
     !real, allocatable :: asnowm(:)
     !real, allocatable :: ar1m(:)
     !real, allocatable :: rainfm(:)
     !real, allocatable :: rhm(:)
     !real, allocatable :: runsrfm(:)
     !real, allocatable :: snowfm(:)
     !real, allocatable :: windm(:)
     !real, allocatable :: tprec10d(:)
     !real, allocatable :: tprec60d(:)
     !real, allocatable :: t2m10d(:)
     !real, allocatable :: sfmcm(:)
     !real, allocatable :: psnsunm(:,:,:)
     !real, allocatable :: psnsham(:,:,:)
     
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
     call formatter%open(filename, pFIO_READ, __RC__)
     meta  = formatter%read(__RC__)
     ntiles = meta%get_dimension('tile', __RC__)
     catch%ntiles = ntiles
     catch%meta  = meta
     catch%time = time
     if (index(cnclm, '40') /=0) then
        catch%VAR_COL = VAR_COL_CLM40
        catch%VAR_PFT = VAR_PFT_CLM40
     endif
     if (index(cnclm, '45') /=0) then
        catch%VAR_COL = VAR_COL_CLM45
        catch%VAR_PFT = VAR_PFT_CLM45
        catch%isCLM45 = .true.
     endif

     if (myid == 0) then
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
     catch%ntiles = meta%get_dimension('tile', __RC__)
     catch%time = time
     catch%meta = meta
     if (index(cnclm, '40') /=0) then
        catch%VAR_COL = VAR_COL_CLM40
        catch%VAR_PFT = VAR_PFT_CLM40
     endif
     if (index(cnclm, '45') /=0) then
        catch%VAR_COL = VAR_COL_CLM45
        catch%VAR_PFT = VAR_PFT_CLM45
        catch%isCLM45 = .true.
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

     myVariable => meta%get_variable("ITY")
     dname => myVariable%get_ith_dimension(2)
     dim1 = meta%get_dimension(dname)
     do j=1,dim1
        call MAPL_VarWrite(formatter,"ITY",this%cnity(:,j),offset1=j)
        call MAPL_VarWrite(formatter,"FVG",this%fvg(:,j),offset1=j)
        call MAPL_VarWrite(formatter,"TG",this%tg(:,j),offset1=j)
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

     allocate (var(dim1),source = 0.)

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

   subroutine allocate_cn(this,rc)
     class(CatchmentCNRst), intent(inout) :: this
     integer, optional, intent(out):: rc
     integer :: status
     integer  :: ncol,npft, ntiles

     ntiles = this%ntiles
     ncol = nzone* this%VAR_COL 
     npft = nzone*nveg*this%VAR_PFT

     call this%CatchmentRst%allocate_catch(__RC__)

     ! W.Jiang notes : some varaiables are not allocated because they are set to zero directly during write
     allocate(this%cnity(ntiles,nveg))
     allocate(this%fvg(ntiles,nveg))
     allocate(this%tg(ntiles,nveg))
     allocate(this%tgwm(ntiles,nzone))
     allocate(this%rzmm(ntiles,nzone))
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
   end subroutine allocate_cn

   SUBROUTINE add_bcs_to_cnrst (this, SURFLAY, OutBcsDir,rc)
    class(CatchmentCNRst), intent(inout) :: this
    real, intent (in)                    :: SURFLAY
    character(*), intent (in)            :: OutBcsDir
    integer, optional, intent(out) :: rc
    real, allocatable :: CLMC_pf1(:), CLMC_pf2(:), CLMC_sf1(:), CLMC_sf2(:)
    real, allocatable :: CLMC_pt1(:), CLMC_pt2(:), CLMC_st1(:), CLMC_st2(:)    
    real, allocatable :: NDEP(:), BVISDR(:), BVISDF(:), BNIRDR(:), BNIRDF(:) 
    real, allocatable :: T2(:), var1(:), hdm(:), fc(:), gdp(:), peatf(:)
    integer, allocatable :: ity(:), abm (:)
    integer       :: STATUS, ntiles, unit27, unit28, unit29, unit30
    integer       :: idum, i,j,n, ib, nv
    real          :: rdum, zdep1, zdep2, zdep3, zmet, term1, term2, bare,fvg(4)
    logical       :: NEWLAND
    logical       :: file_exists

    type(NetCDF4_Fileformatter) :: CatchCNFmt
    character*256        :: Iam = "add_bcs"

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
    allocate (peatf(ntiles), abm(ntiles), var1(ntiles))

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

    if (this%isCLM45 ) then

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

     this%cnity = reshape([CLMC_pt1,CLMC_pt2,CLMC_st1,CLMC_st2],[ntiles,4])
     this%fvg   = reshape([CLMC_pf1,CLMC_pf2,CLMC_sf1,CLMC_sf2],[ntiles,4])

     this%ndep = ndep
     this%t2   = t2
     this%BGALBVR = BVISDR
     this%BGALBVF = BVISDF
     this%BGALBNR = BNIRDR
     this%BGALBNF = BNIRDF
 
     if (this%isCLM45) then
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
         CLMC_pt1, CLMC_pt2,CLMC_st1,CLMC_st2, var_dum2, var_dum3
     integer                :: AGCM_YY,AGCM_MM,AGCM_DD,AGCM_HR=0,AGCM_DATE
     real,    allocatable, dimension(:,:) :: fveg_offl,  ityp_offl, tg_tmp
     real, allocatable :: var_off_col (:,:,:), var_off_pft (:,:,:,:)
     integer :: status, in_ntiles, out_ntiles, numprocs
     logical :: root_proc
     integer :: mpierr, n, i, k, tag, req, st, ed, myid, L, iv, nv,nz, var_col, var_pft
     character(*), parameter :: Iam = "CatchmentCN::Re_tile"


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

     do n = 1, in_ntiles
        tid_offl(n) = n
     enddo

     if (root_proc) then

        allocate (DAYX   (out_ntiles))

        READ(this%time(1:8),'(I8)') AGCM_DATE
        AGCM_YY = AGCM_DATE / 10000
        AGCM_MM = (AGCM_DATE - AGCM_YY*10000) / 100
        AGCM_DD = (AGCM_DATE - AGCM_YY*10000 - AGCM_MM*100)

        call compute_dayx (                                     &
                out_NTILES, AGCM_YY, AGCM_MM, AGCM_DD, AGCM_HR,        &
                this%LATG, DAYX)

        ! save the old vaues dimension (in_ntiles, nv)
        ityp_offl = this%cnity
        fveg_offl = this%fvg

        do n = 1, in_ntiles
           do nv = 1,nveg
              if(ityp_offl(n,nv)<0 .or. ityp_offl(n,nv)>npft)    stop 'ityp'
              if(fveg_offl(n,nv)<0..or. fveg_offl(n,nv)>1.00001) stop 'fveg'
           end do

           if (nint(this%tile_id(n)) /= n) stop ("cannot assign ity_offl to cnity and fvg_offl to fvg")

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
           call MPI_send(this%cnity(st,1),l, MPI_REAL, i, tag, MPI_COMM_WORLD, mpierr)
           call MPI_send(this%cnity(st,2),l, MPI_REAL, i, tag+1, MPI_COMM_WORLD, mpierr)
           call MPI_send(this%cnity(st,3),l, MPI_REAL, i, tag+2, MPI_COMM_WORLD, mpierr)
           call MPI_send(this%cnity(st,4),l, MPI_REAL, i, tag+3, MPI_COMM_WORLD, mpierr)
           call MPI_send(this%fvg(st,1),l,   MPI_REAL, i, tag+4, MPI_COMM_WORLD, mpierr)
           call MPI_send(this%fvg(st,2),l,   MPI_REAL, i, tag+5, MPI_COMM_WORLD, mpierr)
           call MPI_send(this%fvg(st,3),l,   MPI_REAL, i, tag+6, MPI_COMM_WORLD, mpierr)
           call MPI_send(this%fvg(st,4),l,   MPI_REAL, i, tag+7, MPI_COMM_WORLD, mpierr)
        enddo
        st  = low_ind(1)
        l   = nt_local(1)
        ed  = st + l -1
        CLMC_pt1 = this%cnity(st:ed,1)
        CLMC_pt2 = this%cnity(st:ed,2)
        CLMC_st1 = this%cnity(st:ed,3)
        CLMC_st2 = this%cnity(st:ed,4)
        CLMC_pf1 = this%fvg(st:ed,1)
        CLMC_pf2 = this%fvg(st:ed,2)
        CLMC_sf1 = this%fvg(st:ed,3)
        CLMC_sf2 = this%fvg(st:ed,4)
     else
        tag = myid*numprocs
        call MPI_RECV(CLMC_pt1,nt_local(myid+1) , MPI_REAL, 0, tag,  MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
        call MPI_RECV(CLMC_pt2,nt_local(myid+1) , MPI_REAL, 0, tag+1, MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
        call MPI_RECV(CLMC_st1,nt_local(myid+1) , MPI_REAL, 0, tag+2, MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
        call MPI_RECV(CLMC_st2,nt_local(myid+1) , MPI_REAL, 0, tag+3, MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
        call MPI_RECV(CLMC_pf1,nt_local(myid+1) , MPI_REAL, 0, tag+4, MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
        call MPI_RECV(CLMC_pf2,nt_local(myid+1) , MPI_REAL, 0, tag+5, MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
        call MPI_RECV(CLMC_sf1,nt_local(myid+1) , MPI_REAL, 0, tag+6, MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
        call MPI_RECV(CLMC_sf2,nt_local(myid+1) , MPI_REAL, 0, tag+7, MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
    endif

    call MPI_Barrier(MPI_COMM_WORLD, STATUS)
 
    if(root_proc) print*, "GetIDs...."

    call GetIds(this%lonc,this%latc,this%lonn,this%latt,id_loc_cn, tid_offl, &
             CLMC_pf1, CLMC_pf2, CLMC_sf1, CLMC_sf2, CLMC_pt1, CLMC_pt2,CLMC_st1,CLMC_st2, &
             fveg_offl, ityp_offl)

    call MPI_Barrier(MPI_COMM_WORLD, STATUS)

    if(root_proc) allocate (id_glb_cn  (out_ntiles,nveg))

    allocate (id_loc (out_ntiles))
    deallocate (CLMC_pf1, CLMC_pf2, CLMC_sf1, CLMC_sf2)
    deallocate (CLMC_pt1, CLMC_pt2, CLMC_st1, CLMC_st2)

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
        allocate (var_dum2    (1:in_ntiles))

        this%tile_id = [(i*1.0, i=1, out_ntiles)]

        allocate (tg_tmp(out_ntiles, 4),source = 0.)
        do i = 1, 3
          tg_tmp(:,i) = this%tg(this%id_glb(:),i)
        enddo        
        this%tg = tg_tmp

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

        call regrid_carbon (out_NTILES, in_ntiles,id_glb_cn, &
                DAYX, var_off_col,var_off_pft, ityp_offl, fveg_offl)
        deallocate (var_off_col,var_off_pft)
     endif
     call MPI_Barrier(MPI_COMM_WORLD, STATUS)

    _RETURN(_SUCCESS)

  contains
     SUBROUTINE regrid_carbon (NTILES, in_ntiles, id_glb, &
        DAYX, var_off_col, var_off_pft, ityp_offl, fveg_offl)

     ! write out regridded carbon variables
     implicit none
     integer, intent (in) :: NTILES, in_ntiles,id_glb (ntiles,nveg)
     real, intent (in)    :: DAYX (NTILES), var_off_col(in_ntiles,NZONE,var_col), var_off_pft(in_ntiles,NZONE, NVEG, var_pft)
     real, intent (in),  dimension(in_ntiles,nveg) :: fveg_offl,  ityp_offl
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

     CLMC_pt1 = this%cnity(:,1)
     CLMC_pt2 = this%cnity(:,2)
     CLMC_st1 = this%cnity(:,3)
     CLMC_st2 = this%cnity(:,4)
     CLMC_pf1 = this%fvg(:,1)
     CLMC_pf2 = this%fvg(:,2)
     CLMC_sf1 = this%fvg(:,3)
     CLMC_sf2 = this%fvg(:,4)

     allocate (var_col_out (1: NTILES, 1 : nzone,1 : var_col))
     allocate (var_pft_out (1: NTILES, 1 : nzone,1 : nveg, 1 : var_pft))

     var_col_out = 0.
     var_pft_out = NaN

     OUT_TILE : DO N = 1, NTILES

        ! if(mod (n,1000) == 0) print *, myid +1, n, Id_glb(n,:)

        NVLOOP2 : do nv = 1, nveg

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
           if (fveg_new > fmin) then

              offl_cell    = Id_glb(n,nv)

              if(ityp_new      == ityp_offl (offl_cell,nv) .and. fveg_offl (offl_cell,nv)> fmin) then
                 iv = nv                                     ! same type fraction (primary of secondary)                          
              else if(ityp_new == ityp_offl (offl_cell,nx) .and. fveg_offl (offl_cell,nx)> fmin) then
                 iv = nx                                     ! not same fraction
              else if(iclass(ityp_new)==iclass(ityp_offl(offl_cell,nv)) .and. fveg_offl (offl_cell,nv)> fmin) then
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

                 if (nv == 1) ityp_new = CLMC_pt1(n)
                 if (nv == 1) fveg_new = CLMC_pf1(n)
                 if (nv == 2) ityp_new = CLMC_pt2(n)
                 if (nv == 2) fveg_new = CLMC_pf2(n)
                 if (nv == 3) ityp_new = CLMC_st1(n)
                 if (nv == 3) fveg_new = CLMC_sf1(n)
                 if (nv == 4) ityp_new = CLMC_st2(n)
                 if (nv == 4) fveg_new = CLMC_sf2(n)

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

    subroutine compute_dayx (                               &
       NTILES, AGCM_YY, AGCM_MM, AGCM_DD, AGCM_HR,        &
       LATT, DAYX)

      implicit none

      integer, intent (in) :: NTILES,AGCM_YY,AGCM_MM,AGCM_DD,AGCM_HR
      real, dimension (NTILES), intent (in)  :: LATT
      real, dimension (NTILES), intent (out) :: DAYX
      integer, parameter :: DT = 900
      integer, parameter :: ncycle = 1461 ! number of days in a 4-year leap cycle (365*4 + 1)   
      real, dimension(ncycle) :: zc, zs
      integer :: dofyr, sec,YEARS_PER_CYCLE, DAYS_PER_CYCLE, year, iday, idayp1, nn, n
      real    :: fac, YEARLEN, zsin, zcos, declin

      dofyr = AGCM_DD
      if(AGCM_MM >  1) dofyr = dofyr + 31
      if(AGCM_MM >  2) then
         dofyr = dofyr + 28
         if(mod(AGCM_YY,4) == 0) dofyr = dofyr + 1
      endif
      if(AGCM_MM >  3) dofyr = dofyr + 31
      if(AGCM_MM >  4) dofyr = dofyr + 30
      if(AGCM_MM >  5) dofyr = dofyr + 31
      if(AGCM_MM >  6) dofyr = dofyr + 30
      if(AGCM_MM >  7) dofyr = dofyr + 31
      if(AGCM_MM >  8) dofyr = dofyr + 31
      if(AGCM_MM >  9) dofyr = dofyr + 30
      if(AGCM_MM > 10) dofyr = dofyr + 31
      if(AGCM_MM > 11) dofyr = dofyr + 30

      sec = AGCM_HR * 3600 - DT ! subtract DT to get time of previous physics step
      fac = real(sec) / 86400.


      call orbit_create(zs,zc,ncycle) ! GEOS5 leap cycle routine

      YEARLEN = 365.25

      !  Compute length of leap cycle
      !------------------------------

      if(YEARLEN-int(YEARLEN) > 0.) then
         YEARS_PER_CYCLE = nint(1./(YEARLEN-int(YEARLEN)))
      else
         YEARS_PER_CYCLE = 1
      endif

      DAYS_PER_CYCLE=nint(YEARLEN*YEARS_PER_CYCLE)

      ! declination & daylength
      ! -----------------------

      YEAR = mod(AGCM_YY-1,YEARS_PER_CYCLE)

      IDAY = YEAR*int(YEARLEN)+dofyr
      IDAYP1 = mod(IDAY,DAYS_PER_CYCLE) + 1

      ZSin = ZS(IDAYP1)*FAC + ZS(IDAY)*(1.-FAC) !   sine of solar declination
      ZCos = ZC(IDAYP1)*FAC + ZC(IDAY)*(1.-FAC) ! cosine of solar declination

      nn = 0
      do n = 1,days_per_cycle
         nn = nn + 1
         if(nn > 365) nn = nn - 365
         !     print *, 'cycle:',n,nn,asin(ZS(n))
      end do
      declin = asin(ZSin)

      ! compute daylength on input tile space (accounts for any change in physics time step)  
      !  do n = 1,ntiles_cn
      !     fac = -(sin((latc(n)/zoom)*(MAPL_PI/180.))*zsin)/(cos((latc(n)/zoom)*(MAPL_PI/180.))*zcos)
      !     fac = min(1.,max(-1.,fac))
      !     dayl(n) = (86400./MAPL_PI) * acos(fac)   ! daylength (seconds)
      !  end do

      ! compute daylength on output tile space (accounts for lat shift due to split & change in time step)

      do n = 1,ntiles
         fac = -(sin(latt(n)*(MAPL_PI/180.))*zsin)/(cos(latt(n)*(MAPL_PI/180.))*zcos)
         fac = min(1.,max(-1.,fac))
         dayx(n) = (86400./MAPL_PI) * acos(fac)   ! daylength (seconds)
      end do

      ! print *,'DAYX : ', minval(dayx),maxval(dayx), minval(latt), maxval(latt), zsin, zcos, dofyr, iday, idayp1, declin

    end subroutine compute_dayx

  ! *****************************************************************************

    subroutine orbit_create(zs,zc,ncycle)
      implicit none

      integer, intent(in) :: ncycle
      real, intent(out), dimension(ncycle) :: zs, zc

      integer :: YEARS_PER_CYCLE, DAYS_PER_CYCLE
      integer :: K, KP !, KM
      real*8  :: T1, T2, T3, T4, FUN, Y, SOB, OMG, PRH, TT
      real*8  :: YEARLEN

       !  STATEMENT FUNCTION

      FUN(Y) = OMG*(1.0-ECCENTRICITY*cos(Y-PRH))**2

      YEARLEN = 365.25

       !  Factors involving the orbital parameters
       !------------------------------------------

      OMG  = (2.0*MAPL_PI/YEARLEN) / (sqrt(1.-ECCENTRICITY**2)**3)
      PRH  = PERIHELION*(MAPL_PI/180.)
      SOB  = sin(OBLIQUITY*(MAPL_PI/180.))

       !  Compute length of leap cycle
       !------------------------------

      if(YEARLEN-int(YEARLEN) > 0.) then
          YEARS_PER_CYCLE = nint(1./(YEARLEN-int(YEARLEN)))
      else
          YEARS_PER_CYCLE = 1
      endif


      DAYS_PER_CYCLE=nint(YEARLEN*YEARS_PER_CYCLE)

      if(days_per_cycle /= ncycle) stop 'bad cycle'

       !   ZS:   Sine of declination
       !   ZC:   Cosine of declination

       !  Begin integration at vernal equinox

      KP           = EQUINOX
      TT           = 0.0
      ZS(KP) = sin(TT)*SOB
      ZC(KP) = sqrt(1.0-ZS(KP)**2)

       !  Integrate orbit for entire leap cycle using Runge-Kutta

      do K=2,DAYS_PER_CYCLE
          T1 = FUN(TT       )
          T2 = FUN(TT+T1*0.5)
          T3 = FUN(TT+T2*0.5)
          T4 = FUN(TT+T3    )
          KP  = mod(KP,DAYS_PER_CYCLE) + 1
          TT  = TT + (T1 + 2.0*(T2 + T3) + T4) / 6.0
          ZS(KP) = sin(TT)*SOB
          ZC(KP) = sqrt(1.0-ZS(KP)**2)
      end do
    end subroutine orbit_create

  end subroutine re_tile

end module CatchmentCNRstMod
