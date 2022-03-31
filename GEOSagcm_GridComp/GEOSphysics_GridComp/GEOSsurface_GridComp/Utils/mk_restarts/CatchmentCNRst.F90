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

  type, extends(CatchmentRst) :: CatchmentCNRst
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
     procedure :: allocatecn   
  endtype CatchmentCNRst

  interface CatchmentCNRst
     module procedure CatchmentCNRst_Create
  end interface

contains

  function CatchmentCNRst_create(filename,cnclm,rc) result (catch)
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
     
     call formatter%open(filename, pFIO_READ, __RC__)
     meta   = formatter%read(__RC__)
     ntiles = meta%get_dimension('tile', __RC__)
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
     if(index(cnclm,'45') /=0) then
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

   subroutine write_nc4(this, filename, meta, cnclm, rc)
     class(CatchmentCNRst), intent(inout):: this
     character(*), intent(in) :: filename
     type(FileMetadata), intent(inout) :: meta
     character(*), intent(in) :: cnclm
     integer, optional, intent(out):: rc

     type(Netcdf4_fileformatter) :: formatter
     integer :: status
     character(256) :: Iam = "write_nc4"
     logical :: clm45
     integer :: i,j, dim1,dim2
     real, dimension (:), allocatable :: var
     type(Variable), pointer :: myVariable
     character(len=:), pointer :: dname

     clm45 = .false.
     if (index(cnclm,'45') /=0) clm45 = .true.

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

     if (clm45) then
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

   subroutine allocatecn(this, ntiles, ncol, npft, rc)
     class(CatchmentCNRst), intent(inout) :: this
     integer, intent(in) :: ntiles,ncol,npft
     integer, optional, intent(out):: rc
     integer :: status

     call this%CatchmentRst%allocate(ntiles, __RC__)

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

end module CatchmentCNRstMod
