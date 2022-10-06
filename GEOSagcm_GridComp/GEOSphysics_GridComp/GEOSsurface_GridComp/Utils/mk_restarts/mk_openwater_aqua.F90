#define I_AM_MAIN
#include "MAPL_Generic.h"

program makeOpenWaterRestart

! should be called with newTile oldTile oldOpenWaterRst
! open both tile files. For now make sure we have the same number of latitudes
! allocate integer array(NLATS) to store the index to the "old" tile
! loop over the lats; find the nearest ocean point that has comparable LAT (first match); save the index; special care is needed for Antarctica; for now we are replicating "last" ocean LAT south
! loop over the output tiles; find their LAT; use the index to copy data from the "good-old-restart"
! write the new restart and we are done

  use netcdf

  use MAPL
  use mk_restarts_getidsMod, only: ReadTileFile_IntLatLon
  use PFIO
  use gFTL_StringVector

  implicit none

  character*256 :: Usage=" OutTileFile InTileFile InRestart"
  character*256 :: OutTileFile
  character*256 :: InTileFile
  character*256 :: InRestart
  character*256 :: arg

  integer :: i, rc, jc, iostat, iargc, n, mask,j,k,otiles,nsubtiles,l,itiles,nwords
  integer, pointer  :: Lono(:), Lato(:), Id(:), Pf(:)
  integer, pointer  :: Loni(:), Lati(:)
  real, allocatable :: varIn(:),varOut(:)
  real, allocatable :: var2(:,:)
  real, allocatable :: dummy(:)

  type(Netcdf4_Fileformatter) :: InFmt, OutFmt
  type(FileMetadata) :: InCfg, OutCfg
  integer :: dim1, dim2
  integer              :: ndims
  character*256        :: OutFileName
  character*256        :: program_name
  integer              :: filetype,nVars
  type(StringVariableMap), pointer :: variables
  type(Variable), pointer :: myVariable
  type(StringVariableMapIterator) :: var_iter
  type(StringVector), pointer :: var_dimensions
  character(len=:), pointer :: vname,dname
  integer :: dataType
  integer :: status
  integer :: argc
  integer :: lt
  integer, parameter :: JM=180
!--------------------------------------------------------------------

  call get_command_argument(0, program_name)
  argc = command_argument_count()

  if(argc /= 3) then
     print *, "Wrong Number of arguments: ", i
     print *, trim(program_name)//trim(Usage)
     call exit(1)
  end if

  call get_command_argument(1,OutTileFile)
  call get_command_argument(2,InTileFile)
  call get_command_argument(3,InRestart)

! Read Output Tile File .til file
! to get the index into the pfafsttater table

  call ReadTileFile_IntLatLon(OutTileFile,Pf,Id,lono,lato,zoom=1,mask=0)
  deallocate(Pf,Id,lono)

  call ReadTileFile_IntLatLon(InTileFile ,Pf,Id,loni,lati,zoom=1,mask=0)
  deallocate(Pf,Id,loni)

  nullify(Pf)
  nullify(Id)

  itiles = size(lati)  ! Input  Tile Size
  otiles = size(lato)  ! Output Tile Size

  allocate(Id (JM), _STAT)
  Id = 0

  do i = 1, itiles
     j = (lati(i) + 90)+1
!@     print *,'DEBUG:',i,j,lati(i),itiles
     _ASSERT(j>0 .and. j<=JM, "wrong lat index")
     if (id(j) == 0) id(j) = i
  end do

! rounding at the North Pole
  id(jm) = id(jm-1)

! deal with Antarctica
  if (id(JM) == 0) then
     do i = 1,jm
        print *,i, id(i)
     end do
  end if
  _ASSERT(id(JM) /= 0, "North Pole should be ocean")
  do j=jm,1,-1
     if (id(j) == 0) id(j) = id(j+1)
  end do
  _ASSERT(count(id==0) == 0, "zero index is not allowed")

  call MAPL_NCIOGetFileType(InRestart, filetype,_RC)

  _ASSERT(filetype == 0, "This code accepts only NetCDF input")

  call InFmt%open(InRestart,pFIO_READ,_RC)
  InCfg = InFmt%read(_RC)
  call MAPL_IOChangeRes(InCfg,OutCfg,(/'tile'/),(/otiles/),_RC)

  i = index(InRestart,'/',back=.true.)
  OutFileName = "OutData/"//trim(InRestart(i+1:))
  call OutFmt%create(OutFileName,_RC)
  call OutFmt%write(OutCfg,_RC)

  allocate( varIn(itiles) )
  allocate( varOut(otiles) )
  
  call MAPL_IOCountNonDimVars(InCfg,nVars)
  variables => InCfg%get_variables()
  var_iter = variables%begin()
  do while (var_iter /= variables%end())

     vname => var_iter%key()
     myVariable => var_iter%value()
     var_dimensions => myVariable%get_dimensions()
     dataType = myVariable%get_type()
!     _ASSERT(dataType == pFIO_REAL32,"Only real32 is supported") 
     if (.not.InCfg%is_coordinate_variable(vname)) then

        ndims = var_dimensions%size()

        write(*,*)"Writing ",trim(vname)
        if (ndims == 1) then
           call MAPL_VarRead(InFmt,vname,varIn, __RC__)
           do i=1,otiles
              lt = (lato(i) + 90) + 1  !hope ll is between 1 and 180
              lt = max(min(lt,JM),1)
              if (lt < 1 .or. lt>JM) then
                 print *,'ERROR: wrong index ',lt,i,lato(i)
              end if
              varOut(i) = varIn(id(lt))
           enddo
           call MAPL_VarWrite(OutFmt,vname,varOut)
        else if (ndims == 2) then
           dname => myVariable%get_ith_dimension(2)
           dim1=InCfg%get_dimension(dname)
        
           do j=1,dim1
              call MAPL_VarRead(InFmt,vname,varIn,offset1=j, __RC__)
              do i=1,otiles
                 lt = (lato(i) + 90) + 1  !hope ll is between 1 and 180
                 lt = max(min(lt,JM),1)
                 varOut(i) = varIn(id(lt))
              enddo
              call MAPL_VarWrite(OutFmt,vname,varOut,offset1=j)
           enddo
        else 
           _ASSERT(.false., "Unsupported variable rank")
        end if
     end if
       
     call var_iter%next() 
  enddo

  call OutFmt%close()
  deallocate( varIn, varOut )

end program makeOpenWaterRestart

