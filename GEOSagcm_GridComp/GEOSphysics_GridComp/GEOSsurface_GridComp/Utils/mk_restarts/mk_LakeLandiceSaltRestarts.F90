#define I_AM_MAIN
#include "MAPL_Generic.h"
program mk_LakeLandiceSaltRestarts

  use netcdf

  use MAPL
  use mk_restarts_getidsMod, only: GetIDS,ReadTileFile_IntLatLon
  use PFIO
  use gFTL_StringVector

  implicit none

  character*256 :: Usage="mk_LakeLandiceSaltRestarts OutTileFile InTileFile InRestart mask"
  character*256 :: OutTileFile
  character*256 :: InTileFile
  character*256 :: InRestart
  character*256 :: arg

  integer :: i, rc, jc, iostat, iargc, n, mask,j,k,otiles,nsubtiles,l,itiles,nwords
  integer, pointer  :: Lono(:), Lato(:), Id(:), Pf(:)
  integer, pointer  :: Loni(:), Lati(:)
  real, allocatable :: varIn(:),varOut(:)
  real*8, allocatable :: varIn8(:),varOut8(:)
  real, allocatable :: var2(:,:)
  real, allocatable :: dummy(:)

  integer           :: zoom
#ifndef __GFORTRAN__
  integer              :: ftell
  external             :: ftell
#endif
  integer              :: bpos, epos, ntot
  integer, allocatable :: nrecs(:), mrecs(:)
  type(Netcdf4_Fileformatter) :: InFmt, OutFmt
  type(FileMetadata) :: InCfg, OutCfg
  integer :: dim1, dim2
  integer              :: ndims
  character*256        :: OutFileName
  integer              :: dimSizes(3)
  integer              :: filetype,nVars
  type(StringVariableMap), pointer :: variables
  type(Variable), pointer :: myVariable
  type(StringVariableMapIterator) :: var_iter
  type(StringVector), pointer :: var_dimensions
  character(len=:), pointer :: vname,dname
  integer :: dataType
  integer :: status
!---------------------------------------------------------------------------

  I = iargc()

  if(I /= 5) then
     print *, "Wrong Number of arguments: ", i
     print *, trim(Usage)
     call exit(1)
  end if

  call getarg(1,OutTileFile)
  call getarg(2,InTileFile)
  call getarg(3,InRestart)
  call getarg(4,arg)
  read(arg,*) mask
  call getarg(5,arg)
  read(arg,*) zoom

! Read Output Tile File .til file
! to get the index into the pfafsttater table

  call ReadTileFile_IntLatLon(OutTileFile,Pf,Id,lono,lato,zoom,mask)
  deallocate(Pf,Id)

  call ReadTileFile_IntLatLon(InTileFile ,Pf,Id,loni,lati,zoom,mask)
  deallocate(Pf,Id)

  nullify(Pf)
  nullify(Id)

  itiles = size(loni)  ! Input  Tile Size
  otiles = size(lono)  ! Output Tile Size
  allocate(Id (otiles))

  call GetIds(loni,lati,lono,lato,zoom,Id)

  call MAPL_NCIOGetFileType(InRestart, filetype,rc=rc)

  if (filetype == 0) then

     call InFmt%open(InRestart,pFIO_READ,rc=rc)
     InCfg = InFmt%read(rc=rc)
     call MAPL_IOChangeRes(InCfg,OutCfg,(/'tile'/),(/otiles/),rc=rc)

     i = index(InRestart,'/',back=.true.)
     OutFileName = "OutData/"//trim(InRestart(i+1:))
     call OutFmt%create(OutFileName,rc=rc)
     call OutFmt%write(OutCfg,rc=rc)

     allocate( varIn(itiles) )
     allocate( varOut(otiles) )
     allocate( varIn8(itiles) )
     allocate( varOut8(otiles) )

     call MAPL_IOCountNonDimVars(InCfg,nVars)
     variables => InCfg%get_variables()
     var_iter = variables%ftn_begin()
     do while (var_iter /= variables%ftn_end())
        call var_iter%next() 

        vname => var_iter%first()
        myVariable => var_iter%second()
        var_dimensions => myVariable%get_dimensions()
        dataType = myVariable%get_type()

        if (.not.InCfg%is_coordinate_variable(vname)) then

           ndims = var_dimensions%size()

           write(*,*)"Writing ",trim(vname)
           if (ndims == 1) then
              if (dataType == pFIO_REAL64) then
                 call MAPL_VarRead(InFmt,vname,varIn8, __RC__)
                 do i=1,otiles
                    varOut8(i) = varIn8(id(i))
                 enddo
                 call MAPL_VarWrite(OutFmt,vname,varOut8)
              else
                 call MAPL_VarRead(InFmt,vname,varIn, __RC__)
                 do i=1,otiles
                    varOut(i) = varIn(id(i))
                 enddo
                 call MAPL_VarWrite(OutFmt,vname,varOut)
              endif
           else if (ndims == 2) then
           
              dname => myVariable%get_ith_dimension(2)
              dim1=InCfg%get_dimension(dname)
        
              do j=1,dim1
                 if (dataType == pFIO_REAL64) then
                    call MAPL_VarRead(InFmt,vname,varIn8,offset1=j, __RC__)
                    do i=1,otiles
                       varOut8(i) = varIn8(id(i))
                    enddo
                    call MAPL_VarWrite(OutFmt,vname,varOut8,offset1=j)
                 else
                    call MAPL_VarRead(InFmt,vname,varIn,offset1=j, __RC__)
                    do i=1,otiles
                       varOut(i) = varIn(id(i))
                    enddo
                    call MAPL_VarWrite(OutFmt,vname,varOut,offset1=j)
                 endif
              enddo
           else if (ndims == 3) then
              
              dname => myVariable%get_ith_dimension(2)
              dim1=InCfg%get_dimension(dname)
              dname => myVariable%get_ith_dimension(3)
              dim2=InCfg%get_dimension(dname)

              do k=1,dim2
                 do j=1,dim1
                    if (dataType == pFIO_REAL64) then
                       call MAPL_VarRead(InFmt,vname,varIn8,offset1=j,offset2=k, __RC__)
                       do i=1,otiles
                          varOut8(i) = varIn8(id(i))
                       enddo
                       call MAPL_VarWrite(OutFmt,vname,varOut8,offset1=j,offset2=k)
                    else
                       call MAPL_VarRead(InFmt,vname,varIn,offset1=j,offset2=k, __RC__)
                       do i=1,otiles
                          varOut(i) = varIn(id(i))
                       enddo
                       call MAPL_VarWrite(OutFmt,vname,varOut,offset1=j,offset2=k)
                    endif
                 enddo
              enddo
           
           end if

        end if
       
     enddo

     call OutFmt%close()
     deallocate( varIn, varOut )
     deallocate( varIn8, varOut8 )

  else

     i = index(InRestart,'/',back=.true.)

     open(unit=40,FILE="OutData/"//trim(InRestart(i+1:)),form='unformatted',&
          status='unknown',convert='little_endian')

     open(unit=50,FILE=InRestart,form='unformatted',&
          status='old',convert='little_endian')

! Determine NWORDS for Each Input Record
! --------------------------------------
               rc =  0
             bpos =  0
             ntot =  0
     do while( rc.eq.0 )

       read (50,iostat=rc)
       if( rc.eq.0 ) then
           ntot = ntot + 1
           epos = ftell(50)          ! ending position of file pointer
         nwords = (epos-bpos)/4-2    ! record size (in 4 byte words;
         write(6,100) ntot, nwords
         if( ntot.eq.1 ) then
             allocate( nrecs(ntot) )
                       nrecs(ntot) = nwords
         else
             allocate( mrecs(  ntot-1) )
                       mrecs(1:ntot-1) = nrecs
           deallocate( nrecs )
             allocate( nrecs(ntot) )
                       nrecs(1:ntot-1) = mrecs
                       nrecs(ntot)     = nwords
           deallocate( mrecs )
         endif
         bpos = epos
       endif

     enddo

 100 format(1x,'Record #: ',i3,' Total  Tiles: ',i10)
     print *
     rewind(50)

   ! Read and Write Tile or TileTile Data until EOF
   ! ----------------------------------------------
     do n=1,ntot
        nsubtiles = nrecs(n)/itiles
          allocate(  var2(itiles,nsubtiles) )
          read (50)  var2
          if( nsubtiles.eq.1 ) then
              print *, 'Writing Tile_Only Data ...'
          else
              print *, 'Writing TileTile  Data ..., NSUBTILES = ',nsubtiles
          endif
          write(40)((var2(id(i),j),i=1,otiles),j=1,nsubtiles)
         deallocate( var2 )
     enddo

  end if

contains

end program mk_LakeLandiceSaltRestarts

