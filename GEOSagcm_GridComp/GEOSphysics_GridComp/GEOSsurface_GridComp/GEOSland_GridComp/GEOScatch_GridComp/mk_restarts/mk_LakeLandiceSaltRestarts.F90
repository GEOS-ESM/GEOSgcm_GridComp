program mk_LakeLandiceSaltRestarts

! $Id: 

  use netcdf

  use MAPL_ConstantsMod,only: MAPL_PI,  MAPL_radius
  use MAPL_HashMod
  use MAPL_IOMod

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
  type(MAPL_NCIO)    :: InNCIO, OutNCIO
  integer              :: ndims
  character*256        :: OutFileName
  integer              :: dimSizes(3)
  integer              :: filetype,nVars
  character*256        :: vname

  interface GetIds   
     procedure GetIds_fast_1p
     procedure GetIds_accurate_mpi
  end interface

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

  call ReadTileFile(OutTileFile,Pf,Id,lono,lato,mask)
  deallocate(Pf,Id)

  call ReadTileFile(InTileFile ,Pf,Id,loni,lati,mask)
  deallocate(Pf,Id)

  nullify(Pf)
  nullify(Id)

  itiles = size(loni)  ! Input  Tile Size
  otiles = size(lono)  ! Output Tile Size
  allocate(Id (otiles))

  call GetIds(loni,lati,lono,lato,Id)

  call MAPL_NCIOGetFileType(InRestart, filetype,rc=rc)

  if (filetype == 0) then

     InNCIO = MAPL_NCIOOpen(InRestart,rc=rc)

     call MAPL_NCIOChangeRes(InNCIO,OutNCIO,tileSize=otiles,rc=rc)
     i = index(InRestart,'/',back=.true.)
     OutFileName = "OutData/"//trim(InRestart(i+1:))
     call MAPL_NCIOSet( OutNCIO,filename=OutFileName )
     call MAPL_NCIOCreateFile(OutNCIO)

     allocate( varIn(itiles) )
     allocate( varOut(otiles) )
     allocate( varIn8(itiles) )
     allocate( varOut8(otiles) )

     call MAPL_NCIOGetDimSizes(InNCIO,nVars=nVars)
     do n=1,nVars
 
        call MAPL_NCIOGetVarName(InNCIO,n,vname)

        write(*,*)"Writing ",trim(vname)
        call MAPL_NCIOVarGetDims(InNCIO,vname,nDims,dimSizes)
        if (ndims == 1) then
           if(InNCIO%vars(n)%ncDataType == NF90_DOUBLE) then ! R8 vars only from coupled 
              call MAPL_VarRead(InNCIO,vname,varIn8)
              do i=1,otiles
                 varOut8(i) = varIn8(id(i))
              enddo
              call MAPL_VarWrite(OutNCIO,vname,varOut8)
           else
              call MAPL_VarRead(InNCIO,vname,varIn)
              do i=1,otiles
                 varOut(i) = varIn(id(i))
              enddo
              call MAPL_VarWrite(OutNCIO,vname,varOut)
           endif
        else if (ndims == 2) then
           
           do j=1,dimSizes(2)
              if(InNCIO%vars(n)%ncDataType == NF90_DOUBLE) then ! R8 vars only from coupled 
                call MAPL_VarRead(InNCIO,vname,varIn8,offset1=j)
                do i=1,otiles
                   varOut8(i) = varIn8(id(i))
                enddo
                call MAPL_VarWrite(OutNCIO,vname,varOut8,offset1=j)
              else
                call MAPL_VarRead(InNCIO,vname,varIn,offset1=j)
                do i=1,otiles
                   varOut(i) = varIn(id(i))
                enddo
                call MAPL_VarWrite(OutNCIO,vname,varOut,offset1=j)
              endif
           enddo
        else if (ndims == 3) then
           
           do k=1,dimSizes(3)
              do j=1,dimSizes(2)
                 if(InNCIO%vars(n)%ncDataType == NF90_DOUBLE) then ! R8 vars only from coupled 
                    call MAPL_VarRead(InNCIO,vname,varIn8,offset1=j,offset2=k)
                    do i=1,otiles
                       varOut8(i) = varIn8(id(i))
                    enddo
                    call MAPL_VarWrite(OutNCIO,vname,varOut8,offset1=j,offset2=k)
                 else
                    call MAPL_VarRead(InNCIO,vname,varIn,offset1=j,offset2=k)
                    do i=1,otiles
                       varOut(i) = varIn(id(i))
                    enddo
                    call MAPL_VarWrite(OutNCIO,vname,varOut,offset1=j,offset2=k)
                 endif
              enddo
           enddo
        
        end if
        
     enddo

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

#include "getids.H"

end program mk_LakeLandiceSaltRestarts

