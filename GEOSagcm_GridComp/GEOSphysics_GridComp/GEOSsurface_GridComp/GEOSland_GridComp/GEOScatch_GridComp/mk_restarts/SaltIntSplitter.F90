program SaltIntSplitter

! $Id$

  use netcdf
  use MAPL_HashMod
  use MAPL_IOMod

  implicit none


  character*256 :: Usage="SaltIntSplitter InTileFile InRestart"
  character*256 :: InTileFile
  character*256 :: InRestart
  character*256 :: arg


  integer :: i, rc, jc, iostat, iargc, n, mask,j,k,otiles,nsubtiles,l,itiles,nwords
  integer, pointer  :: Lono(:), Lato(:), Id(:), Pf(:)
  integer, pointer  :: Loni(:), Lati(:)
  real, allocatable :: varIn(:),varOut(:)
  real*8, allocatable :: varInR8(:),varOutR8(:)
  real, allocatable :: var2(:,:)
  real, allocatable :: dummy(:)

  integer, parameter   :: zoom=1
#ifndef __GFORTRAN__
  integer              :: ftell
  external             :: ftell
#endif
  integer              :: bpos, epos, ntot
  integer, allocatable :: nrecs(:), mrecs(:)
  type(MAPL_NCIO)      :: InNCIO, OutNCIO
  type(MAPL_NCIO)      :: WaterNCIO, IceNCIO
  integer              :: ndims
  character*256        :: OutFileName
  character*256        :: WaterFileName
  character*256        :: IceFileName
  integer              :: dimSizes(3)
  integer              :: filetype,nVars
  character*256        :: vname

!---------------------------------------------------------------------------

  I = iargc()

  if(I /= 2) then
     print *, "Wrong Number of arguments: ", i
     print *, trim(Usage)
     call exit(1)
  end if

  call getarg(1,InTileFile)
  call getarg(2,InRestart)

! Read Output Tile File .til file
! to get the index into the pfafsttater table

  call ReadTileFile(InTileFile ,Pf,Id,loni,lati, 0)
  deallocate(Pf,Id)

  nullify(Pf)
  nullify(Id)

  itiles = size(loni)  ! Input  Tile Size

  allocate( varIn(itiles) )
  allocate( varOut(itiles) )
  allocate( varInR8(itiles) )
  allocate( varOutR8(itiles) )

  call MAPL_NCIOGetFileType(InRestart, filetype,rc=rc)

  if (filetype == 0) then

     InNCIO = MAPL_NCIOOpen(InRestart,rc=rc)

     call MAPL_NCIOChangeRes(InNCIO,WaterNCIO,tileSize=itiles,rc=rc)
     call MAPL_NCIOChangeRes(InNCIO,IceNCIO,tileSize=itiles,rc=rc)
     i = index(InRestart,'/',back=.true.)
     WaterFileName = "OutData/openwater_internal_rst"
     IceFileName   = "OutData/seaicethermo_internal_rst"
     call MAPL_NCIOSet( WaterNCIO,filename=WaterFileName, overwriteVars=.true.)
     call MAPL_NCIOSet( IceNCIO,  filename=IceFileName, overwriteVars=.true.)
     do i =1,InNCIO%nDims
         if(trim(InNCIO%dims(i)%name) .eq. "subtile") then
            WaterNCIO%dims(i)%len = 1
            IceNCIO%dims(i)%len = InNCIO%dims(i)%len-1
         endif 
         if(trim(InNCIO%dims(i)%name) .eq. "unknown_dim4") then
            WaterNCIO%dims(i)%len = InNCIO%dims(i)%len-1
            IceNCIO%dims(i)%len = InNCIO%dims(i)%len-1
         endif 
     enddo   
     call MAPL_NCIOCreateFile(WaterNCIO)
     call MAPL_NCIOCreateFile(IceNCIO)


     call MAPL_NCIOGetDimSizes(InNCIO,nVars=nVars)
     do n=1,nVars
 
        call MAPL_NCIOGetVarName(InNCIO,n,vname)

        write(*,*)"Writing ",trim(vname)
        call MAPL_NCIOVarGetDims(InNCIO,vname,nDims,dimSizes)
        if (ndims == 1) then
           call MAPL_VarRead(InNCIO,vname,varIn)
           varOut(:) = varIn(:)    
           if (vname(2:6) == 'SKINI') then ! sea ice vars
             call MAPL_VarWrite(IceNCIO,vname,varOut)
           elseif (vname(1:6) == 'SLMASK') then
             call MAPL_VarWrite(IceNCIO,vname,varOut)
           else ! the rest goes to Openwater
             call MAPL_VarWrite(WaterNCIO,vname,varOut)
           endif
        else if (ndims == 2) then
           ! for AMIP rst, ice=1, water=2 
           !if (dimSizes(2) /= 2) then
           !   write(*,*) "not an AMIP rst"
           !   stop
           !endif   
           !print*,trim(vname), dimSizes(1), dimSizes(2)
           if (InNCIO%vars(n)%ncDataType == NF90_DOUBLE) then ! R8 vars only from coupled 
              if (vname(1:2) == 'FR') then ! FR dim changes from 6 to 5
                 do j=2,dimSizes(2)
                   call MAPL_VarRead(InNCIO,vname,varInR8,offset1=j)
                   call MAPL_VarWrite(IceNCIO,vname,varInR8,offset1=j-1)
                 enddo 
              else
                 do j=1,dimSizes(2)
                   call MAPL_VarRead(InNCIO,vname,varInR8,offset1=j)
                   call MAPL_VarWrite(IceNCIO,vname,varInR8,offset1=j)
                 enddo 
              endif
           else if (dimSizes(2) == 2) then ! AMIP
              call MAPL_VarRead(InNCIO,vname,varIn,offset1=1)
              call MAPL_VarWrite(IceNCIO,vname,varIn,offset1=1)
              call MAPL_VarRead(InNCIO,vname,varIn,offset1=2)
              call MAPL_VarWrite(WaterNCIO,vname,varIn,offset1=1)
           else
              if (vname(1:6) == 'TSKINI') then 
                 do j=1,dimSizes(2)
                   call MAPL_VarRead(InNCIO,vname,varIn,offset1=j)
                   call MAPL_VarWrite(IceNCIO,vname,varIn,offset1=j)
                 enddo 
              else
                 call MAPL_VarRead(InNCIO,vname,varIn,offset1=1)
                 call MAPL_VarWrite(WaterNCIO,vname,varIn,offset1=1)
                 do j=2,dimSizes(2)
                   call MAPL_VarRead(InNCIO,vname,varIn,offset1=j)
                   call MAPL_VarWrite(IceNCIO,vname,varIn,offset1=j-1)
                 enddo 
              endif 
           endif
           ! for coupled rst, water=1, ice=2,num_subtiles 
        else if (ndims == 3) then
           ! only coupled model internals conatin ndims=3 vars
           do k=1,dimSizes(3)
              do j=1,dimSizes(2)
                 if (InNCIO%vars(n)%ncDataType == NF90_DOUBLE) then 
                    call MAPL_VarRead(InNCIO,vname,varInR8,offset1=j,offset2=k)
                    call MAPL_VarWrite(IceNCIO,vname,varInR8,offset1=j,offset2=k)
                 else
                    call MAPL_VarRead(InNCIO,vname,varIn,offset1=j,offset2=k)
                    call MAPL_VarWrite(IceNCIO,vname,varIn,offset1=j,offset2=k)
                 endif
              enddo
           enddo
        end if
        
     enddo


  else

      i = index(InRestart,'/',back=.true.)

      open(unit=40,FILE="OutData/openwater_internal_rst",form='unformatted',&
           status='unknown',convert='little_endian')

      open(unit=41,FILE="OutData/seaicethermo_internal_rst",form='unformatted',&
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
              if(n <= 3) then
                 write(40)((var2(i,j),i=1,itiles),j=1,nsubtiles)
              elseif(n <=6 ) then
                 write(41)((var2(i,j),i=1,itiles),j=1,nsubtiles)
              else
                 write(40)((var2(i,j),i=1,itiles),j=1,nsubtiles)
              endif
        else
             print *, 'Writing TileTile  Data ..., NSUBTILES = ',nsubtiles
             ! for MERRA-2  ice=1, water=2
             write(41)((var2(i,j),i=1,itiles),j=1,1)
             write(40)((var2(i,j),i=1,itiles),j=2,2)
        endif
        deallocate( var2 )
     enddo


     !print*, 'Splitter only supports NETCDF rst for now!!' 
     !stop 1

  end if

  deallocate( varIn, varOut )
  deallocate( varInR8, varOutR8 )

contains

#include "getids.H"

end program SaltIntSplitter

