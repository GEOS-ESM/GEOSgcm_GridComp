#define I_AM_MAIN
#include "MAPL_Generic.h"
program SaltIntSplitter

  use MAPL_ConstantsMod,only: MAPL_PI,  MAPL_radius
  use netcdf
  use MAPL
  use mk_restarts_getidsMod, only: ReadTileFile_RealLatLon
  use gFTL_StringVector
  use gFTL_StringIntegerMap 

  implicit none

  character*256 :: Usage="SaltIntSplitter InTileFile InRestart"
  character*256 :: InTileFile
  character*256 :: InRestart
  character*256 :: arg

  integer :: i, rc, jc, iostat, iargc, n, mask,j,k,otiles,nsubtiles,l,itiles,nwords
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
  type(NetCDF4_Fileformatter) :: InFmt,OutFmt,WaterFmt,IceFmt
  type(FileMetadata) :: InCfg,OutCfg,WaterCfg,IceCfg
  integer              :: subtileSize,ungridSize
  type(StringVariableMap), pointer :: variables
  type(Variable), pointer :: myVariable
  type(StringVariableMapIterator) :: var_iter
  character(len=:), pointer :: var_name,dname
  type(StringVector), pointer :: var_dimensions

  integer              :: ndims,dataType
  character*256        :: OutFileName
  character*256        :: WaterFileName
  character*256        :: IceFileName
  integer              :: dimSizes(3)
  integer              :: filetype,nVars
  character*256        :: Iam = "SaltIntSplitter"
  integer :: status
  type (Variable), pointer :: global
  type (StringIntegerMap), pointer :: dimensions
  type (StringIntegerMap) :: water_dimensions
  type (StringIntegerMapIterator) :: d_iter

!---------------------------------------------------------------------------

  I = iargc()

  if(I /= 2) then
     print *, "Wrong Number of arguments: ", i
     print *, trim(Usage)
     call exit(1)
  end if

  call getarg(1,InTileFile)
  call getarg(2,InRestart)


  call ReadTileFile_RealLatLon(InTileFile, itiles, mask=0)

  allocate( varIn(itiles),   source = 0. )
  allocate( varOut(itiles),  source = 0. )
  allocate( varInR8(itiles), source = 0.d0 )
  allocate( varOutR8(itiles),source = 0.d0 )

  call MAPL_NCIOGetFileType(InRestart, filetype,rc=rc)

  if (filetype == 0) then

     call InFmt%open(InRestart,pFIO_Read,rc=rc)
     InCfg=InFmt%read(rc=rc)
     ! what other dimensions do we have
     if (InCfg%has_dimension('subtile')) then
        subtileSize = InCfg%get_dimension('subtile',rc=rc)
     else
        subtileSize = 0
     end if

     if (InCfg%has_dimension('unknown_dim4')) then
        ungridSize = InCfg%get_dimension('unknown_dim4',rc=rc)
     else
        ungridSize = 0
     end if
     
     dimensions => InCfg%get_dimensions()
     global     => Incfg%get_global_var()
     IceCfg   = FileMetaData(dimensions= dimensions, global=global)

     water_dimensions = dimensions
     d_iter = water_dimensions%find('unknown_dim1')
     if (d_iter /= water_dimensions%end()) call water_dimensions%erase(d_iter)
     d_iter = water_dimensions%find('unknown_dim2')
     if (d_iter /= water_dimensions%end()) call water_dimensions%erase(d_iter)
     d_iter = water_dimensions%find('unknown_dim3')
     if (d_iter /= water_dimensions%end()) call water_dimensions%erase(d_iter)
     d_iter = water_dimensions%find('unknown_dim4')
     if (d_iter /= water_dimensions%end()) call water_dimensions%erase(d_iter)

     WaterCfg = FileMetaData(dimensions= water_dimensions, global=global)

     call WaterCfg%modify_dimension('tile', itiles)
     call IceCfg%modify_dimension('tile', itiles)

     if ((subtileSize/=0) .and. (ungridSize==0)) then
        call WaterCfg%modify_dimension('subtile', 1)
        call IceCfg%modify_dimension('subtile', subtileSize-1)
     endif

     if((subtileSize==0) .and. (ungridSize/=0)) then
        call IceCfg%modify_dimension('unknown_dim4', ungridSize-1)
     endif
     
     if ((subtileSize/=0) .and. (ungridSize/=0)) then
        call IceCfg%modify_dimension('unknown_dim4', ungridSize-1)
        call WaterCfg%modify_dimension('subtile', 1)
        call IceCfg%modify_dimension('subtile', subtileSize-1)
     end if

 !########################################

     variables => InCfg%get_variables()
     var_iter = variables%begin()
     do while (var_iter /= variables%end())
 
        var_name => var_iter%key()
        myVariable => var_iter%value()
        var_dimensions => myVariable%get_dimensions()
        ndims = var_dimensions%size()
        dataType = myVariable%get_type()
        if (.not.InCfg%is_coordinate_variable(var_name)) then
           if (ndims == 1) then
              select case (var_name)
              case ('HSKINI','SSKINI','TSKINI', 'SLMASK') ! sea ice vars
                 call IceCfg%add_variable(var_name, myVariable)
              case default
                 call WaterCfg%add_variable(var_name, myVariable)
              end select
           else if (ndims == 2) then
              dname => myVariable%get_ith_dimension(2)
              dimSizes(2)=InCfg%get_dimension(dname)
              call iceCfg%add_variable(var_name, myVariable)
              if (dataType /= pFIO_REAL64) then ! R8 vars only from coupled 
                 if (dimSizes(2) == 2) then ! AMIP
                   call waterCfg%add_variable(var_name, myVariable)
                 else
                   if (var_name /= 'TSKINI' .and. var_name /= 'TAUAGE') then 
                     call waterCfg%add_variable(var_name, myVariable)
                   endif 
                endif
              endif
              ! for coupled rst, water=1, ice=2,num_subtiles 
           else if (ndims == 3) then
              call iceCfg%add_variable(var_name, myVariable)
           end if
        end if

        if (var_name == 'time') then
           call iceCfg%add_variable(var_name, myVariable)
           call waterCfg%add_variable(var_name, myVariable)
        endif
        call var_iter%next()   
     enddo

!####################

     i = index(InRestart,'/',back=.true.)
     WaterFileName = "OutData/openwater_internal_rst"
     IceFileName   = "OutData/seaicethermo_internal_rst"
     call WaterFmt%create(WaterFileName,rc=rc)
     call WaterFmt%write(WaterCfg,rc=rc)
     call IceFmt%create(IceFileName,rc=rc)
     call IceFmt%write(IceCfg,rc=rc)

     variables => InCfg%get_variables()
     var_iter = variables%begin()
     do while (var_iter /= variables%end())
 
        var_name => var_iter%key()
        myVariable => var_iter%value()
        var_dimensions => myVariable%get_dimensions()
        ndims = var_dimensions%size()
        dataType = myVariable%get_type()
        if (.not.InCfg%is_coordinate_variable(var_name)) then

           write(*,*)"Writing ",trim(var_name),ndims
           
           if (ndims == 1) then
              call MAPL_VarRead(InFmt,var_name,varIn, __RC__)
              varOut(:) = varIn(:)
              select case (var_name)
              case ('HSKINI','SSKINI','TSKINI') ! sea ice vars
                 call MAPL_VarWrite(IceFmt,var_name,varOut)
              case ('SLMASK')
                 call MAPL_VarWrite(IceFmt,var_name,varOut)
              case default
                 call MAPL_VarWrite(WaterFmt,var_name,varOut)
              end select
           else if (ndims == 2) then
              dname => myVariable%get_ith_dimension(2)
              dimSizes(2)=InCfg%get_dimension(dname)
              ! for AMIP rst, ice=1, water=2 
              !if (dimSizes(2) /= 2) then
              !   write(*,*) "not an AMIP rst"
              !   stop
              !endif   
              !print*,trim(var_name), dimSizes(1), dimSizes(2)
              if (dataType == pFIO_REAL64) then ! R8 vars only from coupled 
                 if (var_name(1:2) == 'FR') then ! FR dim changes from 6 to 5
                    do j=2,dimSizes(2)
                      call MAPL_VarRead(InFmt,var_name,varInR8,offset1=j, __RC__)
                      call MAPL_VarWrite(IceFmt,var_name,varInR8,offset1=j-1)
                    enddo 
                 else
                    do j=1,dimSizes(2)
                      call MAPL_VarRead(InFmt,var_name,varInR8,offset1=j, __RC__)
                      call MAPL_VarWrite(IceFmt,var_name,varInR8,offset1=j)
                    enddo 
                 endif
              else if (dimSizes(2) == 2) then ! AMIP
                 call MAPL_VarRead(InFmt,var_name,varIn,offset1=1, __RC__)
                 call MAPL_VarWrite(IceFmt,var_name,varIn,offset1=1)
                 call MAPL_VarRead(InFmt,var_name,varIn,offset1=2, __RC__)
                 call MAPL_VarWrite(WaterFmt,var_name,varIn,offset1=1)
              else
                 if (var_name == 'TSKINI' .or. var_name == 'TAUAGE') then 
                    do j=1,dimSizes(2)
                      call MAPL_VarRead(InFmt,var_name,varIn,offset1=j, __RC__)
                      call MAPL_VarWrite(IceFmt,var_name,varIn,offset1=j)
                    enddo 
                 else
                    call MAPL_VarRead(InFmt,var_name,varIn,offset1=1, __RC__)
                    call MAPL_VarWrite(WaterFmt,var_name,varIn,offset1=1)
                    do j=2,dimSizes(2)
                      call MAPL_VarRead(InFmt,var_name,varIn,offset1=j, __RC__)
                      call MAPL_VarWrite(IceFmt,var_name,varIn,offset1=j-1)
                    enddo 
                 endif 
              endif
              ! for coupled rst, water=1, ice=2,num_subtiles 
           else if (ndims == 3) then
              ! only coupled model internals conatin ndims=3 vars
              dname => myVariable%get_ith_dimension(2)
              dimSizes(2)=InCfg%get_dimension(dname)
              dname => myVariable%get_ith_dimension(3)
              dimSizes(3)=InCfg%get_dimension(dname)
              do k=1,dimSizes(3)
                 do j=1,dimSizes(2)
                    if (dataType == pFIO_REAL64) then 
                       call MAPL_VarRead(InFmt,var_name,varInR8,offset1=j,offset2=k, __RC__)
                       call MAPL_VarWrite(IceFmt,var_name,varInR8,offset1=j,offset2=k)
                    else
                       call MAPL_VarRead(InFmt,var_name,varIn,offset1=j,offset2=k, __RC__)
                       call MAPL_VarWrite(IceFmt,var_name,varIn,offset1=j,offset2=k)
                    endif
                 enddo
              enddo
           end if
        end if

        if (var_name == 'time') then
          call MAPL_VarWrite(IceFmt, 'time',[0.0d0])
          call MAPL_VarWrite(waterFmt,'time',[0.0d0])
        endif
     
        call var_iter%next()   
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
        allocate(  var2(itiles,nsubtiles), source = 0. )
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

end program SaltIntSplitter

