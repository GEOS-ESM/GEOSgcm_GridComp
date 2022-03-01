#define I_AM_MAIN
#include "MAPL_Generic.h"

program SaltImpConverter

  use MAPL_ConstantsMod,only: MAPL_PI,  MAPL_radius
  use netcdf
  use MAPL
  use mk_restarts_getidsMod, only: ReadTileFile_IntLatLon
  use gFTL_StringVector
  implicit none

  character*256 :: Usage="SaltImpConverter InTileFile InImpRestart InIntRestart"
  character*256 :: InTileFile
  character*256 :: InRestart
  character*256 :: InImpRestart
  character*256 :: InIntRestart
  character*256 :: arg

  integer :: i, rc, jc, iostat, iargc, n, mask,j,k,otiles,nsubtiles,l,itiles,nwords
  integer, pointer  :: Lono(:), Lato(:), Id(:), Pf(:)
  integer, pointer  :: Loni(:), Lati(:)
  real, allocatable :: varIn(:),varOut(:)
  real, allocatable :: TW(:),SW(:)
  real*8, allocatable :: varInR8(:),varOutR8(:)

  integer, parameter   :: zoom=1
#ifndef __GFORTRAN__
  integer              :: ftell
  external             :: ftell
#endif
  integer              :: bpos, epos, ntot
  integer              :: foutID, status, TimID, TileID 
  integer, allocatable :: nrecs(:), mrecs(:)
  type(Netcdf4_Fileformatter) :: InImpFmt,OutFmt,InIntFmt
  type(FileMetadata)   :: InImpCfg,OutCfg,InIntCfg
  type(StringVariableMap), pointer :: variables
  type(Variable), pointer :: myVariable
  type(StringVariableMapIterator) :: var_iter
  type(StringVector), pointer :: var_dimensions
  type(Attribute), pointer :: attr
  character(len=:), pointer :: var_name
  integer              :: ndims
  character*256        :: OutFileName
  integer              :: dimSizes(3)
  integer              :: filetype,nVars
  integer              :: varid 
  character*256        :: vname
  character*256        :: longname
  character*256        :: units
  character*256        :: impNames(39)
  character*256        :: Iam = "SaltImpConverter"

  INCLUDE 'netcdf.inc'
!---------------------------------------------------------------------------

  Data impNames /    & 
         'ALW',      &
         'BLW',      &
     'LWDNSRF',      &
        'DRPAR'     ,&
        'DFPAR'     ,&
        'DRNIR'     ,&
        'DFNIR'     ,&
        'DRUVR'     ,&
        'DFUVR'     ,&
        'EVAP',      &
          'SH',      &
        'TAUX',      &
        'TAUY',      &
       'DEVAP',      &
         'DSH',      &
         'SNO',      &
          'TA',      &
          'QA',      &
          'UU',      &
 'UWINDLMTILE',      &
 'VWINDLMTILE',      &
          'DZ',      &
          'PS',      &
          'PCU'     ,&
          'PLS'     ,&
       'THATM',      &
       'QHATM',      &
       'UHATM',      &
       'VHATM',      &
       'CTATM',      &
       'CQATM',      &
       'CMATM',      &
     'FRACICE',      &
          'UW',      &
          'UI',      &
          'VW',      &
          'VI',      &
        'KPAR',      &
    'TS_FOUND'       &
      /

  I = iargc()

  if(I /= 3) then
     print *, "Wrong Number of arguments: ", i
     print *, trim(Usage)
     call exit(1)
  end if

  call getarg(1,InTileFile)
  call getarg(2,InImpRestart)
  call getarg(3,InIntRestart)

  InRestart = trim(InImpRestart)

! Read Output Tile File .til file
! to get the index into the pfafsttater table

  call ReadTileFile_IntLatLon(InTileFile ,Pf,Id,loni,lati,zoom, 0)
  deallocate(Pf,Id)

  nullify(Pf)
  nullify(Id)

  itiles = size(loni)  ! Input  Tile Size

  allocate( varIn(itiles) )
  allocate( varOut(itiles) )
  allocate( varInR8(itiles) )
  allocate( varOutR8(itiles) )
  allocate( TW(itiles) )
  allocate( SW(itiles) )

  call MAPL_NCIOGetFileType(InImpRestart, filetype,rc=rc)

  if (filetype == 0) then

     call InImpFmt%open(InImpRestart,pFIO_READ,rc=rc)
     call InIntFmt%open(InIntRestart,pFIO_READ,rc=rc)
     InIntCfg=InIntFmt%read(rc=rc)
     InImpCfg=InImpFmt%read(rc=rc)

     i = index(InImpRestart,'/',back=.true.)
     OutFileName = "OutData/"//trim(InImpRestart(i+1:))
     status = NF_CREATE (OutFileName, NF_NETCDF4, FOutID)
     status = NF_DEF_DIM(FOutID, 'tile', itiles, TileID)
     status = NF_DEF_DIM(FOutID, 'time'   , 1 , TimID)

     variables => InIntCfg%get_variables()
     var_iter = variables%begin()
     do while (var_iter /= variables%end())
        var_name => var_iter%key()
        if(var_name(1:6) == 'TSKINW') & 
           call MAPL_VarRead(InIntFmt,var_name,TW, __RC__)
        if(var_name(1:6) == 'SSKINW') & 
           call MAPL_VarRead(InIntFmt,var_name,SW, __RC__)
        call var_iter%next()
     enddo   

     variables => InImpCfg%get_variables()
     var_iter = variables%begin()
     do while (var_iter /= variables%end())
 
        var_name => var_iter%key()
        myVariable => var_iter%value()
        var_dimensions => myVariable%get_dimensions()
        ndims = var_dimensions%size()      
        if (ndims == 1) then
           status = NF_DEF_VAR(FOutID, var_name , NF_FLOAT, 1 , TileID , varid)
           attr => myVariable%get_attribute('long_name')
           status = NF_PUT_ATT_TEXT(FOutID, varid, 'long_name', &
                    LEN_TRIM(attr%get_string()),           &
                    trim(attr%get_string())) 
           attr => myVariable%get_attribute('units')
           status = NF_PUT_ATT_TEXT(FOutID, varid, 'units',     &
                    LEN_TRIM(attr%get_string()), trim(attr%get_string()) )  
        else 
           write(*,*)"Import States are all TileOnly:, ",trim(var_name), " is not?"
           stop
        endif
        call var_iter%next()
     enddo

     vname = "SS_FOUND"
     longname = "foundation_salinity_for_interface_layer"
     units = "PSU"
     status = NF_DEF_VAR(FOutID, vname , NF_FLOAT, 1 , TileID , varid)
     status = NF_PUT_ATT_TEXT(FOutID, varid, 'long_name', &
              LEN_TRIM(longname),           &
              trim(longname)) 
     status = NF_PUT_ATT_TEXT(FOutID, varid, 'units',     &
              LEN_TRIM(units), trim(units))        
     status = NF_ENDDEF(FOutID)  

     
     
     variables => InImpCfg%get_variables()
     var_iter = variables%begin()
     do while (var_iter /= variables%end())
        var_name => var_iter%key()
        write(*,*)"Writing ",trim(var_name)
        myVariable => var_iter%value()
        var_dimensions => myVariable%get_dimensions()
        ndims = var_dimensions%size()      
        write(*,*)"Writing ",trim(var_name)
        if (ndims == 1) then
           call MAPL_VarRead(InImpFmt,var_name,varIn, __RC__)
           if(vname(1:8) == 'TS_FOUND') then 
              varOut(:) = TW(:)    
           else
              varOut(:) = varIn(:)    
           endif 
           STATUS = NF_INQ_VARID (FoutID, trim(var_name) ,VarID)
           STATUS = NF_PUT_VARA_REAL(FOutID,VarID, (/1/), (/itiles/), varOut)
        else 
           write(*,*)"Import States are all TileOnly:, ",trim(vname), " is not?"
           stop
        endif
        call var_iter%next()
     enddo
     vname = "SS_FOUND"
     STATUS = NF_INQ_VARID (FoutID, trim(VNAME) ,VarID)
     varOut(:) = SW(:)    
     STATUS = NF_PUT_VARA_REAL(FOutID,VarID, (/1/), (/itiles/), varOut)
     
     status = NF_CLOSE (FoutID)  

  else

      i = index(InImpRestart,'/',back=.true.)
      OutFileName = "OutData/"//trim(InImpRestart(i+1:))

      open(unit=50,FILE=InImpRestart,form='unformatted',&
            status='old',convert='little_endian')

      open(unit=51,FILE=InIntRestart,form='unformatted',&
            status='old',convert='little_endian')

      ! get TW and SW from internal
      ! same for AMIP and coupled 
      read(51) varIn
      read(51)  TW
      read(51)  SW

      call create_salt_import_nc4 (itiles, OutFileName, FOutID)

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

     if(ntot /= size(impNames)) then
        print*, 'ntot ', ntot, ' /= ', size(impNames)
        print*, 'this import restart file is NOT compatible with MERRA-2 tag !!!'
        print*, 'DOUBLE CHECK!!! BYE !!!'
        stop 
     endif

     ! Read and Write Tile or TileTile Data until EOF
     ! ----------------------------------------------
     do n=1,ntot
        nsubtiles = nrecs(n)/itiles
        read (50)  varIn
        if( nsubtiles.eq.1 ) then
           STATUS = NF_INQ_VARID (FoutID, trim(impNames(n)) ,VarID)
           varOut(:) = varIn(:)
           STATUS = NF_PUT_VARA_REAL(FOutID,VarID, (/1/), (/itiles/), varOut)
        else
           write(*,*)"Import States are all TileOnly"
           stop
        endif
     enddo

     vname = "SS_FOUND"
     STATUS = NF_INQ_VARID (FoutID, trim(VNAME) ,VarID)
     varOut(:) = SW(:)    
     STATUS = NF_PUT_VARA_REAL(FOutID,VarID, (/1/), (/itiles/), varOut)

     status = NF_CLOSE (FoutID)  

  end if

  deallocate( varIn, varOut )
  deallocate( varInR8, varOutR8 )
  deallocate( TW, SW )

contains

  SUBROUTINE create_salt_import_nc4 (ntiles, fileName, NCFOutID)

    implicit none

    integer, intent (in)       :: ntiles
    character(*), intent(in)   :: fileName
    integer, intent (inout)    :: NCFOutID    
    integer :: CatchID, TimID, VID, status 

    status = NF_CREATE (filename, NF_NETCDF4, NCFOutID)
    status = NF_DEF_DIM(NCFOutID, 'tile', ntiles, CatchID)
    status = NF_DEF_DIM(NCFOutID, 'time'   , 1 , TimID)
    
    status = NF_DEF_VAR(NCFOutID, 'time'  , NF_DOUBLE, 1 , TimID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
         LEN_TRIM('minutes since  2014-01-01 00:00:00'), trim('minutes since 2014-01-01 00:00:00'))   

    status = NF_DEF_VAR(NCFOutID, 'ALW'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('linearization_of_surface_upwelling_longwave_flux'),&
             trim('linearization_of_surface_upwelling_longwave_flux'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('W m-2'), trim('W m-2'))

    status = NF_DEF_VAR(NCFOutID, 'BLW'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('linearization_of_surface_upwelling_longwave_flux'),&
             trim('linearization_of_surface_upwelling_longwave_flux'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('W m-2 K-1'), trim('W m-2 K-1'))

    status = NF_DEF_VAR(NCFOutID, 'CMATM'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('surface_exchange_coefficient_for_momentum'),&
             trim('surface_exchange_coefficient_for_momentum'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('kg m-2 s-1'), trim('kg m-2 s-1'))

    status = NF_DEF_VAR(NCFOutID, 'CQATM'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('surface_exchange_coefficient_for_moisture'),&
             trim('surface_exchange_coefficient_for_moisture'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('kg m-2 s-1'), trim('kg m-2 s-1'))

    status = NF_DEF_VAR(NCFOutID, 'CTATM'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('surface_exchange_coefficient_for_heat'),&
             trim('surface_exchange_coefficient_for_heat'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('kg m-2 s-1'), trim('kg m-2 s-1'))

    status = NF_DEF_VAR(NCFOutID, 'DEVAP'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('derivative_of_evaporation'),&
             trim('derivative_of_evaporation'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('kg m-2 s-1'), trim('kg m-2 s-1'))

    status = NF_DEF_VAR(NCFOutID, 'DFNIR'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('surface_downwelling_nir_diffuse_flux'),&
             trim('surface_downwelling_nir_diffuse_flux'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('W m-2'), trim('W m-2'))

    status = NF_DEF_VAR(NCFOutID, 'DFPAR'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('surface_downwelling_par_diffuse_flux'),&
             trim('surface_downwelling_par_diffuse_flux'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('W m-2'), trim('W m-2'))

    status = NF_DEF_VAR(NCFOutID, 'DFUVR'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('surface_downwelling_uvr_diffuse_flux'),&
             trim('surface_downwelling_uvr_diffuse_flux'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('W m-2'), trim('W m-2'))

    status = NF_DEF_VAR(NCFOutID, 'DRNIR'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('surface_downwelling_nir_beam_flux'),&
             trim('surface_downwelling_nir_beam_flux'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('W m-2'), trim('W m-2'))

    status = NF_DEF_VAR(NCFOutID, 'DRPAR'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('surface_downwelling_par_beam_flux'),&
             trim('surface_downwelling_par_beam_flux'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('W m-2'), trim('W m-2'))

    status = NF_DEF_VAR(NCFOutID, 'DRUVR'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('surface_downwelling_uvr_beam_flux'),&
             trim('surface_downwelling_uvr_beam_flux'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('W m-2'), trim('W m-2'))

    status = NF_DEF_VAR(NCFOutID, 'DSH'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('derivative_of_upward_sensible_heat_flux'),&
             trim('derivative_of_upward_sensible_heat_flux'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('W m-2'), trim('W m-2'))

    status = NF_DEF_VAR(NCFOutID, 'DZ'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('surface_layer_height'),&
             trim('surface_layer_height'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('m'), trim('m'))

    status = NF_DEF_VAR(NCFOutID, 'EVAP'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('evaporation'),&
             trim('evaporation'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('kg m-2 s-1'), trim('kg m-2 s-1'))

    status = NF_DEF_VAR(NCFOutID, 'FRACICE'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('ice_covered_fraction_of_tile'),&
             trim('ice_covered_fraction_of_tile'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('1'), trim('1'))

    status = NF_DEF_VAR(NCFOutID, 'KPAR'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('PAR_extinction_coefficient'),&
             trim('PAR_extinction_coefficient'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('m-1'), trim('m-1'))

    status = NF_DEF_VAR(NCFOutID, 'LWDNSRF'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('surface_downwelling_longwave_flux'),&
             trim('surface_downwelling_longwave_flux'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('W m-2'), trim('W m-2'))

    status = NF_DEF_VAR(NCFOutID, 'PCU'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('liquid_water_convective_precipitation'),&
             trim('liquid_water_convective_precipitation'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('kg m-2 s-1'), trim('kg m-2 s-1'))

    status = NF_DEF_VAR(NCFOutID, 'PLS'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('liquid_water_large_scale_precipitation'),&
             trim('liquid_water_large_scale_precipitation'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('kg m-2 s-1'), trim('kg m-2 s-1'))

    status = NF_DEF_VAR(NCFOutID, 'PS'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('surface_pressure'),&
             trim('surface_pressure'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('Pa'), trim('Pa'))

    status = NF_DEF_VAR(NCFOutID, 'QA'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('surface_air_specific_humidity'),&
             trim('surface_air_specific_humidity'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('kg kg-1'), trim('kg kg-1'))

    status = NF_DEF_VAR(NCFOutID, 'QHATM'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('effective_surface_specific_humidity'),&
             trim('effective_surface_specific_humidity'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('kg kg-1'), trim('kg kg-1'))

    status = NF_DEF_VAR(NCFOutID, 'SH'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('upward_sensible_heat_flux'),&
             trim('upward_sensible_heat_flux'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('W m-2'), trim('W m-2'))

    status = NF_DEF_VAR(NCFOutID, 'SNO'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('snowfall'),&
             trim('snowfall'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('kg m-2 s-1'), trim('kg m-2 s-1'))

    status = NF_DEF_VAR(NCFOutID, 'TA'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('surface_air_temperature'),&
             trim('surface_air_temperature'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('K'), trim('K'))

    status = NF_DEF_VAR(NCFOutID, 'TAUX'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('eastward_surface_stress'),&
             trim('eastward_surface_stress'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('N m-2'), trim('N m-2'))

    status = NF_DEF_VAR(NCFOutID, 'TAUXBOT'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('eastward_stress_at_base_of_ice'),&
             trim('eastward_stress_at_base_of_ice'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('N m-2'), trim('N m-2'))

    status = NF_DEF_VAR(NCFOutID, 'TAUY'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('northward_surface_stress'),&
             trim('northward_surface_stress'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('N m-2'), trim('N m-2'))

    status = NF_DEF_VAR(NCFOutID, 'TAUYBOT'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('northward_stress_at_base_of_ice'),&
             trim('northward_stress_at_base_of_ice'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('N m-2'), trim('N m-2'))

    status = NF_DEF_VAR(NCFOutID, 'THATM'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('effective_surface_skin_temperature'),&
             trim('effective_surface_skin_temperature'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('K'), trim('K'))

    status = NF_DEF_VAR(NCFOutID, 'TS_FOUND'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('foundation_temperature_for_interface_layer'),&
             trim('foundation_temperature_for_interface_layer'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('K'), trim('K'))

    status = NF_DEF_VAR(NCFOutID, 'SS_FOUND'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('foundation_salinity_for_interface_layer'),&
             trim('foundation_salinity_for_interface_layer'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('PSU'), trim('PSU'))

    status = NF_DEF_VAR(NCFOutID, 'UHATM'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('effective_surface_zonal_velocity'),&
             trim('effective_surface_zonal_velocity'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('m s-1'), trim('m s-1'))

    status = NF_DEF_VAR(NCFOutID, 'UI'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('zonal_velocity_of_surface_ice'),&
             trim('zonal_velocity_of_surface_ice'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('m s-1'), trim('m s-1'))

    status = NF_DEF_VAR(NCFOutID, 'UU'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('surface_wind_speed'),&
             trim('surface_wind_speed'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('m s-1'), trim('m s-1'))

    status = NF_DEF_VAR(NCFOutID, 'UW'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('zonal_velocity_of_surface_water'),&
             trim('zonal_velocity_of_surface_water'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('m s-1'), trim('m s-1'))

    status = NF_DEF_VAR(NCFOutID, 'UWINDLMTILE'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('levellm_uwind'),&
             trim('levellm_uwind'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('m s-1'), trim('m s-1'))

    status = NF_DEF_VAR(NCFOutID, 'VHATM'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('effective_surface_meridional_velocity'),&
             trim('effective_surface_meridional_velocity'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('m s-1'), trim('m s-1'))

    status = NF_DEF_VAR(NCFOutID, 'VI'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('meridional_velocity_of_surface_ice'),&
             trim('meridional_velocity_of_surface_ice'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('m s-1'), trim('m s-1'))

    status = NF_DEF_VAR(NCFOutID, 'VW'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('meridional_velocity_of_surface_water'),&
             trim('meridional_velocity_of_surface_water'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('m s-1'), trim('m s-1'))

    status = NF_DEF_VAR(NCFOutID, 'VWINDLMTILE'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',& 
             LEN_TRIM('levellm_vwind'),&
             trim('levellm_vwind'))
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
             LEN_TRIM('m s-1'), trim('m s-1'))

    
    status = NF_ENDDEF(NCFOutID)  

  END SUBROUTINE create_salt_import_nc4

end program SaltImpConverter

