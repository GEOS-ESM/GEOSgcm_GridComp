#define I_AM_MAIN
#include "MAPL_Generic.h"
program  mk_CatchRestarts

!  $Id: 

  use MAPL
  use mk_restarts_getidsMod, only: GetIDs,ReadTileFile_RealLatLon
  use gFTL2_StringVector

  implicit none
  include 'mpif.h'
  ! initialize to non-MPI values

  integer  :: myid=0, numprocs=1, mpierr, mpistatus(MPI_STATUS_SIZE)  
  logical  :: root_proc=.true.

  character*256 :: Usage="mk_CatchRestarts OutTileFile InTileFile InRestart SURFLAY <OutIsOld>"
  character*256 :: OutTileFile
  character*256 :: InTileFile
  character*256 :: InRestart
  character*256 :: OutType
  character*256 :: arg(6)

  integer :: i, k, iargc, n, ntiles,ntiles_in, nplus, req
  integer, pointer  :: Id(:), tid_in (:)
  real,    pointer  :: loni(:),lono(:), lati(:), lato(:)
  real              :: SURFLAY
  integer, allocatable            :: low_ind(:), upp_ind(:), nt_local (:), Id_loc (:)
  real   , pointer, dimension (:) :: LATT, LONN
  logical :: OutIsOld, havedata
  character*256, parameter        :: DataDir="OutData/clsm/"
  real                            :: min_lon, max_lon, min_lat, max_lat
  logical, allocatable, dimension(:)  :: mask
  integer, allocatable, dimension (:) :: sub_tid
  real   , allocatable, dimension (:) :: sub_lon, sub_lat
  integer :: status

  call init_MPI()

!---------------------------------------------------------------------------
  
  I = iargc()

  if( I<4 .or. I>5 ) then
     print *, "Wrong Number of arguments: ", i
     print *, trim(Usage)
     call exit(1)
  end if
 
  do n=1,I
  call getarg(n,arg(n))
  enddo
  read(arg(1),'(a)') OutTileFile
  read(arg(2),'(a)')  InTileFile
  read(arg(3),'(a)')  InRestart
  read(arg(4),*)  SURFLAY

  if(I==5) then
     call getarg(6,OutType)
     OutIsOld = trim(OutType)=="OutIsOld"
  else
     OutIsOld = .false.
  endif

  if (SURFLAY.ne.20 .and. SURFLAY.ne.50) then
     print *, "You must supply a valid SURFLAY value:"
     print *, "(Ganymed-3 and earlier) SURFLAY=20.0 for Old Soil Params"
     print *, "(Ganymed-4 and later  ) SURFLAY=50.0 for New Soil Params"
     call exit(2)
  end if

  inquire(file=trim(DataDir)//"mosaic_veg_typs_fracs",exist=havedata)

  if (root_proc) then

     ! Read Output/Input  .til files
     call ReadTileFile_RealLatLon(OutTileFile, ntiles, lono, lato)  
     call ReadTileFile_RealLatLon(InTileFile,ntiles_in,loni,lati)
     allocate(Id (ntiles))
 !    allocate(mask   (ntiles_in))
 !    allocate(tid_in (ntiles_in))
 !    do n = 1, NTILES_IN
 !       tid_in (n) = n
 !    end do

  endif

  if (havedata) then
     if (root_proc) call read_and_write_rst (NTILES, SURFLAY, OutIsOld, NTILES, __RC__)
  else

     call MPI_BCAST  (ntiles   , 1, MPI_INTEGER, 0,MPI_COMM_WORLD, mpierr)
     call MPI_BCAST  (ntiles_in, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, mpierr)
     
     allocate(low_ind (   numprocs))
     allocate(upp_ind (   numprocs))
     allocate(nt_local(   numprocs))
     
     low_ind (:)    = 1
     upp_ind (:)    = NTILES       
     nt_local(:)    = NTILES 
     
     if (numprocs > 1) then      
        do i = 1, numprocs - 1
           upp_ind(i)   = low_ind(i) + (ntiles/numprocs) - 1 
           low_ind(i+1) = upp_ind(i) + 1
           nt_local(i)  = upp_ind(i) - low_ind(i) + 1
        end do
        nt_local(numprocs) = upp_ind(numprocs) - low_ind(numprocs) + 1
     endif

     ! Get intile lat/lon
     
!     do i = 2, numprocs
!        if (i -1 == myid) then
!           ! receive ntiles_in in the block
!           call MPI_RECV(ntiles_in, 1, MPI_INTEGER,0,999,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
!            ! ALLOCATE
!           allocate (loni   (1:NTILES_IN))
!           allocate (lati   (1:NTILES_IN))
!           allocate (tid_in (1:NTILES_IN))
!
!           ! RECEIVE LAT/LON IN
!           call MPI_RECV(tid_in, ntiles_in, MPI_INTEGER,0,998,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
!           call MPI_RECV(loni  , ntiles_in, MPI_REAL   ,0,997,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
!           call MPI_RECV(lati  , ntiles_in, MPI_REAL   ,0,996,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
!
!        else if (myid == 0) then
!
!           ! Send local ntiles_in 
!           
!           min_lon = MAX(MINVAL(lono (low_ind(i) : upp_ind(i))) - 5, -180.)
!           max_lon = MIN(MAXVAL(lono (low_ind(i) : upp_ind(i))) + 5,  180.)
!           min_lat = MAX(MINVAL(lato (low_ind(i) : upp_ind(i))) - 5,  -90.)
!           max_lat = MIN(MAXVAL(lato (low_ind(i) : upp_ind(i))) + 5,   90.) 
!           mask = .false.
!           mask =  ((lati >= min_lat .and. lati <= max_lat).and.(loni >= min_lon .and. loni <= max_lon))
!           nplus =  count(mask = mask)
!
!           call MPI_ISend(NPLUS ,1,MPI_INTEGER,i-1,999,MPI_COMM_WORLD,req,mpierr)
!           call MPI_WAIT (req,MPI_STATUS_IGNORE,mpierr)
!
!
!           ! SEND LAT/LON IN
!           allocate (sub_tid (1:nplus))
!           allocate (sub_lon (1:nplus))
!           allocate (sub_lat (1:nplus))
!
!           sub_tid = PACK (tid_in , mask= mask) 
!           sub_lon = PACK (loni   , mask= mask)
!           sub_lat = PACK (lati   , mask= mask)
!
!           call MPI_ISend(sub_tid, nplus,MPI_INTEGER,i-1,998,MPI_COMM_WORLD,req,mpierr)
!           call MPI_WAIT (req,MPI_STATUS_IGNORE,mpierr)     
!           call MPI_ISend(sub_lon, nplus,MPI_REAL   ,i-1,997,MPI_COMM_WORLD,req,mpierr)
!           call MPI_WAIT (req,MPI_STATUS_IGNORE,mpierr) 
!           call MPI_ISend(sub_lat, nplus,MPI_REAL   ,i-1,996,MPI_COMM_WORLD,req,mpierr)
!           call MPI_WAIT (req,MPI_STATUS_IGNORE,mpierr)
!           deallocate (sub_tid,sub_lon,sub_lat)
!        endif           
!     end do

     ! Get out tile lat/lots from root

     allocate (id_loc  (nt_local (myid + 1)))
     allocate (lonn    (nt_local (myid + 1)))
     allocate (latt    (nt_local (myid + 1)))  
     
     do i = 1, numprocs
        if((I == 1).and.(myid == 0)) then
           lonn(:) = lono(low_ind(i) : upp_ind(i))
           latt(:) = lato(low_ind(i) : upp_ind(i))
        else if (I > 1) then
           if(I-1 == myid) then
              ! receiving from root
              call MPI_RECV(lonn,nt_local(i) , MPI_REAL, 0,995,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
              call MPI_RECV(latt,nt_local(i) , MPI_REAL, 0,994,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
           else if (myid == 0) then
              ! root sends
              call MPI_ISend(lono(low_ind(i) : upp_ind(i)),nt_local(i),MPI_REAL,i-1,995,MPI_COMM_WORLD,req,mpierr)
              call MPI_WAIT (req,MPI_STATUS_IGNORE,mpierr)    
              call MPI_ISend(lato(low_ind(i) : upp_ind(i)),nt_local(i),MPI_REAL,i-1,994,MPI_COMM_WORLD,req,mpierr)
              call MPI_WAIT (req,MPI_STATUS_IGNORE,mpierr) 
           endif
        endif
     end do

!     call MPI_SCATTERV (                    &
!          lono,nt_local,low_ind-1,MPI_real, &
!          lonn,size(lonn),MPI_real  , &
!          0,MPI_COMM_WORLD, mpierr )
!
!     call MPI_SCATTERV (                    &
!          lato,nt_local,low_ind-1,MPI_real, &
!          latt,size(latt),MPI_real  , &
!          0,MPI_COMM_WORLD, mpierr )

     if(myid > 0) allocate (loni (1:NTILES_IN))
     if(myid > 0) allocate (lati (1:NTILES_IN))
          
     call MPI_BCAST(loni,ntiles_in,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
     call MPI_BCAST(lati,ntiles_in,MPI_REAL,0,MPI_COMM_WORLD,mpierr)  

     allocate(tid_in (ntiles_in))
     do n = 1, NTILES_IN
        tid_in (n) = n
     end do

     call GetIds(loni,lati,lonn,latt,Id_loc, tid_in)
     call MPI_Barrier(MPI_COMM_WORLD, mpierr)
!     call MPI_GATHERV( &
!          id_loc (:), nt_local(myid+1), MPI_real, &
!          id, nt_local,low_ind-1, MPI_real, &
!          0, MPI_COMM_WORLD, mpierr )

     do i = 1, numprocs
        if((I == 1).and.(myid == 0)) then
           id(low_ind(i) : upp_ind(i)) = Id_loc(:)
        else if (I > 1) then
           if(I-1 == myid) then
              ! send to root
              call MPI_ISend(id_loc,nt_local(i),MPI_INTEGER,0,993,MPI_COMM_WORLD,req,mpierr)
              call MPI_WAIT (req,MPI_STATUS_IGNORE,mpierr)    
           else if (myid == 0) then
              ! root receives
              call MPI_RECV(id(low_ind(i) : upp_ind(i)),nt_local(i) , MPI_INTEGER, i-1,993,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
           endif
        endif
     end do
          
     deallocate (loni,lati,lonn,latt, tid_in)
     call MPI_Barrier(MPI_COMM_WORLD, mpierr)
 
     if (root_proc) call read_and_write_rst (NTILES, SURFLAY, OutIsOld, NTILES_IN, id, __RC__)
  
  endif

  call MPI_BARRIER( MPI_COMM_WORLD, mpierr)
  call MPI_FINALIZE(mpierr)   

contains

  SUBROUTINE read_and_write_rst (NTILES,  SURFLAY, OutIsOld, NTILES_IN, idi, rc)

    implicit none
    real,    intent (in) :: SURFLAY
    logical, intent (in) :: OutIsOld
    integer, intent (in) :: NTILES, NTILES_IN
    integer, pointer, dimension(:), optional, intent (in) :: idi    
    integer, optional, intent(out) :: rc
    logical       :: havedata, NewLand
    character(len=256), parameter :: Names(29) = &
         (/'BF1   ','BF2   ','BF3   ','VGWMAX','CDCR1 ', &
         'CDCR2 ','PSIS  ','BEE   ','POROS ','WPWET ', &
         'COND  ','GNU   ','ARS1  ','ARS2  ','ARS3  ', &
         'ARA1  ','ARA2  ','ARA3  ','ARA4  ','ARW1  ', &
         'ARW2  ','ARW3  ','ARW4  ','TSA1  ','TSA2  ', & 
         'TSB1  ','TSB2  ','ATAU  ','BTAU  '/)

    integer, pointer  :: ity(:) 
    real, allocatable :: BF1(:),   BF2(:),   BF3(:),  VGWMAX(:)
    real, allocatable :: CDCR1(:), CDCR2(:), PSIS(:), BEE(:) 
    real, allocatable :: POROS(:), WPWET(:), COND(:), GNU(:)
    real, allocatable :: ARS1(:),  ARS2(:),  ARS3(:)
    real, allocatable :: ARA1(:),  ARA2(:),  ARA3(:), ARA4(:)
    real, allocatable :: ARW1(:),  ARW2(:),  ARW3(:), ARW4(:)
    real, allocatable :: TSA1(:),  TSA2(:),  TSB1(:), TSB2(:)
    real, allocatable :: ATAU2(:), BTAU2(:), DP2BR(:), rity(:)
    
    real              :: zdep1, zdep2, zdep3, zmet, term1, term2, rdum
    real, allocatable :: var1(:),var2(:,:)
    character*256     :: vname
    character*256     :: OutFileName
    integer           :: i, n, j,k,ncatch,idum
    logical,allocatable :: written(:)
    integer             :: ndims,filetype
    integer             :: dimSizes(3),nVars
    logical             :: file_exists
    integer, pointer    :: Ido(:), idx(:), id(:)
    logical             :: InIsOld
    type(NetCDF4_Fileformatter) :: InFmt,OutFmt,CatchFmt
    type(FileMetadata) :: InCfg,OutCfg
    type(StringVariableMap), pointer :: variables
    type(Variable), pointer :: myVariable
    type(StringVariableMapIterator) :: var_iter
    character(len=:), pointer :: var_name,dname
    type(StringVector), pointer :: var_dimensions
    integer :: dim1, dim2
    character(256) :: Iam = "read_and_write_rst"
    integer :: status

    print *, 'SURFLAY: ',SURFLAY 

    inquire(file=trim(DataDir)//"mosaic_veg_typs_fracs",exist=havedata)
    inquire(file=trim(DataDir)//"CLM_veg_typs_fracs"   ,exist=NewLand )

    print *, 'havedata = ',havedata

    call MAPL_NCIOGetFileType(InRestart, filetype,__RC__)

    if (filetype == 0) then
      
       call InFmt%open(InRestart,pFIO_READ,__RC__)
       InCfg=InFmt%read(__RC__)
       call MAPL_IOChangeRes(InCfg,OutCfg,(/'tile'/),(/ntiles/),__RC__)
       i = index(InRestart,'/',back=.true.)
       OutFileName = "OutData/"//trim(InRestart(i+1:))
       call OutFmt%create(OutFileName,__RC__)
       call OutFmt%write(OutCfg,__RC__)
       call MAPL_IOCountNonDimVars(OutCfg,nvars,__RC__)

       allocate(written(nvars))
       written=.false. 
       
    else
       
       open(unit=50,FILE=InRestart,form='unformatted',&
            status='old',convert='little_endian')
       
       do i=1,58
          read(50,end=2001)
       end do
2001   continue
       InIsOld = I==59
       
       rewind(50)
       
       i = index(InRestart,'/',back=.true.)
       
       open(unit=40,FILE="OutData/"//trim(InRestart(i+1:)),form='unformatted',&
            status='unknown',convert='little_endian')
       
    end if
        
    HAVE: if(havedata) then
       
       print *,'Working from Sariths data pretiled for this resolution'
       
       ! Get number of catchments
       
       open(unit=22, &
            file=trim(DataDir)//"catchment.def",status='old',form='formatted')
       
       read (22, *) ncatch
       
       close(22)
       
       if(ncatch==ntiles) then
          print *, "Read ",Ncatch," land tiles."
          allocate (ido (ntiles))
          do i=1,ncatch
             ido(i) = i
          enddo
       else
          print *, "Number of tiles in data, ",Ncatch," does not match number in til file ", size(Ido)
          call exit(1)
       endif
       
       allocate(ity(ncatch),rity(ncatch))
       allocate (   BF1(ncatch),    BF2 (ncatch),     BF3(ncatch)  )
       allocate (VGWMAX(ncatch),   CDCR1(ncatch),   CDCR2(ncatch)  ) 
       allocate (  PSIS(ncatch),     BEE(ncatch),   POROS(ncatch)  ) 
       allocate ( WPWET(ncatch),    COND(ncatch),     GNU(ncatch)  )
       allocate (  ARS1(ncatch),    ARS2(ncatch),    ARS3(ncatch)  )
       allocate (  ARA1(ncatch),    ARA2(ncatch),    ARA3(ncatch)  )
       allocate (  ARA4(ncatch),    ARW1(ncatch),    ARW2(ncatch)  )
       allocate (  ARW3(ncatch),    ARW4(ncatch),    TSA1(ncatch)  )
       allocate (  TSA2(ncatch),    TSB1(ncatch),    TSB2(ncatch)  )
       allocate ( ATAU2(ncatch),    BTAU2(ncatch),   DP2BR(ncatch)  )
       
       inquire(file = trim(DataDir)//'/catch_params.nc4', exist=file_exists)
       
       if(file_exists) then
          print *,'FILE FORMAT FOR LAND BCS IS NC4'
          call CatchFmt%open(trim(DataDir)//'/catch_params.nc4',pFIO_Read, __RC__) 
          call MAPL_VarRead ( catchFmt ,'OLD_ITY', rity, __RC__)
          call MAPL_VarRead ( catchFmt ,'ARA1', ARA1, __RC__)
          call MAPL_VarRead ( catchFmt ,'ARA2', ARA2, __RC__)
          call MAPL_VarRead ( catchFmt ,'ARA3', ARA3, __RC__)
          call MAPL_VarRead ( catchFmt ,'ARA4', ARA4, __RC__)
          call MAPL_VarRead ( catchFmt ,'ARS1', ARS1, __RC__)
          call MAPL_VarRead ( catchFmt ,'ARS2', ARS2, __RC__)
          call MAPL_VarRead ( catchFmt ,'ARS3', ARS3, __RC__)
          call MAPL_VarRead ( catchFmt ,'ARW1', ARW1, __RC__)
          call MAPL_VarRead ( catchFmt ,'ARW2', ARW2, __RC__)
          call MAPL_VarRead ( catchFmt ,'ARW3', ARW3, __RC__)
          call MAPL_VarRead ( catchFmt ,'ARW4', ARW4, __RC__)

          if( SURFLAY.eq.20.0 ) then
             call MAPL_VarRead ( catchFmt ,'ATAU2', ATAU2, __RC__)
             call MAPL_VarRead ( catchFmt ,'BTAU2', BTAU2, __RC__)
          endif
          
          if( SURFLAY.eq.50.0 ) then
             call MAPL_VarRead ( catchFmt ,'ATAU5', ATAU2, __RC__)
             call MAPL_VarRead ( catchFmt ,'BTAU5', BTAU2, __RC__)
          endif
          
          call MAPL_VarRead ( catchFmt ,'PSIS', PSIS, __RC__)
          call MAPL_VarRead ( catchFmt ,'BEE', BEE, __RC__)
          call MAPL_VarRead ( catchFmt ,'BF1', BF1, __RC__)
          call MAPL_VarRead ( catchFmt ,'BF2', BF2, __RC__)
          call MAPL_VarRead ( catchFmt ,'BF3', BF3, __RC__)
          call MAPL_VarRead ( catchFmt ,'TSA1', TSA1, __RC__)
          call MAPL_VarRead ( catchFmt ,'TSA2', TSA2, __RC__)
          call MAPL_VarRead ( catchFmt ,'TSB1', TSB1, __RC__)
          call MAPL_VarRead ( catchFmt ,'TSB2', TSB2, __RC__)
          call MAPL_VarRead ( catchFmt ,'COND', COND, __RC__)
          call MAPL_VarRead ( catchFmt ,'GNU', GNU, __RC__)
          call MAPL_VarRead ( catchFmt ,'WPWET', WPWET, __RC__)
          call MAPL_VarRead ( catchFmt ,'DP2BR', DP2BR, __RC__)
          call MAPL_VarRead ( catchFmt ,'POROS', POROS, __RC__)
          call catchFmt%close(__RC__)

       else
          open(unit=21, file=trim(DataDir)//"mosaic_veg_typs_fracs",status='old',form='formatted')
          open(unit=22, file=trim(DataDir)//'bf.dat'               ,form='formatted')
          open(unit=23, file=trim(DataDir)//'soil_param.dat'       ,form='formatted')
          open(unit=24, file=trim(DataDir)//'ar.new'               ,form='formatted')
          open(unit=25, file=trim(DataDir)//'ts.dat'               ,form='formatted')
          open(unit=26, file=trim(DataDir)//'tau_param.dat'        ,form='formatted')
          
          do n=1,ncatch
             read (21,*) I, j, ITY(N)
             read (22, *) i,j, GNU(n), BF1(n), BF2(n), BF3(n)
             
             read (23, *) i,j, idum, idum, BEE(n), PSIS(n),&
                  POROS(n), COND(n), WPWET(n), DP2BR(n)
             
             read (24, *) i,j, rdum, ARS1(n), ARS2(n), ARS3(n),          &
                  ARA1(n), ARA2(n), ARA3(n), ARA4(n), &
                  ARW1(n), ARW2(n), ARW3(n), ARW4(n)
             
             read (25, *) i,j, rdum, TSA1(n), TSA2(n), TSB1(n), TSB2(n)
             
             if( SURFLAY.eq.20.0 ) read (26, *) i,j, ATAU2(n), BTAU2(n), rdum, rdum   ! for old soil params
             if( SURFLAY.eq.50.0 ) read (26, *) i,j, rdum , rdum, ATAU2(n), BTAU2(n)  ! for new soil params
          end do
          
          rity = float(ity)
          CLOSE (21, STATUS = 'KEEP')
          CLOSE (22, STATUS = 'KEEP')
          CLOSE (23, STATUS = 'KEEP')
          CLOSE (24, STATUS = 'KEEP')
          CLOSE (25, STATUS = 'KEEP')
          
       endif
       
       do n=1,ncatch
          
          zdep2=1000.
          zdep3=amax1(1000.,DP2BR(n))
          
          if (zdep2 > 0.75*zdep3) then
             zdep2  =  0.75*zdep3              
          end if
          
          zdep1=20.
          zmet=zdep3/1000.
          
          term1=-1.+((PSIS(n)-zmet)/PSIS(n))**((BEE(n)-1.)/BEE(n))
          term2=PSIS(n)*BEE(n)/(BEE(n)-1)
          
          VGWMAX(n) = POROS(n)*zdep2   
          CDCR1(n)  = 1000.*POROS(n)*(zmet-(-term2*term1))   
          CDCR2(n)  = (1.-WPWET(n))*POROS(n)*zdep3
       enddo
       
       
       if (filetype /=0) then
          do i=1,30
             read(50)
          enddo
       end if

       idx => ido
    
    else
       
       print *,'Working from restarts alone'
       
       ncatch = NTILES_IN
       
       allocate (   rity(ncatch))
       allocate (   BF1(ncatch),    BF2 (ncatch),     BF3(ncatch)  )
       allocate (VGWMAX(ncatch),   CDCR1(ncatch),   CDCR2(ncatch)  ) 
       allocate (  PSIS(ncatch),     BEE(ncatch),   POROS(ncatch)  ) 
       allocate ( WPWET(ncatch),    COND(ncatch),     GNU(ncatch)  )
       allocate (  ARS1(ncatch),    ARS2(ncatch),    ARS3(ncatch)  )
       allocate (  ARA1(ncatch),    ARA2(ncatch),    ARA3(ncatch)  )
       allocate (  ARA4(ncatch),    ARW1(ncatch),    ARW2(ncatch)  )
       allocate (  ARW3(ncatch),    ARW4(ncatch),    TSA1(ncatch)  )
       allocate (  TSA2(ncatch),    TSB1(ncatch),    TSB2(ncatch)  )
       allocate ( ATAU2(ncatch),    BTAU2(ncatch),   DP2BR(ncatch)  )
       
       if (filetype == 0) then
          
          call MAPL_VarRead(InFmt,names(1),BF1, __RC__)
          call MAPL_VarRead(InFmt,names(2),BF2, __RC__)
          call MAPL_VarRead(InFmt,names(3),BF3, __RC__)
          call MAPL_VarRead(InFmt,names(4),VGWMAX, __RC__)
          call MAPL_VarRead(InFmt,names(5),CDCR1, __RC__)
          call MAPL_VarRead(InFmt,names(6),CDCR2, __RC__)
          call MAPL_VarRead(InFmt,names(7),PSIS, __RC__)
          call MAPL_VarRead(InFmt,names(8),BEE, __RC__)
          call MAPL_VarRead(InFmt,names(9),POROS, __RC__)
          call MAPL_VarRead(InFmt,names(10),WPWET, __RC__)
          
          call MAPL_VarRead(InFmt,names(11),COND, __RC__)
          call MAPL_VarRead(InFmt,names(12),GNU, __RC__)
          call MAPL_VarRead(InFmt,names(13),ARS1, __RC__)
          call MAPL_VarRead(InFmt,names(14),ARS2, __RC__)
          call MAPL_VarRead(InFmt,names(15),ARS3, __RC__)
          call MAPL_VarRead(InFmt,names(16),ARA1, __RC__)
          call MAPL_VarRead(InFmt,names(17),ARA2, __RC__)
          call MAPL_VarRead(InFmt,names(18),ARA3, __RC__)
          call MAPL_VarRead(InFmt,names(19),ARA4, __RC__)
          call MAPL_VarRead(InFmt,names(20),ARW1, __RC__)
          
          call MAPL_VarRead(InFmt,names(21),ARW2, __RC__)
          call MAPL_VarRead(InFmt,names(22),ARW3, __RC__)
          call MAPL_VarRead(InFmt,names(23),ARW4, __RC__)
          call MAPL_VarRead(InFmt,names(24),TSA1, __RC__)
          call MAPL_VarRead(InFmt,names(25),TSA2, __RC__)
          call MAPL_VarRead(InFmt,names(26),TSB1, __RC__)
          call MAPL_VarRead(InFmt,names(27),TSB2, __RC__)
          call MAPL_VarRead(InFmt,names(28),ATAU2, __RC__)
          call MAPL_VarRead(InFmt,names(29),BTAU2, __RC__)
          call MAPL_VarRead(InFmt,'OLD_ITY',rITY, __RC__)
          
       else
          
          read(50) BF1
          read(50) BF2
          read(50) BF3
          read(50) VGWMAX
          read(50) CDCR1
          read(50) CDCR2
          read(50) PSIS
          read(50) BEE
          read(50) POROS
          read(50) WPWET
          
          read(50) COND
          read(50) GNU
          read(50) ARS1
          read(50) ARS2
          read(50) ARS3
          read(50) ARA1
          read(50) ARA2
          read(50) ARA3
          read(50) ARA4
          read(50) ARW1
          
          read(50) ARW2
          read(50) ARW3
          read(50) ARW4
          read(50) TSA1
          read(50) TSA2
          read(50) TSB1
          read(50) TSB2
          read(50) ATAU2
          read(50) BTAU2
          read(50) rITY
          
       end if
 
       idx => idi

    endif HAVE

    if (filetype == 0) then
       call MAPL_VarWrite(OutFmt,names(1),BF1(Idx))
       call MAPL_VarWrite(OutFmt,names(2),BF2(Idx))
       call MAPL_VarWrite(OutFmt,names(3),BF3(Idx))
       call MAPL_VarWrite(OutFmt,names(4),VGWMAX(Idx))
       call MAPL_VarWrite(OutFmt,names(5),CDCR1(Idx))
       call MAPL_VarWrite(OutFmt,names(6),CDCR2(Idx))
       call MAPL_VarWrite(OutFmt,names(7),PSIS(Idx))
       call MAPL_VarWrite(OutFmt,names(8),BEE(Idx))
       call MAPL_VarWrite(OutFmt,names(9),POROS(Idx))
       call MAPL_VarWrite(OutFmt,names(10),WPWET(Idx))
       call MAPL_VarWrite(OutFmt,names(11),COND(Idx))
       call MAPL_VarWrite(OutFmt,names(12),GNU(Idx))
       call MAPL_VarWrite(OutFmt,names(13),ARS1(Idx))
       call MAPL_VarWrite(OutFmt,names(14),ARS2(Idx))
       call MAPL_VarWrite(OutFmt,names(15),ARS3(Idx))
       call MAPL_VarWrite(OutFmt,names(16),ARA1(Idx))
       call MAPL_VarWrite(OutFmt,names(17),ARA2(Idx))
       call MAPL_VarWrite(OutFmt,names(18),ARA3(Idx))
       call MAPL_VarWrite(OutFmt,names(19),ARA4(Idx))
       call MAPL_VarWrite(OutFmt,names(20),ARW1(Idx))
       call MAPL_VarWrite(OutFmt,names(21),ARW2(Idx))
       call MAPL_VarWrite(OutFmt,names(22),ARW3(Idx))
       call MAPL_VarWrite(OutFmt,names(23),ARW4(Idx))
       call MAPL_VarWrite(OutFmt,names(24),TSA1(Idx))
       call MAPL_VarWrite(OutFmt,names(25),TSA2(Idx))
       call MAPL_VarWrite(OutFmt,names(26),TSB1(Idx))
       call MAPL_VarWrite(OutFmt,names(27),TSB2(Idx))
       call MAPL_VarWrite(OutFmt,names(28),ATAU2(Idx))
       call MAPL_VarWrite(OutFmt,names(29),BTAU2(Idx))
       call MAPL_VarWrite(OutFmt,'OLD_ITY',rity(Idx))


       call MAPL_IOCountNonDimVars(InCfg,nvars)

       variables => InCfg%get_variables()
       var_iter = variables%ftn_begin()
       i = 0
       do while (var_iter /= variables%ftn_end())
          call var_iter%next()

          var_name => var_iter%first()
          i=i+1
          do j=1,29
             if ( trim(var_name) == trim(names(j)) ) written(i) = .true.
          enddo
          if (trim(var_name) == "OLD_ITY" ) written(i) = .true.

       enddo

       variables => InCfg%get_variables()
       var_iter = variables%ftn_begin()
       n=0
       allocate(var1(NTILES_IN))
       do while (var_iter /= variables%ftn_end())
          call var_iter%next()

          var_name => var_iter%first()
          myVariable => var_iter%second()

          if (.not.InCfg%is_coordinate_variable(var_name)) then

             n=n+1
             if (.not.written(n) ) then

                var_dimensions => myVariable%get_dimensions()

                ndims = var_dimensions%size()

                if (ndims == 1) then
                   call MAPL_VarRead(InFmt,var_name,var1, __RC__)
                   call MAPL_VarWrite(OutFmt,var_name,var1(idx))
                else if (ndims == 2) then

                   dname => myVariable%get_ith_dimension(2)
                   dim1=InCfg%get_dimension(dname)
                   do j=1,dim1
                      call MAPL_VarRead(InFmt,var_name,var1,offset1=j, __RC__)
                      call MAPL_VarWrite(OutFmt,var_name,var1(idx),offset1=j)
                   enddo
                else if (ndims == 3) then

                   dname => myVariable%get_ith_dimension(2)
                   dim1=InCfg%get_dimension(dname)
                   dname => myVariable%get_ith_dimension(3)
                   dim2=InCfg%get_dimension(dname)
                   do k=1,dim2
                      do j=1,dim1
                         call MAPL_VarRead(InFmt,var_name,var1,offset1=j,offset2=k, __RC__)
                         call MAPL_VarWrite(OutFmt,var_name,var1(idx),offset1=j,offset2=k)
                      enddo
                   enddo

                end if

             end if
          end if

       enddo

    else
       
       write(40) BF1(Idx)
       write(40) BF2(Idx)
       write(40) BF3(Idx)
       write(40) VGWMAX(Idx)
       write(40) CDCR1(Idx)
       write(40) CDCR2(Idx)
       write(40) PSIS(Idx)
       write(40) BEE(Idx)
       write(40) POROS (Idx)
       write(40) WPWET(Idx)
       write(40) COND(Idx)
       write(40) GNU(Idx)
       write(40) ARS1(Idx)
       write(40) ARS2(Idx)
       write(40) ARS3(Idx)
       write(40) ARA1(Idx)
       write(40) ARA2(Idx)
       write(40) ARA3(Idx)
       write(40) ARA4(Idx)
       write(40) ARW1(Idx)
       write(40) ARW2(Idx)
       write(40) ARW3(Idx)
       write(40) ARW4(Idx)
       write(40) TSA1(Idx)
       write(40) TSA2(Idx)
       write(40) TSB1(Idx)
       write(40) TSB2(Idx)
       write(40) ATAU2(Idx)
       write(40) BTAU2(Idx)
       write(40) rITY(Idx)
       
       
       allocate(var1(NTILES_IN))
       allocate(var2(NTILES_IN,4))
       
       ! TC QC
       
       do n=1,2
          read (50) var2
          write(40) ((var2(idx(i),j),i=1,ntiles),j=1,4)
       end do
       
       !CAPAC CATDEF RZEXC SRFEXC ... SNDZN3
       
       do n=1,20
          read (50) var1
          write(40) var1(Idx)
       enddo
       
       ! CH CM CQ FR
       
       do n=1,4
          read (50) var2
          write(40) ((var2(idx(i),j),i=1,ntiles),j=1,4)
       end do
       
       ! These are the 2 prev/next pairs that dont are not
       ! in the internal in fortuna-2_0 and later. Earlier the 
       ! record are there, but their values are not needed, since
       ! they are initialized on start-up.
       
       if(InIsOld) then
          do n=1,4
             read (50)
          enddo
       endif
       
       if(OutIsOld) then
          var1 = 0.0
          do n=1,4
             write(40) (var1(idx(i)),i = 1, ntiles)
          end do
       endif
       
       ! WW
       
       read (50) var2
       write(40) ((var2(idx(i),j),i=1,ntiles),j=1,4)
    end if
    if (present(rc)) rc =0
    !_RETURN(_SUCCESS)
  END SUBROUTINE read_and_write_rst

  ! *****************************************************************************
  
  subroutine init_MPI()
    
    ! initialize MPI
    
    call MPI_INIT(mpierr)
    
    call MPI_COMM_RANK( MPI_COMM_WORLD, myid, mpierr )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, mpierr )

    if (myid .ne. 0)  root_proc = .false.
        
!    write (*,*) "MPI process ", myid, " of ", numprocs, " is alive"    
!    write (*,*) "MPI process ", myid, ": root_proc=", root_proc

  end subroutine init_MPI

end program mk_CatchRestarts

