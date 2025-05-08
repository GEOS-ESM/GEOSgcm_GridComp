#define I_AM_MAIN
#include "MAPL_Generic.h"

program  mk_CatchCNRestarts

!  Usage : mk_CatchCNRestarts OutTileFile InTileFile InRestart SURFLAY RestartTime
!  Version 1 : Sarith Mahanama 
!              sarith.p.mahanama@nasa.gov (Feb 19, 2016) 
!              The program follows the same nearest neighbor based procedure, as in mk_CatchRestarts.F90, 
!                     to regrid hydrological variables and BCs-based parameters. The algorithm developed 
!                     by Greg Walker (~gkwalker/geos5/convert_offline_cn_restart.f90) to regrid carbon 
!                     variables that looks for a neighbor with a similar vegetation type was modified 
!                     to improve efficiency (in subroutine regrid_carbon_vars). The two main 
!                     modifications in this implementation include: (1) instead looping over the globe, 
!                     it starts from a 10 x 10 window and zoom out until a similar type appears, 
!                     (2) uses MPI enabling parrellel computation.
!  Version 2 : Sarith Mahanama (Oct 12, 2016)
!            (1) updated to read both carbon and hydrological variables more recent SMAP M09 simulation from Fanwei.
!            (2) added subroutine reorder_LDASsa_rst
! The program produces catchcn_internal_rst in nc4 format for any user specified AGCM grid resolution.

! regrid.pl visits this program twice during the regridding process. During the first visit, the program does not use BCs data.
!   It just regrids hydrological variables and BCs-based land parameters in InRestart from InTile space to OutTile 
!   space (InRestart could be either a catchcn_internal_rst or a catch_internal_rst). If InRestart is a 
!   catchcn_internal_rst, carbon variables will be regridded using the same simple nearest neighbor algorithm (getids.H) that 
!   was employed for regridding all other variables. If InRestart is a catch_internal_rst, carbon variables will be 
!   filled with zeros.
 
!   During the second visit, the program uses the catchcn_internal_rst produced from the first visit as InRestart (herein 
!   referred to as InRestart2 which is in OutTile space already). The program reads BCs data from BCSDIR, carbon variables  
!   from an offline simulation on the SMAP_EASEv2_M09 grid which  has been initialized by another 3000-year offline simulation, and 
!   hydrological from  
!      InRestart2 in Version 1, 
!      the same offline simulation on the SMAP_EASEv2_M09 in Version 2.
!   Then, they will be regridded to OutTile space. The regridding carbon variables utilizes a more complicated algorithm which looks 
!   for a M09 grid cell in the  neighborhood with a similar vegetation type seperately for each fractional vegetation type within the 
!   catchment-tile. Note, the model can have upto 4 different types per catchment-tile: primary and secondary types 
!   and 2 split types for each primary and secondary type.  
   
! regrid.pl will then execute Scale_CatchCN.F90 which reads catchcn_internal_rst files created in the above 2 steps,
!   and scale soil moisture variables to be consistent with the new BCs-based land parameters to produce the final
!   catchcn_internal_rst file.

! Output file format: Output catchcn_internal_rst is always a nc4 file.

! Here are available options:  
!  (1) OPT1 (for above first step)
!      Input : (1) catchcn_internal_rst from an existing AGCM run (will always be nc4) 
!              (2) InTile and OutTile are DIFFERENT 
!              (3) NO land BCs
!      OutPut: Every variable (BCs-based land parameters, hydrological variables, and carbon parameters) will be regridded
!              from InTile to OutTile space using the simple nearest neighbor algorithm (getids.H)

!  (2) OPT2 (for above first step)
!      Input : (1) catch_internal_rst from an existing AGCM run (either nc4 or binary) 
!              (2) InTile and OutTile are DIFFERENT 
!              (3) NO land BCs 
!      OutPut: BCs-based land parameters, and hydrological variables will regridded from InTile to OutTile space 
!              using the simple nearest neighbor algorithm (getids.H). All carbon variables are filled with zeros.

!  (3) OPT3 (above second step) : 
!      Input : (1) catchcn_internal_rst (file format is always nc4)
!              (2) InTile and OutTile are the same user defined OutTile
!              (3) land BCs, 
!      Output: BCs-based land parameters will be replaced and carbon variables will be filled with regridded (from the
!              nearest offline cell with the same vegetation type) data to produce catchcn_internal_rst

! ----------------------------------------------------------------------------------------------------------------------------------------------

                                                     ! ====================== !
                                                     !          Process       !
                                                     ! ====================== !
                           
!                                                              HAVEDATA
!                                                                 |
!                          _______________________________________________________________________    
!                         |                                                                       |
!                                                                   
!                   NO (OPT1/OPT2)                                                           YES (OPT3)                     
!                   --------------                                                           ----------
!OutTile   :          /= InTile                                                               == InTile
!regridding:   ID  (InTile to OutTile using getids.H)                       ID (one-to-one i.e. 1:NTILES, no regridding)
!                         |                                                                       |
!                     clsmcn_file                                                                 |
!          _____________________________________                                                  |
!         |                                     |                                                 |
!        YES (OPT1)                             NO (OPT2)                                         |
!InRestart : catchcn_internal_rst          catch_internal_rst                            catchcn_internal_rst
!         |                                     |                                                 |
!         |                                  filetype                                             |
!         |                                     |                                                 |
!         |                      _________________________________                                |
!         |                     |                                 |                               |
!         V                     0                                /= 0                             V
!call : read_catchcn_nc4   read_catch_nc4                   read_catch_bin                   read_bcs_data                 
!                               |                                 |
!                               -----------------------------------                
!                                               |                                  
!                                               V                                  
!1) reads InRestart nVars records  (1) reads InCNRestart/regrids/writes (1:65)           (1) reads BCs
!2) regrids                          (takes hydrological initial conditions              (2) writes 1:37; 66:72 
!3) writes                               from offline SMAP M09)                          (3) reads InRestart2/writes 38, 39,40=38,41:65      
!4) close files                    (2) close files                                       (4) call regrid_carbon_vars (from offline SMAP M09)                                  
!                                                                                              (a) reads from InCNRestart                                 
!                                                                                              (b) regrids each veg type from the nearest InRestart cell  
!                                                                                              (c) writes (73-192,193-1080)                               
!                                                                                              (d) close files                                            
!                             
!
!
!                                     OUTPUT catchcn_internal_rst will always be nc4
! ----------------------------------------------------------------------------------------------------------------------------------------------


! The order of the INTERNAL STATE variables in GEOS_CatchCNGridComp
! -----------------------------------------------------------------
!   1: BF1      
!   2: BF2      
!   3: BF3      
!   4: VGWMAX   
!   5: CDCR1    
!   6: CDCR2    
!   7: PSIS     
!   8: BEE      
!   9: POROS    
!  10: WPWET    
!  11: COND     
!  12: GNU      
!  13: ARS1     
!  14: ARS2     
!  15: ARS3     
!  16: ARA1     
!  17: ARA2     
!  18: ARA3     
!  19: ARA4     
!  20: ARW1     
!  21: ARW2     
!  22: ARW3     
!  23: ARW4     
!  24: TSA1     
!  25: TSA2     
!  26: TSB1     
!  27: TSB2     
!  28: ATAU     
!  29: BTAU     
!  30-33: ITY * NUM_VEG
!  34-37: FVEG * NUM_VEG
!  38: ((TC (n,i),n=1,n_catd),i=1,4)      
!  39: ((QC (n,i),n=1,n_catd),i=1,4)  
!  40: ((TG (n,i),n=1,n_catd),i=1,4)       
!  41: CAPAC    
!  42: CATDEF   
!  43: RZEXC    
!  44: SRFEXC   
!  45: GHTCNT1  
!  46: GHTCNT2  
!  47: GHTCNT3  
!  48: GHTCNT4  
!  49: GHTCNT5  
!  50: GHTCNT6  
!  51: TSURF    
!  52: WESNN1   
!  53: WESNN2   
!  54: WESNN3   
!  55: HTSNNN1  
!  56: HTSNNN2  
!  57: HTSNNN3  
!  58: SNDZN1   
!  59: SNDZN2   
!  60: SNDZN3   
!  61: ((CH (n,i),n=1,n_catd),i=1,4)        
!  62: ((CM (n,i),n=1,n_catd),i=1,4)        
!  63: ((CQ (n,i),n=1,n_catd),i=1,4)        
!  64: ((FR (n,i),n=1,n_catd),i=1,4)        
!  65: ((WW (n,i),n=1,n_catd),i=1,4)        
!  66: cat_id   
!  67: ndep   
!  68: cli_t2m   
!  69: BGALBVR
!  70: BGALBVF
!  71: BGALBNR
!  72: BGALBNF
!  73-192: CNCOL (n,nz*VAR_COL)
!  193-1080: CNPFT (n,nz*nv*VAR_PFT)
! 1081-1083: TGWM  (n,nz)
! 1084: SFMCM
! 1085: BFLOWM
! 1086: TOTWATM
! 1087: TAIRM
! 1088: TPM
! 1089: CNSUM
! 1090: SNDZM
! 1091: ASNOWM
! 1092-1103: PSNSUNM (n,nz*nv)
! 1104-1115: PSNSHAM (n,nz*nv)

  use MAPL
  use ESMF
  use gFTL_StringVector
  use ieee_arithmetic, only: isnan => ieee_is_nan
  use mk_restarts_getidsMod, only: GetIDs, ReadTileFile_RealLatLon
  use clm_varpar_shared , only : nzone => NUM_ZON_CN, nveg => NUM_VEG_CN, &
                                 VAR_COL => VAR_COL_40, VAR_PFT => VAR_PFT_40, &
                                 npft => numpft_CN

  implicit none
  include 'mpif.h'
  INCLUDE 'netcdf.inc'

  ! initialize to non-MPI values

  integer  :: myid=0, numprocs=1, mpierr, mpistatus(MPI_STATUS_SIZE)  
  logical  :: root_proc=.true.

  real, parameter :: nan = O'17760000000'
  real, parameter :: fmin= 1.e-4 ! ignore vegetation fractions at or below this value
  integer, parameter :: OutUnit = 40, InUnit = 50

  ! ===============================================================================================
  ! Below hard-wired ldas restart file is from a global offline simulation on the SMAP M09 grid
  ! after 1000s of years of simulations

  integer, parameter :: ntiles_cn = 1684725
  character(len=300), parameter :: &
       InCNRestart = '/discover/nobackup/projects/gmao/ssd/land/l_data/LandRestarts_for_Regridding/CatchCN/M09/20151231/catchcn_internal_rst', &
       InCNTilFile = '/discover/nobackup/projects/gmao/bcs_shared/legacy_bcs/Icarus-NLv3/Icarus-NLv3_EASE/SMAP_EASEv2_M09/SMAP_EASEv2_M09_3856x1624.til'     

  character(len=256), parameter :: CatNames   (57) = &
       (/'BF1    ','BF2    ','BF3    ','VGWMAX ','CDCR1  ', &
         'CDCR2  ','PSIS   ','BEE    ','POROS  ','WPWET  ', &
         'COND   ','GNU    ','ARS1   ','ARS2   ','ARS3   ', &
         'ARA1   ','ARA2   ','ARA3   ','ARA4   ','ARW1   ', &
         'ARW2   ','ARW3   ','ARW4   ','TSA1   ','TSA2   ', & 
         'TSB1   ','TSB2   ','ATAU   ','BTAU   ','OLD_ITY', &
         'TC     ','QC     ','CAPAC  ','CATDEF ','RZEXC  ', &
         'SRFEXC ','GHTCNT1','GHTCNT2','GHTCNT3','GHTCNT4', &
         'GHTCNT5','GHTCNT6','TSURF  ','WESNN1 ','WESNN2 ', &
         'WESNN3 ','HTSNNN1','HTSNNN2','HTSNNN3','SNDZN1 ', &
         'SNDZN2 ','SNDZN3 ','CH     ','CM     ','CQ     ', &
         'FR     ','WW     '/)

  character(len=256), parameter :: CarbNames (68) =  &
       (/'BF1    ','BF2    ','BF3    ','VGWMAX ','CDCR1  ', &
         'CDCR2  ','PSIS   ','BEE    ','POROS  ','WPWET  ', &
         'COND   ','GNU    ','ARS1   ','ARS2   ','ARS3   ', &
         'ARA1   ','ARA2   ','ARA3   ','ARA4   ','ARW1   ', &
         'ARW2   ','ARW3   ','ARW4   ','TSA1   ','TSA2   ', & 
         'TSB1   ','TSB2   ','ATAU   ','BTAU   ','ITY    ', &
         'FVG    ','TC     ','QC     ','TG     ','CAPAC  ', &
         'CATDEF ','RZEXC  ','SRFEXC ','GHTCNT1','GHTCNT2', &
         'GHTCNT3','GHTCNT4','GHTCNT5','GHTCNT6','TSURF  ', &
         'WESNN1 ','WESNN2 ','WESNN3 ','HTSNNN1','HTSNNN2', &
         'HTSNNN3','SNDZN1 ','SNDZN2 ','SNDZN3 ','CH     ', &
         'CM     ','CQ     ','FR     ','WW     ','TILE_ID', &
         'NDEP   ','CLI_T2M','BGALBVR','BGALBVF','BGALBNR', &
         'BGALBNF','CNCOL  ','CNPFT  '   /)
 
  integer :: AGCM_YY, AGCM_MM, AGCM_DD, AGCM_HR

  character*256 :: DataDir="OutData/clsm/"
  character*256 :: Usage="mk_CatchCNRestarts OutTileFile InTileFile InRestart SURFLAY RestartTime"
  character*256 :: OutTileFile, InTileFile, InRestart, arg(6), OutFileName
  character*10  :: RestartTime

  logical :: clsmcn_file = .true., RegridSMAP = .false.
  logical :: havedata
  integer :: i, i1, iargc, n, k, ncatch,ntiles,ntiles_in, filetype, rc, nVars, req, infos, STATUS
  integer, pointer  :: Id(:), id_loc(:), tid_in(:)
  real,    pointer  :: loni(:),lono(:), lati(:), lato(:) , lonn(:), latt(:)
  real              :: SURFLAY
  type(Netcdf4_Fileformatter) :: InFmt,OutFmt
  type(FileMetadata) :: InCfg,OutCfg
  integer, allocatable :: low_ind(:), upp_ind(:), nt_local (:)
  character(256) :: Iam = "mk_CatchCNRestarts"
  
  call init_MPI()
  call MPI_Info_create(infos, STATUS)                                 ; VERIFY_(STATUS)
  call MPI_Info_set(infos, "romio_cb_read", "automatic", STATUS)      ; VERIFY_(STATUS)
  call MPI_Barrier(MPI_COMM_WORLD, STATUS)
  
  !-----------------------------------------------------
  ! Read command-line arguments, file names (inRestart,
  ! inTile, outTile), determine file format, and BCs 
  ! availability.  
  !-----------------------------------------------------

  call ESMF_Initialize(LogKindFlag=ESMF_LOGKIND_NONE) 
  
  I = iargc()
  
  if( I /=5 ) then
     print *, "Wrong Number of arguments: ", i
     print *, trim(Usage)
     stop
  end if
  
  do n=1,I
     call getarg(n,arg(n))
  enddo

  read(arg(1),'(a)') OutTileFile
  read(arg(2),'(a)')  InTileFile
  read(arg(3),'(a)')  InRestart
  read(arg(4),*)  SURFLAY
  read(arg(5),'(a)') RestartTime

  if (SURFLAY.ne.20 .and. SURFLAY.ne.50) then
     print *, "You must supply a valid SURFLAY value:"
     print *, "(Ganymed-3 and earlier) SURFLAY=20.0 for Old Soil Params"
     print *, "(Ganymed-4 and later  ) SURFLAY=50.0 for New Soil Params"
     call exit(2)
  end if

  ! Are BCs data available? 
  ! -----------------------
 
  inquire(file=trim(DataDir)//"CLM_veg_typs_fracs",exist=havedata)
  
  ! Reading restart time stamp and constructing daylength array
  ! -----------------------------------------------------------
  read (RestartTime (1: 4), '(i4)', IOSTAT = K) AGCM_YY ; VERIFY_(K)
  read (RestartTime (5: 6), '(i2)', IOSTAT = K) AGCM_MM ; VERIFY_(K)
  read (RestartTime (7: 8), '(i2)', IOSTAT = K) AGCM_DD ; VERIFY_(K)
  read (RestartTime (9:10), '(i2)', IOSTAT = K) AGCM_HR ; VERIFY_(K)

  MPI_PROC0 : if (root_proc) then
     
     ! Read Output/Input  .til files
     call ReadTileFile_RealLatLon(OutTileFile, ntiles, xlon=lono, xlat=lato)  
     call ReadTileFile_RealLatLon(InTileFile,ntiles_in,xlon=loni, xlat=lati)
     allocate(Id (ntiles))
     
     ! ------------------------------------------------
     ! create output catchcn_internal_rst in nc4 format
     ! ------------------------------------------------
     
     call InFmt%open('/discover/nobackup/projects/gmao/ssd/land/l_data/LandRestarts_for_Regridding/CatchCN/catchcn_internal_dummy',pFIO_READ, __RC__) 
     InCfg=InFmt%read( __RC__)
     call MAPL_IOCountNonDimVars(InCfg,nvars, __RC__)
     call MAPL_IOChangeRes(InCfg,OutCfg,(/'tile'/),(/ntiles/),__RC__)
     i = index(InRestart,'/',back=.true.)     
     OutFileName = "OutData/"//trim(InRestart(i+1:))
     call OutFmt%create(OutFileName, __RC__)
     call OutFmt%write(OutCfg, __RC__)
     i1= index(InRestart,'/',back=.true.)
     i = index(InRestart,'catchcn',back=.true.)
     
  endif MPI_PROC0
  
  call MPI_Barrier(MPI_COMM_WORLD, mpierr)
  call MPI_BCAST(NTILES   ,     1, MPI_INTEGER  ,  0,MPI_COMM_WORLD,mpierr) ; VERIFY_(mpierr)    
  call MPI_BCAST(NTILES_IN,     1, MPI_INTEGER  ,  0,MPI_COMM_WORLD,mpierr) ; VERIFY_(mpierr)
  
  HAVE_DATA :if(havedata) then

     ! OPT3
     ! ----
     ! Get number of catchments
     ! ------------------------
     
     open(unit=22, &
          file=trim(DataDir)//"catchment.def",status='old',form='formatted')
     
     read(22,*) ncatch
     
     close(22)
     
     if(ncatch /= ntiles) then
        print *, "Number of tiles in BCs data, ",Ncatch," does not match number in OutTile file ", NTILES
        print *, trim(OutTileFile)
        stop
     endif
     
     if(ntiles_in /= ntiles) then
        print *, "HAVEDATA : Number of tiles in InTileFile, ",NTILES_IN," does not match number in OutTileFile ", NTILES
        print *, trim ( InTileFile)
        print *, trim (OutTileFile)
        stop
     endif
     
     allocate (Id(ntiles))  
     
     do i = 1,ntiles
        id (i) = i  ! Just one-to-one mapping
     end do
     RegridSMAP = .true.
   
     !OPT3 (Reading/writing BCs/hydrological variables)
   
     if (root_proc) call read_bcs_data  (ntiles, SURFLAY, OutFmt, InRestart, __RC__)
     
  else

     ! What is the format of the InRestart file?
     ! -----------------------------------------

     call MAPL_NCIOGetFileType(InRestart, filetype, __RC__)
     
     if (filetype /= 0) then
        
        ! OPT2 (filetype =/ 0: a binary file must be a catch_internal_rst)
        ! ----
        clsmcn_file = .false. 

        open(unit=InUnit,FILE=InRestart,form='unformatted',  &
             status='old',convert='little_endian')        
        
     else

        ! filetype = 0 : nc4, could be catch_internal_rst or catchcn_internal_rst
        ! check nVars: if nVars > 57 OPT1 (catchcn_internal_rst) ; else OPT2 (catch_internal_rst)
        ! ---------------------------------------------------------------------------------------

        call InFmt%open(InRestart,pFIO_READ, __RC__)
        InCfg = InFmt%read(__RC__)
        call InFmt%close()

        call MAPL_IOCountNonDimVars(InCfg,nvars)
       
        if(nVars == 57) clsmcn_file = .false. 

     endif

     CATCHCN: if (clsmcn_file) then

        ! OPT1
        ! ----

        ! ----------------------------------------------------
        ! INPUT/OUTPUT Mapping since InTileFile =/ OutTileFile
        ! ----------------------------------------------------
        
        if(myid > 0)    allocate (loni   (1:ntiles_in))
        if(myid > 0)    allocate (lati   (1:ntiles_in))
        
        allocate (tid_in (1:ntiles_in))
        do n = 1, NTILES_IN
           tid_in (n) = n
        end do
        
        call MPI_BCAST(loni,ntiles_in,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_BCAST(lati,ntiles_in,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        
        ! Now mapping (Id)
        ! ----------------
        
        allocate (Id(ntiles))  ! Id contains corresponding InTileID after mapping InTiles on to OutTile
        ! call GetIds(loni,lati,lono,lato,zoom,Id)
        allocate(low_ind (   numprocs))
        allocate(upp_ind (   numprocs))
        allocate(nt_local(   numprocs))
        low_ind (:)    = 1
        upp_ind (:)    = NTILES       
        nt_local(:)    = NTILES 
        
        ! Domain decomposition
        ! --------------------
        
        if (numprocs > 1) then      
           do i = 1, numprocs - 1
              upp_ind(i)   = low_ind(i) + (ntiles/numprocs) - 1 
              low_ind(i+1) = upp_ind(i) + 1
              nt_local(i)  = upp_ind(i) - low_ind(i) + 1
           end do
           nt_local(numprocs) = upp_ind(numprocs) - low_ind(numprocs) + 1
        endif
        
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
        
        call GetIds(loni,lati,lonn,latt,id_loc, tid_in)
        call MPI_Barrier(MPI_COMM_WORLD, mpierr)
        
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
        
        if(root_proc) deallocate (lono, lato,lonn,latt, tid_in)
     
        deallocate (loni,lati)


       if (root_proc)  call read_catchcn_nc4 (NTILES_IN, NTILES, OutFmt, ID, InRestart, __RC__)

     else

        call  regrid_hyd_vars (NTILES, OutFmt) 

        ! OPT2
        ! ----
        !        NC4ORBIN: if(filetype ==0) then
        !
        !           call read_catch_nc4 (NTILES_IN, NTILES, OutFmt, ID, InRestart)
        !
        !        else
        !
        !           call read_catch_bin (NTILES_IN, NTILES, OutFmt, ID)
        !
        !        endif NC4ORBIN

     endif CATCHCN
 
  endif HAVE_DATA

  if (root_proc) then
          
     ! -----------------
     ! BEGIN THE PROCESS
     ! -----------------
     
     print *, "                                         "
     print *, "**********************************************************************"
     print *, "                                         "
     print *, "mk_CatchCNRestarts Configuration"
     print *, "--------------------------------"
     print *, "                                         "  
     print '(A22, i4.4,i2.2,i2.2,i2.2)', " Restart Time         :",AGCM_YY,AGCM_MM,AGCM_DD,AGCM_HR 
     print *, 'SURFLAY              : ',SURFLAY
     print *, 'Have BCs data        : ',havedata
     print *, "# of tiles in InTile : ",ntiles_in
     print *, "# of tiles in OutTile: ",ntiles
     
     if(clsmcn_file) then 
        print *,"InRestart is from    : Catchment-carbon AGCM simulation"
     else
        InRestart = trim(InCNRestart)
        print *,"InRestart is from    : offline SMAP_EASEv2_M09"
     endif
     
     print *, "InRestart filename   : ",trim(InRestart)
     print *, "OutRestart filename  : ",trim(OutFileName)
     print *, "OutRestart file fmt  : nc4"
     print *, "                                         "
     print *, "**********************************************************************"
     print *, "                                         "
     
  endif
  
  call MPI_BCAST(OutFileName ,   256, MPI_CHARACTER,  0,MPI_COMM_WORLD,mpierr)
  call MPI_Barrier(MPI_COMM_WORLD, mpierr)
  
  if (RegridSMAP) then
     ntiles_in = ntiles_cn
     !OPT3 (carbon variables from offline SMAP M09) 
     call regrid_carbon_vars (NTILES,AGCM_YY,AGCM_MM,AGCM_DD,AGCM_HR, OutFileName, OutTileFile) 
     !   call regrid_carbon_vars_omp (NTILES,AGCM_YY,AGCM_MM,AGCM_DD,AGCM_HR, OutFileName, OutTileFile) 
     
  endif
call MPI_BARRIER( MPI_COMM_WORLD, mpierr)
call ESMF_Finalize(endflag=ESMF_END_KEEPMPI)
call MPI_FINALIZE(mpierr)
  
contains

  ! *****************************************************************************
  
  SUBROUTINE read_bcs_data (ntiles, SURFLAY, OutFmt, InRestart, rc)

    ! This subroutine :
    !  1) reads BCs from BCSDIR and hydrological varables from InRestart.
    !     InRestart is a catchcn_internal_rst nc4 file.
    !
    !  2) writes out BCs and hydrological variables in catchcn_internal_rst (1:72). 
    !     output catchcn_internal_rst is nc4.

    implicit none
    real, intent (in)                         :: SURFLAY
    integer, intent (in)                      :: ntiles
    character (*), intent (in)                :: InRestart
    type(Netcdf4_Fileformatter), intent (inout)           :: OutFmt
    integer, optional, intent(out) :: rc

    real, allocatable :: CLMC_pf1(:), CLMC_pf2(:), CLMC_sf1(:), CLMC_sf2(:)
    real, allocatable :: CLMC_pt1(:), CLMC_pt2(:), CLMC_st1(:), CLMC_st2(:)    
    real, allocatable :: BF1(:),   BF2(:),   BF3(:),  VGWMAX(:)
    real, allocatable :: CDCR1(:), CDCR2(:), PSIS(:), BEE(:) 
    real, allocatable :: POROS(:), WPWET(:), COND(:), GNU(:)
    real, allocatable :: ARS1(:),  ARS2(:),  ARS3(:)
    real, allocatable :: ARA1(:),  ARA2(:),  ARA3(:), ARA4(:)
    real, allocatable :: ARW1(:),  ARW2(:),  ARW3(:), ARW4(:)
    real, allocatable :: TSA1(:),  TSA2(:),  TSB1(:), TSB2(:)
    real, allocatable :: ATAU2(:), BTAU2(:), DP2BR(:), rity(:), CanopH(:)
    real, allocatable :: NDEP(:), BVISDR(:), BVISDF(:), BNIRDR(:), BNIRDF(:) 
    real, allocatable :: T2(:), var1(:)
    integer, allocatable :: ity(:)
    character*256                            :: vname
    character*256 :: DataDir="OutData/clsm/"
    integer       :: idum, i,j,n, ib, nv
    real          :: rdum, zdep1, zdep2, zdep3, zmet, term1, term2, bare,fvg(4)
    logical       :: file_exists
    type(Netcdf4_Fileformatter)                          :: InFmt,CatchCNFmt, CatchFmt
    integer :: status
  
    allocate (   BF1(ntiles),    BF2 (ntiles),     BF3(ntiles)  )
    allocate (VGWMAX(ntiles),   CDCR1(ntiles),   CDCR2(ntiles)  ) 
    allocate (  PSIS(ntiles),     BEE(ntiles),   POROS(ntiles)  ) 
    allocate ( WPWET(ntiles),    COND(ntiles),     GNU(ntiles)  )
    allocate (  ARS1(ntiles),    ARS2(ntiles),    ARS3(ntiles)  )
    allocate (  ARA1(ntiles),    ARA2(ntiles),    ARA3(ntiles)  )
    allocate (  ARA4(ntiles),    ARW1(ntiles),    ARW2(ntiles)  )
    allocate (  ARW3(ntiles),    ARW4(ntiles),    TSA1(ntiles)  )
    allocate (  TSA2(ntiles),    TSB1(ntiles),    TSB2(ntiles)  )
    allocate ( ATAU2(ntiles),   BTAU2(ntiles),   DP2BR(ntiles)  )
    allocate (BVISDR(ntiles),  BVISDF(ntiles),  BNIRDR(ntiles)  )
    allocate (BNIRDF(ntiles),      T2(ntiles),    NDEP(ntiles)  )    
    allocate (   ity(ntiles),      rity(ntiles),    CanopH(ntiles))
    allocate (CLMC_pf1(ntiles), CLMC_pf2(ntiles), CLMC_sf1(ntiles))
    allocate (CLMC_sf2(ntiles), CLMC_pt1(ntiles), CLMC_pt2(ntiles))
    allocate (CLMC_st1(ntiles), CLMC_st2(ntiles))

    inquire(file = trim(DataDir)//'/catchcn_params.nc4', exist=file_exists)

    if(file_exists) then

       print *,'FILE FORMAT FOR LAND BCS IS NC4'
       call CatchFmt%open(trim(DataDir)//'/catch_params.nc4',pFIO_READ, __RC__)
       call CatchCNFmt%open(trim(DataDir)//'/catchcn_params.nc4',pFIO_READ, __RC__)
       call MAPL_VarRead ( CatchFmt ,'OLD_ITY', rity, __RC__)
       call MAPL_VarRead ( CatchFmt ,'ARA1', ARA1, __RC__)
       call MAPL_VarRead ( CatchFmt ,'ARA2', ARA2, __RC__)
       call MAPL_VarRead ( CatchFmt ,'ARA3', ARA3, __RC__)
       call MAPL_VarRead ( CatchFmt ,'ARA4', ARA4, __RC__)
       call MAPL_VarRead ( CatchFmt ,'ARS1', ARS1, __RC__)
       call MAPL_VarRead ( CatchFmt ,'ARS2', ARS2, __RC__)
       call MAPL_VarRead ( CatchFmt ,'ARS3', ARS3, __RC__)
       call MAPL_VarRead ( CatchFmt ,'ARW1', ARW1, __RC__)
       call MAPL_VarRead ( CatchFmt ,'ARW2', ARW2, __RC__)
       call MAPL_VarRead ( CatchFmt ,'ARW3', ARW3, __RC__)
       call MAPL_VarRead ( CatchFmt ,'ARW4', ARW4, __RC__)

       if( SURFLAY.eq.20.0 ) then
          call MAPL_VarRead ( CatchFmt ,'ATAU2', ATAU2, __RC__)
          call MAPL_VarRead ( CatchFmt ,'BTAU2', BTAU2, __RC__)
       endif

       if( SURFLAY.eq.50.0 ) then
          call MAPL_VarRead ( CatchFmt ,'ATAU5', ATAU2, __RC__)
          call MAPL_VarRead ( CatchFmt ,'BTAU5', BTAU2, __RC__)
       endif

       call MAPL_VarRead ( CatchFmt ,'PSIS', PSIS, __RC__)
       call MAPL_VarRead ( CatchFmt ,'BEE', BEE, __RC__)
       call MAPL_VarRead ( CatchFmt ,'BF1', BF1, __RC__)
       call MAPL_VarRead ( CatchFmt ,'BF2', BF2, __RC__)
       call MAPL_VarRead ( CatchFmt ,'BF3', BF3, __RC__)
       call MAPL_VarRead ( CatchFmt ,'TSA1', TSA1, __RC__)
       call MAPL_VarRead ( CatchFmt ,'TSA2', TSA2, __RC__)
       call MAPL_VarRead ( CatchFmt ,'TSB1', TSB1, __RC__)
       call MAPL_VarRead ( CatchFmt ,'TSB2', TSB2, __RC__)
       call MAPL_VarRead ( CatchFmt ,'COND', COND, __RC__)
       call MAPL_VarRead ( CatchFmt ,'GNU', GNU, __RC__)
       call MAPL_VarRead ( CatchFmt ,'WPWET', WPWET, __RC__)
       call MAPL_VarRead ( CatchFmt ,'DP2BR', DP2BR, __RC__)
       call MAPL_VarRead ( CatchFmt ,'POROS', POROS, __RC__)
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
       call CatchFmt%close()
       call CatchCNFmt%close()
      
    else

       open(unit=22, &
            file=trim(DataDir)//"mosaic_veg_typs_fracs",status='old',form='formatted')
       
       do N=1,ntiles
          read(22,*) I, j, ITY(N),idum, rdum, rdum, CanopH(N)
       enddo
       
       rity(:) = float(ity)
       
       close(22)
       
       open(unit=22, file=trim(DataDir)//'bf.dat'               ,form='formatted')
       open(unit=23, file=trim(DataDir)//'soil_param.dat'       ,form='formatted')
       open(unit=24, file=trim(DataDir)//'ar.new'               ,form='formatted')
       open(unit=25, file=trim(DataDir)//'ts.dat'               ,form='formatted')
       open(unit=26, file=trim(DataDir)//'tau_param.dat'        ,form='formatted')
       open(unit=27, file=trim(DataDir)//'CLM_veg_typs_fracs'   ,form='formatted')
       open(unit=28, file=trim(DataDir)//'CLM_NDep_SoilAlb_T2m' ,form='formatted')
       
       do n=1,ntiles
          read (22, *) i,j, GNU(n), BF1(n), BF2(n), BF3(n)
          
          read (23, *) i,j, idum, idum, BEE(n), PSIS(n),&
               POROS(n), COND(n), WPWET(n), DP2BR(n)
          
          read (24, *) i,j, rdum, ARS1(n), ARS2(n), ARS3(n),          &
               ARA1(n), ARA2(n), ARA3(n), ARA4(n), &
               ARW1(n), ARW2(n), ARW3(n), ARW4(n)
          
          read (25, *) i,j, rdum, TSA1(n), TSA2(n), TSB1(n), TSB2(n)
          
          if( SURFLAY.eq.20.0 ) read (26, *) i,j, ATAU2(n), BTAU2(n), rdum, rdum   ! for old soil params
          if( SURFLAY.eq.50.0 ) read (26, *) i,j, rdum , rdum, ATAU2(n), BTAU2(n)  ! for new soil params
          
          read (27, *) i,j, CLMC_pt1(n), CLMC_pt2(n), CLMC_st1(n), CLMC_st2(n), &
               CLMC_pf1(n), CLMC_pf2(n), CLMC_sf1(n), CLMC_sf2(n)
          
          read (28, *) NDEP(n), BVISDR(n), BVISDF(n), BNIRDR(n), BNIRDF(n), T2(n) ! MERRA-2 Annual Mean Temp is default.
          
       end do

       CLOSE (22, STATUS = 'KEEP')
       CLOSE (23, STATUS = 'KEEP')
       CLOSE (24, STATUS = 'KEEP')
       CLOSE (25, STATUS = 'KEEP')
       CLOSE (26, STATUS = 'KEEP')
       CLOSE (27, STATUS = 'KEEP')
       CLOSE (28, STATUS = 'KEEP')
       
    endif

   do n=1,ntiles

      BVISDR(n) = amax1(1.e-6, BVISDR(n))
      BVISDF(n) = amax1(1.e-6, BVISDF(n))
      BNIRDR(n) = amax1(1.e-6, BNIRDR(n))
      BNIRDF(n) = amax1(1.e-6, BNIRDF(n))
      
      zdep2=1000.
      zdep3=amax1(1000.,DP2BR(n))
      
        if (zdep2 .gt.0.75*zdep3) then
           zdep2  =  0.75*zdep3              
        end if
        
        zdep1=20.
        zmet=zdep3/1000.
        
        term1=-1.+((PSIS(n)-zmet)/PSIS(n))**((BEE(n)-1.)/BEE(n))
        term2=PSIS(n)*BEE(n)/(BEE(n)-1)
        
        VGWMAX(n) = POROS(n)*zdep2   
        CDCR1(n)  = 1000.*POROS(n)*(zmet-(-term2*term1))   
        CDCR2(n)  = (1.-WPWET(n))*POROS(n)*zdep3

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
     end do
     


     ! Now writing BCs (from BCSDIR) and regridded hydrological variables 1-72
     ! -----------------------------------------------------------------------

     call InFmt%open(InRestart,pFIO_READ, __RC__)
     
     call MAPL_VarWrite(OutFmt,trim(CarbNames(1)),BF1)            !   1 
     call MAPL_VarWrite(OutFmt,trim(CarbNames(2)),BF2)            !   2
     call MAPL_VarWrite(OutFmt,trim(CarbNames(3)),BF3)            !   3
     call MAPL_VarWrite(OutFmt,trim(CarbNames(4)),VGWMAX)         !   4
     call MAPL_VarWrite(OutFmt,trim(CarbNames(5)),CDCR1)          !   5
     call MAPL_VarWrite(OutFmt,trim(CarbNames(6)),CDCR2)          !   6
     call MAPL_VarWrite(OutFmt,trim(CarbNames(7)),PSIS)           !   7
     call MAPL_VarWrite(OutFmt,trim(CarbNames(8)),BEE)            !   8
     call MAPL_VarWrite(OutFmt,trim(CarbNames(9)),POROS)          !   9
     call MAPL_VarWrite(OutFmt,trim(CarbNames(10)),WPWET)         !  10
     call MAPL_VarWrite(OutFmt,trim(CarbNames(11)),COND)          !  11
     call MAPL_VarWrite(OutFmt,trim(CarbNames(12)),GNU)           !  12
     call MAPL_VarWrite(OutFmt,trim(CarbNames(13)),ARS1)          !  13
     call MAPL_VarWrite(OutFmt,trim(CarbNames(14)),ARS2)          !  14
     call MAPL_VarWrite(OutFmt,trim(CarbNames(15)),ARS3)          !  15
     call MAPL_VarWrite(OutFmt,trim(CarbNames(16)),ARA1)          !  16
     call MAPL_VarWrite(OutFmt,trim(CarbNames(17)),ARA2)          !  17
     call MAPL_VarWrite(OutFmt,trim(CarbNames(18)),ARA3)          !  18
     call MAPL_VarWrite(OutFmt,trim(CarbNames(19)),ARA4)          !  19
     call MAPL_VarWrite(OutFmt,trim(CarbNames(20)),ARW1)          !  20
     call MAPL_VarWrite(OutFmt,trim(CarbNames(21)),ARW2)          !  21
     call MAPL_VarWrite(OutFmt,trim(CarbNames(22)),ARW3)          !  22
     call MAPL_VarWrite(OutFmt,trim(CarbNames(23)),ARW4)          !  23
     call MAPL_VarWrite(OutFmt,trim(CarbNames(24)),TSA1)          !  24
     call MAPL_VarWrite(OutFmt,trim(CarbNames(25)),TSA2)          !  25
     call MAPL_VarWrite(OutFmt,trim(CarbNames(26)),TSB1)          !  26
     call MAPL_VarWrite(OutFmt,trim(CarbNames(27)),TSB2)          !  27
     call MAPL_VarWrite(OutFmt,trim(CarbNames(28)),ATAU2)         !  28
     call MAPL_VarWrite(OutFmt,trim(CarbNames(29)),BTAU2)         !  29
     call MAPL_VarWrite(OutFmt,'ITY',CLMC_pt1,offset1=1)     !  30
     call MAPL_VarWrite(OutFmt,'ITY',CLMC_pt2,offset1=2)     !  31
     call MAPL_VarWrite(OutFmt,'ITY',CLMC_st1,offset1=3)     !  32
     call MAPL_VarWrite(OutFmt,'ITY',CLMC_st2,offset1=4)     !  33
     call MAPL_VarWrite(OutFmt,'FVG',CLMC_pf1,offset1=1)     !  34
     call MAPL_VarWrite(OutFmt,'FVG',CLMC_pf2,offset1=2)     !  35
     call MAPL_VarWrite(OutFmt,'FVG',CLMC_sf1,offset1=3)     !  36
     call MAPL_VarWrite(OutFmt,'FVG',CLMC_sf2,offset1=4)     !  37
   
     allocate(var1(ntiles))                            

     ! TC QC TG 
     
     do n = 38,40 
        if(n == 38) vname = 'TC'
        if(n == 39) vname = 'QC'
        if(n == 40) vname = 'TG'
        do j = 1,4
           call MAPL_VarRead ( InFmt,vname,var1 ,offset1=j, __RC__)
           call MAPL_VarWrite(OutFmt,vname,var1 ,offset1=j)  ! 38-40
        end do
     end do
     
     ! CAPAC CATDEF RZEXC SRFEXC ... SNDZN3
     
     do n=41,60
        call MAPL_VarRead ( InFmt,trim(CarbNames(n-6)),var1, __RC__)
        call MAPL_VarWrite(OutFmt,trim(CarbNames(n-6)),var1)         !  41-60
     enddo
     
     ! CH CM CQ FR WW
     var1 = 0.

     do n=61,65
        if((n >= 61).AND.(n <= 63)) var1 = 1.e-3
        if(n == 64) var1 = 0.25
        if(n == 65) var1 = 0.1
        do j = 1,4
           
           call MAPL_VarRead ( InFmt,trim(CarbNames(n-6)),var1 ,offset1=j, __RC__)
           call MAPL_VarWrite(OutFmt,trim(CarbNames(n-6)),var1 ,offset1=j)  ! 61-65
        end do
     end do
     
     do i=1,ntiles
        var1(i) = real(i)
     end do
     
     call MAPL_VarWrite(OutFmt,'TILE_ID',var1  )        !  66 : cat_id   
     call MAPL_VarWrite(OutFmt,'NDEP'   ,NDEP  )        !  67 : ndep       
     call MAPL_VarWrite(OutFmt,'CLI_T2M',T2    )        !  68 : cli_t2m    
     call MAPL_VarWrite(OutFmt,'BGALBVR',BVISDR)        !  69 : BGALBVR    
     call MAPL_VarWrite(OutFmt,'BGALBVF',BVISDF)        !  70 : BGALBVF    
     call MAPL_VarWrite(OutFmt,'BGALBNR',BNIRDR)        !  71 : BGALBNR    
     call MAPL_VarWrite(OutFmt,'BGALBNF',BNIRDF)        !  72 : BGALBNF    
     
     deallocate (var1)
     call InFmt%close()
     call OutFmt%close()

! Vegdyn Boundary Condition
! -------------------------
!
!     open(20,file=trim("OutData/vegdyn_internal_rst"), &
!          status="unknown", &
!          form="unformatted",convert="little_endian")
!     write(20) rity  
!     write(20) CanopH
!     close(20)
!     print *, "Wrote vegdyn_internal_restart"

     deallocate (   BF1,     BF2,     BF3  )
     deallocate (VGWMAX,   CDCR1,   CDCR2  ) 
     deallocate (  PSIS,     BEE,   POROS  ) 
     deallocate ( WPWET,    COND,     GNU  )
     deallocate (  ARS1,    ARS2,    ARS3  )
     deallocate (  ARA1,    ARA2,    ARA3  )
     deallocate (  ARA4,    ARW1,    ARW2  )
     deallocate (  ARW3,    ARW4,    TSA1  )
     deallocate (  TSA2,    TSB1,    TSB2  )
     deallocate ( ATAU2,   BTAU2,   DP2BR  )
     deallocate (BVISDR,  BVISDF,  BNIRDR  )
     deallocate (BNIRDF,      T2,    NDEP  )    
     deallocate (   ity,    rity,    CanopH)
     deallocate (CLMC_pf1, CLMC_pf2, CLMC_sf1)
     deallocate (CLMC_sf2, CLMC_pt1, CLMC_pt2)
     deallocate (CLMC_st1,CLMC_st2)
     if (present(rc)) rc = 0
     !_RETURN(_SUCCESS)
  END SUBROUTINE read_bcs_data
  
   ! *****************************************************************************
  
  SUBROUTINE read_catchcn_nc4 (NTILES_IN, NTILES, OutFmt, IDX, InRestart, rc)

    implicit none

    ! Reads catchcn_internal_rst nc4 file, regrids every single variable and writes 
    ! out catchcn_internal_rst in nc4 format.
    ! This subroutine is called when BCs data are not available. 

    integer, intent (in)                     :: NTILES_IN, NTILES
    character(*), intent (in)                :: InRestart
    type(Netcdf4_Fileformatter), intent (inout) :: OutFmt
    integer, dimension (NTILES), intent (in) :: IDX    
    integer, optional, intent(out) :: rc
    type(Netcdf4_Fileformatter)              :: InFmt
    type(FileMetadata)                      :: InCfg
    integer                                  :: n,i,j, ndims, nVars,dim1,dim2
    character(len=:), pointer                :: vname
    real, allocatable                        :: var1 (:), var2 (:)
    integer, allocatable                     :: TILE_ID (:)
    type(StringVariableMap), pointer :: variables
    type(Variable), pointer :: var
    type(StringVariableMapIterator) :: var_iter
    type(StringVector), pointer :: var_dimensions
    character(len=:), pointer :: dname
    integer :: status

    call InFmt%open(InRestart,pFIO_READ, __RC__)
    InCfg = InFmt%read(__RC__)
    
    allocate (var1 (1:NTILES_IN))
    allocate (var2 (1:NTILES_IN))
    allocate (TILE_ID (1:NTILES_IN)) 

    call MAPL_VarRead ( InFmt,'TILE_ID',var1, __RC__)
    do n = 1, NTILES_IN      
       tile_id (NINT (var1(n))) = n
    end do

    variables => InCfg%get_variables()
    var_iter = variables%begin()
    do while (var_iter /= variables%end())

       vname => var_iter%key()     
       var => var_iter%value()
       var_dimensions => var%get_dimensions()
      
       ndims = var_dimensions%size()

       if (ndims == 1) then
          call MAPL_VarRead ( InFmt,vname,var1, __RC__)
          var2 = var1 (tile_id)
          call MAPL_VarWrite(OutFmt,vname,var2(idx))
          
       else if (ndims == 2) then
         
          dname => var%get_ith_dimension(2)
          dim1=InCfg%get_dimension(dname)
          do j=1,dim1
             call MAPL_VarRead ( InFmt,vname,var1     ,offset1=j, __RC__)
             var2 = var1 (tile_id)
             call MAPL_VarWrite(OutFmt,vname,var2(idx),offset1=j)
          enddo
          
       else if (ndims == 3) then
          
          dname => var%get_ith_dimension(2)
          dim1=InCfg%get_dimension(dname)
          dname => var%get_ith_dimension(3)
          dim2=InCfg%get_dimension(dname)
          do i=1,dim2
             do j=1,dim1
                call MAPL_VarRead ( InFmt,vname,var1     ,offset1=j,offset2=i, __RC__)
                var2 = var1 (tile_id)
                call MAPL_VarWrite(OutFmt,vname,var2(idx),offset1=j,offset2=i)
             enddo
          enddo
          
       end if

       call var_iter%next()
    enddo
 
    deallocate (var1, var2, tile_id)
    call InFmt%close()
    call OutFmt%close()
    if (present(rc)) rc = 0
    !_RETURN(_SUCCESS)
  END SUBROUTINE read_catchcn_nc4

  ! *****************************************************************************
  
  SUBROUTINE regrid_carbon_vars (                               &
       NTILES,AGCM_YY,AGCM_MM,AGCM_DD,AGCM_HR, OutFileName, OutTileFile) 

    implicit none
    character (*), intent (in)           :: OutTileFile, OutFileName
    integer, intent (in)                 :: NTILES,AGCM_YY,AGCM_MM,AGCM_DD,AGCM_HR 
    real, allocatable, dimension (:)     :: CLMC_pf1, CLMC_pf2, CLMC_sf1, CLMC_sf2, &
         CLMC_pt1, CLMC_pt2,CLMC_st1,CLMC_st2

    ! ===============================================================================================

    integer :: iclass(npft) = (/1,1,2,3,3,4,5,5,6,7,8,9,10,11,12,11,12,11,12/)
    integer, allocatable, dimension(:,:) :: Id_glb, Id_loc
    integer, allocatable, dimension(:)   :: tid_offl, id_vec
    real,    allocatable, dimension(:,:) :: fveg_offl,  ityp_offl
    real    :: fveg_new, sub_dist
    integer :: n,i,j, k, nv, nx, nz, iv, offl_cell, ityp_new, STATUS,NCFID, req
    integer :: outid, local_id
    integer, allocatable, dimension (:) :: sub_tid, sub_ityp1, sub_ityp2,icl_ityp1
    real   , pointer, dimension (:) :: sub_lon, sub_lat, rev_dist, sub_fevg1, sub_fevg2,&
         lonc, latc, LATT, LONN, DAYX, long, latg, var_dum, TILE_ID, var_dum2
    real, allocatable :: var_off_col (:,:,:), var_off_pft (:,:,:,:) 
    real, allocatable :: var_col_out (:,:,:), var_pft_out (:,:,:,:) 
    integer, allocatable :: low_ind(:), upp_ind(:), nt_local (:)
    integer :: AGCM_YYY, AGCM_MMM, AGCM_DDD, AGCM_HRR, AGCM_MI, AGCM_S, dofyr
    type(MAPL_SunOrbit)         :: ORBIT
    type(ESMF_Time)             :: CURRENT_TIME
    type(ESMF_TimeInterval)     :: timeStep
    type(ESMF_Clock)            :: CLOCK
    type(ESMF_Config)           :: CF


    allocate (tid_offl  (ntiles_cn))
    allocate (ityp_offl (ntiles_cn,nveg))
    allocate (fveg_offl (ntiles_cn,nveg))

    allocate(low_ind (   numprocs))
    allocate(upp_ind (   numprocs))
    allocate(nt_local(   numprocs))

    low_ind (:)    = 1
    upp_ind (:)    = NTILES       
    nt_local(:)    = NTILES 

    ! Domain decomposition
    ! --------------------

    if (numprocs > 1) then      
       do i = 1, numprocs - 1
          upp_ind(i)   = low_ind(i) + (ntiles/numprocs) - 1 
          low_ind(i+1) = upp_ind(i) + 1
          nt_local(i)  = upp_ind(i) - low_ind(i) + 1
       end do
       nt_local(numprocs) = upp_ind(numprocs) - low_ind(numprocs) + 1
    endif

    allocate (id_loc  (nt_local (myid + 1),4))
    allocate (lonn    (nt_local (myid + 1)))
    allocate (latt    (nt_local (myid + 1)))
    allocate (CLMC_pf1(nt_local (myid + 1)))
    allocate (CLMC_pf2(nt_local (myid + 1)))
    allocate (CLMC_sf1(nt_local (myid + 1)))
    allocate (CLMC_sf2(nt_local (myid + 1)))
    allocate (CLMC_pt1(nt_local (myid + 1)))
    allocate (CLMC_pt2(nt_local (myid + 1)))
    allocate (CLMC_st1(nt_local (myid + 1)))
    allocate (CLMC_st2(nt_local (myid + 1)))
    allocate (lonc   (1:ntiles_cn))
    allocate (latc   (1:ntiles_cn))

    if (root_proc) then
       
       ! --------------------------------------------
       ! Read exact lonn, latt from output .til file 
       ! --------------------------------------------

       allocate (long   (ntiles))
       allocate (latg   (ntiles))
       allocate (DAYX   (NTILES))

       call ReadTileFile_RealLatLon (OutTileFile, i, xlon=long, xlat=latg)

       !-----------------------
       ! COMPUTE DAYX
       !-----------------------

        AGCM_YYY = AGCM_YY
        AGCM_MMM = AGCM_MM
        AGCM_DDD = AGCM_DD
        AGCM_HRR = AGCM_HR
        AGCM_MI = 0
        AGCM_S = 0


        call ESMF_CalendarSetDefault ( ESMF_CALKIND_GREGORIAN, rc=status )

       ! get current date & time
       ! -----------------------
        call ESMF_TimeSet  ( CURRENT_TIME, YY = AGCM_YYY,       &
                                           MM = AGCM_MMM,       &
                                           DD = AGCM_DDD,       &
                                           H  = AGCM_HRR,       &
                                           M  = AGCM_MI,       &
                                           S  = AGCM_S ,       &
                                           rc=status )
         VERIFY_(STATUS)

         call ESMF_TimeIntervalSet(TimeStep,  S=450, RC=status)
         clock = ESMF_ClockCreate(TimeStep, startTime = CURRENT_TIME, RC=status)
         VERIFY_(STATUS)
         call ESMF_ClockSet ( clock, CurrTime=CURRENT_TIME, rc=status )

         CF = ESMF_ConfigCreate(RC=STATUS)
         VERIFY_(status)

         ORBIT = MAPL_SunOrbitCreateFromConfig(CF, CLOCK, .false., RC=status)
         VERIFY_(status)
   
     ! compute current daylight duration
        !----------------------------------
        call MAPL_SunGetDaylightDuration(ORBIT,latg,dayx,currTime=CURRENT_TIME,RC=STATUS)
        VERIFY_(STATUS)

       ! ---------------------------------------------
       ! Read exact lonc, latc from offline .til File 
       ! ---------------------------------------------

       call ReadTileFile_RealLatLon(InCNTilFile,i,xlon=lonc,xlat=latc)

    endif

!    call MPI_SCATTERV (                    &
!         long,nt_local,low_ind-1,MPI_real, &
!         lonn,size(lonn),MPI_real  , &
!         0,MPI_COMM_WORLD, mpierr )
!
!    call MPI_SCATTERV (                    &
!         latg,nt_local,low_ind-1,MPI_real, &
!         latt,nt_local(myid+1),MPI_real  , &
!         0,MPI_COMM_WORLD, mpierr )

    do i = 1, numprocs
       if((I == 1).and.(myid == 0)) then
          lonn(:) = long(low_ind(i) : upp_ind(i))
          latt(:) = latg(low_ind(i) : upp_ind(i))
       else if (I > 1) then
          if(I-1 == myid) then
             ! receiving from root
             call MPI_RECV(lonn,nt_local(i) , MPI_REAL, 0,995,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
             call MPI_RECV(latt,nt_local(i) , MPI_REAL, 0,994,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
          else if (myid == 0) then
             ! root sends
             call MPI_ISend(long(low_ind(i) : upp_ind(i)),nt_local(i),MPI_REAL,i-1,995,MPI_COMM_WORLD,req,mpierr)
             call MPI_WAIT (req,MPI_STATUS_IGNORE,mpierr)    
             call MPI_ISend(latg(low_ind(i) : upp_ind(i)),nt_local(i),MPI_REAL,i-1,994,MPI_COMM_WORLD,req,mpierr)
             call MPI_WAIT (req,MPI_STATUS_IGNORE,mpierr) 
          endif
       endif
    end do
    
    if(root_proc) deallocate (long, latg)
 
    call MPI_BCAST(lonc,ntiles_cn,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_BCAST(latc,ntiles_cn,MPI_REAL,0,MPI_COMM_WORLD,mpierr)

    ! Open GKW/Fzeng SMAP M09 catchcn_internal_rst and output catchcn_internal_rst
    ! ----------------------------------------------------------------------------
    !    call MPI_Info_create(info, STATUS)
    !    call MPI_Info_set(info, "romio_cb_read", "automatic", STATUS)   
    !    STATUS = NF_OPEN_PAR   (trim(InCNRestart),IOR(NF_NOWRITE,NF_MPIIO),MPI_COMM_WORLD, info,NCFID)
    !    STATUS = NF_OPEN_PAR   (trim(OutFileName),IOR(NF_WRITE  ,NF_MPIIO),MPI_COMM_WORLD, info,OUTID)
    
    STATUS = NF_OPEN_PAR   (trim(OutFileName),IOR(NF_NOWRITE,NF_MPIIO),MPI_COMM_WORLD, infos,OUTID) ; VERIFY_(STATUS)
 !   if(root_proc) then
 !      STATUS = NF_OPEN (trim(OutFileName),NF_WRITE,OUTID)
 !      
 !   else
 !      STATUS = NF_OPEN (trim(OutFileName),NF_NOWRITE,OUTID)
 !   endif
 !   
    IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS, 'OUTPUT RESTART FAILED')

    STATUS = NF_GET_VARA_REAL(OUTID,VarID(OUTID,'ITY'), (/low_ind(myid+1),1/), (/nt_local(myid+1),1/),CLMC_pt1)
    STATUS = NF_GET_VARA_REAL(OUTID,VarID(OUTID,'ITY'), (/low_ind(myid+1),2/), (/nt_local(myid+1),1/),CLMC_pt2)
    STATUS = NF_GET_VARA_REAL(OUTID,VarID(OUTID,'ITY'), (/low_ind(myid+1),3/), (/nt_local(myid+1),1/),CLMC_st1)
    STATUS = NF_GET_VARA_REAL(OUTID,VarID(OUTID,'ITY'), (/low_ind(myid+1),4/), (/nt_local(myid+1),1/),CLMC_st2)
    STATUS = NF_GET_VARA_REAL(OUTID,VarID(OUTID,'FVG'), (/low_ind(myid+1),1/), (/nt_local(myid+1),1/),CLMC_pf1)
    STATUS = NF_GET_VARA_REAL(OUTID,VarID(OUTID,'FVG'), (/low_ind(myid+1),2/), (/nt_local(myid+1),1/),CLMC_pf2)
    STATUS = NF_GET_VARA_REAL(OUTID,VarID(OUTID,'FVG'), (/low_ind(myid+1),3/), (/nt_local(myid+1),1/),CLMC_sf1)
    STATUS = NF_GET_VARA_REAL(OUTID,VarID(OUTID,'FVG'), (/low_ind(myid+1),4/), (/nt_local(myid+1),1/),CLMC_sf2)
    
    if (root_proc) then
       
       STATUS = NF_OPEN (trim(InCNRestart),NF_NOWRITE,NCFID)
       IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS, 'OFFLINE RESTART FAILED')
       allocate (TILE_ID  (1:ntiles_cn))

       STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'TILE_ID'   ), (/1/), (/NTILES_cn/),TILE_ID)       

       do n = 1,ntiles_cn
  
          K = NINT (TILE_ID (n))

          STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'ITY'), (/n,1/), (/1,4/),ityp_offl(k,:))
          STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'FVG'), (/n,1/), (/1,4/),fveg_offl(k,:))
          
          tid_offl (n) = n
          
          do nv = 1,nveg
             if(ityp_offl(k,nv)<0 .or. ityp_offl(k,nv)>npft)    stop 'ityp'
             if(fveg_offl(k,nv)<0..or. fveg_offl(k,nv)>1.00001) stop 'fveg'             
          end do

          if((ityp_offl(k,3) == 0).and.(ityp_offl(k,4) == 0)) then
             if(ityp_offl(k,1) /= 0) then
                ityp_offl(k,3) = ityp_offl(k,1)
             else
                ityp_offl(k,3) = ityp_offl(k,2)
             endif
          endif
          
          if((ityp_offl(k,1) == 0).and.(ityp_offl(k,2) /= 0)) ityp_offl(k,1) = ityp_offl(k,2)
          if((ityp_offl(k,2) == 0).and.(ityp_offl(k,1) /= 0)) ityp_offl(k,2) = ityp_offl(k,1)
          if((ityp_offl(k,3) == 0).and.(ityp_offl(k,4) /= 0)) ityp_offl(k,3) = ityp_offl(k,4)
          if((ityp_offl(k,4) == 0).and.(ityp_offl(k,3) /= 0)) ityp_offl(k,4) = ityp_offl(k,3)
          
       end do
       
    endif
    
    call MPI_BCAST(tid_offl ,size(tid_offl ),MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_BCAST(ityp_offl,size(ityp_offl),MPI_REAL   ,0,MPI_COMM_WORLD,mpierr)
    call MPI_BCAST(fveg_offl,size(fveg_offl),MPI_REAL   ,0,MPI_COMM_WORLD,mpierr)    

    ! --------------------------------------------------------------------------------
    ! Here we create transfer index array to map offline restarts to output tile space
    ! --------------------------------------------------------------------------------   

    call GetIds(lonc,latc,lonn,latt,id_loc, tid_offl, &
         CLMC_pf1, CLMC_pf2, CLMC_sf1, CLMC_sf2, CLMC_pt1, CLMC_pt2,CLMC_st1,CLMC_st2, &
         fveg_offl, ityp_offl)
    deallocate (CLMC_pf1, CLMC_pf2, CLMC_sf1, CLMC_sf2, CLMC_pt1, CLMC_pt2,CLMC_st1,CLMC_st2,lonc,latc,lonn,latt)
 
     ! update id_glb in root

     if(root_proc)  then
        allocate (id_glb  (ntiles,4))
        allocate (id_vec  (ntiles))
     endif
           
     do nv = 1, nveg
        call MPI_Barrier(MPI_COMM_WORLD, STATUS)
!        call MPI_GATHERV( &
!                   id_loc (:,nv), nt_local(myid+1)  , MPI_real, &
!                   id_vec, nt_local,low_ind-1, MPI_real, &
!                   0, MPI_COMM_WORLD, mpierr )

        do i = 1, numprocs
           if((I == 1).and.(myid == 0)) then
              id_vec(low_ind(i) : upp_ind(i)) = Id_loc(:,nv)
           else if (I > 1) then
              if(I-1 == myid) then
                 ! send to root
                 call MPI_ISend(id_loc(:,nv),nt_local(i),MPI_INTEGER,0,993,MPI_COMM_WORLD,req,mpierr)
                 call MPI_WAIT (req,MPI_STATUS_IGNORE,mpierr)    
              else if (myid == 0) then
                 ! root receives
                 call MPI_RECV(id_vec(low_ind(i) : upp_ind(i)),nt_local(i) , MPI_INTEGER, i-1,993,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
              endif
           endif
        end do
       
        if(root_proc) id_glb (:,nv) = id_vec
        
     end do
     
     call MPI_Barrier(MPI_COMM_WORLD, STATUS)
     STATUS = NF_CLOSE (OutID)
! write out regridded carbon variables

     if(root_proc) then
        
        STATUS = NF_OPEN (trim(OutFileName),NF_WRITE,OUTID) ; VERIFY_(STATUS)
        allocate (CLMC_pf1(NTILES))
        allocate (CLMC_pf2(NTILES))
        allocate (CLMC_sf1(NTILES))
        allocate (CLMC_sf2(NTILES))
        allocate (CLMC_pt1(NTILES))
        allocate (CLMC_pt2(NTILES))
        allocate (CLMC_st1(NTILES))
        allocate (CLMC_st2(NTILES))
        allocate (VAR_DUM (NTILES))
        allocate (var_dum2 (1:ntiles_cn))

        STATUS = NF_GET_VARA_REAL(OUTID,VarID(OUTID,'ITY'), (/1,1/), (/NTILES,1/),CLMC_pt1)
        STATUS = NF_GET_VARA_REAL(OUTID,VarID(OUTID,'ITY'), (/1,2/), (/NTILES,1/),CLMC_pt2)
        STATUS = NF_GET_VARA_REAL(OUTID,VarID(OUTID,'ITY'), (/1,3/), (/NTILES,1/),CLMC_st1)
        STATUS = NF_GET_VARA_REAL(OUTID,VarID(OUTID,'ITY'), (/1,4/), (/NTILES,1/),CLMC_st2)
        STATUS = NF_GET_VARA_REAL(OUTID,VarID(OUTID,'FVG'), (/1,1/), (/NTILES,1/),CLMC_pf1)
        STATUS = NF_GET_VARA_REAL(OUTID,VarID(OUTID,'FVG'), (/1,2/), (/NTILES,1/),CLMC_pf2)
        STATUS = NF_GET_VARA_REAL(OUTID,VarID(OUTID,'FVG'), (/1,3/), (/NTILES,1/),CLMC_sf1)
        STATUS = NF_GET_VARA_REAL(OUTID,VarID(OUTID,'FVG'), (/1,4/), (/NTILES,1/),CLMC_sf2)

        allocate (var_off_col (1: NTILES_CN, 1 : nzone,1 : var_col))
        allocate (var_off_pft (1: NTILES_CN, 1 : nzone,1 : nveg, 1 : var_pft))
        
        allocate (var_col_out (1: NTILES, 1 : nzone,1 : var_col))
        allocate (var_pft_out (1: NTILES, 1 : nzone,1 : nveg, 1 : var_pft)) 

        i = 1
        do nv = 1,VAR_COL
           do nz = 1,nzone
              STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'CNCOL'), (/1,i/), (/NTILES_CN,1 /),VAR_DUM2)
              do k = 1, NTILES_CN
                 var_off_col(TILE_ID(K), nz,nv) = VAR_DUM2(K)
              end do
              i = i + 1
           end do
        end do
        
        i = 1
        do iv = 1,VAR_PFT
           do nv = 1,nveg
              do nz = 1,nzone
                 STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'CNPFT'), (/1,i/), (/NTILES_CN,1 /),VAR_DUM2)
                 do k = 1, NTILES_CN
                    var_off_pft(TILE_ID(K), nz,nv,iv) = VAR_DUM2(K)
                 end do
                 i = i + 1
              end do
           end do
        end do

        var_col_out = 0.
        var_pft_out = NaN
        
        where(isnan(var_off_pft))  var_off_pft = 0.   
        where(var_off_pft /= var_off_pft)  var_off_pft = 0.  

        OUT_TILE : DO N = 1, NTILES

           !if(mod (n,1000) == 0) print *, myid +1, n, Id_glb(n,:)

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
           
           ! column vars
           ! -----------
           !  1 clm3%g%l%c%ccs%col_ctrunc	    
           !  2 clm3%g%l%c%ccs%cwdc			    
           !  3 clm3%g%l%c%ccs%litr1c 	    
           !  4 clm3%g%l%c%ccs%litr2c 	    
           !  5 clm3%g%l%c%ccs%litr3c 	    
           !  6 clm3%g%l%c%ccs%pcs_a%totvegc      
           !  7 clm3%g%l%c%ccs%prod100c	    
           !  8 clm3%g%l%c%ccs%prod10c	    
           !  9 clm3%g%l%c%ccs%seedc  	    
           ! 10 clm3%g%l%c%ccs%soil1c 	    
           ! 11 clm3%g%l%c%ccs%soil2c 	    
           ! 12 clm3%g%l%c%ccs%soil3c 	    
           ! 13 clm3%g%l%c%ccs%soil4c 	    
           ! 14 clm3%g%l%c%ccs%totcolc	    
           ! 15 clm3%g%l%c%ccs%totlitc	    
           ! 16 clm3%g%l%c%cns%col_ntrunc	    
           ! 17 clm3%g%l%c%cns%cwdn			    
           ! 18 clm3%g%l%c%cns%litr1n 	    
           ! 19 clm3%g%l%c%cns%litr2n 	    
           ! 20 clm3%g%l%c%cns%litr3n 	    
           ! 21 clm3%g%l%c%cns%prod100n	    
           ! 22 clm3%g%l%c%cns%prod10n	    
           ! 23 clm3%g%l%c%cns%seedn  	    
           ! 24 clm3%g%l%c%cns%sminn  	    
           ! 25 clm3%g%l%c%cns%soil1n 	    
           ! 26 clm3%g%l%c%cns%soil2n 	    
           ! 27 clm3%g%l%c%cns%soil3n 	    
           ! 28 clm3%g%l%c%cns%soil4n 	    
           ! 29 clm3%g%l%c%cns%totcoln	    
           ! 30 clm3%g%l%c%cps%ann_farea_burned   
           ! 31 clm3%g%l%c%cps%annsum_counter     
           ! 32 clm3%g%l%c%cps%cannavg_t2m		    
           ! 33 clm3%g%l%c%cps%cannsum_npp		    
           ! 34 clm3%g%l%c%cps%farea_burned		    
           ! 35 clm3%g%l%c%cps%fire_prob	    
           ! 36 clm3%g%l%c%cps%fireseasonl		    
           ! 37 clm3%g%l%c%cps%fpg			    
           ! 38 clm3%g%l%c%cps%fpi			    
           ! 39 clm3%g%l%c%cps%me		    
           ! 40 clm3%g%l%c%cps%mean_fire_prob     
           
           ! PFT vars
           ! --------
           !  1 clm3%g%l%c%p%pcs%cpool	       
           !  2 clm3%g%l%c%p%pcs%deadcrootc	       
           !  3 clm3%g%l%c%p%pcs%deadcrootc_storage	
           !  4 clm3%g%l%c%p%pcs%deadcrootc_xfer  
           !  5 clm3%g%l%c%p%pcs%deadstemc 	       
           !  6 clm3%g%l%c%p%pcs%deadstemc_storage 	
           !  7 clm3%g%l%c%p%pcs%deadstemc_xfer   
           !  8 clm3%g%l%c%p%pcs%frootc	       
           !  9 clm3%g%l%c%p%pcs%frootc_storage   
           ! 10 clm3%g%l%c%p%pcs%frootc_xfer      
           ! 11 clm3%g%l%c%p%pcs%gresp_storage    
           ! 12 clm3%g%l%c%p%pcs%gresp_xfer	       
           ! 13 clm3%g%l%c%p%pcs%leafc	       
           ! 14 clm3%g%l%c%p%pcs%leafc_storage	
           ! 15 clm3%g%l%c%p%pcs%leafc_xfer	       
           ! 16 clm3%g%l%c%p%pcs%livecrootc		
           ! 17 clm3%g%l%c%p%pcs%livecrootc_storage	
           ! 18 clm3%g%l%c%p%pcs%livecrootc_xfer  
           ! 19 clm3%g%l%c%p%pcs%livestemc 	       
           ! 20 clm3%g%l%c%p%pcs%livestemc_storage 	
           ! 21 clm3%g%l%c%p%pcs%livestemc_xfer   
           ! 22 clm3%g%l%c%p%pcs%pft_ctrunc	       
           ! 23 clm3%g%l%c%p%pcs%xsmrpool         
           ! 24 clm3%g%l%c%p%pepv%annavg_t2m      
           ! 25 clm3%g%l%c%p%pepv%annmax_retransn 
           ! 26 clm3%g%l%c%p%pepv%annsum_npp      
           ! 27 clm3%g%l%c%p%pepv%annsum_potential_gpp 
           ! 28 clm3%g%l%c%p%pepv%dayl	       
           ! 29 clm3%g%l%c%p%pepv%days_active     
           ! 30 clm3%g%l%c%p%pepv%dormant_flag    
           ! 31 clm3%g%l%c%p%pepv%offset_counter  
           ! 32 clm3%g%l%c%p%pepv%offset_fdd      
           ! 33 clm3%g%l%c%p%pepv%offset_flag     
           ! 34 clm3%g%l%c%p%pepv%offset_swi      
           ! 35 clm3%g%l%c%p%pepv%onset_counter   
           ! 36 clm3%g%l%c%p%pepv%onset_fdd	       
           ! 37 clm3%g%l%c%p%pepv%onset_flag      
           ! 38 clm3%g%l%c%p%pepv%onset_gdd	       
           ! 39 clm3%g%l%c%p%pepv%onset_gddflag   
           ! 40 clm3%g%l%c%p%pepv%onset_swi	       
           ! 41 clm3%g%l%c%p%pepv%prev_frootc_to_litter
           ! 42 clm3%g%l%c%p%pepv%prev_leafc_to_litter 
           ! 43 clm3%g%l%c%p%pepv%tempavg_t2m     
           ! 44 clm3%g%l%c%p%pepv%tempmax_retransn 	
           ! 45 clm3%g%l%c%p%pepv%tempsum_npp     
           ! 46 clm3%g%l%c%p%pepv%tempsum_potential_gpp
           ! 47 clm3%g%l%c%p%pepv%xsmrpool_recover 	
           ! 48 clm3%g%l%c%p%pns%deadcrootn	       
           ! 49 clm3%g%l%c%p%pns%deadcrootn_storage	
           ! 50 clm3%g%l%c%p%pns%deadcrootn_xfer  
           ! 51 clm3%g%l%c%p%pns%deadstemn 	       
           ! 52 clm3%g%l%c%p%pns%deadstemn_storage 	
           ! 53 clm3%g%l%c%p%pns%deadstemn_xfer   
           ! 54 clm3%g%l%c%p%pns%frootn	       
           ! 55 clm3%g%l%c%p%pns%frootn_storage   
           ! 56 clm3%g%l%c%p%pns%frootn_xfer      
           ! 57 clm3%g%l%c%p%pns%leafn	       
           ! 58 clm3%g%l%c%p%pns%leafn_storage    
           ! 59 clm3%g%l%c%p%pns%leafn_xfer	       
           ! 60 clm3%g%l%c%p%pns%livecrootn	       
           ! 61 clm3%g%l%c%p%pns%livecrootn_storage	
           ! 62 clm3%g%l%c%p%pns%livecrootn_xfer  
           ! 63 clm3%g%l%c%p%pns%livestemn 	       
           ! 64 clm3%g%l%c%p%pns%livestemn_storage 	
           ! 65 clm3%g%l%c%p%pns%livestemn_xfer   
           ! 66 clm3%g%l%c%p%pns%npool	    
           ! 67 clm3%g%l%c%p%pns%pft_ntrunc		    
           ! 68 clm3%g%l%c%p%pns%retransn	    
           ! 69 clm3%g%l%c%p%pps%elai 	    
           ! 70 clm3%g%l%c%p%pps%esai 	    
           ! 71 clm3%g%l%c%p%pps%hbot 	      
           ! 72 clm3%g%l%c%p%pps%htop 	      
           ! 73 clm3%g%l%c%p%pps%tlai 	    
           ! 74 clm3%g%l%c%p%pps%tsai 	    
           
        end do OUT_TILE
        
        i = 1
        do nv = 1,VAR_COL
           do nz = 1,nzone
              STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'CNCOL'), (/1,i/), (/NTILES,1 /),var_col_out(:, nz,nv))
              i = i + 1
           end do
        end do
        
        i = 1
        do iv = 1,VAR_PFT
           do nv = 1,nveg
              do nz = 1,nzone
                 STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'CNPFT'), (/1,i/), (/NTILES,1 /),var_pft_out(:, nz,nv,iv))
                 i = i + 1
              end do
           end do
        end do

        VAR_DUM = 0.

        do nz = 1,nzone
           STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'TGWM'), (/1,nz/), (/NTILES,1 /),VAR_DUM(:))
           STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'RZMM'), (/1,nz/), (/NTILES,1 /),VAR_DUM(:))           
        end do

        STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'SFMCM'), (/1/), (/NTILES/),VAR_DUM(:))
        STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'BFLOWM'), (/1/), (/NTILES/),VAR_DUM(:))
        STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'TOTWATM'), (/1/), (/NTILES/),VAR_DUM(:))
        STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'TAIRM'), (/1/), (/NTILES/),VAR_DUM(:))
        STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'TPM'), (/1/), (/NTILES/),VAR_DUM(:))
        STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'CNSUM'), (/1/), (/NTILES/),VAR_DUM(:)) 
        STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'SNDZM'), (/1/), (/NTILES/),VAR_DUM(:))
        STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'ASNOWM'), (/1/), (/NTILES/),VAR_DUM(:))

        do nv = 1,nzone
           do nz = 1,nveg
              STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'PSNSUNM'), (/1,nz,nv/), (/NTILES,1,1/),VAR_DUM(:)) 
              STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'PSNSHAM'), (/1,nz,nv/), (/NTILES,1,1/),VAR_DUM(:)) 
           end do
        end do

        STATUS = NF_CLOSE (NCFID)
        STATUS = NF_CLOSE (OutID)
        
        deallocate (var_off_col,var_off_pft,var_col_out,var_pft_out)  
        deallocate (CLMC_pf1, CLMC_pf2, CLMC_sf1, CLMC_sf2)
        deallocate (CLMC_pt1, CLMC_pt2, CLMC_st1, CLMC_st2)
 
     endif

     call MPI_Barrier(MPI_COMM_WORLD, STATUS)

     
  END SUBROUTINE regrid_carbon_vars

 ! ***************************************************************************** 

  SUBROUTINE NCDF_reshape_getOput (NCFID,CID,col,pft, get_var) 

    implicit none

    integer, intent (in) :: NCFID,CID
    logical, intent (in) :: get_var
    real, intent (inout) :: col (nzone * VAR_COL)
    real, intent (inout) :: pft (nzone * nveg * var_PFT)
    integer              :: STATUS

    if (get_var) then
       STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'CNCOL'), (/CID,1/), (/1,nzone * VAR_COL       /),col)
       STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'CNPFT'), (/CID,1/), (/1,nzone * nveg * var_PFT/),pft)
    else
       STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'CNCOL'), (/CID,1/), (/1,nzone * VAR_COL       /),col)
       STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'CNPFT'), (/CID,1/), (/1,nzone * nveg * var_PFT/),pft)
    endif

    IF ((STATUS .NE. NF_NOERR).and.(get_var)) then
       print *,CID
       CALL HANDLE_ERR(STATUS, 'Out : NCDF_reshape_getOput')  
    ENDIF

    IF ((STATUS .NE. NF_NOERR).and.(.not.get_var)) then
       print *,CID
       CALL HANDLE_ERR(STATUS, 'In : NCDF_reshape_getOput')  
    ENDIF
  END SUBROUTINE NCDF_reshape_getOput

 ! ***************************************************************************** 

  SUBROUTINE NCDF_whole_getOput (NCFID,NTILES,col,pft, get_var) 

    implicit none

    integer, intent (in) :: NCFID,NTILES
    logical, intent (in) :: get_var
    real, intent (inout) :: col (NTILES, nzone * VAR_COL)
    real, intent (inout) :: pft (NTILES, nzone * nveg * var_PFT)
    integer              :: STATUS, J

    if (get_var) then
       DO J = 1,nzone * VAR_COL
          STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'CNCOL'), (/1,J/), (/NTILES,1 /),col(:,j))
       END DO
       DO J = 1, nzone * nveg * var_PFT
          STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'CNPFT'), (/1,J/), (/NTILES,1/),pft(:,J))
       END DO
    else
       DO J = 1,nzone * VAR_COL
          STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'CNCOL'), (/1,J/), (/NTILES,1 /),col(:,J))
       END DO
       DO J = 1, nzone * nveg * var_PFT
          STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'CNPFT'), (/1,J/), (/NTILES,1/) ,pft(:,J))
       END DO
    endif

    IF ((STATUS .NE. NF_NOERR).and.(get_var)) CALL HANDLE_ERR(STATUS, 'Out : NCDF_whole_getOput')  
    IF ((STATUS .NE. NF_NOERR).and.(.not.get_var)) CALL HANDLE_ERR(STATUS, 'In : NCDF_whole_getOput')  

  END SUBROUTINE NCDF_whole_getOput
  
  ! -----------------------------------------------------------------------

   SUBROUTINE HANDLE_ERR(STATUS, Line)

     INTEGER,      INTENT (IN) :: STATUS
     CHARACTER(*), INTENT (IN) :: Line

     IF (STATUS .NE. NF_NOERR) THEN
        PRINT *, trim(Line),': ',NF_STRERROR(STATUS)
        STOP 'Stopped'
     ENDIF

   END SUBROUTINE HANDLE_ERR

  ! *****************************************************************************

   integer function VarID (NCFID, VNAME) 
     
     integer, intent (in)      :: NCFID
     character(*), intent (in) :: VNAME
     integer                   :: status

     STATUS = NF_INQ_VARID (NCFID, trim(VNAME) ,VarID)
     IF (STATUS .NE. NF_NOERR) &
          CALL HANDLE_ERR(STATUS, trim(VNAME))  
     
   end function VarID

  ! *****************************************************************************
  
  SUBROUTINE regrid_hyd_vars (NTILES,  OutFMT) 

    implicit none
    integer, intent (in)           :: NTILES

    ! ===============================================================================================

    integer, allocatable, dimension(:)   :: Id_glb, Id_loc
    integer, allocatable, dimension(:)   :: ld_reorder, tid_offl
    real   , allocatable, dimension(:)   :: tmp_var
    integer :: n,i,j, nv, nx, offl_cell, STATUS,NCFID, req
    integer :: outid, local_id
    real   , pointer, dimension (:) :: lonc, latc, LATT, LONN, long, latg
    integer, allocatable :: low_ind(:), upp_ind(:), nt_local (:)
    type(Netcdf4_Fileformatter) :: InFmt, OutFmt

    allocate (tid_offl  (ntiles_cn))
    allocate (tmp_var   (ntiles_cn))

    allocate(low_ind (   numprocs))
    allocate(upp_ind (   numprocs))
    allocate(nt_local(   numprocs))

    low_ind (:)    = 1
    upp_ind (:)    = NTILES       
    nt_local(:)    = NTILES 

    ! Domain decomposition
    ! --------------------

    if (numprocs > 1) then      
       do i = 1, numprocs - 1
          upp_ind(i)   = low_ind(i) + (ntiles/numprocs) - 1 
          low_ind(i+1) = upp_ind(i) + 1
          nt_local(i)  = upp_ind(i) - low_ind(i) + 1
       end do
       nt_local(numprocs) = upp_ind(numprocs) - low_ind(numprocs) + 1
    endif

    allocate (id_loc (nt_local (myid + 1)))
    allocate (lonn   (nt_local (myid + 1)))
    allocate (latt   (nt_local (myid + 1)))
    allocate (lonc   (1:ntiles_cn))
    allocate (latc   (1:ntiles_cn))

    if (root_proc) then

       allocate (long   (ntiles))
       allocate (latg   (ntiles))
       allocate (ld_reorder(ntiles_cn)) 

       call ReadTileFile_RealLatLon (OutTileFile, i, xlon=long, xlat=latg)

       ! ---------------------------------------------
       ! Read exact lonc, latc from offline .til File 
       ! ---------------------------------------------

       call ReadTileFile_RealLatLon(trim(InCNTilFile), i,xlon=lonc,xlat=latc)

       STATUS = NF_OPEN (trim(InCNRestart),NF_NOWRITE,NCFID)
       STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'TILE_ID'   ), (/1/), (/NTILES_CN/),tmp_var)
       STATUS = NF_CLOSE (NCFID)
 
       do n = 1, ntiles_cn
          ld_reorder ( NINT(tmp_var(n))) = n
          tid_offl(n)    = n
       end do

       deallocate (tmp_var)

    endif

    call MPI_Barrier(MPI_COMM_WORLD, STATUS)

    do i = 1, numprocs
       if((I == 1).and.(myid == 0)) then
          lonn(:) = long(low_ind(i) : upp_ind(i))
          latt(:) = latg(low_ind(i) : upp_ind(i))
       else if (I > 1) then
          if(I-1 == myid) then
             ! receiving from root
             call MPI_RECV(lonn,nt_local(i) , MPI_REAL, 0,995,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
             call MPI_RECV(latt,nt_local(i) , MPI_REAL, 0,994,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
          else if (myid == 0) then
             ! root sends
             call MPI_ISend(long(low_ind(i) : upp_ind(i)),nt_local(i),MPI_REAL,i-1,995,MPI_COMM_WORLD,req,mpierr)
             call MPI_WAIT (req,MPI_STATUS_IGNORE,mpierr)    
             call MPI_ISend(latg(low_ind(i) : upp_ind(i)),nt_local(i),MPI_REAL,i-1,994,MPI_COMM_WORLD,req,mpierr)
             call MPI_WAIT (req,MPI_STATUS_IGNORE,mpierr) 
          endif
       endif
    end do

!    call MPI_SCATTERV (                    &
!         long,nt_local,low_ind-1,MPI_real, &
!         lonn,size(lonn),MPI_real  , &
!         0,MPI_COMM_WORLD, mpierr )
!
!    call MPI_SCATTERV (                    &
!         latg,nt_local,low_ind-1,MPI_real, &
!         latt,nt_local(myid+1),MPI_real  , &
!         0,MPI_COMM_WORLD, mpierr )

    if(root_proc) deallocate (long, latg)
     
     call MPI_BCAST(lonc,ntiles_cn,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
     call MPI_BCAST(latc,ntiles_cn,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
     call MPI_BCAST(tid_offl,size(tid_offl  ),MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)

    ! --------------------------------------------------------------------------------
    ! Here we create transfer index array to map offline restarts to output tile space
    ! --------------------------------------------------------------------------------   

     call GetIds(lonc,latc,lonn,latt,id_loc, tid_offl)      

     ! Loop through NTILES (# of tiles in output array) find the nearest neighbor from Qing.  

     if(root_proc)  allocate (id_glb  (ntiles))

     call MPI_Barrier(MPI_COMM_WORLD, STATUS)

     do i = 1, numprocs
        if((I == 1).and.(myid == 0)) then
           id_glb(low_ind(i) : upp_ind(i)) = Id_loc(:)
        else if (I > 1) then
           if(I-1 == myid) then
              ! send to root
              call MPI_ISend(id_loc,nt_local(i),MPI_INTEGER,0,993,MPI_COMM_WORLD,req,mpierr)
              call MPI_WAIT (req,MPI_STATUS_IGNORE,mpierr)    
           else if (myid == 0) then
              ! root receives
              call MPI_RECV(id_glb(low_ind(i) : upp_ind(i)),nt_local(i) , MPI_INTEGER, i-1,993,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
           endif
        endif
     end do     

!     call MPI_GATHERV( &
!                   id_loc, nt_local(myid+1)  , MPI_real, &
!                   id_glb, nt_local,low_ind-1, MPI_real, &
!                   0, MPI_COMM_WORLD, mpierr )
        
    if (root_proc) call put_land_vars  (NTILES, id_glb, ld_reorder, OutFmt)

    call MPI_Barrier(MPI_COMM_WORLD, STATUS)

   END SUBROUTINE regrid_hyd_vars

  ! *****************************************************************************
   SUBROUTINE put_land_vars (NTILES, id_glb, ld_reorder, OutFmt)

     implicit none
     
     integer, intent (in)       :: NTILES
     integer, intent (in)       :: id_glb(NTILES), ld_reorder (NTILES_CN)
     integer                    :: i,k,n
     real   , dimension (:), allocatable :: var_get, var_put
     type(Netcdf4_Fileformatter) :: OutFmt
     integer         :: nVars, STATUS, NCFID

     allocate (var_get (NTILES_CN))
     allocate (var_put (NTILES))
   
     ! Read catparam
     ! -------------

     STATUS = NF_OPEN (trim(InCNRestart),NF_NOWRITE,NCFID)

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'POROS'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'POROS',var_put)

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'COND'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'COND',var_put)

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'PSIS'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'PSIS',var_put)

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'BEE'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'BEE',var_put) 
     
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'WPWET'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'WPWET',var_put) 
    
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'GNU'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'GNU',var_put)
 
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'VGWMAX'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'VGWMAX',var_put) 
 
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'BF1'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'BF1',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'BF2'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'BF2',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'BF3'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'BF3',var_put) 
     
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'CDCR1'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'CDCR1',var_put) 
     
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'CDCR2'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'CDCR2',var_put) 
     
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'ARS1'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'ARS1',var_put) 
     
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'ARS2'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'ARS2',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'ARS3'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'ARS3',var_put) 
     
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'ARA1'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'ARA1',var_put) 
     
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'ARA2'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'ARA2',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'ARA3'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'ARA3',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'ARA4'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'ARA4',var_put)

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'ARW1'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'ARW1',var_put) 
     
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'ARW2'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'ARW2',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'ARW3'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'ARW3',var_put) 
     
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'ARW4'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'ARW4',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'TSA1'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'TSA1',var_put) 
     
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'TSA2'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'TSA2',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'TSB1'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'TSB1',var_put) 
     
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'TSB2'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'TSB2',var_put) 
     
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'ATAU'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'ATAU',var_put) 
     
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'BTAU'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'BTAU',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'ITY'   ), (/1,1/), (/NTILES_CN,1/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'ITY',var_put, offset1=1) 
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'ITY'   ), (/1,2/), (/NTILES_CN,1/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'ITY',var_put, offset1=2) 
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'ITY'   ), (/1,3/), (/NTILES_CN,1/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'ITY',var_put, offset1=3) 
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'ITY'   ), (/1,4/), (/NTILES_CN,1/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'ITY',var_put, offset1=4)

      STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'FVG'   ), (/1,1/), (/NTILES_CN,1/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'FVG',var_put, offset1=1) 
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'FVG'   ), (/1,2/), (/NTILES_CN,1/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'FVG',var_put, offset1=2) 
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'FVG'   ), (/1,3/), (/NTILES_CN,1/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'FVG',var_put, offset1=3) 
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'FVG'   ), (/1,4/), (/NTILES_CN,1/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'FVG',var_put, offset1=4)

     ! read restart and regrid
     ! -----------------------

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'TG'   ), (/1,1/), (/NTILES_CN,1/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'TG',var_put, offset1=1)  ! if you see offset1=1 it is a 2-D var

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'TG'   ), (/1,2/), (/NTILES_CN,1/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'TG',var_put, offset1=2) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'TG'   ), (/1,3/), (/NTILES_CN,1/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'TG',var_put, offset1=3) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'TC'   ), (/1,1/), (/NTILES_CN,1/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'TC',var_put, offset1=1) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'TC'   ), (/1,2/), (/NTILES_CN,1/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'TC',var_put, offset1=2) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'TC'   ), (/1,3/), (/NTILES_CN,1/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'TC',var_put, offset1=3) 
     
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'QC'   ), (/1,1/), (/NTILES_CN,1/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'QC',var_put, offset1=1) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'QC'   ), (/1,2/), (/NTILES_CN,1/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'QC',var_put, offset1=2)

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'QC'   ), (/1,3/), (/NTILES_CN,1/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'QC',var_put, offset1=3)

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'CAPAC'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'CAPAC',var_put) 
     
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'CATDEF'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'CATDEF',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'RZEXC'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'RZEXC',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'SRFEXC'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'SRFEXC',var_put) 
     
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'GHTCNT1'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'GHTCNT1',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'GHTCNT2'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'GHTCNT2',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'GHTCNT3'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'GHTCNT3',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'GHTCNT4'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'GHTCNT4',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'GHTCNT5'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'GHTCNT5',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'GHTCNT6'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'GHTCNT6',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'WESNN1'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'WESNN1',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'WESNN2'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'WESNN2',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'WESNN3'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'WESNN3',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'HTSNNN1'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'HTSNNN1',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'HTSNNN2'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'HTSNNN2',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'HTSNNN3'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'HTSNNN3',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'SNDZN1'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'SNDZN1',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'SNDZN2'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'SNDZN2',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'SNDZN3'   ), (/1/), (/NTILES_CN/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(OutFmt,'SNDZN3',var_put) 

     STATUS = NF_CLOSE ( NCFID)

     deallocate (var_get, var_put)

   END SUBROUTINE put_land_vars

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

  ! *****************************************************************************
     
end program mk_CatchCNRestarts

