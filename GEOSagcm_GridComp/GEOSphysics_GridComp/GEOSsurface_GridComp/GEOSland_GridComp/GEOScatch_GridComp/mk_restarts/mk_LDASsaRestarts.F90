PROGRAM mk_LDASsaRestarts

! USAGE/HELP (NOTICE mpirun -np 1) 
!  mpirun -np 1 bin/mk_LDASsaRestarts -h
!
! (1) to create an initial catch(cn)_internal_rst file ready for an offline experiment :
! --------------------------------------------------------------------------------------
! (1.1) mpirun -np 1 bin/mk_LDASsaRestarts -a SPONSORCODE -b BCSDIR -m MODEL -s SURFLAY(20/50) -t TILFILE
! where MODEL : catch or catchcn
! (1.2) sbatch mkLDAS.j
!    
! (2) to reorder an LDASsa restart file to the order of the BCs for use in an GCM experiment :
! --------------------------------------------------------------------------------------------
! mpirun -np 1 bin/mk_LDASsaRestarts -b BCSDIR  -d YYYYMMDD -e EXPNAME -l EXPDIR -m MODEL -s SURFLAY(20/50) -r Y -t TILFILE

  use MAPL_ConstantsMod,only: MAPL_PI,  MAPL_radius, MAPL_TICE
  use MAPL_HashMod
  use MAPL_IOMod
  use ieee_arithmetic, only: isnan => ieee_is_nan
  USE STIEGLITZSNOW,   ONLY :                 &
       StieglitzSnow_calc_tpsnow 
  implicit none
  include 'mpif.h'
  INCLUDE 'netcdf.inc'
  
  ! initialize to non-MPI values
  
  integer  :: myid=0, numprocs=1, mpierr
  logical  :: master_proc=.true.
  
  ! Carbon model specifics
  ! ----------------------

  character*256 :: Usage="mk_LDASsaRestarts -a SPONSORCODE -b BCSDIR -d YYYYMMDD  -e EXPNAME -j JOBFILE -l EXPDIR -m MODEL -r REORDER -s SURFLAY -t TILFILE"
  character*256 :: BCSDIR, SPONSORCODE, EXPNAME, EXPDIR, MODEL, TILFILE, YYYYMMDD, SFL
  character*400 :: CMD

  real, parameter :: ECCENTRICITY  = 0.0167
  real, parameter :: PERIHELION    = 102.0
  real, parameter :: OBLIQUITY     = 23.45
  integer, parameter :: EQUINOX    = 80
  
  integer, parameter :: nveg    = 4
  integer, parameter :: nzone   = 3
  integer, parameter :: VAR_COL = 40 ! number of CN column restart variables
  integer, parameter :: VAR_PFT = 74 ! number of CN PFT variables per column
  integer, parameter :: npft    = 19  
  integer, parameter :: npft_clm45    = 27  
  integer, parameter :: VAR_COL_CLM45 = 39 ! number of CN column restart variables
  integer, parameter :: VAR_PFT_CLM45 = 104 ! number of CN PFT variables per column

  real,    parameter :: nan = O'17760000000'
  real,    parameter :: fmin= 1.e-4 ! ignore vegetation fractions at or below this value
  integer, parameter :: OutUnit = 40, InUnit = 50
  character*256      :: arg, tmpstring, ESMADIR
  character*1        :: opt,  REORDER='N', JOBFILE ='N'
  integer            :: ntiles, rc, nxt
  character(len=300) :: OutFileName

  ! ===============================================================================================
  ! Below hard-wired ldas restart file is from a global offline simulation on the SMAP M09 grid
  ! after 1000s of years of simulations

  integer, parameter :: ntiles_cn = 1684725, ntiles_cat = 1653157
  character(len=300), parameter :: &
       InCNRestart = '/gpfsm/dnb42/projects/p16/ssd/land/l_data/LandRestarts_for_Regridding/CatchCN/M09/20151231/catchcn_internal_rst', &
       InCNTilFile = '/discover/nobackup/ltakacs/bcs/Heracles-NL/SMAP_EASEv2_M09/SMAP_EASEv2_M09_3856x1624.til',                        &
       InCatRestart= '/gpfsm/dnb42/projects/p16/ssd/land/l_data/LandRestarts_for_Regridding/Catch/M09/20170101/catch_internal_rst', &
       InCatTilFile= '/discover/nobackup/projects/gmao/ssd/land/l_data/geos5/bcs/CLSM_params/mkCatchParam_SMAP_L4SM_v002/' &
                      //'SMAP_EASEv2_M09/SMAP_EASEv2_M09_3856x1624.til'            

  real        :: SURFLAY = 50
  integer     :: STATUS

  character(len=256), parameter :: CatNames   (57) = &
       (/'BF1    ',  'BF2    ',  'BF3    ',  'VGWMAX ',  'CDCR1  ', &
         'CDCR2  ',  'PSIS   ',  'BEE    ',  'POROS  ',  'WPWET  ', &
         'COND   ',  'GNU    ',  'ARS1   ',  'ARS2   ',  'ARS3   ', &
         'ARA1   ',  'ARA2   ',  'ARA3   ',  'ARA4   ',  'ARW1   ', &
         'ARW2   ',  'ARW3   ',  'ARW4   ',  'TSA1   ',  'TSA2   ', &
         'TSB1   ',  'TSB2   ',  'ATAU   ',  'BTAU   ',  'OLD_ITY', &
         'TC     ',  'QC     ',  'CAPAC  ',  'CATDEF ',  'RZEXC  ', &
         'SRFEXC ',  'GHTCNT1',  'GHTCNT2',  'GHTCNT3',  'GHTCNT4', &
         'GHTCNT5',  'GHTCNT6',  'TSURF  ',  'WESNN1 ',  'WESNN2 ', &
         'WESNN3 ',  'HTSNNN1',  'HTSNNN2',  'HTSNNN3',  'SNDZN1 ', &
         'SNDZN2 ',  'SNDZN3 ',  'CH     ',  'CM     ',  'CQ     ', &
         'FR     ',  'WW     '/)

  character(len=256), parameter :: CarbNames (68) =  &
       (/'BF1    ',  'BF2    ',  'BF3    ',  'VGWMAX ',  'CDCR1  ', &
         'CDCR2  ',  'PSIS   ',  'BEE    ',  'POROS  ',  'WPWET  ', &
         'COND   ',  'GNU    ',  'ARS1   ',  'ARS2   ',  'ARS3   ', &
         'ARA1   ',  'ARA2   ',  'ARA3   ',  'ARA4   ',  'ARW1   ', &
         'ARW2   ',  'ARW3   ',  'ARW4   ',  'TSA1   ',  'TSA2   ', &
         'TSB1   ',  'TSB2   ',  'ATAU   ',  'BTAU   ',  'ITY    ', &
         'FVG    ',  'TC     ',  'QC     ',  'TG     ',  'CAPAC  ', &
         'CATDEF ',  'RZEXC  ',  'SRFEXC ',  'GHTCNT1',  'GHTCNT2', &
         'GHTCNT3',  'GHTCNT4',  'GHTCNT5',  'GHTCNT6',  'TSURF  ', &
         'WESNN1 ',  'WESNN2 ',  'WESNN3 ',  'HTSNNN1',  'HTSNNN2', &
         'HTSNNN3',  'SNDZN1 ',  'SNDZN2 ',  'SNDZN3 ',  'CH     ', &
         'CM     ',  'CQ     ',  'FR     ',  'WW     ',  'TILE_ID', &
         'NDEP   ',  'CLI_T2M',  'BGALBVR',  'BGALBVF',  'BGALBNR', &
         'BGALBNF',  'CNCOL  ',  'CNPFT  '   /)

  CHARACTER( * ), PARAMETER :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
  CHARACTER( * ), PARAMETER :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  logical, parameter  :: clm45 = .false.
  logical             :: second_visit

  call init_MPI()

  ! process commands
  ! ----------------

  CALL get_command (cmd)
  call getenv ("ESMADIR"        ,ESMADIR        )
  nxt = 1
  
  call getarg(nxt,arg)

  do while(arg(1:1)=='-')

     opt=arg(2:2)
     if(len(trim(arg))==2) then
        nxt = nxt + 1
        call getarg(nxt,arg)
     else
        arg = arg(3:)
     end if

     select case (opt)
     case ('a')
        SPONSORCODE = trim(arg)
     case ('b')
        BCSDIR = trim(arg)
     case ('d')
        YYYYMMDD = trim(arg)
     case ('e')
        EXPNAME = trim(arg)
     case ('h')
        print *,'   '
        print *,'(1) to create an initial catch(cn)_internal_rst file ready for an offline experiment :'
        print *,'--------------------------------------------------------------------------------------'
        print *,'(1.1) mpirun -np 1 bin/mk_LDASsaRestarts -a SPONSORCODE -b BCSDIR -m MODEL -s SURFLAY(20/50)'
        print *,'where MODEL : catch or catchcn'
        print *,'(1.2) sbatch mkLDAS.j'
        print *,'   '
        print *,'(2) to reorder an LDASsa restart file to the order of the BCs for use in an GCM experiment :'
        print *,'--------------------------------------------------------------------------------------------'
        print *,'mpirun -np 1 bin/mk_LDASsaRestarts -b BCSDIR  -d YYYYMMDD -e EXPNAME -l EXPDIR -m MODEL -s SURFLAY(20/50) -r Y -t TILFILE'
        stop
     case ('j')
        JOBFILE = trim(arg)
     case ('l')
        EXPDIR = trim(arg) 
     case ('m')
        MODEL = StrUpCase(trim(arg))
     case ('r')
        REORDER = trim(arg)
     case ('s')
        SFL = trim(arg)
        read(arg,*)  SURFLAY
     case ('t')
        TILFILE =  trim(arg)
     case default
        print *, trim(Usage)
        call exit(1)
     end select
     nxt = nxt + 1
     call getarg(nxt,arg)
  end do

  if (trim(model) == 'CATCHCN') then
     if((INDEX(BCSDIR, 'Heracles') == 0).AND.(INDEX(BCSDIR, 'Icarus') == 0)) then
        print *,'Land BCs in : ',trim(BCSDIR)
        print *,'do not support ',trim (model)
        stop
     endif
  endif


  if(trim(REORDER) == 'Y') then

     ! This call is to reorder a LDASsa restart file

     call reorder_LDASsa_restarts (SURFLAY, BCSDIR, YYYYMMDD, EXPNAME, EXPDIR, MODEL)

     stop

  else

     if(JOBFILE == 'N') then

        call system('mkdir -p OutData1/  OutData2/')
        tmpstring =  'ln -s '//trim(BCSDIR)//'/'//trim(TILFILE)//' OutData1/OutTileFile'
        call system(tmpstring)
        tmpstring =  'ln -s '//trim(BCSDIR)//'/'//trim(TILFILE)//' OutData2/OutTileFile'
        call system(tmpstring)
        tmpstring =  'ln -s '//trim(BCSDIR)//'/clsm OutData2/clsm'
        call system(tmpstring)

        open (10, file ='mkLDASsa.j', form = 'formatted', status ='unknown', action = 'write')
        write(10,'(a)')'#!/bin/csh -fx'
        write(10,'(a)')' ' 
        write(10,'(a)')'#SBATCH --account='//trim(SPONSORCODE)
        write(10,'(a)')'#SBATCH --time=1:00:00'
        write(10,'(a)')'#SBATCH --ntasks=56'
        write(10,'(a)')'#SBATCH --job-name=mkLDAS'
        write(10,'(a)')'#SBATCH --constraint=hasw'
        write(10,'(a)')'#SBATCH --output=mkLDAS.o'
        write(10,'(a)')'#SBATCH --error=mkLDAS.e'
        write(10,'(a)')' ' 
        write(10,'(a)')'limit stacksize unlimited'
        write(10,'(a)')'source bin/g5_modules'
        !tmpstring = "set BINDIR=`ls -l bin | cut -d'>' -f2`"
        !write(10,'(a)')trim(tmpstring)
        !tmpstring = "setenv ESMADIR `echo $BINDIR | sed 's/Linux\/bin//g'`"
        write(10,'(a)')'setenv ESMADIR '//trim(ESMADIR)
        write(10,'(a)')'setenv MKL_CBWR SSE4_2 # ensure zero-diff across archs'
        write(10,'(a)')'setenv MV2_ON_DEMAND_THRESHOLD 8192 # MVAPICH2'
        write(10,'(a)')' ' 
        write(10,'(a)')'mpirun -np 56 '//trim(cmd)//' -j Y'

        if(trim(model) == 'CATCHCN') then
           write(10,'(a)')'bin/Scale_CatchCN OutData1/catchcn_internal_rst OutData2/catchcn_internal_rst catchcn_internal_rst '//trim(SFL)
           write(10,'(a)')'bin/Scale_CatchCN OutData1/catchcn_internal_clm45 OutData2/catchcn_internal_clm45 catchcn_internal_clm45 '//trim(SFL)
        else
           write(10,'(a)')'bin/Scale_Catch OutData1/catch_internal_rst OutData2/catch_internal_rst catch_internal_rst '//trim(SFL)
        endif

        close (10, status ='keep')
        call system('chmod 755 mkLDASsa.j')
        stop
     endif
  endif
 
  if (master_proc) then
     
     ! read in ntiles 
     ! ----------------------------
     
     open  (10,file = trim(BCSDIR)//'/clsm/catchment.def', form = 'formatted', status ='old', action = 'read')
     read  (10,*) ntiles
     close (10, status ='keep')

  endif

  call MPI_BCAST(NTILES     ,     1, MPI_INTEGER  ,  0,MPI_COMM_WORLD,mpierr)
  
  ! Regridding
  if(trim(MODEL) == 'CATCH'    )inquire(file='OutData1/catch_internal_rst',exist=second_visit )
  if(trim(MODEL) == 'CATCHCN'  )inquire(file='OutData1/catchcn_internal_rst',exist=second_visit )
  if(.not. second_visit) then
     call regrid_hyd_vars (NTILES, trim(MODEL)) 
     call MPI_Barrier(MPI_COMM_WORLD, STATUS)
     stop
  endif
  if (master_proc) then 
     if(trim(MODEL) == 'CATCH'  ) call read_bcs_data (NTILES, SURFLAY, trim(MODEL),'OutData2/clsm/','OutData2/catch_internal_rst'  )
     if(trim(MODEL) == 'CATCHCN') call read_bcs_data (NTILES, SURFLAY, trim(MODEL),'OutData2/clsm/','OutData2/catchcn_internal_rst')
  endif

  call MPI_Barrier(MPI_COMM_WORLD, STATUS)

  if(trim(MODEL) == 'CATCHCN') then 

     call  regrid_carbon_vars (NTILES)     
     if(CLM45) call  regrid_carbon_vars_clm45 (NTILES)

  endif

  call MPI_FINALIZE(mpierr)
     
contains

  ! *****************************************************************************

  SUBROUTINE  reorder_LDASsa_restarts (SURFLAY, BCSDIR, YYYYMMDD, EXPNAME, EXPDIR, MODEL)

    implicit none

    real, intent (in)         :: SURFLAY
    character(*), intent (in) :: BCSDIR, YYYYMMDD, EXPNAME, EXPDIR, MODEL
    character(256)            :: tile_coord, vname
    character(300)            :: rst_file
    type(MAPL_NCIO)           :: InNCIO, OutNCIO
    integer                   :: NTILES, nVars, i,j,k,n, ndims,dimSizes(3)
    integer, allocatable      :: LDAS2BCS (:), g2d(:), tile_id(:)
    real, allocatable         :: var1(:), var2(:),wesn1(:), htsn1(:)
 
   open (10,file =trim(BCSDIR)//"clsm/catchment.def",status='old',form='formatted')
   read (10,*) ntiles
   close (10, status = 'keep')  

   ! read NTILES from BCs and tile_coord from LDASsa experiment

   tile_coord = trim(EXPDIR)//'rc_out/'//trim(expname)//'.ldas_tilecoord.bin'
   open (10,file =trim(tile_coord),status='old',form='unformatted',convert='big_endian')
   read (10) i
   if (i /= ntiles) then 
      print *,'NTILES BCs/LDASsa mismatch:', i,ntiles
      stop
   endif

   if(trim(MODEL) == 'CATCH') &
        InNCIO = MAPL_NCIOOpen('/discover/nobackup/rreichle/l_data/LandRestarts_for_Regridding/Catch/catch_internal_rst' , rc=rc) 
   if(trim(MODEL) == 'CATCHCN') &
        InNCIO = MAPL_NCIOOpen('/discover/nobackup/rreichle/l_data/LandRestarts_for_Regridding/CatchCN/catchcn_internal_dummy' , rc=rc) 

   call MAPL_NCIOGetDimSizes(InNCIO,nVars=nVars)
   call MAPL_NCIOChangeRes(InNCIO,OutNCIO,tileSize=ntiles,rc=rc)

   if(trim(MODEL) == 'CATCH'  )  call MAPL_NCIOSet( OutNCIO,filename='catch_internal_rst.'//trim(YYYYMMDD))
   if(trim(MODEL) == 'CATCHCN')  call MAPL_NCIOSet( OutNCIO,filename='catchcn_internal_rst.'//trim(YYYYMMDD))
   call MAPL_NCIOCreateFile(OutNCIO)   
   call MAPL_NCIOClose ( InNCIO)


   allocate (tile_id  (1:ntiles))
   allocate (LDAS2BCS (1:ntiles))
   allocate (g2d      (1:ntiles))

   read  (10) LDAS2BCS
   close (10, status = 'keep')

   ! ==========================
   ! READ/WRITE LDASsa RESTARTS
   ! ==========================

   allocate(var1(ntiles))
   allocate(var2(ntiles))
   allocate(wesn1 (ntiles))
   allocate(htsn1 (ntiles))  
   ! CH CM CQ FR WW
   ! WW
   var1 = 0.1
   do j = 1,4
      call MAPL_VarWrite(OutNCIO,'WW',var1 ,offset1=j)  
   end do
   ! FR
   var1 = 0.25
   do j = 1,4
      call MAPL_VarWrite(OutNCIO,'FR',var1 ,offset1=j)  
   end do
   ! CH CM CQ 
   var1 = 0.001
   do j = 1,4
      call MAPL_VarWrite(OutNCIO,'CH',var1 ,offset1=j)  
      call MAPL_VarWrite(OutNCIO,'CM',var1 ,offset1=j)  
      call MAPL_VarWrite(OutNCIO,'CQ',var1 ,offset1=j)  
   end do
   
   tile_id = LDAS2BCS
   do n = 1, NTILES 
      G2D(tile_id(n)) = n
   end do
   
   if(trim(MODEL) == 'CATCH') then   
      rst_file = trim(EXPDIR)//'rs/ens0000/Y'//YYYYMMDD(1:4)//'/M'//YYYYMMDD(5:6)//'/'  &
           //trim(ExpName)//'.ens0000.catch_ldas_rst.'// &
           YYYYMMDD(1:8)//'_0000z.bin'

      open(10, file=trim(rst_file), form='unformatted', status='old', &
           convert='big_endian', action='read')
      
      var1 = real(tile_id)
      var2 = var1 (tile_id)
      do n = 1,  NTILES 
         var2(n) = var1(g2d(n))
      end do
      call MAPL_VarWrite(OutNCIO,'TILE_ID' ,var2)

      read(10) var1
      var2 = var1 (tile_id)
      do n = 1,  NTILES 
         var2(n) = var1(g2d(n))
      end do
      call MAPL_VarWrite(OutNCIO,'TC' ,var2, offset1=1)
      read(10) var1
      var2 = var1 (tile_id)
      do n = 1,  NTILES 
         var2(n) = var1(g2d(n))
      end do
      call MAPL_VarWrite(OutNCIO,'TC' ,var2, offset1=2)
      read(10) var1
      var2 = var1 (tile_id)
      do n = 1,  NTILES 
         var2(n) = var1(g2d(n))
      end do
      call MAPL_VarWrite(OutNCIO,'TC' ,var2, offset1=3)
      
      read(10) var1
      var2 = var1 (tile_id)
      do n = 1,  NTILES 
         var2(n) = var1(g2d(n))
      end do
      call MAPL_VarWrite(OutNCIO,'QC' ,var2, offset1=1)
      read(10) var1
      var2 = var1 (tile_id)
      do n = 1,  NTILES 
         var2(n) = var1(g2d(n))
      end do
      call MAPL_VarWrite(OutNCIO,'QC' ,var2, offset1=2)
      read(10) var1
      var2 = var1 (tile_id)
      do n = 1,  NTILES 
         var2(n) = var1(g2d(n))
      end do
      call MAPL_VarWrite(OutNCIO,'QC' ,var2, offset1=3)
      call MAPL_VarWrite(OutNCIO,'QC' ,var2, offset1=4)
      
      read(10) var1
      var2 = var1 (tile_id)
      do n = 1,  NTILES 
         var2(n) = var1(g2d(n))
      end do
      call MAPL_VarWrite(OutNCIO,'CAPAC' ,var2)
      read(10) var1
      var2 = var1 (tile_id)
      do n = 1,  NTILES 
         var2(n) = var1(g2d(n))
      end do
      call MAPL_VarWrite(OutNCIO,'CATDEF' ,var2)
      read(10) var1
      var2 = var1 (tile_id)
      do n = 1,  NTILES 
         var2(n) = var1(g2d(n))
      end do
      call MAPL_VarWrite(OutNCIO,'RZEXC' ,var2)
      read(10) var1
      var2 = var1 (tile_id)
      do n = 1,  NTILES 
         var2(n) = var1(g2d(n))
      end do
      call MAPL_VarWrite(OutNCIO,'SRFEXC' ,var2)
      read(10) var1
      var2 = var1 (tile_id)
      do n = 1,  NTILES 
         var2(n) = var1(g2d(n))
      end do
      call MAPL_VarWrite(OutNCIO,'GHTCNT1' ,var2)     
      read(10) var1
      var2 = var1 (tile_id) 
      do n = 1,  NTILES 
         var2(n) = var1(g2d(n))
      end do
      call MAPL_VarWrite(OutNCIO,'GHTCNT2' ,var2) 
      read(10) var1
      var2 = var1 (tile_id) 
      do n = 1,  NTILES 
         var2(n) = var1(g2d(n))
      end do
      call MAPL_VarWrite(OutNCIO,'GHTCNT3' ,var2) 
      read(10) var1
      var2 = var1 (tile_id)   
      do n = 1,  NTILES 
         var2(n) = var1(g2d(n))
      end do
      call MAPL_VarWrite(OutNCIO,'GHTCNT4' ,var2) 
      read(10) var1
      var2 = var1 (tile_id)  
      do n = 1,  NTILES 
         var2(n) = var1(g2d(n))
      end do
      call MAPL_VarWrite(OutNCIO,'GHTCNT5' ,var2)   
      read(10) var1
      var2 = var1 (tile_id)  
      do n = 1,  NTILES 
         var2(n) = var1(g2d(n))
      end do
      call MAPL_VarWrite(OutNCIO,'GHTCNT6' ,var2)   
      read(10) var1
      var2 = var1 (tile_id)   
      
      do n = 1,  NTILES 
         var2(n) = var1(g2d(n))
      end do
      wesn1 = var2
      call MAPL_VarWrite(OutNCIO,'WESNN1' ,var2)
      read(10) var1
      var2 = var1 (tile_id) 
      do n = 1,  NTILES 
         var2(n) = var1(g2d(n))
      end do
      call MAPL_VarWrite(OutNCIO,'WESNN2' ,var2)
      read(10) var1
      var2 = var1 (tile_id) 
      do n = 1,  NTILES 
         var2(n) = var1(g2d(n))
      end do
      call MAPL_VarWrite(OutNCIO,'WESNN3' ,var2)
      read(10) var1
      var2 = var1 (tile_id)  
      do n = 1,  NTILES 
         var2(n) = var1(g2d(n))
      end do
      htsn1 = var2
      call MAPL_VarWrite(OutNCIO,'HTSNNN1' ,var2)
      read(10) var1
      var2 = var1 (tile_id)   
      do n = 1,  NTILES 
         var2(n) = var1(g2d(n))
      end do
      call MAPL_VarWrite(OutNCIO,'HTSNNN2' ,var2)
      read(10) var1
      var2 = var1 (tile_id)   
      do n = 1,  NTILES 
         var2(n) = var1(g2d(n))
      end do
      call MAPL_VarWrite(OutNCIO,'HTSNNN3' ,var2)
      read(10) var1
      var2 = var1 (tile_id)   
      do n = 1,  NTILES 
         var2(n) = var1(g2d(n))
      end do
      call MAPL_VarWrite(OutNCIO,'SNDZN1' ,var2)
      read(10) var1
      var2 = var1 (tile_id)   
      do n = 1,  NTILES 
         var2(n) = var1(g2d(n))
      end do
      call MAPL_VarWrite(OutNCIO,'SNDZN2' ,var2)
      read(10) var1
      var2 = var1 (tile_id)   
      do n = 1,  NTILES 
         var2(n) = var1(g2d(n))
      end do
      call MAPL_VarWrite(OutNCIO,'SNDZN3' ,var2)
      call STIEGLITZSNOW_CALC_TPSNOW(NTILES, HTSN1(:), WESN1(:), var2, var1)
      var2 = var2 + 273.16
      call MAPL_VarWrite(OutNCIO,'TC' ,var2, offset1=4)
      deallocate (var1, var2)
      call MAPL_NCIOClose (OutNCIO)

      call read_bcs_data (ntiles, SURFLAY, trim(MODEL), trim(BCSDIR)//'/clsm/','catch_internal_rst.'//trim(YYYYMMDD))

   else
      
       rst_file = trim(EXPDIR)//'rs/ens0000/Y'//YYYYMMDD(1:4)//'/M'//YYYYMMDD(5:6)//'/'//trim(ExpName)//&
            '.ens0000.catchcn_ldas_rst.'//trim(YYYYMMDD)//'_0000z'

       InNCIO = MAPL_NCIOOpen(trim(rst_file), rc=rc) 
       call MAPL_NCIOGetDimSizes(InNCIO,nVars=nVars)
       call MAPL_VarRead ( InNCIO,'TILE_ID',var1)
       if(sum (nint(var1) - LDAS2BCS) /= 0) then
          print *, 'Tile order mismatch ', sum(var1)/ntiles, sum(LDAS2BCS)/ntiles 
          stop
       endif

       do k=1,nVars
          
          call MAPL_NCIOGetVarName(InNCIO,k,vname)
          
          call MAPL_NCIOVarGetDims(InNCIO,vname,nDims,dimSizes)
          if (ndims == 1) then
             call MAPL_VarRead ( InNCIO,vname,var1)
             var2 = var1 (tile_id)   
             do n = 1,  NTILES 
                var2(n) = var1(g2d(n))
             end do
             if(trim(vname) == 'SFMCM'  ) var2 = 0.
             if(trim(vname) == 'BFLOWM' ) var2 = 0.
             if(trim(vname) == 'TOTWATM') var2 = 0.
             if(trim(vname) == 'TAIRM'  ) var2 = 0.
             if(trim(vname) == 'TPM'    ) var2 = 0.
             if(trim(vname) == 'CNSUM'  ) var2 = 0.
             if(trim(vname) == 'SNDZM'  ) var2 = 0.
             if(trim(vname) == 'ASNOWM' ) var2 = 0.
             if(trim(vname) == 'TSURF'  ) var2 = 0.
             
             call MAPL_VarWrite(OutNCIO,vname,var2)
          
       else if (ndims == 2) then
          
          do j=1,dimSizes(2)
             call MAPL_VarRead ( InNCIO,vname,var1 ,offset1=j)
             var2 = var1 (tile_id)
             do n = 1,  NTILES 
                var2(n) = var1(g2d(n))
             end do
             if(trim(vname) == 'TGWM'  ) var2 = 0.
             if(trim(vname) == 'RZMM'  ) var2 = 0.
             if(trim(vname) == 'WW'    ) var2 = 0.1
             if(trim(vname) == 'FR'    ) var2 = 0.25
             if(trim(vname) == 'CQ'    ) var2 = 0.001
             if(trim(vname) == 'CN'    ) var2 = 0.001 
             if(trim(vname) == 'CM'    ) var2 = 0.001 
             if(trim(vname) == 'CH'    ) var2 = 0.001 


             call MAPL_VarWrite(OutNCIO,vname,var2 ,offset1=j)
          enddo
          
       else if (ndims == 3) then
          
          do i=1,dimSizes(3)
             do j=1,dimSizes(2)
                call MAPL_VarRead ( InNCIO,vname,var1 ,offset1=j,offset2=i)
                var2 = var1 (tile_id)
                do n = 1,  NTILES 
                   var2(n) = var1(g2d(n))
                end do
                if(trim(vname) == 'PSNSUNM'  ) var2 = 0.
                if(trim(vname) == 'PSNSHAM'  ) var2 = 0.
                call MAPL_VarWrite(OutNCIO,vname,var2 ,offset1=j,offset2=i)
             enddo
          enddo
          
       end if
    enddo
 
    deallocate (var1, var2, tile_id)
  
    call MAPL_NCIOClose      (InNCIO)
    call MAPL_NCIOClose      (OutNCIO)
    
    call read_bcs_data (ntiles, SURFLAY, trim(MODEL), trim(BCSDIR)//'/clsm/','catchcn_internal_rst.'//trim(YYYYMMDD))    
    
   endif

  END SUBROUTINE reorder_LDASsa_restarts
  
  ! *****************************************************************************
  
  SUBROUTINE regrid_hyd_vars (NTILES, model) 

    implicit none
    integer, intent (in)           :: NTILES
    character(*), intent (in)      :: model

    ! ===============================================================================================

    integer, allocatable, dimension(:)   :: Id_glb, Id_loc
    integer, allocatable, dimension(:)   :: ld_reorder, tid_offl
    real   , allocatable, dimension(:)   :: tmp_var
    logical, allocatable, dimension(:)   :: mask
    real    :: dw, min_lon, max_lon, min_lat, max_lat
    integer :: n,i,nplus, STATUS,NCFID
    integer :: local_id, ntiles_smap
    integer, allocatable, dimension (:) :: sub_tid
    real   , allocatable, dimension (:) :: sub_lon, sub_lat, rev_dist, lonc, latc, LATT, LONN, long, latg
    integer, allocatable :: low_ind(:), upp_ind(:), nt_local (:)

    logical :: all_found

    if(trim(MODEL) == 'CATCHCN') ntiles_smap = ntiles_cn
    if(trim(MODEL) == 'CATCH'  ) ntiles_smap = ntiles_cat

    allocate (tid_offl  (ntiles_smap))
    allocate (tmp_var   (ntiles_smap))
    allocate (mask      (ntiles_smap))

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
    allocate (lonc   (1:ntiles_smap))
    allocate (latc   (1:ntiles_smap))

    if (master_proc) then

       allocate (long   (ntiles))
       allocate (latg   (ntiles))
       allocate (ld_reorder(ntiles_smap)) 

       call ReadCNTilFile ('OutData1/OutTileFile', ntiles, long, latg)

       ! ---------------------------------------------
       ! Read exact lonc, latc from offline .til File 
       ! ---------------------------------------------

       if(trim(MODEL) == 'CATCHCN') call ReadCNTilFile(trim(InCNTilFile ),ntiles_smap,lonc,latc)
       if(trim(MODEL) == 'CATCH'  ) call ReadCNTilFile(trim(InCatTilFile),ntiles_smap,lonc,latc)

       if(trim(MODEL) == 'CATCHCN') STATUS = NF_OPEN (trim(InCNRestart ),NF_NOWRITE,NCFID)
       if(trim(MODEL) == 'CATCH'  ) STATUS = NF_OPEN (trim(InCatRestart),NF_NOWRITE,NCFID)
       STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'TILE_ID'   ), (/1/), (/NTILES_SMAP/),tmp_var)
       STATUS = NF_CLOSE (NCFID)

       do n = 1, ntiles_smap
          ld_reorder ( NINT(tmp_var(n))) = n
          tid_offl(n)    = n
       end do

       deallocate (tmp_var)

    endif

    call MPI_Barrier(MPI_COMM_WORLD, STATUS)

    call MPI_SCATTERV (                    &
         long,nt_local,low_ind-1,MPI_real, &
         lonn,size(lonn),MPI_real  , &
         0,MPI_COMM_WORLD, mpierr )

    call MPI_SCATTERV (                    &
         latg,nt_local,low_ind-1,MPI_real, &
         latt,nt_local(myid+1),MPI_real  , &
         0,MPI_COMM_WORLD, mpierr )

    if(master_proc) deallocate (long, latg)
     
    call MPI_BCAST(lonc,ntiles_smap,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_BCAST(latc,ntiles_smap,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_BCAST(tid_offl,size(tid_offl  ),MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)

    ! --------------------------------------------------------------------------------
    ! Here we create transfer index array to map offline restarts to output tile space
    ! --------------------------------------------------------------------------------   
      
    ! Loop through NTILES (# of tiles in output array) find the nearest neighbor from Qing.  

    Id_loc = -9999
    
    TILES : do n = low_ind (myid + 1), upp_ind (myid + 1)

!       if(MOD(n,10000) == 0) print *,'In HYD', myid,n
       local_id = n - low_ind (myid + 1) + 1
                      
       dw = 0.5 !1.0 ! Start with a 10x10 window, then zoom out by increasing the size by 2-deg until 4 similar tiles are found for 4 PFT types        

       ZOOMOUT : do  

          all_found = .false.   

          ! Min/Max lon/lat of the working window
          ! -------------------------------------
          
          min_lon = MAX(lonn (local_id) - dw, -180.)
          max_lon = MIN(lonn (local_id) + dw,  180.)
          min_lat = MAX(latt (local_id) - dw,  -90.)
          max_lat = MIN(latt (local_id) + dw,   90.) 
          mask = .false.
          mask =  ((latc >= min_lat .and. latc <= max_lat).and.(lonc >= min_lon .and. lonc <= max_lon))
          nplus =  count(mask = mask)
          
          if(nplus < 0) then
             dw = dw + 0.5
             CYCLE
          endif
          
          allocate (sub_tid (1:nplus))
          allocate (sub_lon (1:nplus))
          allocate (sub_lat (1:nplus))
          allocate (rev_dist  (1:nplus))
          
          sub_tid = PACK (tid_offl, mask= mask) 
          sub_lon = PACK (lonc    , mask= mask)
          sub_lat = PACK (latc    , mask= mask)
          
          ! compute distance from the tile
          
          sub_lat = sub_lat * MAPL_PI/180.
          sub_lon = sub_lon * MAPL_PI/180.
          
!          do i = 1,nplus
!             sub_dist(i) = haversine(to_radian(latt(local_id)), to_radian(lonn(local_id)), &
!                  sub_lat(i), sub_lon(i)) 
!             !                   sub_dist(i) = dist
!          end do
             
          SEEK : if(Id_loc(local_id) < 0) then
                             
                rev_dist  = 1.e20
                
                do i = 1,nplus
                      
                   rev_dist(i) = haversine(to_radian(latt(local_id)), to_radian(lonn(local_id)), &
                        sub_lat(i), sub_lon(i))
                                      
                end do

                FOUND : if(minval (rev_dist) < 1.e19) then 
                   
                   Id_loc(local_id) = sub_tid(minloc(rev_dist,1)) ! Cell ID of the nearest neightbor in the offline 
                                                                     ! restart file that has the same veg type
                    all_found = .true.                  
                endif FOUND
                
             endif SEEK
          
          deallocate (sub_tid, sub_lon, sub_lat, rev_dist)
          
          if(all_found) GO TO 100
           
           ! if not increase the window size
           dw = dw + 0.5
           
        end do ZOOMOUT
        
100   continue !  if(mod (n,1000) == 0) print *, myid +1, n, Id_loc(local_id,:)

     END DO TILES
 
     ! update id_glb in root

     if(master_proc)  allocate (id_glb  (ntiles))

     call MPI_Barrier(MPI_COMM_WORLD, STATUS)
     call MPI_GATHERV( &
                   id_loc, nt_local(myid+1)  , MPI_real, &
                   id_glb, nt_local,low_ind-1, MPI_real, &
                   0, MPI_COMM_WORLD, mpierr )
        
    if (master_proc) call put_land_vars  (NTILES, ntiles_smap, id_glb, ld_reorder, model)

    call MPI_Barrier(MPI_COMM_WORLD, STATUS)

   END SUBROUTINE regrid_hyd_vars


  ! *****************************************************************************
  
  SUBROUTINE read_bcs_data (ntiles, SURFLAY,MODEL, DataDir, InRestart)

    ! This subroutine :
    !  1) reads BCs from BCSDIR and hydrological varables from InRestart.
    !     InRestart is a catchcn_internal_rst nc4 file.
    !
    !  2) writes out BCs and hydrological variables in catchcn_internal_rst (1:72). 
    !     output catchcn_internal_rst is nc4.

    implicit none
    real, intent (in)                         :: SURFLAY
    integer, intent (in)                      :: ntiles
    character(*), intent (in)                 :: MODEL, DataDir, InRestart
    real, allocatable :: CLMC_pf1(:), CLMC_pf2(:), CLMC_sf1(:), CLMC_sf2(:)
    real, allocatable :: CLMC_pt1(:), CLMC_pt2(:), CLMC_st1(:), CLMC_st2(:)    
    real, allocatable :: CLMC45_pf1(:), CLMC45_pf2(:), CLMC45_sf1(:), CLMC45_sf2(:)
    real, allocatable :: CLMC45_pt1(:), CLMC45_pt2(:), CLMC45_st1(:), CLMC45_st2(:)    
    real, allocatable :: BF1(:),   BF2(:),   BF3(:),  VGWMAX(:)
    real, allocatable :: CDCR1(:), CDCR2(:), PSIS(:), BEE(:) 
    real, allocatable :: POROS(:), WPWET(:), COND(:), GNU(:)
    real, allocatable :: ARS1(:),  ARS2(:),  ARS3(:)
    real, allocatable :: ARA1(:),  ARA2(:),  ARA3(:), ARA4(:)
    real, allocatable :: ARW1(:),  ARW2(:),  ARW3(:), ARW4(:)
    real, allocatable :: TSA1(:),  TSA2(:),  TSB1(:), TSB2(:)
    real, allocatable :: ATAU2(:), BTAU2(:), DP2BR(:), CanopH(:)
    real, allocatable :: NDEP(:), BVISDR(:), BVISDF(:), BNIRDR(:), BNIRDF(:) 
    real, allocatable :: T2(:), var1(:), hdm(:), fc(:), gdp(:), peatf(:)
    integer, allocatable :: ity(:), abm (:)
    integer       :: NCFID, STATUS
    integer       :: idum, i,j,n, ib, nv
    real          :: rdum, zdep1, zdep2, zdep3, zmet, term1, term2, bare,fvg(4)
    logical       :: NEWLAND
   
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
    allocate (   ity(ntiles),  CanopH(ntiles))
    allocate (CLMC_pf1(ntiles), CLMC_pf2(ntiles), CLMC_sf1(ntiles))
    allocate (CLMC_sf2(ntiles), CLMC_pt1(ntiles), CLMC_pt2(ntiles))
    allocate (CLMC45_pf1(ntiles), CLMC45_pf2(ntiles), CLMC45_sf1(ntiles))
    allocate (CLMC45_sf2(ntiles), CLMC45_pt1(ntiles), CLMC45_pt2(ntiles))
    allocate (CLMC_st1(ntiles), CLMC_st2(ntiles))
    allocate (CLMC45_st1(ntiles), CLMC45_st2(ntiles))
    allocate (hdm(ntiles), fc(ntiles), gdp(ntiles))
    allocate (peatf(ntiles), abm(ntiles), var1(ntiles))

    inquire(file=trim(DataDir)//"CLM_veg_typs_fracs"   ,exist=NewLand )

    open(unit=21, file=trim(DataDir)//'mosaic_veg_typs_fracs',form='formatted') 
    open(unit=22, file=trim(DataDir)//'bf.dat'               ,form='formatted')
    open(unit=23, file=trim(DataDir)//'soil_param.dat'       ,form='formatted')
    open(unit=24, file=trim(DataDir)//'ar.new'               ,form='formatted')
    open(unit=25, file=trim(DataDir)//'ts.dat'               ,form='formatted')
    open(unit=26, file=trim(DataDir)//'tau_param.dat'        ,form='formatted')

    if(NewLand) then
       open(unit=27, file=trim(DataDir)//'CLM_veg_typs_fracs'   ,form='formatted')
       open(unit=28, file=trim(DataDir)//'CLM_NDep_SoilAlb_T2m' ,form='formatted')
       if(clm45) then
          open(unit=29, file=trim(DataDir)//'CLM4.5_veg_typs_fracs',form='formatted')
          open(unit=30, file=trim(DataDir)//'CLM4.5_abm_peatf_gdp_hdm_fc' ,form='formatted')
       endif
    endif

    do n=1,ntiles
       var1 (n) = real (n)

       if (NewLand) then
          read(21,*) I, j, ITY(N),idum, rdum, rdum, CanopH(N)
       else
          read(21,*) I, j, ITY(N),idum, rdum, rdum
       endif

       read (22, *) i,j, GNU(n), BF1(n), BF2(n), BF3(n)
       
       read (23, *) i,j, idum, idum, BEE(n), PSIS(n),&
            POROS(n), COND(n), WPWET(n), DP2BR(n)
       
       read (24, *) i,j, rdum, ARS1(n), ARS2(n), ARS3(n),          &
            ARA1(n), ARA2(n), ARA3(n), ARA4(n), &
            ARW1(n), ARW2(n), ARW3(n), ARW4(n)
       
       read (25, *) i,j, rdum, TSA1(n), TSA2(n), TSB1(n), TSB2(n)
       
       if( SURFLAY.eq.20.0 ) read (26, *) i,j, ATAU2(n), BTAU2(n), rdum, rdum   ! for old soil params
       if( SURFLAY.eq.50.0 ) read (26, *) i,j, rdum , rdum, ATAU2(n), BTAU2(n)  ! for new soil params

       if (NewLand) then
          read (27, *) i,j, CLMC_pt1(n), CLMC_pt2(n), CLMC_st1(n), CLMC_st2(n), &
               CLMC_pf1(n), CLMC_pf2(n), CLMC_sf1(n), CLMC_sf2(n)
          
          read (28, *) NDEP(n), BVISDR(n), BVISDF(n), BNIRDR(n), BNIRDF(n), T2(n) ! MERRA-2 Annual Mean Temp is default.
          if(clm45) then
             read (29, *) i,j, CLMC45_pt1(n), CLMC45_pt2(n), CLMC45_st1(n), CLMC45_st2(n), &
                  CLMC45_pf1(n), CLMC45_pf2(n), CLMC45_sf1(n), CLMC45_sf2(n)
             
             read (30,'(2I8, i3, f8.4, f8.2, f10.2, f8.4)' ) i, j, abm(n), peatf(n), &
                  gdp(n), hdm(n), fc(n)
          endif
       endif

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

       if(trim(model) == 'CATCHCN') then

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
          
          if(CLM45) then
             ! CLM 45
             
             CLMC45_pf1(n) = CLMC45_pf1(n) / 100.
             CLMC45_pf2(n) = CLMC45_pf2(n) / 100.
             CLMC45_sf1(n) = CLMC45_sf1(n) / 100.
             CLMC45_sf2(n) = CLMC45_sf2(n) / 100.
             
             fvg(1) = CLMC45_pf1(n)
             fvg(2) = CLMC45_pf2(n)
             fvg(3) = CLMC45_sf1(n)
             fvg(4) = CLMC45_sf2(n)
             
             BARE = 1.      
             
             DO NV = 1, NVEG
                BARE = BARE - FVG(NV)! subtract vegetated fractions 
             END DO
             
             if (BARE /= 0.) THEN
                IB = MAXLOC(FVG(:),1)
                FVG (IB) = FVG(IB) + BARE ! This also corrects all cases sum ne 0.
             ENDIF
             
             CLMC45_pf1(n) = fvg(1)
             CLMC45_pf2(n) = fvg(2)
             CLMC45_sf1(n) = fvg(3)
             CLMC45_sf2(n) = fvg(4)
          endif
       endif
    enddo

    if(trim(model) == 'CATCHCN') then
       
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
         
          if (clm45) then
             ! CLM45
             if(CLMC45_pf1(n) <= 1.e-4) then
                CLMC45_pf2(n) = CLMC45_pf2(n) + CLMC45_pf1(n)
                CLMC45_pf1(n) = 0.
             endif
             
             if(CLMC45_pf2(n) <= 1.e-4) then
                CLMC45_pf1(n) = CLMC45_pf1(n) + CLMC45_pf2(n)
                CLMC45_pf2(n) = 0.
             endif
             
             if(CLMC45_sf1(n) <= 1.e-4) then
                if(CLMC45_sf2(n) > 1.e-4) then
                   CLMC45_sf2(n) = CLMC45_sf2(n) + CLMC45_sf1(n)
                else if(CLMC45_pf2(n) > 1.e-4) then
                   CLMC45_pf2(n) = CLMC45_pf2(n) + CLMC45_sf1(n)
                else if(CLMC45_pf1(n) > 1.e-4) then
                   CLMC45_pf1(n) = CLMC45_pf1(n) + CLMC45_sf1(n)
                else
                   stop 'fveg3'
                endif
                CLMC45_sf1(n) = 0.
             endif
             
             if(CLMC45_sf2(n) <= 1.e-4) then
                if(CLMC45_sf1(n) > 1.e-4) then
                   CLMC45_sf1(n) = CLMC45_sf1(n) + CLMC45_sf2(n)
                else if(CLMC45_pf2(n) > 1.e-4) then
                   CLMC45_pf2(n) = CLMC45_pf2(n) + CLMC45_sf2(n)
                else if(CLMC45_pf1(n) > 1.e-4) then
                   CLMC45_pf1(n) = CLMC45_pf1(n) + CLMC45_sf2(n)
                else
                   stop 'fveg4'
                endif
                CLMC45_sf2(n) = 0.
             endif
          endif
       end do
    endif
    CLOSE (22, STATUS = 'KEEP')
    CLOSE (23, STATUS = 'KEEP')
    CLOSE (24, STATUS = 'KEEP')
    CLOSE (25, STATUS = 'KEEP')
    CLOSE (26, STATUS = 'KEEP')
    
    if(NewLand) then
        CLOSE (27, STATUS = 'KEEP')
        CLOSE (28, STATUS = 'KEEP')
        if(clm45) then
           CLOSE (29, STATUS = 'KEEP')
           CLOSE (30, STATUS = 'KEEP')
        endif
     endif
     
     ! Vegdyn Boundary Condition
     ! -------------------------
     
     ! open(20,file=trim("vegdyn_internal_rst"), &
     !     status="unknown", &
     !     form="unformatted",convert="little_endian")
     ! write(20) real(ity)
     ! if(NewLand) write(20) CanopH
     ! close(20)
     ! print *, "Wrote vegdyn_internal_restart"

     ! Now writing BCs (from BCSDIR) and regridded hydrological variables 1-72
     ! -----------------------------------------------------------------------

     STATUS = NF_OPEN (trim(InRestart),NF_WRITE,NCFID)
     STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'BF1'), (/1/), (/NTILES/),BF1)
     STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'BF2'), (/1/), (/NTILES/),BF2)
     STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'BF3'), (/1/), (/NTILES/),BF3)
     STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'VGWMAX'), (/1/), (/NTILES/),VGWMAX)
     STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'CDCR1'), (/1/), (/NTILES/),CDCR1)
     STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'CDCR2'), (/1/), (/NTILES/),CDCR2)
     STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'PSIS'), (/1/), (/NTILES/),PSIS)
     STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'BEE'), (/1/), (/NTILES/),BEE)
     STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'POROS'), (/1/), (/NTILES/),POROS)
     STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'WPWET'), (/1/), (/NTILES/),WPWET)
     STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'COND'), (/1/), (/NTILES/),COND)
     STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'GNU'), (/1/), (/NTILES/),GNU)
     STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'ARS1'), (/1/), (/NTILES/),ARS1)
     STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'ARS2'), (/1/), (/NTILES/),ARS2)
     STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'ARS3'), (/1/), (/NTILES/),ARS3)
     STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'ARA1'), (/1/), (/NTILES/),ARA1)
     STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'ARA2'), (/1/), (/NTILES/),ARA2)
     STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'ARA3'), (/1/), (/NTILES/),ARA3)
     STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'ARA4'), (/1/), (/NTILES/),ARA4)
     STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'ARW1'), (/1/), (/NTILES/),ARW1)
     STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'ARW2'), (/1/), (/NTILES/),ARW2)
     STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'ARW3'), (/1/), (/NTILES/),ARW3)
     STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'ARW4'), (/1/), (/NTILES/),ARW4)
     STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'TSA1'), (/1/), (/NTILES/),TSA1)
     STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'TSA2'), (/1/), (/NTILES/),TSA2)
     STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'TSB1'), (/1/), (/NTILES/),TSB1)
     STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'TSB2'), (/1/), (/NTILES/),TSB2)
     STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'ATAU'), (/1/), (/NTILES/),ATAU2)
     STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'BTAU'), (/1/), (/NTILES/),BTAU2)
     STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'TILE_ID'), (/1/), (/NTILES/),VAR1)

     if(trim(model) == 'CATCHCN') then

        STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'ITY'), (/1,1/), (/NTILES,1/),CLMC_pt1)
        STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'ITY'), (/1,2/), (/NTILES,1/),CLMC_pt2)
        STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'ITY'), (/1,3/), (/NTILES,1/),CLMC_st1)
        STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'ITY'), (/1,4/), (/NTILES,1/),CLMC_st2)
        
        STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'FVG'), (/1,1/), (/NTILES,1/),CLMC_pf1)
        STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'FVG'), (/1,2/), (/NTILES,1/),CLMC_pf2)
        STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'FVG'), (/1,3/), (/NTILES,1/),CLMC_sf1)
        STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'FVG'), (/1,4/), (/NTILES,1/),CLMC_sf2)
        

        STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'NDEP'), (/1/), (/NTILES/),NDEP)
        STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'CLI_T2M'), (/1/), (/NTILES/),T2)
        STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'BGALBVR'), (/1/), (/NTILES/),BVISDR)
        STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'BGALBVF'), (/1/), (/NTILES/),BVISDF)
        STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'BGALBNR'), (/1/), (/NTILES/),BNIRDR)
        STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'BGALBNF'), (/1/), (/NTILES/),BNIRDF)
        
        STATUS = NF_CLOSE ( NCFID)

        if(CLM45) then
           ! CLM45
           IF(index(InRestart,'OutData2') /= 0) then
               STATUS = NF_OPEN ('OutData2/catchcn_internal_clm45',NF_WRITE,NCFID)
           ELSE
              STATUS = NF_OPEN ('catchcn_internal_clm45',NF_WRITE,NCFID)
           ENDIF

           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'BF1'), (/1/), (/NTILES/),BF1)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'BF2'), (/1/), (/NTILES/),BF2)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'BF3'), (/1/), (/NTILES/),BF3)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'VGWMAX'), (/1/), (/NTILES/),VGWMAX)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'CDCR1'), (/1/), (/NTILES/),CDCR1)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'CDCR2'), (/1/), (/NTILES/),CDCR2)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'PSIS'), (/1/), (/NTILES/),PSIS)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'BEE'), (/1/), (/NTILES/),BEE)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'POROS'), (/1/), (/NTILES/),POROS)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'WPWET'), (/1/), (/NTILES/),WPWET)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'COND'), (/1/), (/NTILES/),COND)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'GNU'), (/1/), (/NTILES/),GNU)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'ARS1'), (/1/), (/NTILES/),ARS1)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'ARS2'), (/1/), (/NTILES/),ARS2)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'ARS3'), (/1/), (/NTILES/),ARS3)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'ARA1'), (/1/), (/NTILES/),ARA1)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'ARA2'), (/1/), (/NTILES/),ARA2)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'ARA3'), (/1/), (/NTILES/),ARA3)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'ARA4'), (/1/), (/NTILES/),ARA4)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'ARW1'), (/1/), (/NTILES/),ARW1)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'ARW2'), (/1/), (/NTILES/),ARW2)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'ARW3'), (/1/), (/NTILES/),ARW3)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'ARW4'), (/1/), (/NTILES/),ARW4)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'TSA1'), (/1/), (/NTILES/),TSA1)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'TSA2'), (/1/), (/NTILES/),TSA2)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'TSB1'), (/1/), (/NTILES/),TSB1)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'TSB2'), (/1/), (/NTILES/),TSB2)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'ATAU'), (/1/), (/NTILES/),ATAU2)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'BTAU'), (/1/), (/NTILES/),BTAU2)
           
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'ITY'), (/1,1/), (/NTILES,1/),CLMC45_pt1)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'ITY'), (/1,2/), (/NTILES,1/),CLMC45_pt2)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'ITY'), (/1,3/), (/NTILES,1/),CLMC45_st1)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'ITY'), (/1,4/), (/NTILES,1/),CLMC45_st2)
           
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'FVG'), (/1,1/), (/NTILES,1/),CLMC45_pf1)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'FVG'), (/1,2/), (/NTILES,1/),CLMC45_pf2)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'FVG'), (/1,3/), (/NTILES,1/),CLMC45_sf1)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'FVG'), (/1,4/), (/NTILES,1/),CLMC45_sf2)
           
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'TILE_ID'), (/1/), (/NTILES/),VAR1)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'NDEP'), (/1/), (/NTILES/),NDEP)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'CLI_T2M'), (/1/), (/NTILES/),T2)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'BGALBVR'), (/1/), (/NTILES/),BVISDR)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'BGALBVF'), (/1/), (/NTILES/),BVISDF)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'BGALBNR'), (/1/), (/NTILES/),BNIRDR)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'BGALBNF'), (/1/), (/NTILES/),BNIRDF)

           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'ABM'), (/1/), (/NTILES/),real(ABM))
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'FIELDCAP'), (/1/), (/NTILES/),FC)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'HDM'), (/1/), (/NTILES/),HDM)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'GDP'), (/1/), (/NTILES/),GDP)
           STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'PEATF'), (/1/), (/NTILES/),PEATF)
           
           STATUS = NF_CLOSE ( NCFID)
           
        endif
     else
        STATUS = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'OLD_ITY'), (/1/), (/NTILES/),real(ITY))
        STATUS = NF_CLOSE ( NCFID)
     endif

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
     deallocate (   ity,  CanopH)
     deallocate (CLMC_pf1, CLMC_pf2, CLMC_sf1)
     deallocate (CLMC_sf2, CLMC_pt1, CLMC_pt2)
     deallocate (CLMC_st1,CLMC_st2)

  END SUBROUTINE read_bcs_data

  ! *****************************************************************************
  
  SUBROUTINE regrid_carbon_vars (NTILES)

    implicit none

    integer, intent (in)                 :: NTILES
    character*300          :: OutTileFile = 'OutData1/OutTileFile', OutFileName='OutData2/catchcn_internal_rst'
    integer                :: AGCM_YY=2015,AGCM_MM=1,AGCM_DD=1,AGCM_HR=0 
    real, allocatable, dimension (:)     :: CLMC_pf1, CLMC_pf2, CLMC_sf1, CLMC_sf2, &
         CLMC_pt1, CLMC_pt2,CLMC_st1,CLMC_st2

    ! ===============================================================================================

    integer :: iclass(npft) = (/1,1,2,3,3,4,5,5,6,7,8,9,10,11,12,11,12,11,12/)
    integer, allocatable, dimension(:,:) :: Id_glb, Id_loc
    integer, allocatable, dimension(:)   :: tid_offl, id_vec
    logical, allocatable, dimension(:)   :: mask
    real,    allocatable, dimension(:,:) :: fveg_offl,  ityp_offl
    real    :: dw, min_lon, max_lon, min_lat, max_lat, fveg_new, sub_dist
    integer :: n,i,j, k, nplus, nv, nx, nz, iv, offl_cell, ityp_new, STATUS,NCFID
    integer :: outid, local_id
    integer, allocatable, dimension (:) :: sub_tid, sub_ityp1, sub_ityp2,icl_ityp1
    real   , allocatable, dimension (:) :: sub_lon, sub_lat, rev_dist, sub_fevg1, sub_fevg2,&
         lonc, latc, LATT, LONN, DAYX, long, latg, var_dum, TILE_ID, var_dum2
    real, allocatable :: var_off_col (:,:,:), var_off_pft (:,:,:,:) 
    real, allocatable :: var_col_out (:,:,:), var_pft_out (:,:,:,:) 
    integer, allocatable :: low_ind(:), upp_ind(:), nt_local (:)

    logical :: all_found
  
    allocate (tid_offl  (ntiles_cn))
    allocate (mask      (ntiles_cn))
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

    if (master_proc) then
       
       ! --------------------------------------------
       ! Read exact lonn, latt from output .til file 
       ! --------------------------------------------

       allocate (long   (ntiles))
       allocate (latg   (ntiles))
       allocate (DAYX   (NTILES))

       call ReadCNTilFile (OutTileFile, ntiles, long, latg)

       ! Compute DAYX
       ! ------------

       call compute_dayx (                                     &
            NTILES, AGCM_YY, AGCM_MM, AGCM_DD, AGCM_HR,        &
            LATG, DAYX)   

       ! ---------------------------------------------
       ! Read exact lonc, latc from offline .til File 
       ! ---------------------------------------------

       call ReadCNTilFile(trim(InCNTilFile),ntiles_cn,lonc,latc)

    endif

    call MPI_SCATTERV (                    &
         long,nt_local,low_ind-1,MPI_real, &
         lonn,size(lonn),MPI_real  , &
         0,MPI_COMM_WORLD, mpierr )

    call MPI_SCATTERV (                    &
         latg,nt_local,low_ind-1,MPI_real, &
         latt,nt_local(myid+1),MPI_real  , &
         0,MPI_COMM_WORLD, mpierr )

    if(master_proc) deallocate (long, latg)
 
    call MPI_BCAST(lonc,ntiles_cn,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_BCAST(latc,ntiles_cn,MPI_REAL,0,MPI_COMM_WORLD,mpierr)

    ! Open GKW/Fzeng SMAP M09 catchcn_internal_rst and output catchcn_internal_rst
    ! ----------------------------------------------------------------------------
    
    ! NF_OPEN_PAR is no longer needed since IO is done by the root processor.
    !    call MPI_Info_create(info, STATUS)
    !    call MPI_Info_set(info, "romio_cb_read", "automatic", STATUS)   
    !    STATUS = NF_OPEN_PAR   (trim(InCNRestart),IOR(NF_NOWRITE,NF_MPIIO),MPI_COMM_WORLD, info,NCFID)
    !    STATUS = NF_OPEN_PAR   (trim(OutFileName),IOR(NF_WRITE  ,NF_MPIIO),MPI_COMM_WORLD, info,OUTID)
    
    STATUS = NF_OPEN (trim(InCNRestart),NF_NOWRITE,NCFID)
    IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS, 'INPUT RESTART FAILED')
    STATUS = NF_OPEN (trim(OutFileName),NF_WRITE,OUTID)
    IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS, 'OUTPUT RESTART FAILED')
    
    STATUS = NF_GET_VARA_REAL(OUTID,VarID(OUTID,'ITY'), (/low_ind(myid+1),1/), (/nt_local(myid+1),1/),CLMC_pt1)
    STATUS = NF_GET_VARA_REAL(OUTID,VarID(OUTID,'ITY'), (/low_ind(myid+1),2/), (/nt_local(myid+1),1/),CLMC_pt2)
    STATUS = NF_GET_VARA_REAL(OUTID,VarID(OUTID,'ITY'), (/low_ind(myid+1),3/), (/nt_local(myid+1),1/),CLMC_st1)
    STATUS = NF_GET_VARA_REAL(OUTID,VarID(OUTID,'ITY'), (/low_ind(myid+1),4/), (/nt_local(myid+1),1/),CLMC_st2)
    STATUS = NF_GET_VARA_REAL(OUTID,VarID(OUTID,'FVG'), (/low_ind(myid+1),1/), (/nt_local(myid+1),1/),CLMC_pf1)
    STATUS = NF_GET_VARA_REAL(OUTID,VarID(OUTID,'FVG'), (/low_ind(myid+1),2/), (/nt_local(myid+1),1/),CLMC_pf2)
    STATUS = NF_GET_VARA_REAL(OUTID,VarID(OUTID,'FVG'), (/low_ind(myid+1),3/), (/nt_local(myid+1),1/),CLMC_sf1)
    STATUS = NF_GET_VARA_REAL(OUTID,VarID(OUTID,'FVG'), (/low_ind(myid+1),4/), (/nt_local(myid+1),1/),CLMC_sf2)

    if (master_proc) then

       allocate (TILE_ID  (1:ntiles_cn))

       STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'TILE_ID'   ), (/1/), (/NTILES_CN/),TILE_ID)

       do n = 1,ntiles_cn
 
          K = NINT (TILE_ID (n))

          STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'ITY'), (/n,1/), (/1,4/),ityp_offl(K,:))
          STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'FVG'), (/n,1/), (/1,4/),fveg_offl(K,:))
          
          tid_offl (n) = n
          
          do nv = 1,nveg
             if(ityp_offl(K,nv)<0 .or. ityp_offl(K,nv)>npft)    stop 'ityp'
             if(fveg_offl(K,nv)<0..or. fveg_offl(K,nv)>1.00001) stop 'fveg'             
          end do

          if((ityp_offl(K,3) == 0).and.(ityp_offl(K,4) == 0)) then
             if(ityp_offl(K,1) /= 0) then
                ityp_offl(K,3) = ityp_offl(K,1)
             else
                ityp_offl(K,3) = ityp_offl(K,2)
             endif
          endif
          
          if((ityp_offl(K,1) == 0).and.(ityp_offl(K,2) /= 0)) ityp_offl(K,1) = ityp_offl(K,2)
          if((ityp_offl(K,2) == 0).and.(ityp_offl(K,1) /= 0)) ityp_offl(K,2) = ityp_offl(K,1)
          if((ityp_offl(K,3) == 0).and.(ityp_offl(K,4) /= 0)) ityp_offl(K,3) = ityp_offl(K,4)
          if((ityp_offl(K,4) == 0).and.(ityp_offl(K,3) /= 0)) ityp_offl(K,4) = ityp_offl(K,3)
          
       end do

    endif
    
    call MPI_BCAST(tid_offl ,size(tid_offl ),MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_BCAST(ityp_offl,size(ityp_offl),MPI_REAL   ,0,MPI_COMM_WORLD,mpierr)
    call MPI_BCAST(fveg_offl,size(fveg_offl),MPI_REAL   ,0,MPI_COMM_WORLD,mpierr)    

    ! --------------------------------------------------------------------------------
    ! Here we create transfer index array to map offline restarts to output tile space
    ! --------------------------------------------------------------------------------   
      
    ! Loop through NTILES (# of tiles in output array) and find, for each PFT type seperately, 
    !      a tile with the same PFT type in the neighborhood from the SMAP_M09 offline array.  

    Id_loc = -9999

    TILES : do n = low_ind (myid + 1), upp_ind (myid + 1)

!       if(MOD(n,10000) == 0) print *,'In Carbon', myid,n

       local_id = n - low_ind (myid + 1) + 1
                      
       dw = 0.5 !1.0 ! Start with a 1x1 window, then zoom out by increasing the size by 2-deg until 4 similar tiles are found for 4 PFT types        

       ZOOMOUT : do  
                
          ! Min/Max lon/lat of the working window
          ! -------------------------------------
          
          min_lon = MAX(lonn (local_id) - dw, -180.)
          max_lon = MIN(lonn (local_id) + dw,  180.)
          min_lat = MAX(latt (local_id) - dw,  -90.)
          max_lat = MIN(latt (local_id) + dw,   90.) 
          mask = .false.
          mask =  ((latc >= min_lat .and. latc <= max_lat).and.(lonc >= min_lon .and. lonc <= max_lon))
          nplus =  count(mask = mask)
          
          if(nplus < 0) then
             dw = dw + 1.0
             CYCLE
          endif
          
          allocate (sub_tid (1:nplus))
          allocate (sub_lon (1:nplus))
          allocate (sub_lat (1:nplus))
          allocate (rev_dist  (1:nplus))
          allocate (sub_ityp1 (1:nplus))
          allocate (sub_fevg1 (1:nplus))
          allocate (sub_ityp2 (1:nplus))
          allocate (sub_fevg2 (1:nplus))
          allocate (icl_ityp1 (1:nplus))
          
          sub_tid = PACK (tid_offl, mask= mask) 
          sub_lon = PACK (lonc    , mask= mask)
          sub_lat = PACK (latc    , mask= mask)
          
          ! compute distance from the tile
          
          sub_lat = sub_lat * MAPL_PI/180.
          sub_lon = sub_lon * MAPL_PI/180.
          
!          do i = 1,nplus
!             sub_dist(i) = haversine(to_radian(latt(local_id)), to_radian(lonn(local_id)), &
!                  sub_lat(i), sub_lon(i)) 
!             !                   sub_dist(i) = dist
!          end do
          
          ! loop through 4 vegetation types
           
          NVLOOP : do nv = 1, nveg
             
             if (nv == 1) ityp_new = CLMC_pt1(local_id)
             if (nv == 1) fveg_new = CLMC_pf1(local_id)
             if (nv == 2) ityp_new = CLMC_pt2(local_id)
             if (nv == 2) fveg_new = CLMC_pf2(local_id)
             if (nv == 3) ityp_new = CLMC_st1(local_id)
             if (nv == 3) fveg_new = CLMC_sf1(local_id)
             if (nv == 4) ityp_new = CLMC_st2(local_id) 
             if (nv == 4) fveg_new = CLMC_sf2(local_id)
             
             SEEK : if((Id_loc(local_id,nv) < 0).and.(fveg_new > fmin)) then
                
                if(nv <= 2) then ! index for secondary PFT index if primary or primary if secondary
                   nx = nv + 2
                else
                   nx = nv - 2
                endif
                
                sub_ityp1 = ityp_offl (sub_tid,nv)
                sub_fevg1 = fveg_offl (sub_tid,nv)
                sub_ityp2 = ityp_offl (sub_tid,nx)
                sub_fevg2 = fveg_offl (sub_tid,nx)
                
                rev_dist  = 1.e20
                icl_ityp1 = iclass(sub_ityp1)
                
                do i = 1,nplus
                   if((sub_ityp1(i)>fmin .and. (ityp_new ==sub_ityp1(i) .or.   &
                        iclass(ityp_new) ==iclass(sub_ityp1(i)))) .or.             &
                        (sub_fevg2(i)>fmin .and. (ityp_new ==sub_ityp2(i) .or. &
                        iclass(ityp_new)==iclass(sub_ityp2(i))))) then

                      sub_dist = haversine(to_radian(latt(local_id)), to_radian(lonn(local_id)), &
                           sub_lat(i), sub_lon(i))
                      
                      if(ityp_new == sub_ityp1(i) .and. sub_fevg1(i) >fmin) then
                         rev_dist(i) = 1.*sub_dist     ! give priority to same (primary if primary, secondary if secondary)   
                         ! gkw: these weights are tunable
                      else if(ityp_new ==sub_ityp2(i) .and. sub_fevg2(i)>fmin) then
                         rev_dist(i) = 2.*sub_dist     ! lower priority if not same (secondary if primary, primary if secondary)
                      else if(iclass(ityp_new)==iclass(sub_ityp1(i)) .and. sub_fevg1(i)>fmin) then
                         rev_dist(i) = 3.*sub_dist     ! even lower priority if same of some other PFT in same class
                      else if(sub_fevg2(i)>fmin) then
                         rev_dist(i) = 4.*sub_dist     ! even lower priority if not same of some other PFT in same class
                      else
                         rev_dist(i) = 1.e20
                      endif
                   endif                   
                end do

!                where ((sub_ityp1 == ityp_new) .and. (sub_fevg1 > fmin))  
!                   rev_dist = 1. * rev_dist    ! give priority to same (primary if primary, secondary if secondary)  
!                elsewhere ((sub_ityp2 == ityp_new) .and. (sub_fevg2 > fmin))  
!                   rev_dist = 2. * rev_dist    ! lower priority if not same (secondary if primary, primary if secondary)       
!                elsewhere ((icl_ityp1 == iclass(ityp_new)) .and. (sub_fevg1 > fmin))  
!                   rev_dist = 3. * rev_dist    ! even lower priority if same of some other PFT in same class
!                elsewhere (sub_fevg2 > fmin)
!                   rev_dist = 4. * rev_dist    ! even lower priority if not same of some other PFT in same class
!                elsewhere
!                   rev_dist = 1.e20
!                endwhere

                FOUND : if(minval (rev_dist) < 1.e19) then 
                   
                   Id_loc(local_id,nv) = sub_tid(minloc(rev_dist,1)) ! Cell ID of the nearest neightbor in the offline 
                                                                     ! restart file that has the same veg type
                                      
                endif FOUND
                
             endif SEEK
          end do NVLOOP
          
          deallocate (sub_tid, sub_lon, sub_lat, icl_ityp1)
          deallocate (sub_ityp1, sub_fevg1, sub_ityp2, sub_fevg2, rev_dist)  
          
          ! if similar types have been found for ITYPs with fveg > fmin within the window, exit
           
          all_found = .true.
          
          if((all_found).and.((CLMC_pf1(local_id) > fmin).and.(Id_loc(local_id,1) < 0))) all_found = .false.
          if((all_found).and.((CLMC_pf2(local_id) > fmin).and.(Id_loc(local_id,2) < 0))) all_found = .false.
          if((all_found).and.((CLMC_sf1(local_id) > fmin).and.(Id_loc(local_id,3) < 0))) all_found = .false.
          if((all_found).and.((CLMC_sf2(local_id) > fmin).and.(Id_loc(local_id,4) < 0))) all_found = .false.
          
          if(all_found) GO TO 100
           
           ! if not increase the window size
           dw = dw + 1.
           
        end do ZOOMOUT
        
100   continue !  if(mod (n,1000) == 0) print *, myid +1, n, Id_loc(local_id,:)

     END DO TILES

     deallocate (CLMC_pf1, CLMC_pf2, CLMC_sf1, CLMC_sf2)
     deallocate (CLMC_pt1, CLMC_pt2, CLMC_st1, CLMC_st2)
 
     ! update id_glb in root

     if(master_proc)  then
        allocate (id_glb  (ntiles,4))
        allocate (id_vec  (ntiles))
     endif
           
     do nv = 1, nveg
        call MPI_Barrier(MPI_COMM_WORLD, STATUS)
        call MPI_GATHERV( &
                   id_loc (:,nv), nt_local(myid+1)  , MPI_real, &
                   id_vec, nt_local,low_ind-1, MPI_real, &
                   0, MPI_COMM_WORLD, mpierr )
        
        if(master_proc) id_glb (:,nv) = id_vec
        
     end do

! write out regridded carbon variables

     if(master_proc) then

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
              STATUS = NF_PUT_VARA_REAL(OutID,VarID(OutID,'CNCOL'), (/1,i/), (/NTILES,1 /),var_col_out(:, nz,nv))
              i = i + 1
           end do
        end do
        
        i = 1
        do iv = 1,VAR_PFT
           do nv = 1,nveg
              do nz = 1,nzone
                 STATUS = NF_PUT_VARA_REAL(OutID,VarID(OutID,'CNPFT'), (/1,i/), (/NTILES,1 /),var_pft_out(:, nz,nv,iv))
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
  
  SUBROUTINE regrid_carbon_vars_clm45 (NTILES)

    implicit none

    integer, intent (in)                 :: NTILES
    character*300          :: OutTileFile = 'OutData1/OutTileFile', OutFileName='OutData2/catchcn_internal_clm45'
    integer                :: AGCM_YY=2015,AGCM_MM=1,AGCM_DD=1,AGCM_HR=0 
    real, allocatable, dimension (:)     :: CLMC_pf1, CLMC_pf2, CLMC_sf1, CLMC_sf2, &
         CLMC_pt1, CLMC_pt2,CLMC_st1,CLMC_st2

    ! ===============================================================================================

    integer :: iclass(npft_clm45) = &
         (/1,1,2,3,3,4,5,5,6,7,8,9,10,11,12,11,12,11,11,11,11,11,11,11,11,11,11/)
    integer, allocatable, dimension(:,:) :: Id_glb, Id_loc
    integer, allocatable, dimension(:)   :: tid_offl, id_vec
    logical, allocatable, dimension(:)   :: mask
    real,    allocatable, dimension(:,:) :: fveg_offl,  ityp_offl
    real    :: dw, min_lon, max_lon, min_lat, max_lat, fveg_new, sub_dist
    integer :: n,i,j, k, nplus, nv, nx, nz, iv, offl_cell, ityp_new, STATUS,NCFID
    integer :: outid, local_id
    integer, allocatable, dimension (:) :: sub_tid, sub_ityp1, sub_ityp2,icl_ityp1
    real   , allocatable, dimension (:) :: sub_lon, sub_lat, rev_dist, sub_fevg1, sub_fevg2,&
         lonc, latc, LATT, LONN, DAYX, long, latg, var_dum, TILE_ID, var_dum2
    real, allocatable :: var_off_col (:,:,:), var_off_pft (:,:,:,:) 
    real, allocatable :: var_col_out (:,:,:), var_pft_out (:,:,:,:) 
    integer, allocatable :: low_ind(:), upp_ind(:), nt_local (:)

    logical :: all_found
  
    allocate (tid_offl  (ntiles_cn))
    allocate (mask      (ntiles_cn))
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

    if (master_proc) then
       
       ! --------------------------------------------
       ! Read exact lonn, latt from output .til file 
       ! --------------------------------------------

       allocate (long   (ntiles))
       allocate (latg   (ntiles))
       allocate (DAYX   (NTILES))

       call ReadCNTilFile (OutTileFile, ntiles, long, latg)

       ! Compute DAYX
       ! ------------

       call compute_dayx (                                     &
            NTILES, AGCM_YY, AGCM_MM, AGCM_DD, AGCM_HR,        &
            LATG, DAYX)   

       ! ---------------------------------------------
       ! Read exact lonc, latc from offline .til File 
       ! ---------------------------------------------

       call ReadCNTilFile(trim(InCNTilFile),ntiles_cn,lonc,latc)

    endif

    call MPI_SCATTERV (                    &
         long,nt_local,low_ind-1,MPI_real, &
         lonn,size(lonn),MPI_real  , &
         0,MPI_COMM_WORLD, mpierr )

    call MPI_SCATTERV (                    &
         latg,nt_local,low_ind-1,MPI_real, &
         latt,nt_local(myid+1),MPI_real  , &
         0,MPI_COMM_WORLD, mpierr )

    if(master_proc) deallocate (long, latg)
 
    call MPI_BCAST(lonc,ntiles_cn,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_BCAST(latc,ntiles_cn,MPI_REAL,0,MPI_COMM_WORLD,mpierr)

    ! Open GKW/Fzeng SMAP M09 catchcn_internal_rst and output catchcn_internal_rst
    ! ----------------------------------------------------------------------------
    
    ! NF_OPEN_PAR is no longer needed since IO is done by the root processor.
    !    call MPI_Info_create(info, STATUS)
    !    call MPI_Info_set(info, "romio_cb_read", "automatic", STATUS)   
    !    STATUS = NF_OPEN_PAR   (trim(InCNRestart),IOR(NF_NOWRITE,NF_MPIIO),MPI_COMM_WORLD, info,NCFID)
    !    STATUS = NF_OPEN_PAR   (trim(OutFileName),IOR(NF_WRITE  ,NF_MPIIO),MPI_COMM_WORLD, info,OUTID)
    
    STATUS = NF_OPEN (trim(InCNRestart),NF_NOWRITE,NCFID)
    IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS, 'INPUT RESTART FAILED')
    STATUS = NF_OPEN (trim(OutFileName),NF_WRITE,OUTID)
    IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS, 'OUTPUT RESTART FAILED')
    
    STATUS = NF_GET_VARA_REAL(OUTID,VarID(OUTID,'ITY'), (/low_ind(myid+1),1/), (/nt_local(myid+1),1/),CLMC_pt1)
    STATUS = NF_GET_VARA_REAL(OUTID,VarID(OUTID,'ITY'), (/low_ind(myid+1),2/), (/nt_local(myid+1),1/),CLMC_pt2)
    STATUS = NF_GET_VARA_REAL(OUTID,VarID(OUTID,'ITY'), (/low_ind(myid+1),3/), (/nt_local(myid+1),1/),CLMC_st1)
    STATUS = NF_GET_VARA_REAL(OUTID,VarID(OUTID,'ITY'), (/low_ind(myid+1),4/), (/nt_local(myid+1),1/),CLMC_st2)
    STATUS = NF_GET_VARA_REAL(OUTID,VarID(OUTID,'FVG'), (/low_ind(myid+1),1/), (/nt_local(myid+1),1/),CLMC_pf1)
    STATUS = NF_GET_VARA_REAL(OUTID,VarID(OUTID,'FVG'), (/low_ind(myid+1),2/), (/nt_local(myid+1),1/),CLMC_pf2)
    STATUS = NF_GET_VARA_REAL(OUTID,VarID(OUTID,'FVG'), (/low_ind(myid+1),3/), (/nt_local(myid+1),1/),CLMC_sf1)
    STATUS = NF_GET_VARA_REAL(OUTID,VarID(OUTID,'FVG'), (/low_ind(myid+1),4/), (/nt_local(myid+1),1/),CLMC_sf2)

    if (master_proc) then

       allocate (TILE_ID  (1:ntiles_cn))

       STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'TILE_ID'   ), (/1/), (/NTILES_CN/),TILE_ID)

       do n = 1,ntiles_cn
 
          K = NINT (TILE_ID (n))

          STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'ITY'), (/n,1/), (/1,4/),ityp_offl(K,:))
          STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'FVG'), (/n,1/), (/1,4/),fveg_offl(K,:))
          
          tid_offl (n) = n
          
          do nv = 1,nveg
             if(ityp_offl(K,nv)<0 .or. ityp_offl(K,nv)>npft)    stop 'ityp'
             if(fveg_offl(K,nv)<0..or. fveg_offl(K,nv)>1.00001) stop 'fveg'             
          end do

          if((ityp_offl(K,3) == 0).and.(ityp_offl(K,4) == 0)) then
             if(ityp_offl(K,1) /= 0) then
                ityp_offl(K,3) = ityp_offl(K,1)
             else
                ityp_offl(K,3) = ityp_offl(K,2)
             endif
          endif
          
          if((ityp_offl(K,1) == 0).and.(ityp_offl(K,2) /= 0)) ityp_offl(K,1) = ityp_offl(K,2)
          if((ityp_offl(K,2) == 0).and.(ityp_offl(K,1) /= 0)) ityp_offl(K,2) = ityp_offl(K,1)
          if((ityp_offl(K,3) == 0).and.(ityp_offl(K,4) /= 0)) ityp_offl(K,3) = ityp_offl(K,4)
          if((ityp_offl(K,4) == 0).and.(ityp_offl(K,3) /= 0)) ityp_offl(K,4) = ityp_offl(K,3)
          
       end do

    endif
    
    call MPI_BCAST(tid_offl ,size(tid_offl ),MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_BCAST(ityp_offl,size(ityp_offl),MPI_REAL   ,0,MPI_COMM_WORLD,mpierr)
    call MPI_BCAST(fveg_offl,size(fveg_offl),MPI_REAL   ,0,MPI_COMM_WORLD,mpierr)    

    ! --------------------------------------------------------------------------------
    ! Here we create transfer index array to map offline restarts to output tile space
    ! --------------------------------------------------------------------------------   
      
    ! Loop through NTILES (# of tiles in output array) and find, for each PFT type seperately, 
    !      a tile with the same PFT type in the neighborhood from the SMAP_M09 offline array.  

    Id_loc = -9999

    TILES : do n = low_ind (myid + 1), upp_ind (myid + 1)

!       if(MOD(n,10000) == 0) print *,'In Carbon CLM4.5', myid,n

       local_id = n - low_ind (myid + 1) + 1
                      
       dw = 0.5 !1.0 ! Start with a 1x1 window, then zoom out by increasing the size by 2-deg until 4 similar tiles are found for 4 PFT types        

       ZOOMOUT : do  
                
          ! Min/Max lon/lat of the working window
          ! -------------------------------------
          
          min_lon = MAX(lonn (local_id) - dw, -180.)
          max_lon = MIN(lonn (local_id) + dw,  180.)
          min_lat = MAX(latt (local_id) - dw,  -90.)
          max_lat = MIN(latt (local_id) + dw,   90.) 
          mask = .false.
          mask =  ((latc >= min_lat .and. latc <= max_lat).and.(lonc >= min_lon .and. lonc <= max_lon))
          nplus =  count(mask = mask)
          
          if(nplus < 0) then
             dw = dw + 1.0
             CYCLE
          endif
          
          allocate (sub_tid (1:nplus))
          allocate (sub_lon (1:nplus))
          allocate (sub_lat (1:nplus))
          allocate (rev_dist  (1:nplus))
          allocate (sub_ityp1 (1:nplus))
          allocate (sub_fevg1 (1:nplus))
          allocate (sub_ityp2 (1:nplus))
          allocate (sub_fevg2 (1:nplus))
          allocate (icl_ityp1 (1:nplus))
          
          sub_tid = PACK (tid_offl, mask= mask) 
          sub_lon = PACK (lonc    , mask= mask)
          sub_lat = PACK (latc    , mask= mask)
          
          ! compute distance from the tile
          
          sub_lat = sub_lat * MAPL_PI/180.
          sub_lon = sub_lon * MAPL_PI/180.
          
!          do i = 1,nplus
!             sub_dist(i) = haversine(to_radian(latt(local_id)), to_radian(lonn(local_id)), &
!                  sub_lat(i), sub_lon(i)) 
!             !                   sub_dist(i) = dist
!          end do
          
          ! loop through 4 vegetation types
           
          NVLOOP : do nv = 1, nveg
             
             if (nv == 1) ityp_new = CLMC_pt1(local_id)
             if (nv == 1) fveg_new = CLMC_pf1(local_id)
             if (nv == 2) ityp_new = CLMC_pt2(local_id)
             if (nv == 2) fveg_new = CLMC_pf2(local_id)
             if (nv == 3) ityp_new = CLMC_st1(local_id)
             if (nv == 3) fveg_new = CLMC_sf1(local_id)
             if (nv == 4) ityp_new = CLMC_st2(local_id) 
             if (nv == 4) fveg_new = CLMC_sf2(local_id)
             
             SEEK : if((Id_loc(local_id,nv) < 0).and.(fveg_new > fmin)) then
                
                if(nv <= 2) then ! index for secondary PFT index if primary or primary if secondary
                   nx = nv + 2
                else
                   nx = nv - 2
                endif
                
                sub_ityp1 = ityp_offl (sub_tid,nv)
                sub_fevg1 = fveg_offl (sub_tid,nv)
                sub_ityp2 = ityp_offl (sub_tid,nx)
                sub_fevg2 = fveg_offl (sub_tid,nx)
                
                rev_dist  = 1.e20
                icl_ityp1 = iclass(sub_ityp1)
                
                do i = 1,nplus
                   if((sub_ityp1(i)>fmin .and. (ityp_new ==sub_ityp1(i) .or.   &
                        iclass(ityp_new) ==iclass(sub_ityp1(i)))) .or.             &
                        (sub_fevg2(i)>fmin .and. (ityp_new ==sub_ityp2(i) .or. &
                        iclass(ityp_new)==iclass(sub_ityp2(i))))) then

                      sub_dist = haversine(to_radian(latt(local_id)), to_radian(lonn(local_id)), &
                           sub_lat(i), sub_lon(i))
                      
                      if(ityp_new == sub_ityp1(i) .and. sub_fevg1(i) >fmin) then
                         rev_dist(i) = 1.*sub_dist     ! give priority to same (primary if primary, secondary if secondary)   
                         ! gkw: these weights are tunable
                      else if(ityp_new ==sub_ityp2(i) .and. sub_fevg2(i)>fmin) then
                         rev_dist(i) = 2.*sub_dist     ! lower priority if not same (secondary if primary, primary if secondary)
                      else if(iclass(ityp_new)==iclass(sub_ityp1(i)) .and. sub_fevg1(i)>fmin) then
                         rev_dist(i) = 3.*sub_dist     ! even lower priority if same of some other PFT in same class
                      else if(sub_fevg2(i)>fmin) then
                         rev_dist(i) = 4.*sub_dist     ! even lower priority if not same of some other PFT in same class
                      else
                         rev_dist(i) = 1.e20
                      endif
                   endif                   
                end do

!                where ((sub_ityp1 == ityp_new) .and. (sub_fevg1 > fmin))  
!                   rev_dist = 1. * rev_dist    ! give priority to same (primary if primary, secondary if secondary)  
!                elsewhere ((sub_ityp2 == ityp_new) .and. (sub_fevg2 > fmin))  
!                   rev_dist = 2. * rev_dist    ! lower priority if not same (secondary if primary, primary if secondary)       
!                elsewhere ((icl_ityp1 == iclass(ityp_new)) .and. (sub_fevg1 > fmin))  
!                   rev_dist = 3. * rev_dist    ! even lower priority if same of some other PFT in same class
!                elsewhere (sub_fevg2 > fmin)
!                   rev_dist = 4. * rev_dist    ! even lower priority if not same of some other PFT in same class
!                elsewhere
!                   rev_dist = 1.e20
!                endwhere

                FOUND : if(minval (rev_dist) < 1.e19) then 
                   
                   Id_loc(local_id,nv) = sub_tid(minloc(rev_dist,1)) ! Cell ID of the nearest neightbor in the offline 
                                                                     ! restart file that has the same veg type
                                      
                endif FOUND
                
             endif SEEK
          end do NVLOOP
          
          deallocate (sub_tid, sub_lon, sub_lat, icl_ityp1)
          deallocate (sub_ityp1, sub_fevg1, sub_ityp2, sub_fevg2, rev_dist)  
          
          ! if similar types have been found for ITYPs with fveg > fmin within the window, exit
           
          all_found = .true.
          
          if((all_found).and.((CLMC_pf1(local_id) > fmin).and.(Id_loc(local_id,1) < 0))) all_found = .false.
          if((all_found).and.((CLMC_pf2(local_id) > fmin).and.(Id_loc(local_id,2) < 0))) all_found = .false.
          if((all_found).and.((CLMC_sf1(local_id) > fmin).and.(Id_loc(local_id,3) < 0))) all_found = .false.
          if((all_found).and.((CLMC_sf2(local_id) > fmin).and.(Id_loc(local_id,4) < 0))) all_found = .false.
          
          if(all_found) GO TO 100
           
           ! if not increase the window size
           dw = dw + 1.
           
        end do ZOOMOUT
        
100   continue !  if(mod (n,1000) == 0) print *, myid +1, n, Id_loc(local_id,:)

     END DO TILES

     deallocate (CLMC_pf1, CLMC_pf2, CLMC_sf1, CLMC_sf2)
     deallocate (CLMC_pt1, CLMC_pt2, CLMC_st1, CLMC_st2)
 
     ! update id_glb in root

     if(master_proc)  then
        allocate (id_glb  (ntiles,4))
        allocate (id_vec  (ntiles))
     endif
           
     do nv = 1, nveg
        call MPI_Barrier(MPI_COMM_WORLD, STATUS)
        call MPI_GATHERV( &
                   id_loc (:,nv), nt_local(myid+1)  , MPI_real, &
                   id_vec, nt_local,low_ind-1, MPI_real, &
                   0, MPI_COMM_WORLD, mpierr )
        
        if(master_proc) id_glb (:,nv) = id_vec
        
     end do

! write out regridded carbon variables

     if(master_proc) then

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
           !  1 ccs%col_ctrunc_vr   (:,1)     
           !  2 ccs%decomp_cpools_vr(:,1,4)   ! cwdc
           !  3 ccs%decomp_cpools_vr(:,1,1)   ! litr1c
           !  4 ccs%decomp_cpools_vr(:,1,2)   ! litr2c
           !  5 ccs%decomp_cpools_vr(:,1,3)   ! litr3c
           !  6 ccs%totvegc_col         
           !  7 ccs%prod100c            
           !  8 ccs%prod10c             
           !  9 ccs%seedc               
           ! 10 ccs%decomp_cpools_vr(:,1,5)   ! soil1c
           ! 11 ccs%decomp_cpools_vr(:,1,6)   ! soil2c
           ! 12 ccs%decomp_cpools_vr(:,1,7)   ! soil3c
           ! 13 ccs%decomp_cpools_vr(:,1,8)   ! soil4c
           ! 14 ccs%totcolc             
           ! 15 ccs%totlitc             
           ! 16 cns%col_ntrunc_vr   (:,1)  
           ! 17 cns%decomp_npools_vr(:,1,4)   ! cwdn
           ! 18 cns%decomp_npools_vr(:,1,1)   ! litr1n
           ! 19 cns%decomp_npools_vr(:,1,2)   ! litr2n
           ! 20 cns%decomp_npools_vr(:,1,3)   ! litr3n
           ! 21 cns%prod100n            
           ! 22 cns%prod10n             
           ! 23 cns%seedn               
           ! 24 cns%sminn_vr        (:,1)  
           ! 25 cns%decomp_npools_vr(:,1,5)   ! soil1n
           ! 26 cns%decomp_npools_vr(:,1,6)   ! soil2n
           ! 27 cns%decomp_npools_vr(:,1,7)   ! soil3n
           ! 28 cns%decomp_npools_vr(:,1,8)   ! soil4n
           ! 29 cns%totcoln             
           ! 30 cps%altmax              
           ! 31 cps%annsum_counter      
           ! 32 cps%cannavg_t2m         
           ! 33 cps%cannsum_npp         
           ! 34 cps%farea_burned        
           ! 35 cps%altmax_lastyear     
           ! 36 cps%altmax_indx         
           ! 37 cps%fpg                 
           ! 38 cps%fpi_vr          (:,1)  
           ! 39 cps%altmax_lastyear_indx        
           
           
           ! PFT vars
           ! --------
           !  1 pcs%cpool                 
           !  2 pcs%deadcrootc            
           !  3 pcs%deadcrootc_storage    
           !  4 pcs%deadcrootc_xfer       
           !  5 pcs%deadstemc             
           !  6 pcs%deadstemc_storage     
           !  7 pcs%deadstemc_xfer        
           !  8 pcs%frootc                
           !  9 pcs%frootc_storage        
           ! 10 pcs%frootc_xfer           
           ! 11 pcs%gresp_storage         
           ! 12 pcs%gresp_xfer            
           ! 13 pcs%leafc                 
           ! 14 pcs%leafc_storage         
           ! 15 pcs%leafc_xfer            
           ! 16 pcs%livecrootc            
           ! 17 pcs%livecrootc_storage    
           ! 18 pcs%livecrootc_xfer       
           ! 19 pcs%livestemc             
           ! 20 pcs%livestemc_storage     
           ! 21 pcs%livestemc_xfer        
           ! 22 pcs%pft_ctrunc            
           ! 23 pcs%xsmrpool              
           ! 24 pepv%annavg_t2m           
           ! 25 pepv%annmax_retransn      
           ! 26 pepv%annsum_npp           
           ! 27 pepv%annsum_potential_gpp 
           ! 28 pepv%dayl                 
           ! 29 pepv%days_active          
           ! 30 pepv%dormant_flag         
           ! 31 pepv%offset_counter       
           ! 32 pepv%offset_fdd           
           ! 33 pepv%offset_flag          
           ! 34 pepv%offset_swi           
           ! 35 pepv%onset_counter        
           ! 36 pepv%onset_fdd            
           ! 37 pepv%onset_flag           
           ! 38 pepv%onset_gdd            
           ! 39 pepv%onset_gddflag        
           ! 40 pepv%onset_swi            
           ! 41 pepv%prev_frootc_to_litter
           ! 42 pepv%prev_leafc_to_litter 
           ! 43 pepv%tempavg_t2m          
           ! 44 pepv%tempmax_retransn     
           ! 45 pepv%tempsum_npp          
           ! 46 pepv%tempsum_potential_gpp
           ! 47 pepv%xsmrpool_recover     
           ! 48 pns%deadcrootn            
           ! 49 pns%deadcrootn_storage    
           ! 50 pns%deadcrootn_xfer       
           ! 51 pns%deadstemn             
           ! 52 pns%deadstemn_storage     
           ! 53 pns%deadstemn_xfer        
           ! 54 pns%frootn                
           ! 55 pns%frootn_storage        
           ! 56 pns%frootn_xfer           
           ! 57 pns%leafn                 
           ! 58 pns%leafn_storage         
           ! 59 pns%leafn_xfer            
           ! 60 pns%livecrootn            
           ! 61 pns%livecrootn_storage    
           ! 62 pns%livecrootn_xfer       
           ! 63 pns%livestemn             
           ! 64 pns%livestemn_storage     
           ! 65 pns%livestemn_xfer        
           ! 66 pns%npool                 
           ! 67 pns%pft_ntrunc            
           ! 68 pns%retransn              
           ! 69 pps%elai                  
           ! 70 pps%esai                  
           ! 71 pps%hbot                  
           ! 72 pps%htop                  
           ! 73 pps%tlai                  
           ! 74 pps%tsai                  
           ! 75 pps%gddplant              
           ! 76 pps%gddtsoi               
           ! 77 pps%peaklai               
           ! 78 pps%idop                  
           ! 79 pps%aleaf                 
           ! 80 pps%aleafi                
           ! 81 pps%astem                 
           ! 82 pps%astemi                
           ! 83 pps%htmx                  
           ! 84 pps%hdidx                 
           ! 85 pps%vf                    
           ! 86 pps%cumvd                 
           ! 87 pps%croplive              
           ! 88 pps%cropplant             
           ! 89 pps%harvdate              
           ! 90 pps%gdd1020               
           ! 91 pps%gdd820                
           ! 92 pps%gdd020                
           ! 93 pps%gddmaturity           
           ! 94 pps%huileaf               
           ! 95 pps%huigrain              
           ! 96 pcs%grainc                
           ! 97 pcs%grainc_storage        
           ! 98 pcs%grainc_xfer           
           ! 99 pns%grainn                
           !100 pns%grainn_storage        
           !101 pns%grainn_xfer           
           !102 pepv%fert_counter         
           !103 pnf%fert                  
           !104 pepv%grain_flag           
           
        end do OUT_TILE
        
        i = 1
        VAR_DUM = 0.
        do nv = 1,VAR_COL_CLM45
           do nz = 1,nzone
              STATUS = NF_PUT_VARA_REAL(OutID,VarID(OutID,'CNCOL'), (/1,i/), (/NTILES,1 /),var_col_out(:, nz,nv))
              if((nv == 30).OR.(nv == 35).OR.(nv == 36).OR.(nv == 39)) &
                 STATUS = NF_PUT_VARA_REAL(OutID,VarID(OutID,'CNCOL'), (/1,i/), (/NTILES,1 /),var_dum)  
              i = i + 1
           end do
        end do
        
        i = 1
        do iv = 1,VAR_PFT_CLM45
           do nv = 1,nveg
              do nz = 1,nzone
                 if(iv <= 74) then
                    STATUS = NF_PUT_VARA_REAL(OutID,VarID(OutID,'CNPFT'), (/1,i/), (/NTILES,1 /),var_pft_out(:, nz,nv,iv))
                 else
                    if((iv == 78) .OR. (iv == 89)) then    ! idop and harvdate
                       var_dum = 999
                       STATUS = NF_PUT_VARA_REAL(OutID,VarID(OutID,'CNPFT'), (/1,i/), (/NTILES,1 /),var_dum)
                    else
                       var_dum = 0.
                       STATUS = NF_PUT_VARA_REAL(OutID,VarID(OutID,'CNPFT'), (/1,i/), (/NTILES,1 /),var_dum)
                    endif
                 endif
                 i = i + 1
              end do
           end do
        end do

        VAR_DUM = 0.

        do nz = 1,nzone
           STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'TGWM'), (/1,nz/), (/NTILES,1 /),VAR_DUM(:))
           STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'RZMM'), (/1,nz/), (/NTILES,1 /),VAR_DUM(:)) 
           STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'SFMM'), (/1,nz/), (/NTILES,1 /),VAR_DUM(:))          
        end do

        STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'BFLOWM'), (/1/), (/NTILES/),VAR_DUM(:))
        STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'TOTWATM'), (/1/), (/NTILES/),VAR_DUM(:))
        STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'TAIRM'), (/1/), (/NTILES/),VAR_DUM(:))
        STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'TPM'), (/1/), (/NTILES/),VAR_DUM(:))
        STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'CNSUM'), (/1/), (/NTILES/),VAR_DUM(:)) 
        STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'SNDZM'), (/1/), (/NTILES/),VAR_DUM(:))
        STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'ASNOWM'), (/1/), (/NTILES/),VAR_DUM(:))
        STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'AR1M'), (/1/), (/NTILES/),VAR_DUM(:))
        STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'RAINFM'), (/1/), (/NTILES/),VAR_DUM(:))
        STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'RHM'), (/1/), (/NTILES/),VAR_DUM(:))
        STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'RUNSRFM'), (/1/), (/NTILES/),VAR_DUM(:))
        STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'SNOWFM'), (/1/), (/NTILES/),VAR_DUM(:))
        STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'WINDM'), (/1/), (/NTILES/),VAR_DUM(:))
        STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'TPREC10D'), (/1/), (/NTILES/),VAR_DUM(:))
        STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'TPREC60D'), (/1/), (/NTILES/),VAR_DUM(:))
        STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'T2M1H'), (/1/), (/NTILES/),VAR_DUM(:))
        STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'GDD0'), (/1/), (/NTILES/),VAR_DUM(:))
        STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'GDD8'), (/1/), (/NTILES/),VAR_DUM(:))
        STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'GDD10'), (/1/), (/NTILES/),VAR_DUM(:))
 
        do nv = 1,nzone
           do nz = 1,nveg
              STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'PSNSUNM'), (/1,nz,nv/), (/NTILES,1,1/),VAR_DUM(:)) 
              STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'PSNSHAM'), (/1,nz,nv/), (/NTILES,1,1/),VAR_DUM(:)) 
           end do
        end do

        VAR_DUM = MAPL_TICE + 20.
        STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'T2M10D'), (/1/), (/NTILES/),VAR_DUM(:))

        VAR_DUM = -1.e36
        STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'T2MMAXI'), (/1/), (/NTILES/),VAR_DUM(:))

        VAR_DUM = 1.e36
        STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'T2MMAX'), (/1/), (/NTILES/),VAR_DUM(:))
        STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'T2MMINI'), (/1/), (/NTILES/),VAR_DUM(:))
        STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'T2MMIN'), (/1/), (/NTILES/),VAR_DUM(:))

        VAR_DUM = MAPL_TICE
        STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'T2MIN5D'), (/1/), (/NTILES/),VAR_DUM(:))
        STATUS = NF_PUT_VARA_REAL(OutID,VarID(OUTID,'T2MIN10D'), (/1/), (/NTILES/),VAR_DUM(:))

        STATUS = NF_CLOSE (NCFID)
        STATUS = NF_CLOSE (OutID)

        deallocate (var_off_col,var_off_pft,var_col_out,var_pft_out)  
        deallocate (CLMC_pf1, CLMC_pf2, CLMC_sf1, CLMC_sf2)
        deallocate (CLMC_pt1, CLMC_pt2, CLMC_st1, CLMC_st2)
 
     endif

     call MPI_Barrier(MPI_COMM_WORLD, STATUS)

   END SUBROUTINE regrid_carbon_vars_clm45


  ! *****************************************************************************

   SUBROUTINE put_land_vars (NTILES, ntiles_smap, id_glb, ld_reorder, model)

     implicit none
     character(*), intent (in)  :: model   
     integer, intent (in)       :: NTILES, ntiles_smap
     integer, intent (in)       :: id_glb(NTILES), ld_reorder (ntiles_smap)
     integer                    :: k
     real   , dimension (:), allocatable :: var_get, var_put
     type(MAPL_NCIO) :: InNCIO2, InNCIO,  InCLM45Dummy, InNCIO3 
     integer         :: nVars, STATUS, NCFID,nVars45

     allocate (var_get (NTILES_SMAP))
     allocate (var_put (NTILES))
   
     ! create output catchcn_internal_rst
     if(trim(model) == 'CATCHCN') InNCIO = MAPL_NCIOOpen(trim(InCNRestart ) ,rc=rc) 
     if(trim(model) == 'CATCH'  ) InNCIO = MAPL_NCIOOpen(trim(InCatRestart) ,rc=rc) 
     call MAPL_NCIOGetDimSizes(InNCIO,nVars=nVars)
     call MAPL_NCIOChangeRes(InNCIO,InNCIO2,tileSize=ntiles,rc=rc)

     if(trim(model) == 'CATCHCN') OutFileName = "OutData1/catchcn_internal_rst"
     if(trim(model) == 'CATCH'  ) OutFileName = "OutData1/catch_internal_rst"
     call MAPL_NCIOSet( InNCIO2,filename=OutFileName )
     call MAPL_NCIOCreateFile(InNCIO2) 
     call MAPL_NCIOClose ( InNCIO)
    
     if(CLM45) then
        ! create output catchcn_internal_clm45
        InCLM45Dummy = MAPL_NCIOOpen('/discover/nobackup/rreichle/l_data/LandRestarts_for_Regridding/CatchCN/catchcn_internal_clm45', rc=rc)
        call MAPL_NCIOGetDimSizes(InCLM45Dummy, nVars=nVars45)
        call MAPL_NCIOChangeRes(InCLM45Dummy, InNCIO3, tileSize=ntiles,rc=rc)
        
        OutFileName = "OutData1/catchcn_internal_clm45"
        call MAPL_NCIOSet(InNCIO3 ,filename=OutFileName )
        call MAPL_NCIOCreateFile(InNCIO3) 
        call MAPL_NCIOClose (InCLM45Dummy)
     endif

     ! Read catparam
     ! -------------

     if(trim(model) == 'CATCHCN') STATUS = NF_OPEN (trim(InCNRestart ),NF_NOWRITE,NCFID)
     if(trim(model) == 'CATCH')   STATUS = NF_OPEN (trim(InCatRestart),NF_NOWRITE,NCFID)
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'POROS'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'POROS',var_put)
     if(CLM45) call MAPL_VarWrite(InNCIO3,'POROS',var_put)

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'COND'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'COND',var_put)
     if(CLM45) call MAPL_VarWrite(InNCIO3,'COND',var_put)

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'PSIS'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'PSIS',var_put)
     if(CLM45) call MAPL_VarWrite(InNCIO3,'PSIS',var_put)

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'BEE'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'BEE',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'BEE',var_put) 
     
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'WPWET'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'WPWET',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'WPWET',var_put) 
    
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'GNU'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'GNU',var_put)
     if(CLM45) call MAPL_VarWrite(InNCIO3,'GNU',var_put)
 
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'VGWMAX'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'VGWMAX',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'VGWMAX',var_put) 
 
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'BF1'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'BF1',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'BF1',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'BF2'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'BF2',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'BF2',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'BF3'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'BF3',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'BF3',var_put) 
     
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'CDCR1'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'CDCR1',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'CDCR1',var_put) 
     
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'CDCR2'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'CDCR2',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'CDCR2',var_put) 
     
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'ARS1'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'ARS1',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'ARS1',var_put) 
     
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'ARS2'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'ARS2',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'ARS2',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'ARS3'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'ARS3',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'ARS3',var_put) 
     
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'ARA1'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'ARA1',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'ARA1',var_put) 
     
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'ARA2'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'ARA2',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'ARA2',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'ARA3'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'ARA3',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'ARA3',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'ARA4'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'ARA4',var_put)
     if(CLM45) call MAPL_VarWrite(InNCIO3,'ARA4',var_put)

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'ARW1'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'ARW1',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'ARW1',var_put) 
     
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'ARW2'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'ARW2',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'ARW2',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'ARW3'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'ARW3',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'ARW3',var_put) 
     
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'ARW4'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'ARW4',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'ARW4',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'TSA1'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'TSA1',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'TSA1',var_put) 
     
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'TSA2'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'TSA2',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'TSA2',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'TSB1'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'TSB1',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'TSB1',var_put) 
     
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'TSB2'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'TSB2',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'TSB2',var_put) 
     
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'ATAU'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'ATAU',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'ATAU',var_put) 
     
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'BTAU'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'BTAU',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'BTAU',var_put) 

     if(trim(model) == 'CATCHCN') then

        STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'ITY'   ), (/1,1/), (/NTILES_SMAP,1/),var_get)
        do k = 1, NTILES
           VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
        end do
        call MAPL_VarWrite(InNCIO2,'ITY',var_put, offset1=1) 
        if(CLM45) call MAPL_VarWrite(InNCIO3,'ITY',var_put, offset1=1) 
        STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'ITY'   ), (/1,2/), (/NTILES_SMAP,1/),var_get)
        do k = 1, NTILES
           VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
        end do
        call MAPL_VarWrite(InNCIO2,'ITY',var_put, offset1=2) 
        if(CLM45) call MAPL_VarWrite(InNCIO3,'ITY',var_put, offset1=2) 
        STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'ITY'   ), (/1,3/), (/NTILES_SMAP,1/),var_get)
        do k = 1, NTILES
           VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
        end do
        call MAPL_VarWrite(InNCIO2,'ITY',var_put, offset1=3) 
        if(CLM45) call MAPL_VarWrite(InNCIO3,'ITY',var_put, offset1=3) 
        STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'ITY'   ), (/1,4/), (/NTILES_SMAP,1/),var_get)
        do k = 1, NTILES
           VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
        end do
        call MAPL_VarWrite(InNCIO2,'ITY',var_put, offset1=4)
        if(CLM45) call MAPL_VarWrite(InNCIO3,'ITY',var_put, offset1=4)
        
        STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'FVG'   ), (/1,1/), (/NTILES_SMAP,1/),var_get)
        do k = 1, NTILES
           VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
        end do
        call MAPL_VarWrite(InNCIO2,'FVG',var_put, offset1=1) 
        if(CLM45) call MAPL_VarWrite(InNCIO3,'FVG',var_put, offset1=1) 
        STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'FVG'   ), (/1,2/), (/NTILES_SMAP,1/),var_get)
        do k = 1, NTILES
           VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
        end do
        call MAPL_VarWrite(InNCIO2,'FVG',var_put, offset1=2) 
        if(CLM45) call MAPL_VarWrite(InNCIO3,'FVG',var_put, offset1=2) 
        STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'FVG'   ), (/1,3/), (/NTILES_SMAP,1/),var_get)
        do k = 1, NTILES
           VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
        end do
        call MAPL_VarWrite(InNCIO2,'FVG',var_put, offset1=3) 
        if(CLM45) call MAPL_VarWrite(InNCIO3,'FVG',var_put, offset1=3) 
        STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'FVG'   ), (/1,4/), (/NTILES_SMAP,1/),var_get)
        do k = 1, NTILES
           VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
        end do
        call MAPL_VarWrite(InNCIO2,'FVG',var_put, offset1=4)
        if(CLM45) call MAPL_VarWrite(InNCIO3,'FVG',var_put, offset1=4)


     ! read restart and regrid
     ! -----------------------

        STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'TG'   ), (/1,1/), (/NTILES_SMAP,1/),var_get)
        do k = 1, NTILES
           VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
        end do
        call MAPL_VarWrite(InNCIO2,'TG',var_put, offset1=1)  ! if you see offset1=1 it is a 2-D var
        if(CLM45) call MAPL_VarWrite(InNCIO3,'TG',var_put, offset1=1)  ! if you see offset1=1 it is a 2-D var
        
        STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'TG'   ), (/1,2/), (/NTILES_SMAP,1/),var_get)
        do k = 1, NTILES
           VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
        end do
        call MAPL_VarWrite(InNCIO2,'TG',var_put, offset1=2) 
        if(CLM45) call MAPL_VarWrite(InNCIO3,'TG',var_put, offset1=2) 
        
        STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'TG'   ), (/1,3/), (/NTILES_SMAP,1/),var_get)
        do k = 1, NTILES
           VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
        end do
        call MAPL_VarWrite(InNCIO2,'TG',var_put, offset1=3) 
        if(CLM45) call MAPL_VarWrite(InNCIO3,'TG',var_put, offset1=3) 

     endif
        
        STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'TC'   ), (/1,1/), (/NTILES_SMAP,1/),var_get)
        do k = 1, NTILES
           VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
        end do
       
     call MAPL_VarWrite(InNCIO2,'TC',var_put, offset1=1) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'TC',var_put, offset1=1) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'TC'   ), (/1,2/), (/NTILES_SMAP,1/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'TC',var_put, offset1=2) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'TC',var_put, offset1=2) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'TC'   ), (/1,3/), (/NTILES_SMAP,1/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'TC',var_put, offset1=3) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'TC',var_put, offset1=3) 
     
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'QC'   ), (/1,1/), (/NTILES_SMAP,1/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'QC',var_put, offset1=1) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'QC',var_put, offset1=1) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'QC'   ), (/1,2/), (/NTILES_SMAP,1/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'QC',var_put, offset1=2)
     if(CLM45) call MAPL_VarWrite(InNCIO3,'QC',var_put, offset1=2)

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'QC'   ), (/1,3/), (/NTILES_SMAP,1/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'QC',var_put, offset1=3)
     if(CLM45) call MAPL_VarWrite(InNCIO3,'QC',var_put, offset1=3)

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'CAPAC'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'CAPAC',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'CAPAC',var_put) 
     
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'CATDEF'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'CATDEF',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'CATDEF',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'RZEXC'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'RZEXC',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'RZEXC',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'SRFEXC'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'SRFEXC',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'SRFEXC',var_put) 
     
     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'GHTCNT1'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'GHTCNT1',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'GHTCNT1',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'GHTCNT2'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'GHTCNT2',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'GHTCNT2',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'GHTCNT3'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'GHTCNT3',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'GHTCNT3',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'GHTCNT4'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'GHTCNT4',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'GHTCNT4',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'GHTCNT5'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'GHTCNT5',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'GHTCNT5',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'GHTCNT6'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'GHTCNT6',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'GHTCNT6',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'WESNN1'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'WESNN1',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'WESNN1',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'WESNN2'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'WESNN2',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'WESNN2',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'WESNN3'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'WESNN3',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'WESNN3',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'HTSNNN1'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'HTSNNN1',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'HTSNNN1',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'HTSNNN2'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'HTSNNN2',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'HTSNNN2',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'HTSNNN3'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'HTSNNN3',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'HTSNNN3',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'SNDZN1'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'SNDZN1',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'SNDZN1',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'SNDZN2'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'SNDZN2',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'SNDZN2',var_put) 

     STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'SNDZN3'   ), (/1/), (/NTILES_SMAP/),var_get)
     do k = 1, NTILES
        VAR_PUT(k) = var_get(ld_reorder(id_glb(k)))
     end do
     call MAPL_VarWrite(InNCIO2,'SNDZN3',var_put) 
     if(CLM45) call MAPL_VarWrite(InNCIO3,'SNDZN3',var_put) 

     call MAPL_NCIOClose ( InNCIO2)
     if(CLM45) call MAPL_NCIOClose ( InNCIO3)
     STATUS = NF_CLOSE ( NCFID)

     deallocate (var_get, var_put)
  
     if(trim(MODEL) == 'CATCHCN') then
        call system('/bin/cp OutData1/catchcn_internal_rst OutData2/catchcn_internal_rst')
        call system('/bin/cp OutData1/catchcn_internal_clm45 OutData2/catchcn_internal_clm45')
     endif

     if(trim(MODEL) == 'CATCH') then
        call system('/bin/cp OutData1/catch_internal_rst OutData2/catch_internal_rst')
     endif

   END SUBROUTINE put_land_vars

  ! *****************************************************************************
  
  subroutine ReadCNTilFile (InCNTileFile, nt, xlon, xlat)
    
     implicit none
     character(*), intent (in) ::  InCNTileFile
     integer , intent (in) :: nt
     real, dimension (nt), intent(inout) :: xlon, xlat
     integer :: n,icnt,ityp,status
     real    :: xval,yval, pf
     
   open(11,file=InCNTileFile, &
        form='formatted',action='read',status='old')

   do n = 1,8 ! skip header
      read(11,*)
   end do
   
   icnt = 0
   ityp = 100

   do while (ityp == 100) ! loop over land tiles
      read(11,*, iostat = status) ityp,pf,xval,yval
      if (status /= 0) exit
      if((ityp == 100).and.(status == 0)) then
         icnt = icnt + 1
         xlon(icnt) = xval
         xlat(icnt) = yval
      endif
   end do

   close(11)
    

   end subroutine ReadCNTilFile

  ! *****************************************************************************
  
  subroutine init_MPI()
    
    ! initialize MPI
    
    call MPI_INIT(mpierr)
    
    call MPI_COMM_RANK( MPI_COMM_WORLD, myid, mpierr )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, mpierr )

    if (myid .ne. 0)  master_proc = .false.
    
!    call init_MPI_types()
    
    write (*,*) "MPI process ", myid, " of ", numprocs, " is alive"    
    write (*,*) "MPI process ", myid, ": master_proc=", master_proc

  end subroutine init_MPI
  
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

  ! *****************************************************************************

   function to_radian(degree) result(rad)

     ! degrees to radians
     real,intent(in) :: degree
     real :: rad

     rad = degree*MAPL_PI/180.

   end function to_radian

   ! *****************************************************************************
   
   real function haversine(deglat1,deglon1,deglat2,deglon2)
     ! great circle distance -- adapted from Matlab 
     real,intent(in) :: deglat1,deglon1,deglat2,deglon2
     real :: a,c, dlat,dlon,lat1,lat2
     real,parameter :: radius = MAPL_radius
     
!     dlat = to_radian(deglat2-deglat1)
!     dlon = to_radian(deglon2-deglon1)
     !     lat1 = to_radian(deglat1)
!     lat2 = to_radian(deglat2)
     dlat = deglat2-deglat1
     dlon = deglon2-deglon1
     lat1 = deglat1
     lat2 = deglat2     
     a = (sin(dlat/2))**2 + cos(lat1)*cos(lat2)*(sin(dlon/2))**2
     if(a>=0. .and. a<=1.) then
        c = 2*atan2(sqrt(a),sqrt(1-a))
        haversine = radius*c / 1000.
     else
        haversine = 1.e20
     endif
   end function

  ! ----------------------------------------------------------------------

   integer function VarID (NCFID, VNAME) 
     
     integer, intent (in)      :: NCFID
     character(*), intent (in) :: VNAME
     integer                   :: status

     STATUS = NF_INQ_VARID (NCFID, trim(VNAME) ,VarID)
     IF (STATUS .NE. NF_NOERR) &
          CALL HANDLE_ERR(STATUS, trim(VNAME))  
     
   end function VarID
   ! -----------------------------------------------------------------------------
 
   
   FUNCTION StrUpCase ( Input_String ) RESULT ( Output_String ) 
     ! -- Argument and result 
     CHARACTER( * ), INTENT( IN ) :: Input_String 
     CHARACTER( LEN( Input_String ) ) :: Output_String 
     ! -- Local variables 
     INTEGER :: i, n 


     ! -- Copy input string 
     Output_String = Input_String 
     ! -- Loop over string elements 
     DO i = 1, LEN( Output_String ) 
       ! -- Find location of letter in lower case constant string 
       n = INDEX( LOWER_CASE, Output_String( i:i ) ) 
       ! -- If current substring is a lower case letter, make it upper case 
       IF ( n /= 0 ) Output_String( i:i ) = UPPER_CASE( n:n ) 
     END DO 
   END FUNCTION StrUpCase 

   ! -----------------------------------------------------------------------------

   FUNCTION StrLowCase ( Input_String ) RESULT ( Output_String ) 
     ! -- Argument and result 
     CHARACTER( * ), INTENT( IN ) :: Input_String 
     CHARACTER( LEN( Input_String ) ) :: Output_String 
     ! -- Local variables 
     INTEGER :: i, n 

     ! -- Copy input string 
     Output_String = Input_String 
     ! -- Loop over string elements 
     DO i = 1, LEN( Output_String ) 
       ! -- Find location of letter in upper case constant string 
       n = INDEX( UPPER_CASE, Output_String( i:i ) ) 
       ! -- If current substring is an upper case letter, make it lower case 
       IF ( n /= 0 ) Output_String( i:i ) = LOWER_CASE( n:n ) 
     END DO 
   END FUNCTION StrLowCase 

   ! -----------------------------------------------------------------------------
   
   FUNCTION StrExtName ( Input_String ) RESULT ( Output_String ) 
     ! -- Argument and result 
     CHARACTER( * ), INTENT( IN ) :: Input_String 
     CHARACTER( LEN( Input_String ) ) :: Output_String 
     ! -- Local variables 
     INTEGER :: i, n1, n2, n3, n4, n5, n, k 

     ! -- Copy input string 
     ! Output_String = Input_String 
     ! -- Loop over string elements 

     k = 1

     DO i = 1, LEN( Input_String ) 

       ! -- Find location of letter in upper case constant string 
       n1 = INDEX( UPPER_CASE, Input_String( i:i ) )
       n2 = INDEX( LOWER_CASE, Input_String( i:i ) )
       n3 = INDEX( '.', Input_String( i:i ) )
       n4 = INDEX( '-', Input_String( i:i ) )
       n5 = INDEX( '_', Input_String( i:i ) )

       n = 0
       Output_String(i:i) = ''

       if (n1 /= 0) n = n1
       if (n2 /= 0) n = n2
       if (n3 /= 0) n = n3
       if (n4 /= 0) n = n4
       if (n5 /= 0) n = n5
            
       ! -- If current substring is acceptable 
       IF ( n /= 0 ) then
          Output_String( k:k ) = Input_String( i:i )
          k = k + 1
       endif

     END DO 

   END FUNCTION StrExtName

   ! -----------------------------------------------------------------------------

 END PROGRAM mk_LDASsaRestarts
