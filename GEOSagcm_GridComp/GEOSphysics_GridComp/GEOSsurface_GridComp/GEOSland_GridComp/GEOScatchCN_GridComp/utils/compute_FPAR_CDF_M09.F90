#include "Raster.h"

PROGRAM comp_FPAR_CDF

  use math_routines
  use MAPL_SortMod
  use date_time_util,                   ONLY:     &
       date_time_type, augment_date_time
  use ieee_arithmetic, only: isnan => ieee_is_nan  
  IMPLICIT NONE
  
  INCLUDE 'netcdf.inc'
  INCLUDE 'mpif.h'

  integer              :: comm_rank, comm_size, error, info, DOY,NCInID,NCOutID, NCOutID2, n, iargc, STATUS
  integer              :: stat(MPI_STATUS_SIZE), NTILES, NCatch, NDATA1, NDATA2,N_VIS_DATA, N_NIR_DATA
  character*400        :: arg(2)
  REAL, PARAMETER      :: TILESIZE = 10., TILEINT = 2. 
  INTEGER, PARAMETER   :: NPFT = 19, NOCTAD = 46, NSETS = 2 ,MXCNT = 500000 
  logical, parameter   :: onebin = .true.
  INTEGER, PARAMETER   :: YearB = 2003, YearE = 2016
  character*400, PARAMETER ::                &
       BCSDIR = 'SMAP_EASEv2_M09/', &
       EXPDIR = '/discover/nobackup/fzeng/Catchment/SMAP_EASEv2_M09/e0004s_wet2/output/SMAP_EASEv2_M09_GLOBAL/', &
       EXNAME = 'e0004s_wet2',                                                             &
       OUTFIL = 'global_alb_mu_std/FPAR_CDF_Params-M09',                                       &
       LOGFIL = 'global_alb_mu_std/FPAR_CDF_Params-M09_log.',                                  & 
       GFILE  = 'SMAP_EASEv2_M09_3856x1624'

  logical              :: file_exists, put_aux = .true.
!  real, dimension (4)  :: limits = (/20., -130., 60., -60./)
  real, dimension (4)  :: limits = (/-90., -180., 90., 180./)
  real, dimension (NBINS)     :: MODIS_BINS, CLM4_BINS
  INTEGER              :: I_INDEX(10),JM, IM, NC, NT, year, day, maxcat, ThisTile, req
  type(date_time_type), dimension (YearE - YearB + 1) :: octad_time
  CHARACTER*8          :: YYYYMMDD
  CHARACTER*6          :: YYYYMM
  CHARACTER*4          :: YYYY
  CHARACTER*2          :: MM, DD,MMR, TSLICE
  INTEGER              :: i, j,k, pf, i1,i2,i3,i4, f1,f2,f3,f4, ND
  INTEGER, DIMENSION (:,:), allocatable :: veg_index, catchs_all
  integer, allocatable, dimension (:) :: ldas2bcs, Catchs, NCATCH_ALL
  real                                :: yr,mn,dy,dum,yr1,mn1,dy1, lw, up, db, mean, std, skew, &
       minv, maxv,  minv1, maxv1, minv2, maxv2, VISmean, FPARmean,  FPARstd, NIRmean, NIRstd, VISstd, r2, var1, var2
  real, dimension (:,:), allocatable  :: modis_fpar, clm4_fpar, modis_visdf, modis_nirdf
  real, dimension (:), allocatable    :: modis_cdf, clm4_cdf, data_read, data_save, modis_hist, clm4_hist
  real, dimension (:), allocatable    :: vis_this, nir_this, fpar_this, est_this
  real (kind=8), dimension(NBINS)     :: dbins,dcdf 
  real (kind=8), dimension(3)         :: modis_param, clm4_param
  real , dimension(NBINS*2 + 13 + 6)  :: tmp_real
  character*300                       :: tmpstring

  call MPI_Init(error)
  call MPI_COMM_Size(MPI_COMM_WORLD,comm_size,error)
  call MPI_COMM_Rank(MPI_COMM_WORLD,comm_rank,error)
  
  call MPI_Info_create(info, error)
  call MPI_Info_set(info, "romio_cb_read", "automatic", error)
  
  write (TSLICE ,'(i2.2)') comm_rank + 1 
  open (99,file=trim(logfil)//TSLICE, form ='formatted', action='write', status= 'unknown')

  ! STEP 1 check/create CDF params files
  
  inquire(file=trim(OUTFIL)//'.nc4', exist = file_exists)
  
  if(.not.file_exists) then
     if(comm_rank == 0) then 
        call create_CDF_ParamFile
        WRITE (99,*)'CREATED CDF PARAM FILE : ', trim(OUTFIL)//'.nc4'
     endif  
  endif

  call MPI_BARRIER( MPI_COMM_WORLD, error)

  ! READ # OF 10x10 MODIS TILES AND SMAP_TILE_IDs THAT CONTRIBUTE TO EACH MODIS TILE

  STATUS = NF_OPEN (trim(OUTFIL)//'_aux.nc4', NF_NOWRITE,NCInID) ; VERIFY_(STATUS)
  STATUS = NF_INQ_DIMID    (NCInID, 'tile10D', K); VERIFY_(STATUS)
  CALL HANDLE_ERR(STATUS, 'INQ_DIM')  
  STATUS = NF_INQ_DIMLEN (NCInID, K, NTILES); VERIFY_(STATUS)
  CALL HANDLE_ERR(STATUS, 'DIMLEN_NTILES') 
  WRITE (99,*)'NOF 10D CELLS : ', NTILES

  allocate (NCatch_all (1:NTILES))
  allocate (Catchs_all (1:16000,1:NTILES)) 

  STATUS = NF_GET_VARA_INT (NCInID, VarID(NCInID,'nSMAP'),  (/1/), (/NTILES/), NCatch_ALL); VERIFY_(STATUS)
  STATUS = NF_GET_VARA_INT (NCInID, VarID(NCInID,'SMAPID'),  (/1, 1/), (/16000,NTILES/), Catchs_ALL); VERIFY_(STATUS)
  status = NF_CLOSE (NCInID)

  call MPI_BARRIER( MPI_COMM_WORLD, error)

  if(comm_rank == 0) then    

    ! ROOT PROCESSOR OPENS FILES TO UPDATE CDF PARAMS

     STATUS = NF_OPEN (trim(OUTFIL)//'.nc4'    ,NF_WRITE,NCOutID )
     VERIFY_(STATUS)         
     STATUS = NF_OPEN (trim(OUTFIL)//'_aux.nc4',NF_WRITE,NCOutID2)
     VERIFY_(STATUS)

  endif

  ! read maxcat, LDASsa tile order, veg types and define binvals

  open (10,file=trim(BCSDIR)//'clsm/catchment.def', status='old',action='read',  &
       form='formatted')

  read (10,*) maxcat

  close (10, status ='keep')

  allocate (veg_index   (1: NPFT, 1: MAXCAT))
  allocate (ldas2bcs    (1: maxcat))
  allocate (data_read   (1: maxcat))
  allocate (data_save   (1: maxcat))
  allocate (modis_fpar  (1: maxcat, yearE - yearB + 1))
  allocate (clm4_fpar   (1: maxcat, yearE - yearB + 1))
  allocate (modis_visdf (1: maxcat, yearE - yearB))
  allocate (modis_nirdf (1: maxcat, yearE - yearB))

  clm4_fpar  = 0.
  modis_fpar = 0.
  data_save  = 0.
  veg_index  = -9999

  open (10,file =trim(EXPDIR)//'rc_out/'//trim(EXNAME)//'.ldas_tilecoord.bin',status='old',form='unformatted',convert='big_endian')
  read (10) i
  read  (10)  LDAS2BCS
  close (10, status = 'keep')
  
  open (10,file=trim(BCSDIR)//'clsm/CLM_veg_typs_fracs', status='old',action='read',  &
       form='formatted')

  do i = 1, maxcat

     read (10,*) pf, pf, i1,i2,i3,i4, f1,f2,f3,f4

     if(f1 >= f2) then 
        veg_index (i1,i) = i
     else
        veg_index (i2,i) = i
     endif

  end do

  close (10, status ='keep')

!  OCTAD_LOOP : DO NT = 1, NOCTAD   

  NT = comm_rank + 1
  ND = 8
  
  ! BEGIN READING MODIS FPAR
  
  MODIS_LOOP : DO year = YearB, YearE
     
     write (YYYY ,'(i4.4)')  year
     
     WRITE (99,*)'FPAR, VISDF and NIRDF FILES : '
     WRITE (99,*) YYYY//'//fpar.dat'
     
     open (10, file = trim(BCSDIR)//'MODIS6/'//YYYY//'//fpar.dat', form = 'unformatted', action = 'read')
     
     if(year < yearE) then
        open (11, file = trim(BCSDIR)//'MODIS6/'//YYYY//'//visdf.dat', form = 'unformatted', action = 'read')
        open (12, file = trim(BCSDIR)//'MODIS6/'//YYYY//'//nirdf.dat', form = 'unformatted', action = 'read')
        WRITE (99,*) YYYY//'//visdf.dat'  
        WRITE (99,*) YYYY//'//nirdf.dat'
        WRITE (99,*) ' '
     endif
     
     do k = 0, nt ! your processor rank
        
        read (10) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
        read (10) modis_fpar (:, year - yearB + 1)
        
        if ((k > 0) .and.(year < yearE)) then
           read (11) modis_visdf (:, year - yearB + 1)
           read (12) modis_nirdf (:, year - yearB + 1)
        endif

        IF (k == NT) WRITE (99,*) 'PROCESSING TIME SLICE : ', yr,mn,dy,yr1,mn1,dy1

     end do
     
     close (10, status = 'keep')

     if(year < yearE) then
        close (11, status = 'keep')
        close (12, status = 'keep')
     endif
     
  END DO MODIS_LOOP

  ! END READING MODIS FPAR and BEGIN READING CLM4 FPAR 
  
  WRITE (99,*) ' '
  WRITE (99,*) 'READING CLM4 FPAR: '
  
  CLM4_LOOP : DO year = YearB, YearE
     
     octad_time(year - yearB + 1)%year  = year - 1
     octad_time(year - yearB + 1)%month = 12
     octad_time(year - yearB + 1)%day   = 31
     octad_time(year - yearB + 1)%hour  = 0
     octad_time(year - yearB + 1)%min   = 0
     octad_time(year - yearB + 1)%sec   = 0
     
     do k = 1, nt
        
        if((K == NT).and.(NT == 46)) ND = 5
        
        DO day = 1,ND
           call augment_date_time( 86400, octad_time(year - yearB + 1))
           
           if(k == nt) then
              write (YYYY, '(i4.4)') octad_time(year - yearB + 1)%year
              write (MM  , '(i2.2)') octad_time(year - yearB + 1)%month
              write (DD  , '(i2.2)') octad_time(year - yearB + 1)%day
              YYYYMMDD = YYYY//MM//DD
              WRITE (99,*) trim(EXNAME)//'.ens_avg.ldas_tile_daily_out.'//YYYYMMDD//'.bin'
              open (10, file = trim(EXPDIR)//'cat/ens_avg/Y'//YYYY//'/M'//MM//'/'// &
                   trim(EXNAME)//'.ens_avg.ldas_tile_daily_out.'//YYYYMMDD//'.bin', &
                   form = 'unformatted', convert='big_endian', action = 'read')   
              
              do n = 1,3
                 read (10) data_read
                 if(n == 2) data_save = data_read
              end do
              
              clm4_fpar (:,year - yearB + 1) = clm4_fpar (:,year - yearB + 1) + data_save / (data_read + 1.e-20)/ real (ND)
              
              close (10, status = 'keep')
              
           endif
        end do
        
        ND = 8
        
     end do
     
     ! reoder to the order of BCs
     
     data_read = clm4_fpar (:,year - yearB + 1)
     
     do n = 1, maxcat
        clm4_fpar (LDAS2BCS(n),year - yearB + 1) = data_read (n) 
     end do
     
  END DO CLM4_LOOP
  
  ! Now compute CDFs
  ! ----------------
  
  ! loop through tiles
  
  allocate (modis_hist (1:MXCNT))
  allocate (clm4_hist  (1:MXCNT))
  allocate (vis_this   (1:MXCNT))
  allocate (nir_this   (1:MXCNT))
  allocate (fpar_this  (1:MXCNT))
  allocate (est_this   (1:MXCNT))
  allocate (modis_cdf  (1:NBINS))
  allocate (clm4_cdf   (1:NBINS))
  
  TILE_LOOP : DO ThisTile = 1,NTILES

     NCATCH = NCATCH_ALL (ThisTile)
     
     if(NCatch < 1) then
        WRITE (99,*) 'nSMAP problem ', NCatch, ThisTile
     endif
     
     allocate (Catchs (1: NCatch))
     !        STATUS = NF_GET_VARA_INT (NCInID, VarID(NCInID,'SMAPID'),  (/1, ThisTile/), (/NCatch,1/), Catchs)
     Catchs (1: NCatch) = Catchs_ALL (1: NCatch, ThisTile) 
     
     PFT_LOOP: DO k =1, NPFT 
        
        modis_hist = 0.
        clm4_hist  = 0.
        modis_cdf  = 0.
        clm4_cdf   = 0.
        NDATA1     = 1
        NDATA2     = 1
        N_VIS_DATA = 1
        N_NIR_DATA = 1
        MODIS_BINS = 0.
        CLM4_BINS  = 0.
        
        DO year = YearB, YearE
           DO N = 1, NCatch                    
              if((veg_index(k,Catchs(n)) > 0).and.(modis_fpar (veg_index(k,Catchs(n)), year - yearB + 1) > 0. )) then
                 modis_hist (NDATA1) = modis_fpar (veg_index(k,Catchs(n)), year - yearB + 1) 
                 NDATA1 = NDATA1 + 1
              endif
              
              if((veg_index(k,Catchs(n)) > 0).and.(clm4_fpar (veg_index(k,Catchs(n)),year - yearB + 1)  > 0. )) then                    
                 clm4_hist (NDATA2) = clm4_fpar (veg_index(k,Catchs(n)),year - yearB + 1)
                 NDATA2 = NDATA2 + 1              
              endif

              if (year < yearE) then 
                 if((veg_index(k,Catchs(n)) > 0).and.(modis_fpar (veg_index(k,Catchs(n)), year - yearB + 1) >= 0. ).and.( modis_visdf(veg_index(k,Catchs(n)), year - yearB + 1) > 0.) &
                      .and.( modis_nirdf(veg_index(k,Catchs(n)), year - yearB + 1) > 0.)) then
                    vis_this (N_VIS_DATA) = modis_visdf (veg_index(k,Catchs(n)), year - yearB + 1) 
                    fpar_this(N_VIS_DATA) = modis_fpar  (veg_index(k,Catchs(n)), year - yearB + 1)
                    nir_this (N_VIS_DATA) = modis_nirdf (veg_index(k,Catchs(n)), year - yearB + 1)  
                    N_VIS_DATA = N_VIS_DATA  + 1
                 endif
                 
!                 if((veg_index(k,Catchs(n)) > 0).and.(modis_fpar (veg_index(k,Catchs(n)), year - yearB + 1) >= 0. ).and.( modis_nirdf(veg_index(k,Catchs(n)), year - yearB + 1) > 0.)) then
!                    nir_this (N_NIR_DATA) = modis_nirdf (veg_index(k,Catchs(n)), year - yearB + 1) 
!                    N_NIR_DATA = N_NIR_DATA  + 1
!                 endif
              endif
              
              if((ndata1 > MXCNT).or.(ndata2 > MXCNT).or.(N_VIS_DATA > MXCNT)) then 
                 WRITE (99,*) 'NDATA1 or NDATA2 exceeded ',ndata1, ndata2, N_VIS_DATA, N_NIR_DATA
                 stop
              endif
              
           END DO
        END DO
        
        NDATA1 = NDATA1 - 1
        NDATA2 = NDATA2 - 1
        N_NIR_DATA = N_NIR_DATA - 1
        N_VIS_DATA = N_VIS_DATA - 1
                
        MINV1   = -9999.
        MAXV1   = -9999.
        MINV2   = -9999.
        MAXV2   = -9999.
        
        ! curve fitting
        
        modis_param = -9999.
        clm4_param  = -9999.
        FPARmean    = -9999.
        FPARstd     = -9999.

        if((NDATA1 > 10).and.(NDATA2 > 10)) then
              
           WRITE (99,*) '# of SMAP DATA CELLS ',ThisTile, K, NDATA1,NDATA2, N_NIR_DATA,N_VIS_DATA
           
           maxv = MAXVAL ((/MAXVAL(modis_hist (1: NDATA1)),MAXVAL(clm4_hist (1: NDATA2))/)) 
           minv = MINVAL ((/MINVAL(modis_hist (1: NDATA1)),MINVAL(clm4_hist (1: NDATA2))/)) 
           
           if (maxv > minv) then
              
              modis_hist (1: NDATA1) = (modis_hist (1: NDATA1) - minv)/maxv
              call prob_den_func (ndata1, modis_hist (1: NDATA1) , modis_cdf, MODIS_bins, lwval= 0.,upval = 1.)    
              
              dbins = MODIS_bins        
              dcdf  = modis_cdf
              modis_param (1) = .5
              modis_param (2) = .5
              modis_param (3) = 0.9           
              call optimiz (nbins,dbins,dcdf,modis_param)
              
              FPARmean =  SUM (clm4_hist (1:NDATA2)) / real (NDATA2)
              var1 = 0.
              
              do i = 1,NDATA2
                 var1 = var1 + (clm4_hist(i) - FPARmean)*(clm4_hist(i) - FPARmean)
              end do
              
              FPARstd = sqrt (var1/real(NDATA2 - 1))
              
              clm4_hist (1: NDATA2) = (clm4_hist (1: NDATA2) - minv)/maxv
              call prob_den_func (ndata2, clm4_hist (1: NDATA2) , clm4_cdf, CLM4_bins, lwval= 0.,upval = 1.)    
              
              dbins = CLM4_bins       
              dcdf  = clm4_cdf
              clm4_param (1) = .5
              clm4_param (2) = .5
              clm4_param (3) = 0.9           
              call optimiz (NBINS,dbins,dcdf,clm4_param)
           endif
        endif
     
        ! albedo parameters

        VISmean  = -9999.
        NIRmean  = -9999.
        VISstd   = -9999.
        NIRstd   = -9999.
        
        if (N_VIS_DATA > 12)  then 
           NIRmean = SUM (nir_this (1:N_VIS_DATA)) / real (N_VIS_DATA)
           VISmean = SUM (vis_this (1:N_VIS_DATA)) / real (N_VIS_DATA)

           var1 = 0.
           var2 = 0.
           
           do i = 1,N_VIS_DATA
              var1 = var1 + (vis_this(i) - VISmean)*(vis_this(i) - VISmean)
              var2 = var2 + (nir_this(i) - NIRmean)*(nir_this(i) - NIRmean)
           end do

           VISstd = sqrt (var1/real(N_VIS_DATA - 1))
           NIRstd = sqrt (var2/real(N_VIS_DATA - 1))  

        endif
        
        do i = 1, comm_size
           
           tmp_real = -9999.
           
           if((I == 1).and.(comm_rank == 0)) then
              
              if(put_aux) then 
                 STATUS = NF_PUT_VARA_REAL(NCOutID2,VARID(NCOutID2,'CDF'    ),(/1,ThisTile,k,NT,1/), (/NBINS,1,1,1,1/), modis_cdf );VERIFY_(STATUS)
                 STATUS = NF_PUT_VARA_REAL(NCOutID2,VARID(NCOutID2,'CDF'    ),(/1,ThisTile,k,NT,2/), (/NBINS,1,1,1,1/), clm4_cdf  );VERIFY_(STATUS)
              endif
              
              STATUS = NF_PUT_VARA_REAL(NCOutID, VARID(NCOutID ,'Kappa' ),(/ThisTile,k,NT,1/), (/1,1,1,1/), REAL (modis_param(3))); VERIFY_(STATUS)
              STATUS = NF_PUT_VARA_REAL(NCOutID, VARID(NCOutID ,'Lambda'),(/ThisTile,k,NT,1/), (/1,1,1,1/), REAL (modis_param(1))); VERIFY_(STATUS)
              STATUS = NF_PUT_VARA_REAL(NCOutID, VARID(NCOutID ,'Mu'    ),(/ThisTile,k,NT,1/), (/1,1,1,1/), REAL (modis_param(2))) ; VERIFY_(STATUS)            
              STATUS = NF_PUT_VARA_REAL(NCOutID, VARID(NCOutID ,'Kappa' ),(/ThisTile,k,NT,2/), (/1,1,1,1/), REAL (clm4_param(3))); VERIFY_(STATUS)
              STATUS = NF_PUT_VARA_REAL(NCOutID, VARID(NCOutID ,'Lambda'),(/ThisTile,k,NT,2/), (/1,1,1,1/), REAL (clm4_param(1))); VERIFY_(STATUS)
              STATUS = NF_PUT_VARA_REAL(NCOutID, VARID(NCOutID ,'Mu'    ),(/ThisTile,k,NT,2/), (/1,1,1,1/), REAL (clm4_param(2))) ; VERIFY_(STATUS)             
              STATUS = NF_PUT_VARA_REAL(NCOutID, VARID(NCOutID ,'MinVal'),(/ThisTile,k,NT/)  , (/1,1,1/), MINV); VERIFY_(STATUS);VERIFY_(STATUS)
              STATUS = NF_PUT_VARA_REAL(NCOutID, VARID(NCOutID ,'MaxVal'),(/ThisTile,k,NT/)  , (/1,1,1/), MAXV); VERIFY_(STATUS);VERIFY_(STATUS)
              STATUS = NF_PUT_VARA_REAL(NCOutID, VARID(NCOutID ,'MODISVISmean' ),(/ThisTile,k,NT/), (/1,1,1/),VISmean ); VERIFY_(STATUS);VERIFY_(STATUS)
              STATUS = NF_PUT_VARA_REAL(NCOutID, VARID(NCOutID ,'MODISNIRmean' ),(/ThisTile,k,NT/), (/1,1,1/),NIRmean ); VERIFY_(STATUS);VERIFY_(STATUS)
              STATUS = NF_PUT_VARA_REAL(NCOutID, VARID(NCOutID ,'MODISVISstd'  ),(/ThisTile,k,NT/), (/1,1,1/),VISstd  ); VERIFY_(STATUS);VERIFY_(STATUS)
              STATUS = NF_PUT_VARA_REAL(NCOutID, VARID(NCOutID ,'MODISNIRstd'  ),(/ThisTile,k,NT/), (/1,1,1/),NIRstd  ); VERIFY_(STATUS);VERIFY_(STATUS)
              STATUS = NF_PUT_VARA_REAL(NCOutID, VARID(NCOutID ,'MODELFPARmean'),(/ThisTile,k,NT/), (/1,1,1/),FPARmean); VERIFY_(STATUS);VERIFY_(STATUS)
              STATUS = NF_PUT_VARA_REAL(NCOutID, VARID(NCOutID ,'MODELFPARstd' ),(/ThisTile,k,NT/), (/1,1,1/),FPARstd ); VERIFY_(STATUS);VERIFY_(STATUS)
              WRITE (99,*) 'Writing Out OCTAD', ThisTile,k,NT

           else if (I > 1) then
              
              if(I-1 == comm_rank) then
                 !                 print *, comm_rank, 'sending', ThisTile, k,nt
                 tmp_real (1) = real (ThisTile)
                 tmp_real (2) = real(k)
                 tmp_real (3) = real(NT)
                 tmp_real (4) = REAL (modis_param(3))
                 tmp_real (5) = REAL (modis_param(1))
                 tmp_real (6) = REAL (modis_param(2))
                 tmp_real (7) = REAL (clm4_param(3))
                 tmp_real (8) = REAL (clm4_param(1))
                 tmp_real (9) = REAL (clm4_param(2))
                 tmp_real(10) = MINV 
                 tmp_real(11) = MINV 
                 tmp_real(12) = MAXV 
                 tmp_real(13) = MAXV 
                 tmp_real(14) = VISmean 
                 tmp_real(15) = NIRmean 
                 tmp_real(16) = VISstd  
                 tmp_real(17) = NIRstd  
                 tmp_real(18) = FPARmean
                 tmp_real(19) = FPARstd 
                 
                 NC = 19
                 n = NC + NBINS
                 tmp_real (NC+1: N) = modis_cdf(:)
                 NC = n
                 n = NC + NBINS
                 tmp_real (NC+1: N) = clm4_cdf(:)
                 NC = n
                 
                 call MPI_ISend(tmp_real ,2*NBINS + 19,MPI_real,0,999,MPI_COMM_WORLD,req,status)
                 call MPI_WAIT (req,MPI_STATUS_IGNORE,status)

              else if (comm_rank == 0) then

                 call MPI_RECV(tmp_real,2*NBINS + 19,MPI_real,I-1,999,MPI_COMM_WORLD,MPI_STATUS_IGNORE,status)

                 IM = NINT (tmp_real (1))
                 JM = NINT (tmp_real (2))
                 I1 = NINT (tmp_real (3))
                                 
                 STATUS = NF_PUT_VARA_REAL(NCOutID, VARID(NCOutID,'Kappa' ),(/IM,JM,I1,1/), (/1,1,1,1/), tmp_real (4)); VERIFY_(STATUS);VERIFY_(STATUS)
                 STATUS = NF_PUT_VARA_REAL(NCOutID, VARID(NCOutID,'Lambda'),(/IM,JM,I1,1/), (/1,1,1,1/), tmp_real (5)); VERIFY_(STATUS);VERIFY_(STATUS)
                 STATUS = NF_PUT_VARA_REAL(NCOutID, VARID(NCOutID,'Mu'    ),(/IM,JM,I1,1/), (/1,1,1,1/), tmp_real (6)) ; VERIFY_(STATUS) ;VERIFY_(STATUS)   
                
                 STATUS = NF_PUT_VARA_REAL(NCOutID, VARID(NCOutID,'Kappa' ),(/IM,JM,I1,2/), (/1,1,1,1/), tmp_real (7)); VERIFY_(STATUS);VERIFY_(STATUS)
                 STATUS = NF_PUT_VARA_REAL(NCOutID, VARID(NCOutID,'Lambda'),(/IM,JM,I1,2/), (/1,1,1,1/), tmp_real (8)); VERIFY_(STATUS);VERIFY_(STATUS)
                 STATUS = NF_PUT_VARA_REAL(NCOutID, VARID(NCOutID,'Mu'    ),(/IM,JM,I1,2/), (/1,1,1,1/), tmp_real (9)); VERIFY_(STATUS);VERIFY_(STATUS)     
                
                 STATUS = NF_PUT_VARA_REAL(NCOutID, VARID(NCOutID,'MinVal'),(/IM,JM,I1/)  , (/1,1,1/)  , tmp_real(10)) ; VERIFY_(STATUS);VERIFY_(STATUS)    
                 STATUS = NF_PUT_VARA_REAL(NCOutID, VARID(NCOutID,'MaxVal'),(/IM,JM,I1/)  , (/1,1,1/)  , tmp_real(12)) ; VERIFY_(STATUS) ;VERIFY_(STATUS)   
                
                 STATUS = NF_PUT_VARA_REAL(NCOutID, VARID(NCOutID ,'MODISVISmean' ),(/IM,JM,I1/), (/1,1,1/), tmp_real(14)); VERIFY_(STATUS);VERIFY_(STATUS)
                 STATUS = NF_PUT_VARA_REAL(NCOutID, VARID(NCOutID ,'MODISNIRmean' ),(/IM,JM,I1/), (/1,1,1/), tmp_real(15)); VERIFY_(STATUS);VERIFY_(STATUS)
                 STATUS = NF_PUT_VARA_REAL(NCOutID, VARID(NCOutID ,'MODISVISstd'  ),(/IM,JM,I1/), (/1,1,1/), tmp_real(16)); VERIFY_(STATUS);VERIFY_(STATUS)
                 STATUS = NF_PUT_VARA_REAL(NCOutID, VARID(NCOutID ,'MODISNIRstd'  ),(/IM,JM,I1/), (/1,1,1/), tmp_real(17)); VERIFY_(STATUS);VERIFY_(STATUS)
                 STATUS = NF_PUT_VARA_REAL(NCOutID, VARID(NCOutID ,'MODELFPARmean'),(/IM,JM,I1/), (/1,1,1/), tmp_real(18)); VERIFY_(STATUS);VERIFY_(STATUS)
                 STATUS = NF_PUT_VARA_REAL(NCOutID, VARID(NCOutID ,'MODELFPARstd' ),(/IM,JM,I1/), (/1,1,1/), tmp_real(19)); VERIFY_(STATUS);VERIFY_(STATUS)
                 
                 NC = 19
                 n = NC + NBINS
                 if(put_aux) STATUS = NF_PUT_VARA_REAL(NCOutID2,VARID(NCOutID2,'CDF'     ),(/1,IM,JM,I1,1/), (/NBINS,1,1,1,1/), tmp_real (nc+1: N)); VERIFY_(STATUS)
                 NC = n
                 n = NC + NBINS
                 if(put_aux) STATUS = NF_PUT_VARA_REAL(NCOutID2,VARID(NCOutID2,'CDF'     ),(/1,IM,JM,I1,2/), (/NBINS,1,1,1,1/), tmp_real (nc+1: N)); VERIFY_(STATUS)
                 NC = n
!                 n = NC + NBINS
!                 if(put_aux) STATUS = NF_PUT_VARA_REAL(NCOutID2,VARID(NCOutID2,'BINS'    ),(/1,IM,JM,I1/)  , (/NBINS,1,1,1/)  , tmp_real (nc+1: N)); VERIFY_(STATUS)
!                 NC = n
!                 n = NC + NBINS
!                 if(put_aux) STATUS = NF_PUT_VARA_REAL(NCOutID2,VARID(NCOutID2,'BINS'    ),(/1,IM,JM,I1,2/), (/NBINS,1,1,1,1/), tmp_real (nc+1: N)); VERIFY_(STATUS)

                 WRITE (99,*) ' RECEIVED Writing Out OCTAD', I-1,IM,JM,I1
              endif
           endif
        end do
           
     END DO PFT_LOOP

     deallocate (Catchs)

  END DO TILE_LOOP

  ! END DO OCTAD_LOOP
  if(comm_rank == 0) then 
     status = NF_CLOSE (NCOutID )
     status = NF_CLOSE (NCOutID2)
  endif
 close (99, status = 'keep') 

 call MPI_BARRIER( MPI_COMM_WORLD, error)
 call MPI_Finalize(STATUS)

STOP

CONTAINS 

!________________________________________________________________________


SUBROUTINE create_CDF_ParamFile

  implicit none

  integer :: NCFOutID, NCFOutID2, status, pid, tid, did, lid, bid, vid,CID
  integer :: i,j, k, n, maxcat,ii,jj, tile_count, nplus, nb, NB_max
  integer :: nc_rst = 43200, nr_rst = 21600, DIJ, MXT = 16000
  integer :: PID2, TID2, DID2, LID2,BID2, CID2
  real    :: dxy, lw, up, db
  real,    allocatable,           dimension (:) :: abins
  integer, allocatable, target, dimension (:,:) :: tile_id
  integer, pointer, dimension (:,:)             :: tile_id_box
  integer, allocatable, dimension (:)           :: tile_id_vec, bcs2ldas
  integer, allocatable, dimension (:) :: density, loc_int
  integer, allocatable, dimension (:) :: loc_val
  logical, allocatable, dimension (:) :: unq_mask
  character (22)             :: time_stamp
  integer, dimension(8)      :: date_time_values


  status = NF_CREATE (trim(OUTFIL)//'.nc4'    , NF_NETCDF4, NCFOutID );VERIFY_(STATUS)
  status = NF_CREATE (trim(OUTFIL)//'_aux.nc4', NF_NETCDF4, NCFOutID2);VERIFY_(STATUS)        
  ! Define Dimensions

  status = NF_DEF_DIM(NCFOutID, 'pft'        , NPFT  , PID);VERIFY_(STATUS)
  status = NF_DEF_DIM(NCFOutID, 'octad'      , NOCTAD, TID);VERIFY_(STATUS) 
  status = NF_DEF_DIM(NCFOutID, 'data'       , NSETS , DID);VERIFY_(STATUS)
  status = NF_DEF_DIM(NCFOutID, 'tile10D',NF_UNLIMITED,LID);VERIFY_(STATUS)


  status = NF_DEF_DIM(NCFOutID2, 'pft'        , NPFT  , PID2);VERIFY_(STATUS)
  status = NF_DEF_DIM(NCFOutID2, 'octad'      , NOCTAD, TID2);VERIFY_(STATUS)
  status = NF_DEF_DIM(NCFOutID2, 'data'       , NSETS , DID2);VERIFY_(STATUS)
  status = NF_DEF_DIM(NCFOutID2, 'tile10D',NF_UNLIMITED,LID2);VERIFY_(STATUS)
  status = NF_DEF_DIM(NCFOutID2, 'bin'        ,  nbins, BID2);VERIFY_(STATUS)
  status = NF_DEF_DIM(NCFOutID2, 'nCELLS'     ,    MXT, CID2);VERIFY_(STATUS)

  ! Define variables

  status = NF_DEF_VAR(NCFOutID, 'lon'    , NF_FLOAT ,1 ,(/LID/), vid);VERIFY_(STATUS)
  status = NF_DEF_VAR(NCFOutID, 'lat'    , NF_FLOAT ,1 ,(/LID/), vid);VERIFY_(STATUS)
  status = NF_DEF_VAR(NCFOutID, 'Kappa'  , NF_FLOAT ,4 ,(/LID,PID,TID,DID/), vid);VERIFY_(STATUS)
  status = NF_DEF_VAR(NCFOutID, 'Lambda' , NF_FLOAT ,4 ,(/LID,PID,TID,DID/), vid);VERIFY_(STATUS)
  status = NF_DEF_VAR(NCFOutID, 'Mu'     , NF_FLOAT ,4 ,(/LID,PID,TID,DID/), vid);VERIFY_(STATUS)
  status = NF_DEF_VAR(NCFOutID, 'MinVal' , NF_FLOAT ,3 ,(/LID,PID,TID/), vid);VERIFY_(STATUS)
  status = NF_DEF_VAR(NCFOutID, 'MaxVal' , NF_FLOAT ,3 ,(/LID,PID,TID/), vid);VERIFY_(STATUS)
  status = NF_DEF_VAR(NCFOutID, 'MODISVISmean', NF_FLOAT ,3 ,(/LID,PID,TID/), vid);VERIFY_(STATUS)
  status = NF_DEF_VAR(NCFOutID, 'MODISNIRmean', NF_FLOAT ,3 ,(/LID,PID,TID/), vid);VERIFY_(STATUS)
  status = NF_DEF_VAR(NCFOutID, 'MODISVISstd' , NF_FLOAT ,3 ,(/LID,PID,TID/), vid);VERIFY_(STATUS)
  status = NF_DEF_VAR(NCFOutID, 'MODISNIRstd' , NF_FLOAT ,3 ,(/LID,PID,TID/), vid);VERIFY_(STATUS)
  status = NF_DEF_VAR(NCFOutID, 'MODELFPARmean',NF_FLOAT ,3 ,(/LID,PID,TID/), vid);VERIFY_(STATUS)
  status = NF_DEF_VAR(NCFOutID, 'MODELFPARstd' ,NF_FLOAT ,3 ,(/LID,PID,TID/), vid);VERIFY_(STATUS)
  status = NF_DEF_VAR(NCFOutID2, 'lon'   , NF_FLOAT ,1 ,(/LID2/), vid);VERIFY_(STATUS)
  status = NF_DEF_VAR(NCFOutID2, 'lat'   , NF_FLOAT ,1 ,(/LID2/), vid);VERIFY_(STATUS)
  status = NF_DEF_VAR(NCFOutID2, 'nSMAP' , NF_SHORT ,1 ,(/LID2/), vid);VERIFY_(STATUS)
  status = NF_DEF_VAR(NCFOutID2, 'BINS'  , NF_FLOAT ,1 ,(/BID2/), vid);VERIFY_(STATUS)
  status = NF_DEF_VAR(NCFOutID2, 'SMAPID', NF_INT   ,2 ,(/CID2, LID/), vid);VERIFY_(STATUS)
  status = NF_DEF_VAR(NCFOutID2, 'CDF'   , NF_FLOAT ,5 ,(/BID2,LID2,PID2,TID2,DID2/), vid);VERIFY_(STATUS)

!  Global attributes
!
  call date_and_time(VALUES=date_time_values)
          
  write (time_stamp,'(i4.4,a1,i2.2,a1,i2.2,1x,a2,1x,i2.2,a1,i2.2,a1,i2.2)')      &
       date_time_values(1),'-',date_time_values(2),'-',date_time_values(3),'at', &
       date_time_values(5),':',date_time_values(6),':',date_time_values(7)
  
  status = NF_PUT_ATT_TEXT(NCFOutID, NF_GLOBAL, 'CreatedBy', LEN_TRIM("Sarith Mahanama"),  &
         trim("Sarith Mahanama"))
  status = NF_PUT_ATT_TEXT(NCFOutID, NF_GLOBAL, 'Date'   , LEN_TRIM(time_stamp),trim(time_stamp))

  status = NF_PUT_ATT_TEXT(NCFOutID2, NF_GLOBAL, 'CreatedBy', LEN_TRIM("Sarith Mahanama"),  &
         trim("Sarith Mahanama"))
  status = NF_PUT_ATT_TEXT(NCFOutID2, NF_GLOBAL, 'Date'   , LEN_TRIM(time_stamp),trim(time_stamp))

  status = NF_ENDDEF(NCFOutID )
  status = NF_ENDDEF(NCFOutID2)

 allocate (abins (1:nbins)) 
 
  lw = 0.
  up = 1.
  db = (up - lw)/real(nbins)
  do n = 1, nbins
     abins (n) = lw + real(n)*db - db/2.  
  end do

  STATUS = NF_PUT_VARA_REAL(NCFOutID2,VARID(NCFOutID2, 'BINS'),(/1/)  , (/NBINS/)  , abins); VERIFY_(STATUS)
 
  ! create TILESIZE x TILESIZE tiles at TILEINT 

  dij = nint (TILESIZE *  nc_rst/360)
  dxy = 360./nc_rst

  ! read maxcat from catchment.def

  open (10,file=trim(BCSDIR)//'clsm/catchment.def', status='old',action='read',  &
       form='formatted')

  read (10,*) maxcat

  close (10, status ='keep')

  ! read tilecoord for tile order in LDASsa

  open (10,file =trim(EXPDIR)//'rc_out/'//trim(EXNAME)//'.ldas_tilecoord.bin',status='old',form='unformatted',convert='big_endian')
  read (10) i
  if (i /= maxcat) then 
     print *,'NTILES BCs/LDASsa mismatch:', i,maxcat
     stop
  endif
  
  allocate (tile_id_vec (1: maxcat))
  allocate (bcs2ldas    (1: maxcat))
  read  (10) tile_id_vec
  close (10, status = 'keep')

  ! indexing to the LDASsa order

  do i = 1, maxcat
      BCS2LDAS(tile_id_vec(i)) = i
  end do

  ! read tile_id raster and index according to the order of LDASsa

  open (10,file=trim(BCSDIR)//'rst/'//trim(GFILE)//'.rst',status='old',action='read',  &
       form='unformatted',convert='little_endian')

  ALLOCATE (tile_id     (1:nc_rst,1:nr_rst))

  DO j = 1, nr_rst
     read (10) tile_id (:, j)
  END DO

  close (10, status ='keep')

  deallocate (tile_id_vec)

  ! find catchment-tiles that contribute to each 10x10 tile

  tile_count = 0
  Nb_max = 0

  DO J = FLOOR (limits(1) + TILESIZE/2), CEILING (limits (3) - TILESIZE/2), NINT(TILEINT)
     DO I = FLOOR (limits(2) + TILESIZE/2), CEILING (limits (4) - TILESIZE/2), NINT(TILEINT)
        if (associated (tile_id_box)) NULLIFY (tile_id_box)
        jj = (j +  90)*nc_rst/360 -  dij/2
        ii = (i + 180)*nc_rst/360 -  dij/2
        tile_id_box => tile_id (ii + 1 : ii + dij, jj +1 : jj + dij)

        NPLUS = count(tile_id_box >= 1 .and. tile_id_box <= maxcat)

        if(NPLUS > 0)  then

           allocate (loc_int (1:NPLUS))
           allocate (unq_mask(1:NPLUS))
           loc_int = pack(tile_id_box,mask = (tile_id_box >= 1 .and. tile_id_box <= maxcat))
           call MAPL_Sort (loc_int)
           unq_mask = .true.
           do n = 2,NPLUS 
              unq_mask(n) = .not.(loc_int(n) == loc_int(n-1))
           end do
           NB = count(unq_mask)
           tile_count = tile_count + 1
           allocate(loc_val (1:NB))
           loc_val = 1.*pack(loc_int,mask =unq_mask)
           
           IF(NB_MAX < NB) NB_MAX = NB

           if(NB > MXT) then
              print *, 'NB EXCEEDED MXT', NB, tile_count
              stop
           endif

           STATUS = NF_PUT_VARA_INT (NCFOutID2,VARID(NCFOutID2,'nSMAP' ),(/tile_count/), (/1/), NB);VERIFY_(STATUS)
           STATUS = NF_PUT_VARA_REAL(NCFOutID2,VARID(NCFOutID2,'lat'   ),(/tile_count/), (/1/), real(j));VERIFY_(STATUS)
           STATUS = NF_PUT_VARA_REAL(NCFOutID2,VARID(NCFOutID2,'lon'   ),(/tile_count/), (/1/), real(i));VERIFY_(STATUS)
           STATUS = NF_PUT_VARA_INT (NCFOutID2,VARID(NCFOutID2,'SMAPID'),(/1,tile_count/),(/NB, 1/), loc_val);VERIFY_(STATUS)
           STATUS = NF_PUT_VARA_REAL(NCFOutID ,VARID(NCFOutID, 'lat'   ),(/tile_count/), (/1/), real(j));VERIFY_(STATUS)
           STATUS = NF_PUT_VARA_REAL(NCFOutID ,VARID(NCFOutID, 'lon'   ),(/tile_count/), (/1/), real(i));VERIFY_(STATUS)

!           do k = 1,NBINS
!              print *,k, loc_val(k), BCS2LDAS(loc_val(k))
!              STATUS = NF_PUT_VARA_INT (NCFOutID,VARID(NCFOutID,'SMAPID' ),(/k,tile_count/), (/1, 1/), BCS2LDAS(loc_val(k)));VERIFY_(STATUS)
!              print *, k,nbins,BCS2LDAS(loc_val(k))
!           end do
           deallocate (loc_val,loc_int,unq_mask)

        endif
     END DO
  END DO
 
  PRINT *, 'NB_MAX :', NB_MAX

  status = NF_CLOSE (NCFOutID )
  status = NF_CLOSE (NCFOutID2)

  deallocate (abins)

END SUBROUTINE create_CDF_ParamFile

! ----------------------------------------------------------------------

integer function VarID (NCFID, VNAME) 
  
  integer, intent (in)      :: NCFID
  character(*), intent (in) :: VNAME
  integer                   :: status
  
  STATUS = NF_INQ_VARID (NCFID, trim(VNAME) ,VarID)
  IF (STATUS .NE. NF_NOERR) &
       CALL HANDLE_ERR(STATUS, trim(VNAME))  
  
end function VarID

! -----------------------------------------------------------------------

SUBROUTINE HANDLE_ERR(STATUS, Line)
  
  INTEGER,      INTENT (IN) :: STATUS
  CHARACTER(*), INTENT (IN) :: Line
  
  IF (STATUS .NE. NF_NOERR) THEN
     PRINT *, trim(Line),': ',STATUS, NF_STRERROR(STATUS)
     STOP 'Stopped'
  ENDIF
  
END SUBROUTINE HANDLE_ERR

! ----------------------------------------------------------------------

SUBROUTINE HISTOGRAM (NLENS, NBINS, density, loc_val, x, BIN)
  
  implicit none
  
  integer, intent (in) :: NBINS, NLENS
  real,    intent (in) :: x (NLENS)
  integer, intent (out):: density (NBINS)
  real,    intent (inout) :: loc_val (NBINS)
  real,    intent (in), optional :: bin
  real :: xdum(NLENS), xl, xu, min_value
  integer :: n
  
  if(present (bin)) min_value =  real(floor(minval(x)))
  
  DO N = 1, NBINS
     if(present (bin)) then
        xl = (N - 1)*BIN + min_value
        loc_val (n) = xl
        xu = xl + bin
        XDUM = 0.
        where((x >= xl).and.(x < xu))XDUM = 1
     else
        XDUM = 0.
        where(x == loc_val (n)) XDUM = 1
     endif
     density(n) = int(sum(XDUM))       
  END DO
  
END SUBROUTINE HISTOGRAM

!----------------------------------------------------------------------

END PROGRAM comp_FPAR_CDF
