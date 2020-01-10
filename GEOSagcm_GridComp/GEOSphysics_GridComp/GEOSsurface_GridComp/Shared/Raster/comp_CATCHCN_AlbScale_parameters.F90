#define VERIFY_(A)   IF(A/=0)THEN;PRINT *,'ERROR AT LINE ', __LINE__;STOP;ENDIF
#define ASSERT_(A)   if(.not.A)then;print *,'Error:',__FILE__,__LINE__;stop;endif

MODULE comp_CATCHCN_AlbScale_parameters

  use date_time_util,                   ONLY:     &
       date_time_type, augment_date_time

  implicit none
  INCLUDE 'netcdf.inc'

  private 

  public :: albedo4catchcn 

  character*400, PARAMETER ::                                                                          &
                 InBCSDIR = '/gpfsm/dnb02/smahanam/MERRA3/FPAR/SMAP_EASEv2_M09/',                      &
                 EXPDIR   = '/archive/u/smahanam/FPAR-ALB/e0004s_wet2/output/SMAP_EASEv2_M09_GLOBAL/', &
                 EXNAME   = 'e0004s_wet2',                                                             &
                 InGFILE  = 'SMAP_EASEv2_M09_3856x1624'
                                               
  ! character*400 :: GFILE = 'til/CF0180x6C_TM0720xTM0410-Pfafstetter.til'
  real, parameter                         :: MAPL_PI     =  3.14159265358979323846d0
  integer, parameter :: yearB = 2001, yearE = 2015, InNTILES = 1684725, NOCTAD = 46
  integer, parameter :: yearB1= 2002
  contains
    
    SUBROUTINE albedo4catchcn (gfile)
      
      implicit none
      character (*), intent (in)          :: gfile
      integer                             :: NTILES
      integer, dimension (:), allocatable :: id_loc
            
      call preprocess_m09
      open  (10, file = 'clsm/catchment.def', form = 'formatted', status= 'old', action = 'read')
      read  (10, *) NTILES
      close (10, status = 'keep')

      allocate (id_loc (1: NTILES))
  
      call get_id_loc (NTILES, trim (GFILE)//'.til',  id_loc) 
      call regrid_alb (NTILES, id_loc)

    end SUBROUTINE albedo4catchcn

   ! -------------------------------------------------------------------------------

    SUBROUTINE get_id_loc (NT, gfile, id_loc)

      implicit none

      integer, intent (in)                    :: NT
      integer, dimension (NT), intent (inout) :: id_loc
      character(*), intent (in)               :: gfile
      integer                                 :: n, i, nplus, t_count
      real, dimension (:), allocatable        :: lon, lat, m09_lon, m09_lat, tid_m09
      integer, allocatable, dimension (:)     :: sub_tid
      real   , allocatable, dimension (:)     :: sub_lon, sub_lat, rev_dist
      real                                    :: dw, dx, dy, min_lon, max_lon, min_lat, max_lat
      logical                                 :: tile_found
      logical, allocatable, dimension(:)      :: mask
      integer, allocatable :: low_ind(:), upp_ind(:)

! --------- VARIABLES FOR *OPENMP* PARALLEL ENVIRONMENT ------------
!
! NOTE: "!$" is for conditional compilation
!
logical :: running_omp = .false.
!
!$ integer :: omp_get_thread_num, omp_get_num_threads
!
integer :: n_threads=1
!
!
! ----------- OpenMP PARALLEL ENVIRONMENT ----------------------------
!
! FIND OUT WHETHER -omp FLAG HAS BEEN SET DURING COMPILATION
!
!$ running_omp = .true.         ! conditional compilation
!
! ECHO BASIC OMP VARIABLES
!
!$OMP PARALLEL DEFAULT(NONE) SHARED(running_omp,n_threads) 
!
!$OMP SINGLE
!
!$ n_threads = omp_get_num_threads()
!
!$ write (*,*) 'running_omp = ', running_omp
!$ write (*,*)
!$ write (*,*) 'parallel OpenMP with ', n_threads, 'threads'
!$ write (*,*)
!$OMP ENDSINGLE
!
!$OMP CRITICAL
!$ write (*,*) 'thread ', omp_get_thread_num(), ' alive'
!$OMP ENDCRITICAL
!
!$OMP BARRIER
!
!$OMP ENDPARALLEL
    
      allocate (lon      (1: NT))
      allocate (lat      (1: NT))
      allocate (m09_lon  (1: InNTILES))
      allocate (m09_lat  (1: InNTILES))
      allocate (tid_m09  (1: InNTILES))
      
      call ReadCNTilFile (trim(InBCSDIR)//trim(InGFILE)//'.til', InNTILES, m09_lon, m09_lat)
      call ReadCNTilFile (trim(GFILE), NT, lon, lat) 

      Id_loc = -9999
      do n = 1,  InNTILES
         tid_m09(n) = n 
      end do

      ! Domain decomposition
      ! --------------------
      
      allocate(low_ind(n_threads))
      allocate(upp_ind(n_threads))
      low_ind(1)         = 1
      upp_ind(n_threads) = nt
      
      if (running_omp)  then
         do i=1,n_threads-1  
            upp_ind(i)   = low_ind(i) + (NT/n_threads) - 1 
            low_ind(i+1) = upp_ind(i) + 1
            !           print *,i,low_ind(i),upp_ind(i)
         end do
         !        print *,i,low_ind(i),upp_ind(i)
      end if

!$OMP PARALLELDO DEFAULT(NONE)                          &
!$OMP SHARED( n_threads, low_ind, upp_ind, Id_loc,      &
!$OMP         lon, lat, m09_lon, m09_lat, tid_m09)      &
!$OMP PRIVATE(n,i,t_count,min_lon, max_lon, min_lat, max_lat,   &
!$OMP         sub_tid, nplus, sub_lon, sub_lat, rev_dist, dw, &
!$OMP         tile_found, mask)

      DO t_count = 1,n_threads

      allocate (mask     (1: InNTILES))      

      OUT_TILES : do n = low_ind(t_count),upp_ind(t_count)
         if(MOD(n,10000) == 0) print *,'In ID_LOC', t_count,n
         dw = 0.25

         ZOOMOUT : do  

            tile_found = .false. 
            
            ! Min/Max lon/lat of the working window
            ! -------------------------------------
            
            min_lon = MAX(lon (n) - dw, -180.)
            max_lon = MIN(lon (n) + dw,  180.)
            min_lat = MAX(lat (n) - dw,  -90.)
            max_lat = MIN(lat (n) + dw,   90.) 

            mask = .false.
            mask =  ((m09_lat >= min_lat .and. m09_lat <= max_lat).and.(m09_lon >= min_lon .and. m09_lon <= max_lon))
            nplus =  count(mask = mask)

            if(nplus < 0) then
               dw = dw + 0.5
               CYCLE
            endif

            allocate (sub_tid (1:nplus))
            allocate (sub_lon (1:nplus))
            allocate (sub_lat (1:nplus))
            allocate (rev_dist(1:nplus))
            
            sub_tid = PACK (tid_m09  , mask= mask) 
            sub_lon = PACK (m09_lon  , mask= mask)
            sub_lat = PACK (m09_lat  , mask= mask)
            
            ! compute distance from the tile
            
            sub_lat = sub_lat * MAPL_PI/180.
            sub_lon = sub_lon * MAPL_PI/180.
            
            SEEK : if(Id_loc(n) < 0) then
               
               rev_dist  = 1.e20
               
               do i = 1,nplus
                  
                  rev_dist(i) = haversine(to_radian(lat(n)), to_radian(lon(n)), &
                       sub_lat(i), sub_lon(i))
                  
               end do
               
               FOUND : if(minval (rev_dist) < 1.e19) then               
                  Id_loc(n) = sub_tid(minloc(rev_dist,1)) 
                  tile_found = .true.   

                  if(Id_loc(n) ==0) then
                     print *, rev_dist
                     print *,  sub_tid
                     print *, minval (rev_dist)
                     print *, minloc(rev_dist,1)
                     stop
                  endif
               endif FOUND
               
            endif SEEK
            
            deallocate (sub_tid, sub_lon, sub_lat, rev_dist)
            
            if(tile_found) GO TO 100
            
            ! if not increase the window size
            dw = dw + 0.25
            
         end do ZOOMOUT
         
100      continue
    
      END do OUT_TILES

      deallocate (mask)

   end DO ! PARALLEL
   
!$OMP ENDPARALLELDO   
 
    END SUBROUTINE get_id_loc

    ! *****************************************************************************

    function to_radian(degree) result(rad)
      
      ! degrees to radians
      real,intent(in) :: degree
      real :: rad
      
      rad = degree*MAPL_PI/180.
      
    end function to_radian

    ! -------------------------------------------------------------------------------------------------

    SUBROUTINE regrid_alb (NTILES, id_loc)

      implicit none
      integer, intent (in) :: NTILES
      integer, dimension (NTILES), intent (in) :: id_loc
      character*10 :: string
      integer  :: STATUS, ncid, NCOutID, t, time_slice, time_slice_next, yr, mn, dd, yr1, mn1, dd1, n_tslices
      character (len=4), dimension (:), allocatable :: MMDD, MMDD_next
      real, allocatable, dimension (:) :: varin, varout

      n_tslices = NOCTAD

      status = NF_OPEN('data/CATCH/MODIS-Albedo2/MCD43GF_wsa_H11V13.nc',NF_NOWRITE, ncid); VERIFY_(STATUS)
      status = NF_INQ_DIM (ncid,3,string, n_tslices); VERIFY_(STATUS)
      allocate (MMDD      (0: n_tslices + 1))
      allocate (MMDD_next (0: n_tslices + 1))
      status = NF_GET_VARA_text(ncid, 3,(/1,1/),(/4,n_tslices/),MMDD(1:n_tslices)); VERIFY_(STATUS)
      status = NF_CLOSE(ncid); VERIFY_(STATUS)

      mmdd(0) = mmdd(n_tslices)
      mmdd(n_tslices + 1)= mmdd(1)
      mmdd_next(0:n_tslices - 1) =  mmdd(1:n_tslices)
      mmdd_next(n_tslices: n_tslices + 1) = mmdd (1:2)

      allocate (varin (1:InNTILES))
      allocate (varout(1:NTILES))

      STATUS = NF_OPEN('data/CATCH/CATCHCN_fPAR_Alb_stats2.nc4',      NF_NOWRITE,NCOutID) ; VERIFY_(STATUS)
      open (10, file = 'clsm/MODISVISmean.dat'  ,form='unformatted',status='unknown',convert='little_endian')
      open (11, file = 'clsm/MODISNIRmean.dat'  ,form='unformatted',status='unknown',convert='little_endian')
      open (12, file = 'clsm/MODISVISstd.dat'   ,form='unformatted',status='unknown',convert='little_endian')
      open (13, file = 'clsm/MODISNIRstd.dat'   ,form='unformatted',status='unknown',convert='little_endian')
      open (14, file = 'clsm/MODELFPARmean.dat' ,form='unformatted',status='unknown',convert='little_endian')
      open (15, file = 'clsm/MODELFPARstd.dat'  ,form='unformatted',status='unknown',convert='little_endian')
      open (16, file = 'clsm/MODISFPARmean.dat' ,form='unformatted',status='unknown',convert='little_endian')
      open (17, file = 'clsm/MODISFPARstd.dat'  ,form='unformatted',status='unknown',convert='little_endian')

      do t =0,n_tslices+1
         
         time_slice = t
         yr = 1
         yr1= 1
         if(t == 0) then
            time_slice =  n_tslices
            yr         =  1 - 1
         endif
         
         if(t >= n_tslices) then 
            yr1 = 1 + 1
            if(t ==n_tslices + 1) then
               time_slice =  1
               yr = 1 + 1
            endif
         endif
         
         read(mmdd(t),'(i2.2,i2.2)') mn,dd
         read(mmdd_next(t),'(i2.2,i2.2)') mn1,dd1

         STATUS = NF_GET_VARA_REAL(NCOutID,VARID(NCOutID,'MODISVISmean'   ),(/1,time_slice/), (/InNTILES,1/), varin ) ; VERIFY_(STATUS)
         varout = varin (id_loc)
         write(10) float((/yr,mn,dd,0,0,0,yr1,mn1,dd1,0,0,0,ntiles,1/))
         write(11) float((/yr,mn,dd,0,0,0,yr1,mn1,dd1,0,0,0,ntiles,1/))
         write(12) float((/yr,mn,dd,0,0,0,yr1,mn1,dd1,0,0,0,ntiles,1/))
         write(13) float((/yr,mn,dd,0,0,0,yr1,mn1,dd1,0,0,0,ntiles,1/))
         write(14) float((/yr,mn,dd,0,0,0,yr1,mn1,dd1,0,0,0,ntiles,1/))
         write(15) float((/yr,mn,dd,0,0,0,yr1,mn1,dd1,0,0,0,ntiles,1/))
         write(16) float((/yr,mn,dd,0,0,0,yr1,mn1,dd1,0,0,0,ntiles,1/))
         write(17) float((/yr,mn,dd,0,0,0,yr1,mn1,dd1,0,0,0,ntiles,1/))

         STATUS = NF_GET_VARA_REAL(NCOutID,VARID(NCOutID,'MODISVISmean'   ),(/1,time_slice/), (/InNTILES,1/), varin ) ; VERIFY_(STATUS)
         varout = varin (id_loc)
         write (10) varout

         STATUS = NF_GET_VARA_REAL(NCOutID,VARID(NCOutID,'MODISNIRmean'   ),(/1,time_slice/), (/InNTILES,1/), varin ) ; VERIFY_(STATUS)
         varout = varin (id_loc)
         write (11) varout

         STATUS = NF_GET_VARA_REAL(NCOutID,VARID(NCOutID,'MODISVISstd'   ),(/1,time_slice/) , (/InNTILES,1/), varin ) ; VERIFY_(STATUS)
         varout = varin (id_loc)
         write (12) varout

         STATUS = NF_GET_VARA_REAL(NCOutID,VARID(NCOutID,'MODISNIRstd'   ),(/1,time_slice/) , (/InNTILES,1/), varin ) ; VERIFY_(STATUS)
         varout = varin (id_loc)
         write (13) varout

         STATUS = NF_GET_VARA_REAL(NCOutID,VARID(NCOutID,'MODELFPARmean'   ),(/1,time_slice/), (/InNTILES,1/), varin ) ; VERIFY_(STATUS)
         varout = varin (id_loc)
         write (14) varout

         STATUS = NF_GET_VARA_REAL(NCOutID,VARID(NCOutID,'MODELFPARstd'   ),(/1,time_slice/) , (/InNTILES,1/), varin ) ; VERIFY_(STATUS)
         varout = varin (id_loc)
         write (15) varout

         STATUS = NF_GET_VARA_REAL(NCOutID,VARID(NCOutID,'MODISFPARmean'   ),(/1,time_slice/), (/InNTILES,1/), varin ) ; VERIFY_(STATUS)
         varout = varin (id_loc)
         write (16) varout

         STATUS = NF_GET_VARA_REAL(NCOutID,VARID(NCOutID,'MODISFPARstd'   ),(/1,time_slice/) , (/InNTILES,1/), varin ) ; VERIFY_(STATUS)
         varout = varin (id_loc)
         write (17) varout
      end do

      deallocate (varin, varout)

      close (10 , status = 'keep')
      close (11 , status = 'keep')
      close (12 , status = 'keep')
      close (13 , status = 'keep')
      close (14 , status = 'keep')
      close (15 , status = 'keep')
      close (16 , status = 'keep')
      close (17 , status = 'keep')
      
    END SUBROUTINE regrid_alb

    ! -------------------------------------------------------------------------------------------------

    SUBROUTINE preprocess_m09

      implicit none
      logical     :: file_exists
      INTEGER     :: NT, ND, DAY, year, STATUS, NCOutID, k
      CHARACTER*8 :: YYYYMMDD
      CHARACTER*6 :: YYYYMM
      CHARACTER*4 :: YYYY
      CHARACTER*2 :: MM, DD
      real        :: yr,mn,dy,dum,yr1,mn1,dy1

      real,    allocatable, dimension (:,:) :: MODIS_VISDF, MODIS_NIRDF, MODEL_fPAR 
      real,    allocatable, dimension (:)   :: data_read, data_save,  &
               MODISVISmean, MODISNIRmean, MODISVISstd, MODISNIRstd, MODELFPARmean, MODELFPARstd, &
               MODISFPARmean, MODISFPARstd
      integer, allocatable, dimension (:)   :: ldas2bcs
      type(date_time_type), dimension (YearE - YearB + 1) :: octad_time


      inquire(file='data/CATCH/CATCHCN_fPAR_Alb_stats2.nc4', exist=file_exists)
      if(.not. file_exists) call create_stat_file

!      open (99,file='clsm/comp_CATCHCN_AlbScale_parameters.log', form ='formatted', action='write', status= 'unknown')

      STATUS = NF_OPEN ('data/CATCH/CATCHCN_fPAR_Alb_stats2.nc4', NF_WRITE,NCOutID) ; VERIFY_(STATUS)

      allocate (MODIS_VISDF (1:InNTILES, yearE - yearB + 1))
      allocate (MODIS_NIRDF (1:InNTILES, yearE - yearB + 1))
      allocate (MODEL_fPAR  (1:InNTILES, yearE - yearB + 1))
      allocate (ldas2bcs    (1:InNTILES))
      allocate (data_read   (1:InNTILES))
      allocate (data_save   (1:InNTILES))   
      allocate (MODISVISmean (1:InNTILES))   
      allocate (MODISNIRmean (1:InNTILES))   
      allocate (MODISVISstd  (1:InNTILES))   
      allocate (MODISNIRstd  (1:InNTILES))   
      allocate (MODELFPARmean(1:InNTILES))   
      allocate (MODELFPARstd (1:InNTILES))   
      allocate (MODISFPARmean(1:InNTILES))   
      allocate (MODISFPARstd (1:InNTILES))   

      open (10,file =trim(EXPDIR)//'rc_out/'//trim(EXNAME)//'.ldas_tilecoord.bin',status='old',form='unformatted',convert='big_endian')
      read (10) k
      read (10)  LDAS2BCS
      close(10, status = 'keep')

      OPEN_FILES1 : DO year = YearB, YearE

         write (YYYY ,'(i4.4)')  year
         open (10 + year - yearB, file = trim(InBCSDIR)//'MODIS6/'//YYYY//'//visdf.dat', form = 'unformatted', action = 'read')
         open (30 + year - yearB, file = trim(InBCSDIR)//'MODIS6/'//YYYY//'//nirdf.dat', form = 'unformatted', action = 'read')
!         WRITE (99,*)10 + year - yearB, YYYY//'//visdf.dat'  
!         WRITE (99,*)30 + year - yearB, YYYY//'//nirdf.dat'
!         WRITE (99,*) ' '
         WRITE (*,*)10 + year - yearB, YYYY//'//visdf.dat'  
         WRITE (*,*)30 + year - yearB, YYYY//'//nirdf.dat'
         WRITE (*,*) ' '
         octad_time(year - yearB + 1)%year  = year - 1
         octad_time(year - yearB + 1)%month = 12
         octad_time(year - yearB + 1)%day   = 31
         octad_time(year - yearB + 1)%hour  = 0
         octad_time(year - yearB + 1)%min   = 0
         octad_time(year - yearB + 1)%sec   = 0  
         
      END DO OPEN_FILES1

      ND = 8

      OCTAD_LOOP : DO NT = 1, NOCTAD   

         ! BEGIN READING MODIS VISDF/NIRDF

         MODIS_VISDF = 0.
         MODIS_NIRDF = 0.
         MODEL_fPAR  = 0.

         if(NT == NOCTAD) ND = 5
!         WRITE (99,*) NT, ND, yearB, yearE
         print *, NT, ND, yearB, yearE

         READ_YEARS : DO year = YearB, YearE

            read (10 + year - yearB) modis_visdf (:, year - yearB + 1)
            read (30 + year - yearB) modis_nirdf (:, year - yearB + 1)

            DAILY_LOOP : DO day = 1,ND

               call augment_date_time(86400, octad_time(year - yearB + 1))

               write (YYYY, '(i4.4)') octad_time(year - yearB + 1)%year
               write (MM  , '(i2.2)') octad_time(year - yearB + 1)%month
               write (DD  , '(i2.2)') octad_time(year - yearB + 1)%day

               YYYYMMDD = YYYY//MM//DD
!               WRITE (99,*) trim(EXNAME)//'.ens_avg.ldas_tile_daily_out.'//YYYYMMDD//'.bin'
               print *, trim(EXNAME)//'.ens_avg.ldas_tile_daily_out.'//YYYYMMDD//'.bin'
               open (60, file = trim(EXPDIR)//'cat/ens_avg/Y'//YYYY//'/M'//MM//'/'// &
                   trim(EXNAME)//'.ens_avg.ldas_tile_daily_out.'//YYYYMMDD//'.bin',  &
                   form = 'unformatted', convert='big_endian', action = 'read')   

               do k = 1,3
                  read (60) data_read
                  if(k == 2) data_save = data_read
               end do
              
               MODEL_fPAR (:,year - yearB + 1) = MODEL_fPAR (:,year - yearB + 1) + data_save / (data_read + 1.e-20)/ real (ND)              
               close (60, status = 'keep')

            END DO DAILY_LOOP

            ! reoder to the order of BCs
     
            data_read = MODEL_fPAR (:,year - yearB + 1)
     
            do k = 1, InNTILES
               MODEL_fPAR (LDAS2BCS(k),year - yearB + 1) = data_read (k) 
            end do
         END DO READ_YEARS
         
         ! COMPUTE STATS 
         
         CALL compute_stats (MODIS_VISDF, MODIS_NIRDF, MODEL_fPAR, &
              MODISVISmean, MODISNIRmean, MODISVISstd, MODISNIRstd, MODELFPARmean, MODELFPARstd)

         STATUS = NF_PUT_VARA_REAL(NCOutID,VARID(NCOutID,'MODISVISmean'   ),(/1,NT/), (/InNTILES,1/), MODISVISmean ) ; VERIFY_(STATUS)
         STATUS = NF_PUT_VARA_REAL(NCOutID,VARID(NCOutID,'MODISNIRmean'   ),(/1,NT/), (/InNTILES,1/), MODISNIRmean ) ; VERIFY_(STATUS)     
         STATUS = NF_PUT_VARA_REAL(NCOutID,VARID(NCOutID,'MODISVISstd'    ),(/1,NT/), (/InNTILES,1/), MODISVISstd  ) ; VERIFY_(STATUS)
         STATUS = NF_PUT_VARA_REAL(NCOutID,VARID(NCOutID,'MODISNIRstd'    ),(/1,NT/), (/InNTILES,1/), MODISNIRstd  ) ; VERIFY_(STATUS)     
         STATUS = NF_PUT_VARA_REAL(NCOutID,VARID(NCOutID,'MODELFPARmean'  ),(/1,NT/), (/InNTILES,1/), MODELFPARmean) ; VERIFY_(STATUS)
         STATUS = NF_PUT_VARA_REAL(NCOutID,VARID(NCOutID,'MODELFPARstd'   ),(/1,NT/), (/InNTILES,1/), MODELFPARstd ) ; VERIFY_(STATUS)     

      END DO OCTAD_LOOP

      CLOSE_FILES1 : DO year = YearB, YearE

         close (10 + year - yearB, status = 'keep')
         close (30 + year - yearB, status = 'keep')

       END DO CLOSE_FILES1

      ! MODIS FPAR
      ! ----------

      deallocate (MODIS_VISDF, MODIS_NIRDF, MODEL_fPAR)
      allocate (MODIS_VISDF (1:InNTILES, yearE - yearB1 + 1))
      allocate (MODIS_NIRDF (1:InNTILES, yearE - yearB1 + 1))
      allocate (MODEL_fPAR  (1:InNTILES, yearE - yearB1 + 1))

      OPEN_FILES2 : DO year = YearB1, YearE

         write (YYYY ,'(i4.4)')  year
         open (10 + year - yearB1, file = trim(InBCSDIR)//'MODIS6/'//YYYY//'//visdf.dat', form = 'unformatted', action = 'read')
         open (30 + year - yearB1, file = trim(InBCSDIR)//'MODIS6/'//YYYY//'//nirdf.dat', form = 'unformatted', action = 'read')
         open (50 + year - yearB1, file = trim(InBCSDIR)//'MODIS6/'//YYYY//'//fpar.dat', form = 'unformatted', action = 'read')
!         WRITE (99,*) "              "
!         WRITE (99,*) "MODIS FPAR"
!         WRITE (99,*) "=========="
!         WRITE (99,*) "              "

!         WRITE (99,*)10 + year - yearB1, YYYY//'//visdf.dat'  
!         WRITE (99,*)30 + year - yearB1, YYYY//'//nirdf.dat'
!         WRITE (99,*)50 + year - yearB1, YYYY//'//fpar.dat'
!         WRITE (99,*) ' '
         WRITE (*,*)10 + year - yearB1, YYYY//'//visdf.dat'  
         WRITE (*,*)30 + year - yearB1, YYYY//'//nirdf.dat'
         WRITE (*,*)50 + year - yearB1, YYYY//'//fpar.dat'
         WRITE (*,*) ' '
         read (50 + year - yearB1) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
         WRITE (*,*) yr,mn,dy,yr1,mn1,dy1
         read (50 + year - yearB1) MODEL_FPAR (:, year - yearB1 + 1)

         octad_time(year - yearB1 + 1)%year  = year - 1
         octad_time(year - yearB1 + 1)%month = 12
         octad_time(year - yearB1 + 1)%day   = 31
         octad_time(year - yearB1 + 1)%hour  = 0
         octad_time(year - yearB1 + 1)%min   = 0
         octad_time(year - yearB1 + 1)%sec   = 0  
         
      END DO OPEN_FILES2

      ND = 8

      OCTAD_LOOP2 : DO NT = 1, NOCTAD   

         ! BEGIN READING MODIS VISDF/NIRDF/FPAR

         MODIS_VISDF = 0.
         MODIS_NIRDF = 0.
         MODEL_fPAR  = 0.

         if(NT == NOCTAD) ND = 5
!         WRITE (99,*) NT, ND, yearB1, yearE
         print *, NT, ND, yearB1, yearE
         READ_YEARS2 : DO year = YearB1, YearE

            read (10 + year - yearB1) modis_visdf (:, year - yearB1 + 1)
            read (30 + year - yearB1) modis_nirdf (:, year - yearB1 + 1)
            read (50 + year - yearB1) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
            WRITE (*,*) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
            read (50 + year - yearB1) MODEL_FPAR (:, year - yearB1 + 1)
         END DO READ_YEARS2
         
         ! COMPUTE STATS 
         
         CALL compute_stats (MODIS_VISDF, MODIS_NIRDF, MODEL_fPAR, &
              MODISVISmean, MODISNIRmean, MODISVISstd, MODISNIRstd, MODELFPARmean, MODELFPARstd)

         STATUS = NF_PUT_VARA_REAL(NCOutID,VARID(NCOutID,'MODISFPARmean'  ),(/1,NT/), (/InNTILES,1/), MODELFPARmean) ; VERIFY_(STATUS)
         STATUS = NF_PUT_VARA_REAL(NCOutID,VARID(NCOutID,'MODISFPARstd'   ),(/1,NT/), (/InNTILES,1/), MODELFPARstd ) ; VERIFY_(STATUS)     

      END DO OCTAD_LOOP2

      CLOSE_FILES2 : DO year = YearB1, YearE

         close (10 + year - yearB1, status = 'keep')
         close (30 + year - yearB1, status = 'keep')
         close (50 + year - yearB1, status = 'keep')

      END DO CLOSE_FILES2

      STATUS = NF_CLOSE (NCOutID )

    END SUBROUTINE preprocess_m09

    ! -----------------------------------------------------------------

    SUBROUTINE create_stat_file

      implicit none

      integer                :: NCFOutID, STATUS, vid, tid, lid, n, k
      character (22)         :: time_stamp, tmpstr
      integer, dimension(8)  :: date_time_values
      real, dimension (:), allocatable :: lons, lats

      STATUS = NF_CREATE ('data/CATCH/CATCHCN_fPAR_Alb_stats2.nc4',  NF_NETCDF4, NCFOutID );VERIFY_(STATUS)
      STATUS = NF_DEF_DIM(NCFOutID, 'octad', NOCTAD, TID) ; VERIFY_(STATUS) 
      STATUS = NF_DEF_DIM(NCFOutID, 'tiles', InNTILES, LID) ; VERIFY_(STATUS)
      STATUS = NF_DEF_VAR(NCFOutID, 'lon'    , NF_FLOAT ,1 ,(/LID/), vid);VERIFY_(STATUS)
      STATUS = NF_DEF_VAR(NCFOutID, 'lat'    , NF_FLOAT ,1 ,(/LID/), vid);VERIFY_(STATUS)
      status = NF_DEF_VAR(NCFOutID, 'MODISVISmean', NF_FLOAT ,2 ,(/LID,TID/), vid);VERIFY_(STATUS)
      status = NF_DEF_VAR(NCFOutID, 'MODISNIRmean', NF_FLOAT ,2 ,(/LID,TID/), vid);VERIFY_(STATUS)
      status = NF_DEF_VAR(NCFOutID, 'MODISVISstd' , NF_FLOAT ,2 ,(/LID,TID/), vid);VERIFY_(STATUS)
      status = NF_DEF_VAR(NCFOutID, 'MODISNIRstd' , NF_FLOAT ,2 ,(/LID,TID/), vid);VERIFY_(STATUS)
      status = NF_DEF_VAR(NCFOutID, 'MODELFPARmean',NF_FLOAT ,2 ,(/LID,TID/), vid);VERIFY_(STATUS)
      status = NF_DEF_VAR(NCFOutID, 'MODELFPARstd' ,NF_FLOAT ,2 ,(/LID,TID/), vid);VERIFY_(STATUS)
      status = NF_DEF_VAR(NCFOutID, 'MODISFPARmean',NF_FLOAT ,2 ,(/LID,TID/), vid);VERIFY_(STATUS)
      status = NF_DEF_VAR(NCFOutID, 'MODISFPARstd' ,NF_FLOAT ,2 ,(/LID,TID/), vid);VERIFY_(STATUS)
      !  Global attributes
      
      call date_and_time(VALUES=date_time_values)
      
      write (time_stamp,'(i4.4,a1,i2.2,a1,i2.2,1x,a2,1x,i2.2,a1,i2.2,a1,i2.2)')      &
           date_time_values(1),'-',date_time_values(2),'-',date_time_values(3),'at', &
           date_time_values(5),':',date_time_values(6),':',date_time_values(7)
      
      status = NF_PUT_ATT_TEXT(NCFOutID, NF_GLOBAL, 'CreatedBy', LEN_TRIM("Sarith Mahanama"),  &
           trim("Sarith Mahanama"))
      status = NF_PUT_ATT_TEXT(NCFOutID, NF_GLOBAL, 'Date'   , LEN_TRIM(time_stamp),trim(time_stamp))
      
      status = NF_ENDDEF(NCFOutID )      

      ! Read and put lat/lon data

      allocate (lons (1:InNTILES))
      allocate (lats (1:InNTILES))

      open (10, file = trim(InBCSDIR)//trim(InGFILE)//'.til', form = 'formatted', status = 'old')
      
      do n = 1,8
         read (10,*) tmpstr
      end do

      do n = 1, InNTILES
         read (10,*) k, k, lons(n), lats(n)
      end do

      close (10, status = 'keep')

      STATUS = NF_PUT_VARA_REAL(NCFOutID,VARID(NCFOutID,'lon'   ),(/1/), (/InNTILES/), lons) ; VERIFY_(STATUS)
      STATUS = NF_PUT_VARA_REAL(NCFOutID,VARID(NCFOutID,'lat'   ),(/1/), (/InNTILES/), lats) ; VERIFY_(STATUS)      

      STATUS = NF_CLOSE (NCFOutID )

    END SUBROUTINE create_stat_file

    ! ----------------------------------------------------------------------

    SUBROUTINE compute_stats (MODIS_VISDF, MODIS_NIRDF, MODEL_fPAR, &
         MODISVISmean, MODISNIRmean, MODISVISstd, MODISNIRstd, MODELFPARmean, MODELFPARstd)

      implicit none
      real, dimension (:,:), intent (in)    :: MODIS_VISDF, MODIS_NIRDF, MODEL_fPAR
      real, dimension (:)  , intent (inout) :: MODISVISmean, MODISNIRmean, MODISVISstd, & 
                                               MODISNIRstd, MODELFPARmean, MODELFPARstd
      integer :: NX, NY, N, t
      REAL    :: MF, MV, MN, SF, SV, SN, ZV, ZN, CV, CN, CF

      NX = size (MODIS_VISDF,1)
      NY = size (MODIS_VISDF,2)
      print *,'Entered compute_stats', NX, NY
      if (NX /= InNTILES) then 
         print *, 'NX NTILLES MISMAATCH : ', InNTILES, NX, NY 
         STOP 
      ENDIF
      
      DO N = 1, NX

!          MF = SUM (MODEL_fPAR  (N,:)) / REAL (NY)
!          MV = SUM (MODIS_VISDF (N,:)) / REAL (NY) 
!          MN = SUM (MODIS_nirDF (N,:)) / REAL (NY)
         MF = 0.
         MV = 0.
         MN = 0.
         CV = 0.
         CN = 0.
         CF = 0.

         do T = 1, NY
            if ((MODEL_fPAR  (N,T) >= 0.).AND.(MODEL_fPAR  (N,T) <= 1.)) then
               MF = MF + MODEL_fPAR  (N,T)
               CF = CF + 1.
            endif
            if ((MODIS_VISDF (N,T) >= 0.).AND.(MODIS_VISDF (N,T) <= 1.)) then
               MV = MV + MODIS_VISDF (N,T)
               CV = CV + 1.
            endif
            if ((MODIS_NIRDF (N,T) >= 0.).AND.(MODIS_NIRDF (N,T) <= 1.)) then
               MN = MN + MODIS_NIRDF (N,T)
               CN = CN + 1.
            endif
         end do

         IF(CF > 0) MF = MF / CF
         IF(CV > 0) MV = MV / CV
         IF(CN > 0) MN = MN / CN
 
          ! STANDARD DEVIATION

          SF = 1.e-15
          SV = 1.e-15
          SN = 1.e-15
       
          do T = 1, NY
             if ((MODEL_fPAR  (N,T) >= 0.).AND.(MODEL_fPAR  (N,T) <= 1.)) SF = SF + (MODEL_fPAR  (N,t) - MF)*(MODEL_fPAR  (N,t) - MF) 
             if ((MODIS_VISDF (N,T) >= 0.).AND.(MODIS_VISDF (N,T) <= 1.)) SV = SV + (MODIS_VISDF (N,t) - MV)*(MODIS_VISDF (N,t) - MV) 
             if ((MODIS_NIRDF (N,T) >= 0.).AND.(MODIS_NIRDF (N,T) <= 1.)) SN = SN + (MODIS_NIRDF (N,t) - MN)*(MODIS_NIRDF (N,t) - MN) 
          end do

          IF(CF > 0) SF = SQRT (SF / CF) 
          IF(CV > 0) SV = SQRT (SV / CV) 
          IF(CN > 0) SN = SQRT (SN / CN) 

          ! CORRELATION

          ZV = 0.
          ZN = 0.

          DO T = 1, NY
             ZV = ZV + (MODEL_fPAR  (N,t) - MF)*(MODIS_VISDF (N,t) - MV)/SF/SV
             ZN = ZN + (MODEL_fPAR  (N,t) - MF)*(MODIS_NIRDF (N,t) - MN)/SF/SV 
          END DO

          ZV = ZV / REAL (NY)
          ZN = ZN / REAL (NY)

          MODISVISmean (N) = MV  
          MODISNIRmean (N) = MN   
          MODISVISstd  (N) = SV   
          MODISNIRstd  (N) = SN   
          MODELFPARmean(N) = MF   
          MODELFPARstd (N) = SF

          if(ZV < 0.) MODISVISstd  (N) = -1. * MODISVISstd  (N)
          if(Zn < 0.) MODISnirstd  (N) = -1. * MODISnirstd  (N)

      END DO

      print *,'Leaving compute_stats'

    END SUBROUTINE compute_stats

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

   ! *****************************************************************************

   subroutine ReadCNTilFile (InCNTileFile, nt, xlon, xlat)

     implicit none
     character(*), intent (in) ::  InCNTileFile
     integer , intent (in) :: nt
     real, dimension (nt), intent(inout) :: xlon, xlat
     integer :: n,icnt,ityp
     real    :: xval,yval, pf
     
   open(11,file=InCNTileFile, &
        form='formatted',action='read',status='old')

   do n = 1,8 ! skip header
      read(11,*)
   end do
   
   icnt = 0
   ityp = 100

   do while (ityp == 100) ! loop over land tiles
      read(11,*) ityp,pf,xval,yval
      if(ityp == 100) then
         icnt = icnt + 1
         xlon(icnt) = xval
         xlat(icnt) = yval
      endif
   end do

   close(11)
    
   end subroutine ReadCNTilFile
   
   ! *****************************************************************************
   
   real function haversine(deglat1,deglon1,deglat2,deglon2)
     ! great circle distance -- adapted from Matlab 
     real,intent(in) :: deglat1,deglon1,deglat2,deglon2
     real :: a,c, dlat,dlon,lat1,lat2
     real,parameter :: radius = 6371.0E3 
     
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
   
   ! *****************************************************************************
    
 END MODULE comp_CATCHCN_AlbScale_parameters
