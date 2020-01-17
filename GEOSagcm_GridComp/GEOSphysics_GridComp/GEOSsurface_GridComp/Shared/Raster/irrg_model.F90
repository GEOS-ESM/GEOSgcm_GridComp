#define VERIFY_(A)   IF(A/=0)THEN;PRINT *,'ERROR AT LINE ', __LINE__;STOP;ENDIF
#define ASSERT_(A)   if(.not.A)then;print *,'Error:',__FILE__,__LINE__;stop;endif

PROGRAM irrg_model

  use rmTinyCatchParaMod

  implicit none

     integer, parameter                    :: ncols = 86400, nrows_data = 36000, nrows = 43200, nc = 43200, nr = 21600
     real*4, allocatable                   :: var_in (:,:), cnt_pix1 (:), cnt_pix2 (:), cnt_pix3 (:), tot_cnt (:) 
     real, allocatable                     :: lai_min(:), lai_max(:), lai(:)  
     real                                  :: yr,mo,dy,hr,mn,ss,yr1,mo1,dy1,hr1,mn1,ss1,nx,ny
     integer                               :: i,j, n, r, ntiles,status, cellid, vid, NCFID
     character*400                         :: GFILE,arg, LAIFILE
     integer, allocatable, dimension (:,:) :: tile_id
     integer, pointer                      :: iraster  (:,:)
     real, allocatable                     :: CLM4_pf(:), CLM4_sf(:), CLM4_pt(:), CLM4_st(:) 
     real                                  :: CLMC_pf1, CLMC_pf2, CLMC_sf1, CLMC_sf2
     real                                  :: CLMC_pt1, CLMC_pt2, CLMC_st1, CLMC_st2    
     INCLUDE 'netcdf.inc'

     I = iargc()     
     GFILE =''
     call getarg(1,arg)
     read(arg,'(a)') GFILE
     call getarg(2,arg)
     read(arg,'(a)') LAIFILE

     allocate( var_in(ncols,nrows)) 
     var_in = -9999.

!     open ( 10, file = '/gpfsm/dnb43/projects/p03/LS_PARAMETERS/irrigation/global_gripc/irrigtype_salmon2013.flt', &
     open ( 10, file = '/discover/nobackup/rreichle/l_data/LandBCs_files_for_mkCatchParam/V001/irrigtype_salmon2013.flt', &
          form = 'unformatted', access='direct', recl=(ncols))

     !- Read input file::
     print *, " -- Reading in input file -- "
     do j = 1, nrows_data
        r = nrows -j + 1
!        print *,j,r
        read(10,rec=j) var_in(:, r)
        do i = 1, ncols
           if( var_in(i, r) == 0. ) var_in(i, r) = -9999.
           if( var_in(i, r) == 4. ) var_in(i, r) = -9999.
        end do
     end do
     close( 10 )
    
! Reading rst file

   open (10,file='rst/'//trim(gfile)//'.rst',status='old',action='read',  &
        form='unformatted',convert='little_endian')
   allocate (tile_id    (1:nc,1:nr))         
   
   do j=1,nr
      read(10)tile_id(:,j)
   end do
   close (10,status='keep')

! Reading number of catchments

   open (10,file='clsm/catchment.def',status='old',action='read',  &
        form='formatted')
   read (10, *) ntiles
   close (10, status = 'keep')
   
   allocate(iraster(ncols,nrows),stat=STATUS); VERIFY_(STATUS)
   call RegridRaster(tile_id,iraster)

   allocate (cnt_pix1 (1:ntiles))
   allocate (cnt_pix2 (1:ntiles))
   allocate (cnt_pix3 (1:ntiles))
   allocate (tot_cnt  (1:ntiles))
   allocate (CLM4_pf  (1:ntiles))
   allocate (CLM4_sf  (1:ntiles))
   allocate (CLM4_pt  (1:ntiles))
   allocate (CLM4_st  (1:ntiles))
   allocate (lai_min  (ntiles))
   allocate (lai_max  (ntiles))
   allocate (lai      (ntiles))

   cnt_pix1 = 0.
   cnt_pix2 = 0.
   cnt_pix3 = 0.
   tot_cnt  = 0.

   do j = 1,nrows
      do i =  1,ncols
         if((iraster (i,j) >=1).and.(iraster (i,j) <=ntiles)) then
            tot_cnt (iraster (i,j)) = tot_cnt (iraster (i,j)) + 1.
            if (var_in(i,j) == 1) cnt_pix1(iraster (i,j)) = cnt_pix1(iraster (i,j)) + 1.
            if (var_in(i,j) == 2) cnt_pix2(iraster (i,j)) = cnt_pix2(iraster (i,j)) + 1.
            if (var_in(i,j) == 3) cnt_pix3(iraster (i,j)) = cnt_pix3(iraster (i,j)) + 1.            
         endif
      end do
   end do

   cnt_pix1 = cnt_pix1 / tot_cnt
   cnt_pix2 = cnt_pix2 / tot_cnt
   cnt_pix3 = cnt_pix3 / tot_cnt

   ! CLM typs and laimin, laimax
   ! ---------------------------
   
   open(unit=27, file='clsm/CLM_veg_typs_fracs'   ,form='formatted')

   do n=1,ntiles
      read (27, *) i,j, CLMC_pt1, CLMC_pt2, CLMC_st1, CLMC_st2, &
           CLMC_pf1, CLMC_pf2, CLMC_sf1, CLMC_sf2,              & 
           CLM4_pt(n), CLM4_st(n), CLM4_pf(n),CLM4_sf(n) 
   end do

   CLOSE (27, STATUS = 'KEEP')

   lai_max = -9999.
   lai_min = 9999.

   open (31, file =trim(LAIFILE), form = 'unformatted', action = 'read', status = 'old')
   
   READ_LAIFILE : do i = 1,100
      
      read (31, IOSTAT= STATUS)  yr,mo,dy,hr,mn,ss,yr1,mo1,dy1,hr1,mn1,ss1,nx,ny
      if(STATUS /= 0) exit
      read(31)lai
      
      do n = 1, ntiles
         if(lai (n) < lai_min (n)) lai_min (n) = lai (n)
         if(lai (n) > lai_max (n)) lai_max (n) = lai (n)
      end do
      
   end do READ_LAIFILE

   CLOSE (31, STATUS = 'KEEP')
   
!   open (10,file = 'clsm/gripc.data', form = 'unformatted', action = 'write', status = 'unknown')
!   write (10) cnt_pix1
!   write (10) cnt_pix2
!   write (10) cnt_pix3
!   close (10, status = 'keep')

        status = NF_CREATE ('clsm/irrigation_internal_rst', NF_NETCDF4, NCFID)     ; VERIFY_(STATUS)
        status = NF_DEF_DIM(NCFID, 'tile' , NTILES, CellID)                   ; VERIFY_(STATUS)
        status = NF_DEF_VAR(NCFID, 'IRRIGFRAC'    , NF_FLOAT, 1 ,CellID, vid) ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCFID, vid, 'long_name',         &   
           LEN_TRIM('Fraction of irrigated cropland'),           &
           'fraction of irrigated cropland')                                 ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCFID, vid, 'units', 1,'1')                  ; VERIFY_(STATUS)

        status = NF_DEF_VAR(NCFID, 'PADDYFRAC'    , NF_FLOAT, 1 ,CellID, vid) ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCFID, vid, 'long_name',         &   
            LEN_TRIM('fraction of paddy cropland'),               &
            'fraction of paddy cropland')                                     ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCFID, vid, 'units', 1,'1')                  ; VERIFY_(STATUS)

        status = NF_DEF_VAR(NCFID, 'LAIMIN'    , NF_FLOAT, 1 ,CellID, vid)    ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCFID, vid, 'long_name',        &                             
            LEN_TRIM('Minimum LAI'), 'Minimum LAI')                           ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCFID, vid, 'units', 1,'1')                  ; VERIFY_(STATUS) 

        status = NF_DEF_VAR(NCFID, 'LAIMAX'    , NF_FLOAT, 1 ,CellID, vid)    ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCFID, vid, 'long_name',        &                             
            LEN_TRIM('Maximum LAI'), 'Maximum LAI')                           ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCFID, vid, 'units', 1,'1')                  ; VERIFY_(STATUS)

        status = NF_DEF_VAR(NCFID, 'CLMPT'    , NF_FLOAT, 1 ,CellID, vid)     ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCFID, vid, 'long_name',        &                             
            LEN_TRIM('CLM primary type'), 'CLM primary type')                 ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCFID, vid, 'units', 1,'1')                  ; VERIFY_(STATUS)

        status = NF_DEF_VAR(NCFID, 'CLMST'    , NF_FLOAT, 1 ,CellID, vid)     ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCFID, vid, 'long_name',        &                             
            LEN_TRIM('CLM secondary type'), 'CLM secondary type')             ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCFID, vid, 'units', 1,'1')                  ; VERIFY_(STATUS)

        status = NF_DEF_VAR(NCFID, 'CLMPF'    , NF_FLOAT, 1 ,CellID, vid)     ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCFID, vid, 'long_name',        &                             
            LEN_TRIM('CLM primary fraction'), 'CLM primary fraction')         ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCFID, vid, 'units', 1,'1')                  ; VERIFY_(STATUS)

        status = NF_DEF_VAR(NCFID, 'CLMSF'    , NF_FLOAT, 1 ,CellID, vid)     ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCFID, vid, 'long_name',        &                             
            LEN_TRIM('CLM secondary fraction'), 'CLM secondary fraction')     ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCFID, vid, 'units', 1,'1')                  ; VERIFY_(STATUS)

        status  = NF_ENDDEF(NCFID)

       status = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'IRRIGFRAC') ,(/1/),(/NTILES/), cnt_pix2) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'PADDYFRAC') ,(/1/),(/NTILES/), cnt_pix3) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'LAIMIN'   ) ,(/1/),(/NTILES/), lai_min)  ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'LAIMAX'   ) ,(/1/),(/NTILES/), lai_max)  ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'CLMPT'    ) ,(/1/),(/NTILES/), CLM4_pt)  ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'CLMST'    ) ,(/1/),(/NTILES/), CLM4_st)  ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'CLMPF'    ) ,(/1/),(/NTILES/), CLM4_pf/100.); VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCFID,VarID(NCFID,'CLMSF'    ) ,(/1/),(/NTILES/), CLM4_sf/100.); VERIFY_(STATUS)
       STATUS   = NF_CLOSE (NCFID)
       
       contains

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
        PRINT *, trim(Line),': ',NF_STRERROR(STATUS)
        STOP 'Stopped'
     ENDIF

   END SUBROUTINE HANDLE_ERR

   ! -----------------------------------------------------------------------------
 
END PROGRAM irrg_model
