PROGRAM STEP3_process_soildata

! Legacy starting code was written by Sarith Mahanama since then code evolved with updated data sets, file formats and bug fixes.

! Pipeline summary (STEP3):
!  1) Read NGDC/Reynolds 0.5° (4320×2160) fields from `path` (native bin/img inputs).
!  2) Regrid NGDC/Reynolds to 30 arcsec (43200×21600) and write F77-unformatted
!     binaries into `path_regrid` (v2/MERGED-DATA/):
!       PRtestBiljana_{gtext,clay,silt,sand,om,poros}_{top,sub}_30sec.dat
!  3) Merge regridded NGDC/Reynolds (layer1) with HWSDv2 NetCDF fields (layer2),
!     using the HWSD land mask, and write merged 30-arcsec binaries into `path3`
!     (v2/MERGED-DATA/):
!       {clay,silt,sand,oc,gtext}_{top,sub}_30sec_testBiljana.dat
!     (symlinks may be created separately for downstream convenience).
!=================================================
! Build / Run notes (Discover / SLES15)
!
! 1) Load Intel + IntelMPI environment (so mpiifort is available)
! module load comp/gcc/11.4.0
! module load comp/intel/2024.2.0
! module load mpi/openmpi/4.1.6/intel-2024.2.0 
! source /home/wjiang/swdev/github/GEOSgcm/@env/g5_modules.sh
!
! 2) Compile against ESMA Baselibs NetCDF (ifort + intelmpi, SLES15):
!
!    BASE=/discover/swdev/gmao_SIteam/Baselibs/ESMA-Baselibs-7.24.0/x86_64-pc-linux-gnu/ifort_2021.6.0-intelmpi_2021.13.0-SLES15/Linux
!
!    mpiifort -O2 -g STEP3_process_soildata.f90 -o STEP3_process_soildata.x \
!      -I$BASE/include/netcdf \
!      -L$BASE/lib \
!      -lnetcdff -lnetcdf \
!      -lhdf5_hl -lhdf5 -lz -lsz -lbz2 \
!     -lmfhdf -ldf \
!      -ljpeg \
!      -lcurl -lssl -lcrypto \
!      -lbrotlidec -lbrotlicommon -lzstd \
!     -lxml2 \
!      -ltirpc \
!      -ldl -lm
! OR for copy paste below:
!mpiifort -O2 -g STEP3_process_soildata.f90 -o STEP3_process_soildata.x -I$BASE/include/netcdf -L$BASE/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lsz -lbz2 -lmfhdf -ldf -ljpeg -lcurl -lssl -lcrypto -lbrotlidec -lbrotlicommon -lzstd -lxml2 -ltirpc -ldl -lm

!    Note: ifort prints a deprecation remark (#10448); it is harmless.
!
! 3) Run (ensure runtime libs are found):
!
!    export LD_LIBRARY_PATH=$BASE/lib:$LD_LIBRARY_PATH
!    ./STEP3_process_soildata.x
!
!=================================================
!  Debug memory issues run: valgrind ./STEP3_process_soildata.x
!=================================================

use netcdf

implicit none

INTEGER            :: i,j,k,ii,c1,c2,irrecs
INTEGER, PARAMETER :: nc_ngdc = 4320, nr_ngdc = 2160,nc=43200,nr=21600
INTEGER (kind=1), allocatable, dimension (:,:) :: gtext,clay,silt,sand,om,poros

character(len=300) :: path, path1, path2, path3, path_regrid

LOGICAL, parameter :: regrid_ngdc = .true. !.false.

    ! ----------------------- !
    !  for gtext file         !
    integer(kind=1), allocatable :: data_bytes(:)         ! Array to hold all bytes
    integer(kind=1), allocatable :: data_reshaped(:,:)   ! Reshaped array (nr_ngdc, 8 + nc_ngdc)
    integer(kind=1), allocatable :: data_array(:,:)      ! Final data array (nr_ngdc, nc_ngdc)
    integer :: ios
    integer :: actual_size
    integer :: jjj, iii
    integer, parameter :: extra_bytes_per_row = 8
    integer, parameter :: bytes_per_row = extra_bytes_per_row + nc_ngdc
    integer, parameter :: expected_data_size = nc_ngdc * nr_ngdc
    integer, parameter :: expected_total_size = expected_data_size + (nr_ngdc * extra_bytes_per_row)

! The goal is to produce top- and sub-layer soil property fields on the 30-arcsec global grid
! (nc=43200, nr=21600). Inputs arrive at different native resolutions and in separate files
! (gtext, clay, silt, sand, OM/SOC, porosity). When needed, NGDC/Reynolds fields are first
! regridded to the 30-arcsec grid, then merged with HWSD on a land-mask basis.
!
! NGDC/Reynolds soil texture input is on a 0.5° × 0.5° grid (4320 × 2160) and is regridded
! to 30 arcsec (43200 × 21600) when regrid_ngdc = .true.

! if there is regridding to be done (regrid_ngdc = .true.), read in the datasets and regrid them
if(regrid_ngdc) then
   path='/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/soil/soil_properties/v2/NGDC-Reynolds/'
   path_regrid = '/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/soil/soil_properties/v2/MERGED-DATA/' 

   allocate (gtext (1:nc_ngdc,1:nr_ngdc))
   allocate (clay  (1:nc_ngdc,1:nr_ngdc))
   allocate (silt  (1:nc_ngdc,1:nr_ngdc))
   allocate (sand  (1:nc_ngdc,1:nr_ngdc))
   allocate (om    (1:nc_ngdc,1:nr_ngdc))
   allocate (poros (1:nc_ngdc,1:nr_ngdc))

   ! Doing Top Layer first. 

   ! set everything to zero
   gtext =0
   clay  =0
   silt  =0
   sand  =0
   om    =0
   poros =0

   ! calculate number of "records" in i-direction
     irrecs = nc_ngdc / 4
     if (4*irrecs /= nc_ngdc) then
        print *, 'Error: nc_ngdc not divisible by 4:', nc_ngdc
        stop
     end if

    ! ----------------------- !
    ! ----------------------- !
    !   For now, gtext file needs a special reading sequence !
    ! ----------------------- !
   open(unit=10, file=trim(path)//'dtex_tp1.bin', form='unformatted', access='stream', &
        status='old', action='read', iostat=ios)

    !      Get File Size      !
    inquire(unit=10, size=actual_size)
    if (actual_size /= expected_total_size) then
        print *, 'Error: Unexpected data size for', trim('dtex_tp1.bin')
        print *, 'Expected:', expected_total_size, 'bytes, Got:', actual_size, 'bytes.'
        close(10)
        stop
    end if
    print *, 'File size is as expected:', actual_size, 'bytes.'

    !    Allocate Arrays      !
    allocate(data_bytes(actual_size))
    allocate(data_reshaped(nr_ngdc, bytes_per_row))
    allocate(data_array(nr_ngdc, nc_ngdc))

    !      Read Data Bytes    !
    read(10) data_bytes
    close(10)

    !    Reshape Data Bytes   !
    ! Each row has (8 + nc_ngdc) bytes
    ! Assign data_bytes to data_reshaped (row-major to column-major)
    do jjj = 1, nr_ngdc
        data_reshaped(jjj, :) = data_bytes( (jjj-1)*bytes_per_row + 1 : jjj*bytes_per_row )
    end do

    !   Extract Data Array    !
    ! Skip the first 8 bytes of each row
    do jjj = 1, nr_ngdc
        data_array(jjj, :) = data_reshaped(jjj, extra_bytes_per_row + 1 : extra_bytes_per_row + nc_ngdc )
    end do

    ! Correct the shape
    do iii =1,nr_ngdc 
        gtext(iii        ,:)=data_array(nr_ngdc:1:-1,iii        )
        gtext(iii+nr_ngdc,:)=data_array(nr_ngdc:1:-1,iii+nr_ngdc)
    end do

    !    Deallocate Arrays    !
    deallocate(data_bytes)
    deallocate(data_reshaped)
    deallocate(data_array)

    ! ----------------------- !
    !  DONE WITH GTEXT! Other parameters are read in a regular way !
    ! ----------------------- !

   open (unit=11, file=trim(path)//'clay_tp1.img',form='unformatted', status='old',&
        access='direct',recl=1,convert = 'little_endian',action='read')
   open (unit=12, file=trim(path)//'silt_tp1.img',form='unformatted', status='old',&
        access='direct',recl=1,convert = 'little_endian',action='read')
   open (unit=13, file=trim(path)//'sand_tp1.img',form='unformatted', status='old',&
        access='direct',recl=1,convert = 'little_endian',action='read')
   open (unit=14, file=trim(path)//'om_top1.img' ,form='unformatted', status='old',&
        access='direct',recl=1,convert = 'little_endian',action='read')
   open (unit=15, file=trim(path)//'por_top1.img',form='unformatted', status='old',&
        access='direct',recl=1,convert = 'little_endian',action='read')

   ! read the files to get clay, silt, sand, poros(ity), om, and gtext
   k =0
   do j=nr_ngdc,1,-1 ! loop over rows
      do i=1,irrecs ! loop over the "records" in i-direction
         k=k+1  ! which record to read
         c1 = (4*i)-3 ! first element to read from this record
         c2 = (4*i) ! last element to read from this record
     !   read (10, rec=k) (gtext(ii,j), ii=c1,c2)
         read (11, rec=k) (clay (ii,j), ii=c1,c2) 
         read (12, rec=k) (silt (ii,j), ii=c1,c2) 
         read (13, rec=k) (sand (ii,j), ii=c1,c2) 
         read (14, rec=k) (om   (ii,j), ii=c1,c2) 
         read (15, rec=k) (poros(ii,j), ii=c1,c2) 
      end do
   end do

   close (11,status='keep')
   close (12,status='keep')
   close (13,status='keep')
   close (14,status='keep')
   close (15,status='keep')   

   ! data for TOP layer are now read. Start regridding 
   ! To do so, we use RegridRaster fucntion (see the bottom of the script). Here,
   ! on the input we provide info on the size of the input (nc_ngdc,nr_ngdc) and 
   ! output (nc,nr) grids, as well as the data itself, and the output filename. 

   call RegridRaster(nc_ngdc,nr_ngdc,nc,nr,gtext,trim(path_regrid)//'PRtestBiljana_gtext_top')
   call RegridRaster(nc_ngdc,nr_ngdc,nc,nr,clay ,trim(path_regrid)//'PRtestBiljana_clay_top' )
   call RegridRaster(nc_ngdc,nr_ngdc,nc,nr,silt ,trim(path_regrid)//'PRtestBiljana_silt_top' )
   call RegridRaster(nc_ngdc,nr_ngdc,nc,nr,sand ,trim(path_regrid)//'PRtestBiljana_sand_top' )
   call RegridRaster(nc_ngdc,nr_ngdc,nc,nr,om   ,trim(path_regrid)//'PRtestBiljana_om_top'   )
   call RegridRaster(nc_ngdc,nr_ngdc,nc,nr,poros,trim(path_regrid)//'PRtestBiljana_poros_top')

   ! Done regridding. Move to SUB layer and repeat the same steps as in TOP layer

! Sub layer
   gtext =0
   clay  =0
   silt  =0
   sand  =0
   om    =0
   poros =0

   irrecs = nc_ngdc / 4
   if (4*irrecs /= nc_ngdc) then
      print *, 'Error: nc_ngdc not divisible by 4:', nc_ngdc
      stop
   end if
   
    ! ----------------------- !
    !   For now, gtext file needs a special reading sequence !
    ! ----------------------- !
    open(unit=10, file=trim(path)//'dtex_sb1.bin', form='unformatted', access='stream', &
         status='old', action='read', iostat=ios)

    !      Get File Size      !
    inquire(unit=10, size=actual_size)
    if (actual_size /= expected_total_size) then
        print *, 'Error: Unexpected data size for', trim('dtex_sb1.bin')
        print *, 'Expected:', expected_total_size, 'bytes, Got:', actual_size, 'bytes.'
        close(10)
        stop
    end if
    print *, 'File size is as expected:', actual_size, 'bytes.'

    !    Allocate Arrays      !
    allocate(data_bytes(actual_size))
    allocate(data_reshaped(nr_ngdc, bytes_per_row))
    allocate(data_array(nr_ngdc, nc_ngdc))

    !      Read Data Bytes    !
    read(10) data_bytes
    close(10)

    print*, "done reading sub gtext"
    !    Reshape Data Bytes   !
    ! Each row has (8 + nc_ngdc) bytes
    ! Assign data_bytes to data_reshaped (row-major to column-major)
    do jjj = 1, nr_ngdc
        data_reshaped(jjj, :) = data_bytes( (jjj-1)*bytes_per_row + 1 : jjj*bytes_per_row )
    end do

    !   Extract Data Array    !
    ! Skip the first 8 bytes of each row
    do jjj = 1, nr_ngdc
        data_array(jjj, :) = data_reshaped(jjj, extra_bytes_per_row + 1 : extra_bytes_per_row + nc_ngdc )
    end do

    ! Correct the shape
    do iii =1,nr_ngdc 
        gtext(iii        ,:)=data_array(nr_ngdc:1:-1,iii        )
        gtext(iii+nr_ngdc,:)=data_array(nr_ngdc:1:-1,iii+nr_ngdc)
    end do

    !    Deallocate Arrays    !
    deallocate(data_bytes)
    deallocate(data_reshaped)
    deallocate(data_array)

    ! ----------------------- !
    !  DONE WITH GTEXT! Other parameters are read in a regular way !
    ! ----------------------- !

   open (unit=11, file=trim(path)//'clay_sb1.img',form='unformatted', status='old',&
        access='direct',recl=1,convert = 'little_endian',action='read')
   open (unit=12, file=trim(path)//'silt_sb1.img',form='unformatted', status='old',&
        access='direct',recl=1,convert = 'little_endian',action='read')
   open (unit=13, file=trim(path)//'sand_sb1.img',form='unformatted', status='old',&
        access='direct',recl=1,convert = 'little_endian',action='read')
   open (unit=14, file=trim(path)//'om_sub1.img' ,form='unformatted', status='old',&
        access='direct',recl=1,convert = 'little_endian',action='read')
   open (unit=15, file=trim(path)//'por_sub1.img',form='unformatted', status='old',&
        access='direct',recl=1,convert = 'little_endian',action='read')

   k =0
   do j=nr_ngdc,1,-1
      do i=1,irrecs
         k=k+1
         c1 = (4*i)-3
         c2 = (4*i)
  !      read (10, rec=k) (gtext(ii,j), ii=c1,c2)
         read (11, rec=k) (clay (ii,j), ii=c1,c2) 
         read (12, rec=k) (silt (ii,j), ii=c1,c2) 
         read (13, rec=k) (sand (ii,j), ii=c1,c2) 
         read (14, rec=k) (om   (ii,j), ii=c1,c2) 
         read (15, rec=k) (poros(ii,j), ii=c1,c2) 
      end do
   end do

   close (11,status='keep')
   close (12,status='keep')
   close (13,status='keep')
   close (14,status='keep')
   close (15,status='keep')   

   call RegridRaster(nc_ngdc,nr_ngdc,nc,nr,gtext,trim(path_regrid)//'PRtestBiljana_gtext_sub')
   call RegridRaster(nc_ngdc,nr_ngdc,nc,nr,clay ,trim(path_regrid)//'PRtestBiljana_clay_sub' )
   call RegridRaster(nc_ngdc,nr_ngdc,nc,nr,silt ,trim(path_regrid)//'PRtestBiljana_silt_sub' )
   call RegridRaster(nc_ngdc,nr_ngdc,nc,nr,sand ,trim(path_regrid)//'PRtestBiljana_sand_sub' )
   call RegridRaster(nc_ngdc,nr_ngdc,nc,nr,om   ,trim(path_regrid)//'PRtestBiljana_om_sub'   )
   call RegridRaster(nc_ngdc,nr_ngdc,nc,nr,poros,trim(path_regrid)//'PRtestBiljana_poros_sub')

endif

! Regridding is performed only if needed (regrid_ngdc = .true.).
! At this point, all datasets are on the same target grid.
! Originally, NGDC/Reynolds data are on nc_ngdc = 4320 × nr_ngdc = 2160 (0.5° grid).
! After regridding, all fields are on nc = 43200 × nr = 21600 (30-arcsec grid).
!
! We can now proceed with merging the variables.
! The merge is mask-free: HWSDv2 values overwrite NGDC/Reynolds values
! wherever HWSDv2 provides valid (non-missing) data; otherwise NGDC is retained.
!
! <<<<< MERGING BEGINS >>>>>
! STEP1: NGDC-HWSD Merger
! LAYER1: NGDC Reynolds Data
! LAYER2: HWSD DATA

! ------------------------------------------------------------
! STEP3 merge setup:
!
!   Layer 1 (NGDC/Reynolds):
!     - Already regridded to 30 arcsec and written to:
!         path_regrid = .../v2/MERGED-DATA/
!     - Files: PRtestBiljana_*_{top,sub}_30sec.dat
!
!   Layer 2 (HWSDv2):
!     - Read directly from NetCDF files in:
!         path2 = .../v2/input_STEP3_also_output_STEP2/
!
!   Merge logic (mask-free):
!     - If HWSDv2 provides a valid (non-missing) value at a grid cell,
!       it overwrites the NGDC/Reynolds value.
!     - If HWSDv2 is missing, the NGDC/Reynolds value is retained.
!     - No external land mask (e.g., hwsd_mask.bin) is used.
!
!   Outputs:
!     - Merged 30-arcsec binaries written to:
!         path3 = .../v2/MERGED-DATA/
! ------------------------------------------------------------
path2  = '/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/' // &
         'land/soil/soil_properties/v2/input_STEP3_also_output_STEP2/'
path3  = '/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/' // &
         'land/soil/soil_properties/v2/MERGED-DATA/'

path1 = path_regrid

! Ready to Merge!!!
! Merging is done in same manner first for the top- and then for sub-layer. 

print *,'Step 1: Merging NGDC & HWSD '
print *,'        Top Layer '
print *,'============================'
print *,'             '

! in case the regridding was performed (regrid_ngdc = .true.)
! each of the regridded datasets (i.e., gtext,clay,silt,sand,om_su,poros)
! needs to be merged.

if(regrid_ngdc) then

  ! call merge_data function (see bottom of the script) and send the regridded file
  ! (e.g., clay_top_30sec.dat) and the "other" file (e.g., clay_top.bin), together 
  ! with grid size (nc,nr) and their paths/filenames. (Mask-free merge: overwrite where HWSDv2 is valid.)

  call merge_data (1.,nc,nr,             &
       path1,path2,path3,                          &
       'PRtestBiljana_sand_top_30sec.dat',           &
       'sand_top.nc4',                             &
       'sand_top_30sec_testBiljana.dat')
  call merge_data (1.,nc,nr,             &
       path1,path2,path3,                          &
       'PRtestBiljana_clay_top_30sec.dat',           & ! dataset-1 to be merged (this is what the regriding part of the code created)
       'clay_top.nc4',                             & ! dataset-2 to be merged (where did this one come from?)
       'clay_top_30sec_testBiljana.dat')                         ! output (merged data) file
      
  call merge_data (1.,nc,nr,             &
       path1,path2,path3,                          &
       'PRtestBiljana_silt_top_30sec.dat',           &
       'silt_top.nc4',                             &
       'silt_top_30sec_testBiljana.dat')
  
  call merge_data (1.0,nc,nr,           &
       path1,path2,path3,                          &
       'PRtestBiljana_om_top_30sec.dat',             &
       'OC_top.nc4',                               &
       'oc_top_30sec_testBiljana.dat')
endif

call merge_data (2.,nc,nr,             &
     path1,path2,path3,                          &
     'PRtestBiljana_gtext_top_30sec.dat',          &
     'gtext_top.nc4',                             &
     'gtext_top_30sec_testBiljana.dat')

! repeat for Sub-layer
print *,'             '
print *,'        Sub Layer '
print *,'============================'
print *,'             '
if(regrid_ngdc) then

  call merge_data (1.,nc,nr,             &
       path1,path2,path3,                          &
       'PRtestBiljana_clay_sub_30sec.dat',           &
       'clay_sub.nc4',                             &
       'clay_sub_30sec_testBiljana.dat')
      
  call merge_data (1.,nc,nr,             &
       path1,path2,path3,                          &
       'PRtestBiljana_silt_sub_30sec.dat',           &
       'silt_sub.nc4',                             &
       'silt_sub_30sec_testBiljana.dat')
  
  call merge_data (1.,nc,nr,             &
       path1,path2,path3,                          &
       'PRtestBiljana_sand_sub_30sec.dat',           &
       'sand_sub.nc4',                             &
       'sand_sub_30sec_testBiljana.dat')
  
  call merge_data (1.0,nc,nr,           &
       path1,path2,path3,                          &
       'PRtestBiljana_om_sub_30sec.dat',             &
       'OC_sub.nc4',                               &
       'oc_sub_30sec_testBiljana.dat')
endif
call merge_data (2.,nc,nr,             &
     path1,path2,path3,                          &
     'PRtestBiljana_gtext_sub_30sec.dat',          &
     'gtext_sub.nc4',                             &
     'gtext_sub_30sec_testBiljana.dat')

! LAYER1: NGDC-HWSD merged data
! LAYER2: STATSGO

! LAYER1: NGDC-HWSD-STATSGO merged data
! LAYER2: CANADA

! LAYER1: NGDC-HWSD-STATSGO-CANADA merged data
! LAYER2: AUSTRALIA


END PROGRAM STEP3_process_soildata

! -----------------------------------------------------------------------------------

subroutine RegridRaster(nx,ny,nc,nr,Rin,path)

  ! this is regridding routine, that dilevers the input data
  ! on a new grid. It can upscale or downscale. All it needs
  ! is the input and output grid size (nx,ny,nc,nr), the data
  ! to be regridded, and the name of the output file.
  integer, intent (in) :: nx,ny,nc,nr
  integer (kind=1),intent(IN)  :: Rin (nx,ny)
  real,allocatable, dimension (:,:) :: Rout
  CHARACTER(LEN=*) :: path
  REAL    :: xx, yy
  integer :: i,j,ii,jj
  allocate (Rout(1:nc,1:nr))
  Rout=0
  xx = size(Rin ,1)/float(size(Rout,1))
  yy = size(Rin ,2)/float(size(Rout,2))

  print *,'Writing ..:',trim(path) ! this line should really go after the do loop

  do j=1,nr
     jj = int((j-1)*yy) + 1
     jj = max(1, min(jj, ny))
     do i=1,nc
        ii = int((i-1)*xx) + 1
        ii = max(1, min(ii, nx))
        Rout(i,j) = Rin(ii,jj)*1.
     end do
  end do

  ! write out the regridded data array
  open (10,file=trim(path)//'_30sec.dat',form='unformatted',status='unknown')
  do j=1,nr
     write (10) Rout(:,j)
  end do
  close (10,status='keep')
  deallocate (Rout)
end subroutine RegridRaster

! -----------------------------------------------------------------------------------

subroutine merge_data(sf, nc, nr, path1, path2, path3, file1, file2, file3)
  use netcdf
  use ieee_arithmetic  ! For ieee_is_nan function
  implicit none

  ! Inputs
  integer, intent(in) :: nc, nr
  real, intent(in) :: sf
  character(len=300), intent(in) :: path1, path2, path3
  character(len=*), intent(in) :: file1, file2, file3

  ! Local variables
  integer :: i, j, retval, ncid, varid, lat_len, lon_len
  logical :: is_netcdf
  real, dimension(13) :: HWSD_CLASS
  real, allocatable :: data1(:, :)
  real(4), allocatable :: data2(:, :)
  character(len=100) :: variable_name
  character(len=100) :: base_file2
  character(len=128) :: errmsg
  real(4) :: fill_value
  integer :: ndims
  integer, allocatable :: dimids(:)
  character(len=NF90_MAX_NAME) :: dim_name
  integer :: dim_len
  integer(8) :: count_positive
  integer :: iostat

  ! Initialize HWSD_CLASS
  data HWSD_CLASS /12., 11., 12., 8., 9., 5., 4., 10., 6., 7., 3., 2., 1./

  ! Allocate arrays
  allocate(data1(nc, nr))
  allocate(data2(nc, nr))  ! Matching dimensions to the global grid

  ! Extract base file name for file2
  call get_base_file_name(trim(file2), base_file2)
  print *, "Base file name extracted: ", trim(base_file2)//'END'

  ! Assign variable name based on base file name
  if (base_file2 == 'sand_top.nc4' .or. base_file2 == 'sand_sub.nc4') then
     variable_name = 'Sand'
  else if (base_file2 == 'silt_top.nc4' .or. base_file2 == 'silt_sub.nc4') then
     variable_name = 'Silt'
  else if (base_file2 == 'clay_top.nc4' .or. base_file2 == 'clay_sub.nc4') then
     variable_name = 'Clay'
  else if (base_file2 == 'OC_top.nc4'   .or. base_file2 == 'OC_sub.nc4'  ) then
     variable_name = 'SOC'
  else if (base_file2 == 'gtext_top.nc4' .or. base_file2 == 'gtext_sub.nc4') then
     variable_name = 'Coarse'
  else
     print *, "Error: Unknown file name for variable assignment:", trim(base_file2)
     stop
  end if

  print *, "Variable name assigned: <", trim(variable_name), ">"

  ! Read data1 from file1 (binary)
  open(unit=10, file=trim(path1)//trim(file1), form="unformatted", &
       status="old", action="read")
  do j = 1, nr
     read(10) data1(:, j)
  end do
  close(10)
  ! Harmonize units for SOC merge: Reynolds file is OM; convert to OC so output is OC everywhere
  if (variable_name == 'SOC') then
    if (index(trim(file1), '_om_') > 0) then
      where (data1 >= 0.0) data1 = data1 / 1.72
       print *, 'INFO: Converted layer-1 OM->OC for file1=', trim(file1)  
    end if
  end if  

  ! Debugging: Print sample value from data1
  print *, "Sample data1(1, 1): ", data1(1, 1)

  ! Open NetCDF file
  retval = nf90_open(trim(path2)//trim(file2), nf90_nowrite, ncid)
  if (retval /= nf90_noerr) then
     errmsg = nf90_strerror(retval)
     print *, "Error opening NetCDF file: ", trim(errmsg)
     stop
  end if
  print *, "Detected NetCDF format for file2."

  ! Attempt to inquire about the variable using nf90_inq_varid
  retval = nf90_inq_varid(ncid, trim(variable_name), varid)
  if (retval /= nf90_noerr) then
     errmsg = nf90_strerror(retval)
     print *, "Error: Variable <", trim(variable_name), "> not found in NetCDF file."
     print *, "NetCDF error code: ", retval
     print *, "NetCDF error message: ", trim(errmsg)
     stop
  end if
  print *, "Successfully found variable: <", trim(variable_name), ">, varid: ", varid

  ! Inquire about the number of dimensions
  retval = nf90_inquire_variable(ncid, varid, ndims=ndims)
  if (retval /= nf90_noerr) then
     errmsg = nf90_strerror(retval)
     print *, "Error in nf90_inquire_variable for varid=", varid
     print *, "NetCDF error message: ", trim(errmsg)
     stop
  end if

  ! Allocate dimids array based on ndims
  allocate(dimids(ndims))

  ! Inquire about dimension IDs
  retval = nf90_inquire_variable(ncid, varid, dimids=dimids)
  if (retval /= nf90_noerr) then
     errmsg = nf90_strerror(retval)
     print *, "Error in nf90_inquire_variable for varid=", varid
     print *, "NetCDF error message: ", trim(errmsg)
     stop
  end if

  if (ndims /= 2) then
     print *, "Error: Expected variable to have 2 dimensions, but got ", ndims
     stop
  end if

  ! Get dimensions lengths and names
  retval = nf90_inquire_dimension(ncid, dimids(1), name=dim_name, len=dim_len)
  if (retval /= nf90_noerr) then
     errmsg = nf90_strerror(retval)
     print *, "Error in nf90_inquire_dimension for dimid=", dimids(1)
     print *, "NetCDF error message: ", trim(errmsg)
     stop
  end if
  print *, "Dimension 1: name=", trim(dim_name), " length=", dim_len

  ! Assign lengths based on dimension names
  if (trim(dim_name) == 'lat') then
     lat_len = dim_len
  else if (trim(dim_name) == 'lon') then
     lon_len = dim_len
  else
     print *, "Error: Unexpected dimension name: ", trim(dim_name)
     stop
  end if

  retval = nf90_inquire_dimension(ncid, dimids(2), name=dim_name, len=dim_len)
  if (retval /= nf90_noerr) then
     errmsg = nf90_strerror(retval)
     print *, "Error in nf90_inquire_dimension for dimid=", dimids(2)
     print *, "NetCDF error message: ", trim(errmsg)
     stop
  end if
  print *, "Dimension 2: name=", trim(dim_name), " length=", dim_len

  ! Assign lengths based on dimension names
  if (trim(dim_name) == 'lat') then
     lat_len = dim_len
  else if (trim(dim_name) == 'lon') then
     lon_len = dim_len
  else
     print *, "Error: Unexpected dimension name: ", trim(dim_name)
     stop
  end if

  ! Check if dimensions match
  if (lat_len /= nr) then
     print *, "Error: Latitude dimension mismatch. Expected ", nr, " but got ", lat_len
     stop
  end if
  if (lon_len /= nc) then
     print *, "Error: Longitude dimension mismatch. Expected ", nc, " but got ", lon_len
     stop
  end if

  print *, "Dimension checks passed."
  print *, "Attempting to read variable: ", trim(variable_name)

  print *, "Expected dimensions: nc=", nc, " nr=", nr
  print *, "NetCDF dimensions: lon_len=", lon_len, " lat_len=", lat_len

  ! Try to get the 'missing_value' attribute
  retval = nf90_get_att(ncid, varid, 'missing_value', fill_value)
  if (retval /= nf90_noerr) then
     ! If 'missing_value' is not found, default to -9999.0
     fill_value = -9999.0
     print *, "No 'missing_value' attribute found; defaulting fill_value to ", fill_value
  else
     print *, "'missing_value' attribute found; fill_value set to ", fill_value
  end if

  ! Read the variable into data2
  retval = nf90_get_var(ncid, varid, data2)
  if (retval /= nf90_noerr) then
     errmsg = nf90_strerror(retval)
     print *, "Error reading variable from NetCDF file: ", trim(variable_name)
     print *, "NetCDF error message: ", trim(errmsg)
     stop
  end if
  print *, "Successfully read variable: ", trim(variable_name)

  ! Close the NetCDF file
  retval = nf90_close(ncid)
  if (retval /= nf90_noerr) then
     print *, "Error closing NetCDF file."
     stop
  end if

  ! Deallocate dimids array
  deallocate(dimids)

  ! Debugging: Print a sample value from data2
  print *, "Sample data2(1, 1): ", data2(1, 1)

  ! Replace missing values in data2 with NaN
  do j = 1, nr
     do i = 1, nc
        if (data2(i, j) == fill_value) data2(i, j) = ieee_value(0.0, ieee_quiet_nan)
     end do
  end do

  ! Merge data
  count_positive = 0_8

   do j = 1, nr
     do i = 1, nc
       if (.not. ieee_is_nan(data2(i,j))) then
         count_positive = count_positive + 1_8
   
         if (sf == 2.0) then
           data1(i,j) = HWSD_CLASS(int(data2(i,j)))
         else
           data1(i,j) = data2(i,j) * sf
         end if
       end if
     end do
   end do

  print *, "Number of valid data2 values processed: ", count_positive

  

  ! Debugging: Print a sample merged value
  print *, "Sample merged data1(1, 1): ", data1(1, 1)

  ! Write the merged data array to file3
  open(unit=10, file=trim(path3)//trim(file3), form="unformatted", &
       status="unknown", action="write", iostat=iostat)
  if (iostat /= 0) then
     print *, "Error opening file for writing: ", trim(path3)//trim(file3)
     stop
  end if

  do j = 1, nr
     write(10, iostat=iostat) data1(:, j)
     if (iostat /= 0) then
        print *, "Error writing to file at row ", j
        stop
     end if
  end do
  close(10)

  print *, "Merged data written to file: ", trim(file3)

end subroutine merge_data

subroutine get_base_file_name(file_path, base_name)
  implicit none
  character(len=*), intent(in) :: file_path
  character(len=100), intent(out) :: base_name
  integer :: slash_pos

  ! Find the last slash
  slash_pos = len_trim(file_path)
  do while (slash_pos > 0 .and. file_path(slash_pos:slash_pos) /= '/')
     slash_pos = slash_pos - 1
  end do

  ! Extract the base file name
  if (slash_pos > 0) then
     base_name = file_path(slash_pos+1:)
  else
     base_name = file_path
  end if
  base_name = trim(base_name)
end subroutine get_base_file_name

