PROGRAM STEP4_process_soildata

! Legacy starting code was written by Sarith Mahanama since then code evolved with updated data sets, file formats and bug fixes.  

! Differs from process_soildata.f90 in several ways. Some include:
! - Uses different version of data HWSDv1.21 instead of HWSD
! - Does not read in pro_ files (porosity?!)
! - Populates with missing data in the main part of the code rather than in subroutine
! - Current workflow uses merge_min_statsgo/merge_carbon_statsgo/merge_min/merge_carbon.
! - has a number of conditions when merging and/or regriding. 

!=================================================
! Build / Run notes (Discover / SLES15)
!
! This program uses NetCDF-Fortran (`use netcdf`). On Discover, avoid nf-config from
! GEOSpyD (can point to a netcdf.mod built by a different compiler). Use ESMA Baselibs
! NetCDF built with Intel Fortran instead.
!
! Compile:
!   BASE=/discover/swdev/gmao_SIteam/Baselibs/ESMA-Baselibs-7.24.0/x86_64-pc-linux-gnu/ifort_2021.6.0-intelmpi_2021.13.0-SLES15/Linux
!
!   mpiifort -O2 -g STEP4_process_soildata.f90 -o STEP4_process_soildata.x \
!     -I$BASE/include/netcdf \
!     -L$BASE/lib \
!     -lnetcdff -lnetcdf \
!     -lhdf5_hl -lhdf5 -lz -lsz -lbz2 \
!     -lmfhdf -ldf -ljpeg \
!     -lcurl -lssl -lcrypto \
!     -lbrotlidec -lbrotlicommon -lzstd \
!     -lxml2 -ltirpc \
!     -ldl -lm

! one line to use: 
!mpiifort -O2 -g STEP4_process_soildata.f90 -o STEP4_process_soildata.x -I$BASE/include/netcdf -L$BASE/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lsz -lbz2 -lmfhdf -ldf -ljpeg -lcurl -lssl -lcrypto -lbrotlidec -lbrotlicommon -lzstd -lxml2 -ltirpc -ldl -lm

!
! Run:
!   export LD_LIBRARY_PATH=$BASE/lib:$LD_LIBRARY_PATH
!   ./process_soildata.x
!
! Inputs:
!   - STEP3 merged binaries:  v2/MERGED-DATA/*_{top,sub}_30sec*.dat
!   - STATSGO2:               .../STATSGO2/statsgo_*
!   - AFSIS:                  .../Africa-AFSIS/afsis_*
!
! Outputs:
!   - v2/HWSDv2-NGDC-STATSGO/              (+ data_sources.msk)
!   - v2/HWSDv2-NGDC-STATSGO-AFSIS/        (+ data_sources.msk)
!=================================================

implicit none
INTEGER, PARAMETER :: nc_ngdc = 4320, nr_ngdc = 2160,nc=43200,nr=21600
INTEGER (kind=1), allocatable, dimension (:,:) :: gtext,clay,silt,sand,om,poros
INTEGER , allocatable, dimension (:,:) :: src_mask
INTEGER :: i,j,k,ii,c1,c2,irrecs
LOGICAL, parameter :: regrid_ngdc = .false.
CHARACTER*300 :: path
CHARACTER*300 :: path1,path2,path3
logical, parameter :: ver_merge = .true.
integer :: layer_id ! (1=>HWSD; 2=> STATSGO)

! This NGDC part wont be used below, so the flag is set to .false.
! Reading and Regridding NGDC Reynolds data
if(regrid_ngdc) then
   path='/discover/nobackup/smahanam/NGDC-Reynolds/'
   allocate (gtext (1:nc_ngdc,1:nr_ngdc))
   allocate (clay  (1:nc_ngdc,1:nr_ngdc))
   allocate (silt  (1:nc_ngdc,1:nr_ngdc))
   allocate (sand  (1:nc_ngdc,1:nr_ngdc))
   allocate (om    (1:nc_ngdc,1:nr_ngdc))
   allocate (poros (1:nc_ngdc,1:nr_ngdc))

   ! Top Layer
   gtext =0
   clay  =0
   silt  =0
   sand  =0
   om    =0
   poros =0
   irrecs = nint (nc_ngdc / 4.0) 
   open (unit=10, file=trim(path)//'dtex_tp1.img',form='unformatted', status='old',&
        access='direct',recl=1,convert = 'little_endian',action='read')
   open (unit=11, file=trim(path)//'clay_tp1.img',form='unformatted', status='old',&
        access='direct',recl=1,convert = 'little_endian',action='read')
   open (unit=12, file=trim(path)//'silt_tp1.img',form='unformatted', status='old',&
        access='direct',recl=1,convert = 'little_endian',action='read')
   open (unit=13, file=trim(path)//'sand_tp1.img',form='unformatted', status='old',&
        access='direct',recl=1,convert = 'little_endian',action='read')
   open (unit=14, file=trim(path)//'om_top1.img' ,form='unformatted', status='old',&
        access='direct',recl=1,convert = 'little_endian',action='read')
!   open (unit=15, file=trim(path)//'por_top1.img',form='unformatted', status='old',&
!        access='direct',recl=1,convert = 'little_endian',action='read')

   k =0
   do j=nr_ngdc,1,-1
      do i=1,irrecs
         k=k+1
         c1 = (4*i)-3
         c2 = (4*i)
         read (10, rec=k) (gtext(ii,j), ii=c1,c2)
         read (11, rec=k) (clay (ii,j), ii=c1,c2) 
         read (12, rec=k) (silt (ii,j), ii=c1,c2) 
         read (13, rec=k) (sand (ii,j), ii=c1,c2) 
         read (14, rec=k) (om   (ii,j), ii=c1,c2) 
!         read (15, rec=k) (poros(ii,j), ii=c1,c2) 
!         print *, gtext(ii,j),sand (ii,j),clay (ii,j)
         do ii = c1,c2
            if((gtext(ii,j) ==0).or.((sand (ii,j) + clay (ii,j)) == 0.)) then
               gtext(ii,j) = -99
               clay (ii,j) = -99
               silt (ii,j) = -99
               sand (ii,j) = -99
               om   (ii,j) = -99
!            pause
            endif
         end do
      end do
   end do

   close (10,status='keep')
   close (11,status='keep')
   close (12,status='keep')
   close (13,status='keep')
   close (14,status='keep')
   close (15,status='keep')   

   call RegridRaster(nc_ngdc,nr_ngdc,nc,nr,gtext,trim(path)//'gtext_top')
   call RegridRaster(nc_ngdc,nr_ngdc,nc,nr,clay ,trim(path)//'clay_top' )
   call RegridRaster(nc_ngdc,nr_ngdc,nc,nr,silt ,trim(path)//'silt_top' )
   call RegridRaster(nc_ngdc,nr_ngdc,nc,nr,sand ,trim(path)//'sand_top' )
   call RegridRaster(nc_ngdc,nr_ngdc,nc,nr,om   ,trim(path)//'om_top'   )
!   call RegridRaster(nc_ngdc,nr_ngdc,nc,nr,poros,trim(path)//'poros_top')

! Sub layer
   gtext =0
   clay  =0
   silt  =0
   sand  =0
   om    =0
   poros =0
   irrecs = nint (nc_ngdc / 4.0) 
   open (unit=10, file=trim(path)//'dtex_sb1.img',form='unformatted', status='old',&
        access='direct',recl=1,convert = 'little_endian',action='read')
   open (unit=11, file=trim(path)//'clay_sb1.img',form='unformatted', status='old',&
        access='direct',recl=1,convert = 'little_endian',action='read')
   open (unit=12, file=trim(path)//'silt_sb1.img',form='unformatted', status='old',&
        access='direct',recl=1,convert = 'little_endian',action='read')
   open (unit=13, file=trim(path)//'sand_sb1.img',form='unformatted', status='old',&
        access='direct',recl=1,convert = 'little_endian',action='read')
   open (unit=14, file=trim(path)//'om_sub1.img' ,form='unformatted', status='old',&
        access='direct',recl=1,convert = 'little_endian',action='read')
!   open (unit=15, file=trim(path)//'por_sub1.img',form='unformatted', status='old',&
!        access='direct',recl=1,convert = 'little_endian',action='read')

   k =0
   do j=nr_ngdc,1,-1
      do i=1,irrecs
         k=k+1
         c1 = (4*i)-3
         c2 = (4*i)
         read (10, rec=k) (gtext(ii,j), ii=c1,c2)
         read (11, rec=k) (clay (ii,j), ii=c1,c2) 
         read (12, rec=k) (silt (ii,j), ii=c1,c2) 
         read (13, rec=k) (sand (ii,j), ii=c1,c2) 
         read (14, rec=k) (om   (ii,j), ii=c1,c2) 
!         read (15, rec=k)(poros(ii,j), ii=c1,c2) 
         do ii = c1,c2
            if((gtext(ii,j) ==0).or.((sand (ii,j) + clay (ii,j)) == 0.))  then
              gtext(ii,j) = -99
              clay (ii,j) = -99
              silt (ii,j) = -99
              sand (ii,j) = -99
              om   (ii,j) = -99
             endif
         end do
      end do
   end do

   close (10,status='keep')
   close (11,status='keep')
   close (12,status='keep')
   close (13,status='keep')
   close (14,status='keep')
   close (15,status='keep')   

   call RegridRaster(nc_ngdc,nr_ngdc,nc,nr,gtext,trim(path)//'gtext_sub')
   call RegridRaster(nc_ngdc,nr_ngdc,nc,nr,clay ,trim(path)//'clay_sub' )
   call RegridRaster(nc_ngdc,nr_ngdc,nc,nr,silt ,trim(path)//'silt_sub' )
   call RegridRaster(nc_ngdc,nr_ngdc,nc,nr,sand ,trim(path)//'sand_sub' )
   call RegridRaster(nc_ngdc,nr_ngdc,nc,nr,om   ,trim(path)//'om_sub'   )
!   call RegridRaster(nc_ngdc,nr_ngdc,nc,nr,poros,trim(path)//'poros_sub')

endif

!
! <<<<< OVERLAYING BEGINS >>>>>
! >>>>> WE DON'T USE NGDC ANY MORE <<<<<
! STEP1: Overlay HWSD on NGDC
! LAYER1: NGDC Reynolds Data
! LAYER2: HWSDv2
!
! STEP 2: Overlay Statsgo on HWSDv2
! LAYER1: HWSDv2
! LAYER2: STATSGO2

allocate (src_mask(1:nc,1:nr))
src_mask = 0

! NGDC+HWSD_v2 data (have to be created first; and they are)
! combined HWSDv2 and reynolds
path1 ='/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/soil/soil_properties/v2/MERGED-DATA/' 
! STATSGO file 
path2 ='/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/soil/soil_properties/v2/STATSGO2/statsgo_'
! Output file
path3 ='/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/soil/soil_properties/v2/HWSDv2-NGDC-STATSGO/'

print *,'Step 1: Overlaying STATSGO on HWSDv1.21 while vertical merging '
print *,'==============================================================='
print *,'             '

call merge_min_statsgo   (.true., ver_merge,1, nc,nr,src_mask,    &
     path1,path2,path3)

print *,maxval(src_mask),minval(src_mask)

call merge_carbon_statsgo(.true., ver_merge,1, nc,nr,src_mask,    &
     path1,path2,path3)
print *,maxval(src_mask),minval(src_mask)

  open (10,file=trim(path3)//'data_sources.msk',form='unformatted', &
       status='unknown',action='write')
  do j = 1,nr
       write (10) (src_mask (i,j),i =1,nc) 
  end do

  close (10, status ='keep')

! Data generated in the step above
path1 ='/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/soil/soil_properties/v2/HWSDv2-NGDC-STATSGO/'

! AFSIS data 
path2 ='/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/soil/soil_properties/v2/Africa-AFSIS/afsis_'

! output file
path3 ='/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/soil/soil_properties/v2/HWSDv2-NGDC-STATSGO-AFSIS/'

print *,'Step 2: Overlaying AFSIS on HWSDv1.21-STATSGO while vertical merging '
print *,'====================================================================='
print *,'             '


call merge_min   (.false., ver_merge, 2, nc,nr,src_mask,    &
     path1,path2,path3)

print *,maxval(src_mask),minval(src_mask)

call merge_carbon(.false., ver_merge, 2, nc,nr,src_mask,    &
     path1,path2,path3)
print *,maxval(src_mask),minval(src_mask)
  
  open (10,file=trim(path3)//'data_sources.msk',form='unformatted', &
       status='unknown',action='write')
  do j = 1,nr
       write (10) (src_mask (i,j),i =1,nc) 
  end do

  close (10, status ='keep')

END PROGRAM STEP4_process_soildata
!
! -----------------------------------------------------------------------------------
!
subroutine RegridRaster(nx,ny,nc,nr,Rin,path)
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

  print *,'Writing ..:',trim(path)

  do j=1,size(Rout,2)
     jj = (j-1)*yy + 1
     do i=1,size(Rout,1)
        ii = (i-1)*xx + 1
        Rout(i,j) = Rin(ii,jj)
        if (Rin(ii,jj) == -99) Rout(i,j) = -9999.
     end do
  end do

  open (10,file=trim(path)//'_30sec.dat',form='unformatted',status='unknown')
  do j=1,nr
     write (10) Rout(:,j)
  end do
  close (10,status='keep')
  deallocate (Rout)
end subroutine RegridRaster

! -----------------------------------------------------------------------------------
! to accomodate WHSD-v2 input format (in .nc) a separate [copy] of merging subroutines 
! is created and use in Step-1, although the mathodology of merging is same in both
! _statsgo and _afsis
subroutine merge_min_statsgo (self_fill,ver_merge,layer_id,nc,nr,src_mask, &
  path1,path2,path3)

  use netcdf
  use, intrinsic :: ieee_arithmetic  ! Import IEEE module for NaN detection

  implicit none
  integer,  intent(in) :: nc,nr,layer_id
  logical, intent (in) :: self_fill,ver_merge
  character(len=300),intent (in) :: path1,path2,path3
  character(len=100) :: file1,file2,file3,file4
  character(len=100) :: file1_bin,file2_bin,file3_bin,file4_bin
  character(len=100) :: file5,file6,file7,file8
  integer, intent (inout), dimension (nc,nr) :: src_mask

  integer :: i,j,lp1

  real, allocatable, dimension (:) :: data1, data2, data3, data4
  real, allocatable, dimension (:) :: data5, data6, data7, data8
  real :: sum_top, sum_sub

  integer            :: retval,ncid1,ncid2,ncid3,ncid4,varid1,varid2,varid3,varid4
  real(4)            :: fill_value1,fill_value2,fill_value3,fill_value4

  allocate (data1(1:nc))
  allocate (data2(1:nc))
  allocate (data3(1:nc))
  allocate (data4(1:nc))
  allocate (data5(1:nc))
  allocate (data6(1:nc))
  allocate (data7(1:nc))
  allocate (data8(1:nc))  

  file1 = 'clay_top_30sec.dat'
  file2 = 'sand_top_30sec.dat'
  file3 = 'clay_sub_30sec.dat'
  file4 = 'sand_sub_30sec.dat'

  ! note: for HWDS_v2-NGDC files, we are reading binary format
  ! These files need to be generated first!
  file1_bin = 'clay_top_30sec_testBiljana.dat'
  file2_bin = 'sand_top_30sec_testBiljana.dat'
  file3_bin = 'clay_sub_30sec_testBiljana.dat'
  file4_bin = 'sand_sub_30sec_testBiljana.dat'

  print *,'file1 :',trim(path1)//trim(file1)
  print *,'file2 :',trim(path2)//trim(file1)
  print *,'file3 :',trim(path3)//trim(file1)
  print *,'file1 :',trim(path1)//trim(file1_bin)
  print *,'file2 :',trim(path2)//trim(file1_bin)
  print *,'file3 :',trim(path3)//trim(file1_bin)

  lp1 = layer_id + 1

  if(self_fill) print *,'Vertical self-filling of first layer'
  if(ver_merge)  print *,'Vertical MERGING layers'

  ! These are HWSDv2-NGDC data
  open (10,file=trim(path1)//trim(file1_bin),form='unformatted', &
       status='old',action='read')
  open (11,file=trim(path1)//trim(file2_bin),form='unformatted', &
       status='old',action='read')
  open (12,file=trim(path1)//trim(file3_bin),form='unformatted', &
       status='old',action='read')
  open (13,file=trim(path1)//trim(file4_bin),form='unformatted', &
       status='old',action='read')

  ! These are STATSFO data
  open (20,file=trim(path2)//trim(file1),form='unformatted', &
       status='old',action='read')
  open (21,file=trim(path2)//trim(file2),form='unformatted', &
       status='old',action='read')
  open (22,file=trim(path2)//trim(file3),form='unformatted', &
       status='old',action='read')
  open (23,file=trim(path2)//trim(file4),form='unformatted', &
       status='old',action='read')

  open (30,file=trim(path3)//trim(file1),form='unformatted', &
       status='unknown',action='write')
  open (31,file=trim(path3)//trim(file2),form='unformatted', &
       status='unknown',action='write')
  open (32,file=trim(path3)//trim(file3),form='unformatted', &
       status='unknown',action='write')
  open (33,file=trim(path3)//trim(file4),form='unformatted', &
       status='unknown',action='write')

 
  do j=1,nr

       ! These are HWSD-v2-NGDC data
     read (10) (data1 (i),i=1,nc)
     read (11) (data2 (i),i=1,nc)
     read (12) (data3 (i),i=1,nc)
     read (13) (data4 (i),i=1,nc)
    
       ! Replace missing values (if any) with -9999.0
     do i = 1, nc
        if (ieee_is_nan(data1(i))) data1(i) = -9999.0
        if (ieee_is_nan(data2(i))) data2(i) = -9999.0
        if (ieee_is_nan(data3(i))) data3(i) = -9999.0
        if (ieee_is_nan(data4(i))) data4(i) = -9999.0
     end do
    
     ! These are STATSFO data
     read (20) (data5 (i),i=1,nc)
     read (21) (data6 (i),i=1,nc)
     read (22) (data7 (i),i=1,nc)
     read (23) (data8 (i),i=1,nc)

     do i=1,nc

	sum_top = data1(i) + data2(i)
        sum_sub = data3(i) + data4(i)

	if(layer_id == 1) then
	   if(sum_top >= 0.) src_mask(i,j) = 1
	   if(sum_sub >= 0.) src_mask(i,j) = src_mask(i,j) + 10 
	endif

        if(ver_merge) then
           if(self_fill) then
              if((sum_top >=0.).and.(sum_sub < 0.)) then                
                 data3 (i) = data1 (i)
                 data4 (i) = data2 (i)
	         src_mask(i,j) = src_mask(i,j) + 80
              endif
              
              if((sum_top < 0.).and.(sum_sub >= 0.)) then
                 data1 (i) = data3 (i)
                 data2 (i) = data4 (i)
	         src_mask(i,j) = src_mask(i,j) + 8                
              endif
           endif
        endif
        
        sum_top = data5(i) + data6(i)
        sum_sub = data7(i) + data8(i)
          
        if((sum_top >=0.).and.(sum_sub >=0.)) then !SM: remember to switch
           data1 (i) = data5 (i)
           data2 (i) = data6 (i)
           data3 (i) = data7 (i)
           data4 (i) = data8 (i)
	   src_mask(i,j) = 100*(src_mask(i,j)/100) + 10*lp1 + lp1
        endif
        
        if(ver_merge) then 
           if((sum_top >=0.).and.(sum_sub < 0.)) then
              data1 (i) = data5 (i)
              data2 (i) = data6 (i)
              data3 (i) = data5 (i)  
              data4 (i) = data6 (i)  
	      src_mask(i,j) = 100*(src_mask(i,j)/100) + 80 + lp1         
           endif
           
           if((sum_top < 0.).and.(sum_sub >= 0.)) then
              data1 (i) = data7 (i)  
              data2 (i) = data8 (i)  
              data3 (i) = data7 (i)
              data4 (i) = data8 (i)           
              src_mask(i,j) = 100*(src_mask(i,j)/100) + 10*lp1 + 8
           endif
        endif
!	if((data1(i) >= 0.).and.(src_mask(i,j) ==0.)) &
!	   print *,i,j,data1(i),data2(i),data3(i),data4(i),src_mask(i,j)
!	if((data1(i) < 0.).and.(src_mask(i,j) > 0.)) then
!	   print *,i,j,data1(i),data2(i),data3(i),data4(i),src_mask(i,j), &
!	          data5(i),data6(i),data7(i),data8(i)
!	  pause
!	endif
     end do

     write (30) (data1 (i),i=1,nc)
     write (31) (data2 (i),i=1,nc)
     write (32) (data3 (i),i=1,nc)
     write (33) (data4 (i),i=1,nc)
  end do

  deallocate (data1, data2, data3, data4)
  deallocate (data5, data6, data7, data8)

    ! These are HWSD-v2-NGDC files
  close (10,status='keep')
  close (11,status='keep')
  close (12,status='keep')
  close (13,status='keep')

    ! These are STATSFO files
  close (20,status='keep')
  close (21,status='keep')
  close (22,status='keep')
  close (23,status='keep')

    ! These are output WHSD-v2-STATSFO files
  close (30,status='keep')
  close (31,status='keep')
  close (32,status='keep')
  close (33,status='keep')

end subroutine merge_min_statsgo

! -----------------------------------------------------------------------------------

subroutine merge_carbon_statsgo (self_fill,ver_merge, layer_id,nc,nr,src_mask, &
  path1,path2,path3)

  use netcdf
  use ieee_arithmetic  ! For ieee_is_nan function
  implicit none
  integer,  intent(in) :: nc,nr
  logical,  intent(in) :: self_fill, ver_merge
  integer,  intent(in) :: layer_id
  real                 :: SUB2TOP 
  character(len=300),intent (in) :: path1,path2,path3
  character(len=100) :: file1,file2,file1_bin,file2_bin

  integer :: i,j, lp1,min_dig
  real, allocatable, dimension (:) :: data1, data2, data3, data4
  integer, intent (inout), dimension (nc,nr) :: src_mask
  real :: sum_top, sum_sub

 !character(len=100) :: variable_name, file1_whds, file2_whds
  integer            :: retval,ncid1,ncid2,varid1,varid2
  real(4)            :: fill_value1,fill_value2

  allocate (data1(1:nc))
  allocate (data2(1:nc))
  allocate (data3(1:nc))
  allocate (data4(1:nc))

  file1 = 'oc_top_30sec.dat'
  file2 = 'oc_sub_30sec.dat'

  file1_bin = 'oc_top_30sec_testBiljana.dat'
  file2_bin = 'oc_sub_30sec_testBiljana.dat'

  print *,'file1 :',trim(path1)//trim(file1)
  print *,'file2 :',trim(path2)//trim(file1)
  print *,'file2 :',trim(path3)//trim(file1)

  if(self_fill) print *,'Vertical self-filling of first layer'

  ! These are HWSD-v2-NGDC data
  open (10,file=trim(path1)//trim(file1_bin),form='unformatted', &
        status='old',action='read')
  open (11,file=trim(path1)//trim(file2_bin),form='unformatted', &
        status='old',action='read')
 
  ! These are STATSFO data
  open (20,file=trim(path2)//trim(file1),form='unformatted', &
       status='old',action='read')
  open (21,file=trim(path2)//trim(file2),form='unformatted', &
       status='old',action='read')

  open (30,file=trim(path3)//trim(file1),form='unformatted', &
       status='unknown',action='write')
  open (31,file=trim(path3)//trim(file2),form='unformatted', &
       status='unknown',action='write')
	
  lp1 = layer_id + 1 

  do j=1,nr

       ! These are HWSD-v2 data
     read (10) (data1 (i),i=1,nc)
     read (11) (data2 (i),i=1,nc)
    
     ! If needed, convert om into oc (apply factor to match units): [om] = 1.72 * [oc] => [oc] = [om]/1.72
     !data1=data1/1.72
     !data2=data2/1.72

       ! Replace missing values (if any) with -9999.0
     do i = 1, nc
        if (ieee_is_nan(data1(i))) data1(i) = -9999.0
        if (ieee_is_nan(data2(i))) data2(i) = -9999.0
     end do
    
     ! These are STATSFO data
     read (20) (data3 (i),i=1,nc)
     read (21) (data4 (i),i=1,nc)
	
     do i=1,nc

	if (layer_id == 1) then
           if(data1 (i) >=0.) data1(i) = data1(i) 
           if(data2 (i) >=0.) data2(i) = data2(i) 
        endif

        sum_top = data1(i) 
        sum_sub = data2(i) 
	SUB2TOP = 1.

	if(layer_id == 1) then
	   if(sum_top >= 0.) src_mask(i,j) = src_mask (i,j) +  100
	   if(sum_sub >= 0.) src_mask(i,j) = src_mask (i,j) + 1000 
	endif
	
        if(ver_merge) then
           if(self_fill) then
              if((sum_top >=0.).and.(sum_sub < 0.)) then
	         if(data1(i) < 15./1.72) SUB2TOP = 0.3
                 data2 (i) = data1 (i) * SUB2TOP
	         src_mask(i,j) = src_mask (i,j) + 8000 
              endif
              
              if((sum_top < 0.).and.(sum_sub >= 0.)) then
                 if(data2(i) < 15./1.72) SUB2TOP = 0.3
                 data1 (i) = data2 (i) / SUB2TOP
                 src_mask(i,j) = src_mask (i,j) + 800 
              endif              
           endif
        endif

        sum_top = data3(i) 
        sum_sub = data4(i) 
        SUB2TOP = 1.
	min_dig = src_mask (i,j) - 100*(src_mask (i,j)/100)

        if((sum_top >=0.).and.(sum_sub >=0.)) then  !SM: switch to and
           data1 (i) = data3 (i) 
           data2 (i) = data4 (i) 
           src_mask(i,j) =  1000*lp1 + 100*lp1 + min_dig
           if(src_mask(i,j) == 3344) print *,i,j,sum_top,sum_sub,min_dig,lp1,layer_id
        endif

        if(ver_merge) then
           if((sum_top >=0.).and.(sum_sub < 0.)) then
	      if(data3(i) < 15./1.72) SUB2TOP = 0.3
              data1 (i) = data3 (i) 
              data2 (i) = data1 (i) * SUB2TOP 
	      src_mask(i,j) =  8000 + 100*lp1 + min_dig
           endif
           
           if((sum_top < 0.).and.(sum_sub >= 0.)) then
              if(data4(i) < 15./1.72) SUB2TOP = 0.3
              data2 (i) = data4 (i) 
              data1 (i) = data2 (i) / SUB2TOP  
              src_mask(i,j) = 1000*lp1 + 800 + min_dig
           endif
        endif
        if(src_mask(i,j) == 3344) print *,i,j,sum_top,sum_sub,min_dig,lp1,layer_id
     end do

     write (30) (data1 (i),i=1,nc)
     write (31) (data2 (i),i=1,nc)

  end do

  deallocate (data1,data2,data3,data4)

    ! These are WHSD-v2 files
  close (10,status='keep')
  close (11,status='keep')

    ! These are STATSFO files
  close (20,status='keep')
  close (21,status='keep')

    ! These are WHSD-v2-STATSFO files
  close (30,status='keep')
  close (31,status='keep')

end subroutine merge_carbon_statsgo

! -----------------------------------------------------------------------------------

subroutine merge_min (self_fill,ver_merge,layer_id,nc,nr,src_mask, &
  path1,path2,path3)
  implicit none
  integer,  intent(in) :: nc,nr,layer_id
  logical, intent (in) :: self_fill,ver_merge
  character(len=300),intent (in) :: path1,path2,path3
  character(len=100) :: file1,file2,file3,file4
  character(len=100) :: file5,file6,file7,file8
  integer, intent (inout), dimension (nc,nr) :: src_mask

  integer :: i,j,lp1

  real, allocatable, dimension (:) :: data1, data2, data3, data4
  real, allocatable, dimension (:) :: data5, data6, data7, data8
  real :: sum_top, sum_sub

  allocate (data1(1:nc))
  allocate (data2(1:nc))
  allocate (data3(1:nc))
  allocate (data4(1:nc))
  allocate (data5(1:nc))
  allocate (data6(1:nc))
  allocate (data7(1:nc))
  allocate (data8(1:nc))  

  file1 = 'clay_top_30sec.dat'
  file2 = 'sand_top_30sec.dat'
  file3 = 'clay_sub_30sec.dat'
  file4 = 'sand_sub_30sec.dat'

  print *,'file1 :',trim(path1)//trim(file1)
  print *,'file2 :',trim(path2)//trim(file1)
  print *,'file3 :',trim(path3)//trim(file1)

  lp1 = layer_id + 1

  if(self_fill) print *,'Vertical self-filling of first layer'
  if(ver_merge)  print *,'Vertical MERGING layers'
  open (10,file=trim(path1)//trim(file1),form='unformatted', &
       status='old',action='read')
  open (11,file=trim(path1)//trim(file2),form='unformatted', &
       status='old',action='read')
  open (12,file=trim(path1)//trim(file3),form='unformatted', &
       status='old',action='read')
  open (13,file=trim(path1)//trim(file4),form='unformatted', &
       status='old',action='read')

  open (20,file=trim(path2)//trim(file1),form='unformatted', &
       status='old',action='read')
  open (21,file=trim(path2)//trim(file2),form='unformatted', &
       status='old',action='read')
  open (22,file=trim(path2)//trim(file3),form='unformatted', &
       status='old',action='read')
  open (23,file=trim(path2)//trim(file4),form='unformatted', &
       status='old',action='read')

  open (30,file=trim(path3)//trim(file1),form='unformatted', &
       status='unknown',action='write')
  open (31,file=trim(path3)//trim(file2),form='unformatted', &
       status='unknown',action='write')
  open (32,file=trim(path3)//trim(file3),form='unformatted', &
       status='unknown',action='write')
  open (33,file=trim(path3)//trim(file4),form='unformatted', &
       status='unknown',action='write')

 
  do j=1,nr

     read (10) (data1 (i),i=1,nc)
     read (11) (data2 (i),i=1,nc)
     read (12) (data3 (i),i=1,nc)
     read (13) (data4 (i),i=1,nc)

     read (20) (data5 (i),i=1,nc)
     read (21) (data6 (i),i=1,nc)
     read (22) (data7 (i),i=1,nc)
     read (23) (data8 (i),i=1,nc)

     do i=1,nc

	sum_top = data1(i) + data2(i)
        sum_sub = data3(i) + data4(i)

	if(layer_id == 1) then
	   if(sum_top >= 0.) src_mask(i,j) = 1
	   if(sum_sub >= 0.) src_mask(i,j) = src_mask(i,j) + 10 
	endif

        if(ver_merge) then
           if(self_fill) then
              if((sum_top >=0.).and.(sum_sub < 0.)) then                
                 data3 (i) = data1 (i)
                 data4 (i) = data2 (i)
	         src_mask(i,j) = src_mask(i,j) + 80
              endif
              
              if((sum_top < 0.).and.(sum_sub >= 0.)) then
                 data1 (i) = data3 (i)
                 data2 (i) = data4 (i)
	         src_mask(i,j) = src_mask(i,j) + 8                
              endif
           endif
        endif
        
        sum_top = data5(i) + data6(i)
        sum_sub = data7(i) + data8(i)
          
        if((sum_top >=0.).and.(sum_sub >=0.)) then !SM: remember to switch
           data1 (i) = data5 (i)
           data2 (i) = data6 (i)
           data3 (i) = data7 (i)
           data4 (i) = data8 (i)
	   src_mask(i,j) = 100*(src_mask(i,j)/100) + 10*lp1 + lp1
        endif
        
        if(ver_merge) then 
           if((sum_top >=0.).and.(sum_sub < 0.)) then
              data1 (i) = data5 (i)
              data2 (i) = data6 (i)
              data3 (i) = data5 (i)  
              data4 (i) = data6 (i)  
	      src_mask(i,j) = 100*(src_mask(i,j)/100) + 80 + lp1         
           endif
           
           if((sum_top < 0.).and.(sum_sub >= 0.)) then
              data1 (i) = data7 (i)  
              data2 (i) = data8 (i)  
              data3 (i) = data7 (i)
              data4 (i) = data8 (i)           
              src_mask(i,j) = 100*(src_mask(i,j)/100) + 10*lp1 + 8
           endif
        endif
!	if((data1(i) >= 0.).and.(src_mask(i,j) ==0.)) &
!	   print *,i,j,data1(i),data2(i),data3(i),data4(i),src_mask(i,j)
!	if((data1(i) < 0.).and.(src_mask(i,j) > 0.)) then
!	   print *,i,j,data1(i),data2(i),data3(i),data4(i),src_mask(i,j), &
!	          data5(i),data6(i),data7(i),data8(i)
!	  pause
!	endif
     end do

     write (30) (data1 (i),i=1,nc)
     write (31) (data2 (i),i=1,nc)
     write (32) (data3 (i),i=1,nc)
     write (33) (data4 (i),i=1,nc)
  end do

  deallocate (data1, data2, data3, data4)
  deallocate (data5, data6, data7, data8)

  close (10,status='keep')
  close (11,status='keep')
  close (12,status='keep')
  close (13,status='keep')

  close (20,status='keep')
  close (21,status='keep')
  close (22,status='keep')
  close (23,status='keep')

  close (30,status='keep')
  close (31,status='keep')
  close (32,status='keep')
  close (33,status='keep')

end subroutine merge_min

! -----------------------------------------------------------------------------------

subroutine merge_carbon (self_fill,ver_merge, layer_id,nc,nr,src_mask, &
  path1,path2,path3)
  implicit none
  integer,  intent(in) :: nc,nr
  logical,  intent(in) :: self_fill, ver_merge
  integer,  intent(in) :: layer_id
  real                 :: SUB2TOP 
  character(len=300),intent (in) :: path1,path2,path3
  character(len=100) :: file1,file2,file3,file4

  integer :: i,j, lp1,min_dig
   real, allocatable, dimension (:) :: data1, data2, data3, data4
  integer, intent (inout), dimension (nc,nr) :: src_mask
  real :: sum_top, sum_sub

  allocate (data1(1:nc))
  allocate (data2(1:nc))
  allocate (data3(1:nc))
  allocate (data4(1:nc))

  file1 = 'oc_top_30sec.dat'
  file2 = 'oc_sub_30sec.dat'

  print *,'file1 :',trim(path1)//trim(file1)
  print *,'file2 :',trim(path2)//trim(file1)
  print *,'file3 :',trim(path3)//trim(file1)

  if(self_fill) print *,'Vertical self-filling of first layer'

  open (10,file=trim(path1)//trim(file1),form='unformatted', &
       status='old',action='read')
  open (11,file=trim(path1)//trim(file2),form='unformatted', &
       status='old',action='read')

  open (20,file=trim(path2)//trim(file1),form='unformatted', &
       status='old',action='read')
  open (21,file=trim(path2)//trim(file2),form='unformatted', &
       status='old',action='read')

  open (30,file=trim(path3)//trim(file1),form='unformatted', &
       status='unknown',action='write')
  open (31,file=trim(path3)//trim(file2),form='unformatted', &
       status='unknown',action='write')
	
  lp1 = layer_id + 1 

  do j=1,nr

     read (10) (data1 (i),i=1,nc)
     read (11) (data2 (i),i=1,nc)

     read (20) (data3 (i),i=1,nc)
     read (21) (data4 (i),i=1,nc)
	
     do i=1,nc

	if (layer_id == 1) then
           if(data1 (i) >=0.) data1(i) = data1(i) 
           if(data2 (i) >=0.) data2(i) = data2(i) 
        endif

        sum_top = data1(i) 
        sum_sub = data2(i) 
	SUB2TOP = 1.

	if(layer_id == 1) then
	   if(sum_top >= 0.) src_mask(i,j) = src_mask (i,j) +  100
	   if(sum_sub >= 0.) src_mask(i,j) = src_mask (i,j) + 1000 
	endif
	
        if(ver_merge) then
           if(self_fill) then
              if((sum_top >=0.).and.(sum_sub < 0.)) then
	         if(data1(i) < 15./1.72) SUB2TOP = 0.3
                 data2 (i) = data1 (i) * SUB2TOP
	         src_mask(i,j) = src_mask (i,j) + 8000 
              endif
              
              if((sum_top < 0.).and.(sum_sub >= 0.)) then
                 if(data2(i) < 15./1.72) SUB2TOP = 0.3
                 data1 (i) = data2 (i) / SUB2TOP
                 src_mask(i,j) = src_mask (i,j) + 800 
              endif              
           endif
        endif

        sum_top = data3(i) 
        sum_sub = data4(i) 
        SUB2TOP = 1.
	min_dig = src_mask (i,j) - 100*(src_mask (i,j)/100)

        if((sum_top >=0.).and.(sum_sub >=0.)) then  !SM: switch to and
           data1 (i) = data3 (i) 
           data2 (i) = data4 (i) 
           src_mask(i,j) =  1000*lp1 + 100*lp1 + min_dig
           if(src_mask(i,j) == 3344) print *,i,j,sum_top,sum_sub,min_dig,lp1,layer_id
        endif

        if(ver_merge) then
           if((sum_top >=0.).and.(sum_sub < 0.)) then
	      if(data3(i) < 15./1.72) SUB2TOP = 0.3
              data1 (i) = data3 (i) 
              data2 (i) = data1 (i) * SUB2TOP 
	      src_mask(i,j) =  8000 + 100*lp1 + min_dig
           endif
           
           if((sum_top < 0.).and.(sum_sub >= 0.)) then
              if(data4(i) < 15./1.72) SUB2TOP = 0.3
              data2 (i) = data4 (i) 
              data1 (i) = data2 (i) / SUB2TOP  
              src_mask(i,j) = 1000*lp1 + 800 + min_dig
           endif
        endif
        if(src_mask(i,j) == 3344) print *,i,j,sum_top,sum_sub,min_dig,lp1,layer_id
     end do

     write (30) (data1 (i),i=1,nc)
     write (31) (data2 (i),i=1,nc)

  end do

  deallocate (data1,data2,data3,data4)

  close (10,status='keep')
  close (11,status='keep')

  close (20,status='keep')
  close (21,status='keep')

  close (30,status='keep')
  close (31,status='keep')

end subroutine merge_carbon

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
