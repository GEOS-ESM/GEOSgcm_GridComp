!-*- F90 -*-
module CUB2LATLON_mod
  !--------------------------------------------------------------------!
  ! author:  Michael Herzog                                            !
  ! email:   Michael.Herzog@noaa.gov                                   !
  ! date:    May 2006                                                  !
  ! version: 0.1                                                       !
  !                                                                    !
  ! routines for interpolation from cubed sphere to latlon             ! 
  !--------------------------------------------------------------------!
  use fv_arrays_mod, only: REAL4, REAL8

  implicit none

  private
  public :: read_c2l_namelist,                                           &
            init_cubsph_grid,  init_latlon_grid,                         &
            read_c2l_weight,   write_c2l_weight,                         &
            new_get_c2l_weight, get_c2l_weight,    interpolate_data,                         &
            read_netcdf_grid,  read_grid_dimensions,                     &
            do_c2l_interpolation,                                        &
            do_c2l_interpolation_r4

! interface do_c2l_interpolation
!    module procedure do_c2l_interpolation_r8
!    module procedure do_c2l_interpolation_r4
! end interface

  real(REAL8), parameter :: pi = 3.141592653589793

contains
!======================================================================!
  subroutine read_c2l_namelist(ntiles, nlon, nlat, finer_steps,         &
                               read_res, write_res, memphis,            &
                               grid_file, data_file, data_out,          &
                               uname, vname, missing_value, fill_missing)
    !------------------------------------------------------------------!
    ! read namelist files                                              !
    !------------------------------------------------------------------!
    integer, intent(out) :: ntiles, nlon, nlat, finer_steps
    logical, intent(out) :: read_res, write_res, memphis, fill_missing
    character(len=120), intent(out) :: grid_file, data_file, data_out,  &
                                       uname, vname, missing_value
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    integer :: ios, l, nlon_fine, nlat_fine
    character(len=80) :: filename
    logical :: exists

    namelist /c2l_nml/ ntiles, nlon, nlat, finer_steps,                 &
                       write_res, read_res, memphis,                    &
                       grid_file, data_file, data_out,                  &
                       uname, vname, missing_value, fill_missing
    !------------------------------------------------------------------!
    ! read file names with cubed sphere data from namelist file        !
    !------------------------------------------------------------------!
    filename = "c2l.nml"
    inquire(file=filename,exist=exists)

    if (.not. exists) then
       write(6,100) trim(filename)
100    format (/,"namelist file ",a," doesn't exist",/)
       stop
    else
       ntiles=6
       nlon=0
       nlat=0
       read_res=.false.
       write_res=.false.
       memphis=.false.
       finer_steps=0
       grid_file="grid_spec"
       data_file="?"
       data_out="?"
       uname="ucomp"
       vname="vcomp"
       missing_value="missing_value"
       fill_missing=.false.
       !---------------------------------------------------------------!
       ! read main namelist                                            !
       !---------------------------------------------------------------!
       open (10,file=filename)
       read (10,nml=c2l_nml,iostat=ios)
       close(10)
       if (ios > 0) then
          write(6,101) trim(filename), ios
101       format(/,"c2l_nml ERROR: reading ",a,", iostat=",i4,/)
          stop
       endif

       if (nlon*nlat==0) then
          write(6,111) trim(filename)
111       format(/,"nlon, nlat must be specified in :",a,/)
          stop
       endif
       nlon_fine=2**finer_steps*nlon
       nlat_fine=2**finer_steps*(nlat-1)+1
       write(6,106) nlon, nlat, nlon_fine, nlat_fine
106    format(/,"will output data on a latlon grid with nlon =",i4,", nlat =",i4,/, &
                "starting from nlon =",i4,", nlat =",i4," for interpolation",/)

       if (ntiles/=6) write(6,112) ntiles
112    format(/,"WARNING: ntiles not equal 6! ntiles = ",i3)

       if (trim(data_out)=="?") data_out=trim(data_file)

    endif

  end subroutine read_c2l_namelist
  !====================================================================!
  subroutine read_grid_dimensions(grid_file, npx, npy)
    !------------------------------------------------------------------!
    ! write grid data in netcdf format                                 !
    !------------------------------------------------------------------!
#include "netcdf.inc"
    character(len=120), intent(in) :: grid_file
    integer, intent(out) :: npx, npy
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    integer :: status, ncid, lon_id, ndims, dimids(2), l

    character(len=120) :: filename
    logical :: exists

#if !defined(MAPL_MODE)

    l=1
    write(filename,100) trim(grid_file), l
100 format(a,'.tile',i1,'.nc')
    inquire(file=filename, exist=exists)
    if (.not. exists) then
       print 110, trim(filename)
110    format(/,"grid file ",a," doesn't exist",/)
       stop
    endif
    status = nf_open(trim(filename), 0, ncid)
    if (status /= nf_noerr) then
       print 120, trim(filename)
120    format (/,"nf_open: could not open file: ",a,/)
       stop
    endif

    status = nf_inq_varid(ncid, "grid_lon", lon_id)
    status = nf_inq_varndims(ncid, lon_id, ndims)
    if (ndims/=2) then
       print *, " unexpected number of dimension for grid_lon: ", ndims
       stop
    endif
    status = nf_inq_vardimid(ncid, lon_id, dimids)
    status = nf_inq_dimlen(ncid, dimids(1), npx)
    status = nf_inq_dimlen(ncid, dimids(2), npy)

    write(*,200) npx, npy
200 format (" grid dimension, npx/npy = ", 2i5)

    status = nf_close(ncid)
#endif

  end subroutine read_grid_dimensions
  !====================================================================!
  subroutine init_cubsph_grid(npx, npy, ntiles, grid_file, sph_corner)
    !------------------------------------------------------------------!
    ! read in cubed sphere grid from file,                             !
    ! calculate cell center from cell corners                          !
    !                                                                  !
    ! input:                                                           !
    ! npx, npy, ntiles       number of grid points and tiles           !
    !                                                                  !
    ! output:                                                          !
    ! sph_corner             cell corners in spherical coor            !
    !------------------------------------------------------------------!
    use GHOST_CUBSPH_mod, only: B_grid, ghost_cubsph_update

    integer, intent(in) :: npx, npy, ntiles
    real(REAL8), dimension(2,0:npx+1,0:npy+1,ntiles), intent(out) :: sph_corner
    character(len=120), intent(in) :: grid_file
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    integer :: l
    !------------------------------------------------------------------!
    ! read sph_corner from file                                        !
    ! do transposition to accomodate Bill Putman's array definition    ! 
    !------------------------------------------------------------------!
    call read_netcdf_grid(npx, npy, ntiles, grid_file, sph_corner)
    !------------------------------------------------------------------!
    ! do halo update                                                   !
    !------------------------------------------------------------------!
    do l=1,ntiles
       sph_corner(1:2,0    ,0    ,l)=0.
       sph_corner(1:2,npx+1,0    ,l)=0.
       sph_corner(1:2,0    ,npy+1,l)=0.
       sph_corner(1:2,npx+1,npy+1,l)=0.
       call ghost_cubsph_update(sph_corner(1:1,0:npx+1,0:npy+1,:), 0, npx+1, 0, npy+1, 1, &
                                1, ntiles, 1, 1, l, B_grid)
       call ghost_cubsph_update(sph_corner(2:2,0:npx+1,0:npy+1,:), 0, npx+1, 0, npy+1, 1, &
                                1, ntiles, 1, 1, l, B_grid)
    enddo

  end subroutine init_cubsph_grid
  !====================================================================!
  subroutine read_netcdf_grid(npx, npy, ntiles, grid_file, sph_corner)
    !------------------------------------------------------------------!
    ! write grid data in netcdf format                                 !
    !------------------------------------------------------------------!
#include "netcdf.inc"
    integer, intent(in) :: npx, npy, ntiles
    real(REAL8), dimension(2,0:npx+1,0:npy+1,ntiles), intent(out) :: sph_corner
    character(len=120), intent(in) :: grid_file
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    real(REAL4), dimension(npx,npy) :: var_r4
    real(REAL8), dimension(npx,npy) :: var_r8
    integer*2, dimension(npx,npy) :: var_i2
    integer :: l, status, ncid, lon_id, lat_id, ndims, dimids(2), n1, n2, type

    character(len=120) :: filename
    logical :: exists

#if !defined(MAPL_MODE)
    do l=1,ntiles
       write(filename,100) trim(grid_file),l
100    format(a,'.tile',i1,'.nc')
       inquire(file=filename, exist=exists)
       if (.not. exists) then
          print 110, trim(filename)
110       format(/,"grid file ",a," doesn't exist",/)
          stop
       endif
       status = nf_open(trim(filename), 0, ncid)
       if (status /= nf_noerr) then
          print 120, trim(filename)
120       format (/,"nf_open: could not open file: ",a,/)
          stop
       endif
       
       status = nf_inq_varid(ncid, "grid_lon", lon_id)
       status = nf_inq_varndims(ncid, lon_id, ndims)
       if (ndims/=2) then
          print *, " unexpected number of dimension for grid_lon: ", ndims
          stop
       endif
       status = nf_inq_vardimid(ncid, lon_id, dimids)
       status = nf_inq_dimlen(ncid, dimids(1), n1)
       status = nf_inq_dimlen(ncid, dimids(2), n2)
       if (n1/=npx .or. n2/=npy) then
          print *, " unexpected grid dimension, npx/npy = ", n1, n2
          stop
       endif
       
       status = nf_inq_varid(ncid, "grid_lat", lat_id)
       status = nf_inq_varndims(ncid, lat_id, ndims)
       if (ndims/=2) then
          print *, " unexpected number of dimension for grid_lat: ", ndims
          stop
       endif
       status = nf_inq_vardimid(ncid, lat_id, dimids)
       status = nf_inq_dimlen(ncid, dimids(1), n1)
       status = nf_inq_dimlen(ncid, dimids(2), n2)
       if (n1/=npx .or. n2/=npy) then
          print *, " unexpected grid dimension, npx/npy = ", n1, n2
          stop
       endif
       
       status = nf_inq_vartype(ncid, lon_id, type)
       if (type==nf_double) then
          status = nf_get_var_double(ncid, lon_id, var_r8)
          sph_corner(1,1:npx,1:npy,l)=var_r8(1:npx,1:npy)
       elseif (type==nf_real(REAL8)) then
          status = nf_get_var_real(ncid, lon_id, var_r4)
          sph_corner(1,1:npx,1:npy,l)=var_r4(1:npx,1:npy)
       else
          print *, " unrecognized var_type: ", type
       endif
       sph_corner(1,1:npx,1:npy,l)=sph_corner(1,1:npx,1:npy,l)*pi/180.
       
       status = nf_inq_vartype(ncid, lat_id, type)
       if (type==nf_double) then
          status = nf_get_var_double(ncid, lat_id, var_r8)
          sph_corner(2,1:npx,1:npy,l)=var_r8(1:npx,1:npy)
       elseif (type==nf_real(REAL8)) then
          status = nf_get_var_real(ncid, lat_id, var_r4)
          sph_corner(2,1:npx,1:npy,l)=var_r4(1:npx,1:npy)
       else
          print *, " unrecognized var_type: ", type
       endif
       sph_corner(2,1:npx,1:npy,l)=sph_corner(2,1:npx,1:npy,l)*pi/180.
       
       status = nf_close(ncid)
    enddo
    write(*,200) trim(grid_file)
200 format(/," grid information read from file ",a,".tile?.nc",/)
#endif
    
  end subroutine read_netcdf_grid
  !====================================================================!
  subroutine init_latlon_grid(xlon, ylat, nlon, nlat)
    !------------------------------------------------------------------!
    ! initialize latlon grid                                           !
    !------------------------------------------------------------------!
    integer, intent(in)  :: nlon, nlat
    real(REAL8),    intent(out) :: xlon(nlon), ylat(nlat)
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    real(REAL8)    :: dlon, dlat
    integer :: i, j
    !------------------------------------------------------------------!
    ! calculate location of cell centers                               !
    !------------------------------------------------------------------!
    dlon=(pi+pi)/real(nlon)
    dlat=pi/real(nlat-1)
    
    do i=1,nlon
       xlon(i)=(real(i)-1.)*dlon
    enddo
    ylat(1)   =-0.5*pi
    ylat(nlat)= 0.5*pi
    do j=2,nlat-1
       ylat(j)=ylat(1)+(real(j)-1.)*dlat
    enddo
    
  end subroutine init_latlon_grid
  !====================================================================!
  subroutine read_c2l_weight(c2l_index, c2l_weight, nlon, nlat, npx, npy, ntiles, &
                elon_cubsph, elat_cubsph, elon_latlon, elat_latlon, found)
    !------------------------------------------------------------------!
    ! read restart file c2l_res.nc                                     !
    !------------------------------------------------------------------!
#include "netcdf.inc"
    integer, intent(in) :: nlon, nlat, npx, npy, ntiles
    real(REAL8),    dimension(4, nlon, nlat), intent(out) :: c2l_weight
    integer, dimension(3, nlon, nlat), intent(out) :: c2l_index

    real(REAL8), dimension(3, 0:npx, 0:npy ,ntiles), intent(out) :: elon_cubsph, elat_cubsph
    real(REAL8), dimension(3, nlon, nlat)          , intent(out) :: elon_latlon, elat_latlon

    logical, intent(out) :: found
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    real(REAL8), dimension(:,:,:),   allocatable :: var3d
    real(REAL8), dimension(:,:,:,:), allocatable :: var4d

    integer :: ncid_res, status, c2l_weight_id, c2l_index_id, var_id,   &
               dimids(3), nx, ny, nt,start(4), count(4)

    status = nf_open("c2l_res.nc", 0, ncid_res)
    if (status == nf_noerr) status = nf_inq_varid(ncid_res, "c2l_weight", c2l_weight_id)
    if (status == nf_noerr) then
       status = nf_inq_vardimid(ncid_res, c2l_weight_id, dimids)
       status = nf_inq_dimlen(ncid_res, dimids(2), nx)
       status = nf_inq_dimlen(ncid_res, dimids(3), ny)
       if (nx/=nlon .or. ny/=nlat) status = nf_fatal

       status = nf_inq_dimid(ncid_res, "npx", dimids(1))
       status = nf_inq_dimid(ncid_res, "npy", dimids(2))
       status = nf_inq_dimid(ncid_res, "ntiles", dimids(3))
       status = nf_inq_dimlen(ncid_res, dimids(1), nx)
       status = nf_inq_dimlen(ncid_res, dimids(2), ny)
       status = nf_inq_dimlen(ncid_res, dimids(3), nt)
       if (nx/=npx .or. ny/=npy .or. nt/=ntiles) status = nf_fatal
    endif

    start(1:3)=(/1,1,1/)
    count(1:3)=(/4,nlon,nlat/)
    allocate(var3d(4, nlon, nlat))
    if (status == nf_noerr) status = nf_get_vara_double(ncid_res, c2l_weight_id, start, count, var3d)
    if (status == nf_noerr) c2l_weight(:,:,:) = var3d(:,:,:)
    deallocate(var3d)

    start(1:3)=(/1,1,1/)
    count(1:3)=(/3,nlon,nlat/)
    if (status == nf_noerr) status = nf_inq_varid(ncid_res, "c2l_index",  c2l_index_id)
    if (status == nf_noerr) status = nf_get_vara_int(ncid_res, c2l_index_id, start, count, c2l_index)

    start(1:3)=(/1,1,1/)
    count(1:3)=(/3,nlon,nlat/)
    allocate(var3d(3, nlon, nlat))
    if (status == nf_noerr) status = nf_inq_varid(ncid_res, "elon_latlon", var_id)
    if (status == nf_noerr) status = nf_get_vara_double(ncid_res, var_id, start, count, var3d)
    if (status == nf_noerr) elon_latlon(:,:,:) = var3d(:,:,:)
    if (status == nf_noerr) status = nf_inq_varid(ncid_res, "elat_latlon", var_id)
    if (status == nf_noerr) status = nf_get_vara_double(ncid_res, var_id, start, count, var3d)
    if (status == nf_noerr) elat_latlon(:,:,:) = var3d(:,:,:)
    deallocate(var3d)

    start(1:4)=(/1,1,1,1/)
    count(1:4)=(/3,npx+1,npy+1,ntiles/)
    allocate(var4d(3, npx+1, npy+1, ntiles))
    if (status == nf_noerr) status = nf_inq_varid(ncid_res, "elon_cubsph", var_id)
    if (status == nf_noerr) status = nf_get_vara_double(ncid_res, var_id, start, count, var4d)
    if (status == nf_noerr) elon_cubsph(:,0:npx,0:npy,:) = var4d(:,1:npx+1,1:npy+1,:)
    if (status == nf_noerr) status = nf_inq_varid(ncid_res, "elat_cubsph", var_id)
    if (status == nf_noerr) status = nf_get_vara_double(ncid_res, var_id, start, count, var4d)
    if (status == nf_noerr) elat_cubsph(:,0:npx,0:npy,:) = var4d(:,1:npx+1,1:npy+1,:)
    deallocate(var4d)

!!$    allocate(var4d(3, 0:npx, 0:npy, ntiles))
!!$    if (status == nf_noerr) status = nf_inq_varid(ncid_res, "elon_cubsph", var_id)
!!$    if (status == nf_noerr) status = nf_get_var_double(ncid_res, var_id, var4d)
!!$    if (status == nf_noerr) elon_cubsph(:,:,:,:) = var4d(:,:,:,:)
!!$    if (status == nf_noerr) status = nf_inq_varid(ncid_res, "elat_cubsph", var_id)
!!$    if (status == nf_noerr) status = nf_get_var_double(ncid_res, var_id, var4d)
!!$    if (status == nf_noerr) elat_cubsph(:,:,:,:) = var4d(:,:,:,:)
!!$    deallocate(var4d)

    if (status == nf_noerr) then
       status = nf_close(ncid_res)
       found = .true.
#if !defined(MAPL_MODE)
       print*, 'c2l restart file succesfully read'
#endif
    else
       found = .false.
#if !defined(MAPL_MODE)
       print*, 'ERROR reading c2l restart - will recalculate coefficients'
#endif
    endif

#if !defined(MAPL_MODE)
    write(*,100) 
100 format (/," done reading c2l_index and c2l_weight",/)
#endif

  end subroutine read_c2l_weight
  !====================================================================!
  subroutine write_c2l_weight(c2l_index, c2l_weight, nlon, nlat, npx, npy, ntiles, &
                              elon_cubsph, elat_cubsph, elon_latlon, elat_latlon)
    !------------------------------------------------------------------!
    ! write restart file c2l_res.nc                                    !
    !------------------------------------------------------------------!
#include "netcdf.inc"
    integer, intent(in) :: nlon, nlat, npx, npy, ntiles
    real(REAL8),    dimension(4, nlon, nlat), intent(in) :: c2l_weight
    integer, dimension(3, nlon, nlat), intent(in) :: c2l_index
    real(REAL8), dimension(3, 0:npx, 0:npy ,ntiles), intent(in) :: elon_cubsph, elat_cubsph
    real(REAL8), dimension(3, nlon, nlat)          , intent(in) :: elon_latlon, elat_latlon
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    real(REAL8), dimension(:,:,:),   allocatable :: var3d
    real(REAL8), dimension(:,:,:,:), allocatable :: var4d
    integer :: ncid_res, status, nlon_dim, nlat_dim, n3_dim,            &
               n4_dim, dims(4), c2l_index_id, c2l_weight_id,            &
               elon1_id, elat1_id, elon2_id, elat2_id,                  &
               npx_dim, npy_dim, npx1_dim, npy1_dim, ntiles_dim,        &
               start(4), count(4)

    status = nf_create("c2l_res.nc", nf_clobber, ncid_res)

    status = nf_def_dim(ncid_res, "n3" , 3, n3_dim)
    status = nf_def_dim(ncid_res, "n4", 4, n4_dim)
    status = nf_def_dim(ncid_res, "nlon", nlon, nlon_dim)
    status = nf_def_dim(ncid_res, "nlat", nlat, nlat_dim)
    status = nf_def_dim(ncid_res, "npx", npx, npx_dim)
    status = nf_def_dim(ncid_res, "npy", npy, npy_dim)
    status = nf_def_dim(ncid_res, "npx1", npx+1, npx1_dim)
    status = nf_def_dim(ncid_res, "npy1", npy+1, npy1_dim)
    status = nf_def_dim(ncid_res, "ntiles", ntiles, ntiles_dim)

    dims(1) = n3_dim
    dims(2) = nlon_dim
    dims(3) = nlat_dim
    status = nf_def_var(ncid_res, "c2l_index", nf_int, 3, dims, c2l_index_id)

    dims(1) = n4_dim
    dims(2) = nlon_dim
    dims(3) = nlat_dim
    status = nf_def_var(ncid_res, "c2l_weight", nf_double, 3, dims, c2l_weight_id)

    dims(1) = n3_dim
    dims(2) = nlon_dim
    dims(3) = nlat_dim
    status = nf_def_var(ncid_res, "elon_latlon", nf_double, 3, dims, elon1_id)
    status = nf_def_var(ncid_res, "elat_latlon", nf_double, 3, dims, elat1_id)

    dims(1) = n3_dim
    dims(2) = npx1_dim
    dims(3) = npy1_dim
    dims(4) = ntiles
    status = nf_def_var(ncid_res, "elon_cubsph", nf_double, 4, dims, elon2_id)
    status = nf_def_var(ncid_res, "elat_cubsph", nf_double, 4, dims, elat2_id)

    status = nf_put_att_text(ncid_res, c2l_index_id, "long_name", 60,   &
         "index location nearest cubsph grid point to the lower-left")
    status = nf_put_att_text(ncid_res, c2l_weight_id, "long_name", 50,  &
         "weights for linear cubsph_to_latlon interpolation")
    status = nf_put_att_text(ncid_res, elon1_id, "long_name", 30,  &
         "lon unit vector for latlon grid")
    status = nf_put_att_text(ncid_res, elat1_id, "long_name", 30,  &
         "lat unit vector for latlon grid")
    status = nf_put_att_text(ncid_res, elon2_id, "long_name", 30,  &
         "lon unit vector for cubsph grid")
    status = nf_put_att_text(ncid_res, elat2_id, "long_name", 30,  &
         "lat unit vector for cubsph grid")
    status = nf_put_att_text(ncid_res, nf_global, "file_type", 50,      &
         "restart file for  cubsph_to_latlon interpolation")

    status = nf_enddef(ncid_res)

    start(1:3)=(/1,1,1/)
    count(1:3)=(/3,nlon,nlat/)
    status = nf_put_vara_int(ncid_res, c2l_index_id, start, count, c2l_index)

    start(1:3)=(/1,1,1/)
    count(1:3)=(/4,nlon,nlat/)
    allocate(var3d(4, nlon, nlat))
    var3d(:,:,:)=c2l_weight(:,:,:)
    status = nf_put_vara_double(ncid_res, c2l_weight_id, start, count, var3d)
    deallocate(var3d)

    start(1:3)=(/1,1,1/)
    count(1:3)=(/3,nlon,nlat/)
    allocate(var3d(3, nlon, nlat))
    var3d(:,:,:)=elon_latlon(:,:,:)
    status = nf_put_vara_double(ncid_res, elon1_id, start, count, var3d)
    var3d(:,:,:)=elat_latlon(:,:,:)
    status = nf_put_vara_double(ncid_res, elat1_id, start, count, var3d)
    deallocate(var3d)


    start(1:4)=(/1,1,1,1/)
    count(1:4)=(/3,npx+1,npy+1,ntiles/)
    allocate(var4d(3, npx+1, npy+1, ntiles))
    var4d(:,1:npx+1,1:npy+1,:)=elon_cubsph(:,0:npx,0:npy,:)
    status = nf_put_vara_double(ncid_res, elon2_id, start, count, var4d)
    var4d(:,1:npx+1,1:npy+1,:)=elat_cubsph(:,0:npx,0:npy,:)
    status = nf_put_vara_double(ncid_res, elat2_id, start, count, var4d)
    deallocate(var4d)

    status = nf_close(ncid_res)

#ifndef MAPL_MODE
    write(*,100) 
100 format (/," done writing c2l_index and c2l_weight",/)
#endif

  end subroutine write_c2l_weight

  !====================================================================!
  subroutine new_get_c2l_weight(sph_corner, npx, npy, ntiles, xlon, ylat, nlon, nlat,  &
       c2l_index, c2l_weight, elon_cubsph, elat_cubsph, elon_latlon, elat_latlon)
    !------------------------------------------------------------------!
    ! calculate weights for bilinear interpolation                     !
    ! from cubed sphere to latlon grid                                 !
    !                                                                  !
    ! input:                                                           !
    ! sph_corner      cubed sphere corner location in spherical coor   !
    ! npx, npy        number of corners per tile                       !
    ! ntiles          number of tiles                                  !
    ! xlon, ylat      latlon grid coor                                 !
    ! nlon, nlat      latlon grid dimension                            !
    !                                                                  !
    ! output:                                                          !
    ! c2l_index       cubed sphere index for latlon interpolation      !
    ! c2l_weight      weights for cubsph_to_latlon interpolation       !
    ! elon_cubsph     lon unit vector for cubed sphere center          !
    ! elat_cubsph     lat unit vector for cubed sphere center          !
    ! elon_latlon     lon unit vector for latlon grid                  !
    ! elat_latlon     lat unit vector for latlon grid                  !
    !------------------------------------------------------------------!
    use GRID_UTILS_mod, only: great_circle, dist2side, spherical_angle, &
                              xyz2latlon, latlon2xyz, unit_vect_latlon, &
                              vect_cross, normalize_vect

    integer, intent(in) :: npx, npy, ntiles, nlon, nlat
    real(REAL8), dimension(2,0:npx+1, 0:npy+1, ntiles), intent(in) :: sph_corner
    real(REAL8), intent(in) :: xlon(nlon), ylat(nlat)


    real(REAL8), optional, dimension(3, 0:npx, 0:npy ,ntiles), intent(out) :: elon_cubsph, elat_cubsph
    real(REAL8), optional, dimension(3, nlon, nlat)          , intent(out) :: elon_latlon, elat_latlon

    real(REAL8),    dimension(4, nlon, nlat), intent(out) :: c2l_weight
    integer, dimension(3, nlon, nlat), intent(out) :: c2l_index
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    real(REAL8) :: xyz_corner(3,4)
    real(REAL8) :: sph_latlon(2), sph_center(2)
    real(REAL8) :: xyz_center(3, 0:npx  , 0:npy  , ntiles)
    real(REAL8) :: xyz_latlon(3, nlon, nlat)
    real(REAL8) :: abs_center, dcub, dlon, dlat, coslat,                       &
            distance, shortest(ntiles),                                 &
            angle_1, angle_1a, angle_1b,                                &
            angle_2, angle_2a, angle_2b,                                &
            angle_3, angle_3a, angle_3b,                                &
            angle_4, angle_4a, angle_4b,                                &
            dist1, dist2, dist3, dist4, sum

    integer :: i, j, l, n, ic, jc, lc, icc, jcc, index(3,ntiles),          &
               i_min, i_max, j_min, j_max, iter

    logical :: found(nlon,nlat)

    real(REAL8) :: faceEdgeNormals(3,4,ntiles) 
    real(REAL8) :: p1(3), p2(3)

    real(REAL8) :: xNormals(3,npx,ntiles)
    real(REAL8) :: yNormals(3,npy,ntiles)

    integer :: iTile
    integer :: cubedSphereIndex(3)
    real(REAL8) :: q(3)

    !------------------------------------------------------------------!
    ! cubed sphere: cartesian coordinates of cell corners,             !
    !               cell lenghts between corners,                      !
    !               cartesian and spherical coordinates of cell centers! 
    !               calculate latlon unit vector                       !
    !------------------------------------------------------------------!
    do l=1,ntiles
       do j=0,npx
          do i=0,npy
             call latlon2xyz(sph_corner(:,i  ,j  ,l), xyz_corner(:,1))
             call latlon2xyz(sph_corner(:,i+1,j  ,l), xyz_corner(:,2))
             call latlon2xyz(sph_corner(:,i  ,j+1,l), xyz_corner(:,3))
             call latlon2xyz(sph_corner(:,i+1,j+1,l), xyz_corner(:,4))

             xyz_center(:,i,j,l)=0.25*(xyz_corner(:,1)+xyz_corner(:,2)  &
                                      +xyz_corner(:,3)+xyz_corner(:,4))
             abs_center=xyz_center(1,i,j,l)*xyz_center(1,i,j,l)         &
                       +xyz_center(2,i,j,l)*xyz_center(2,i,j,l)         &
                       +xyz_center(3,i,j,l)*xyz_center(3,i,j,l)
             if (abs_center>0.) xyz_center(:,i,j,l)=xyz_center(:,i,j,l)/sqrt(abs_center)
             if (present(elon_cubsph) .and. present(elat_cubsph)) then
                call xyz2latlon(xyz_center(:,i,j,l), sph_center(:))
                call unit_vect_latlon(sph_center(:), elon_cubsph(:,i,j,l), elat_cubsph(:,i,j,l))
             endif
          enddo
       enddo
    enddo

    !------------------------------------------------------------------!
    ! latlon: cartesian coordinates of cell centers                    !
    !         calculate latlon unit vector                             !
    !------------------------------------------------------------------!
    do j=1,nlat
       do i=1,nlon
          sph_latlon(1)=xlon(i)
          sph_latlon(2)=ylat(j)
          call latlon2xyz(sph_latlon, xyz_latlon(:,i,j))
          if (present(elon_cubsph) .and. present(elat_cubsph)) then
             call unit_vect_latlon(sph_latlon, elon_latlon(:,i,j), elat_latlon(:,i,j))
          endif
       enddo
    enddo

  do l=1,ntiles
   !  print*, 'Tile: ', l
     ! Bottom Edge 
      n = 1
      call latlon2xyz(sph_corner(:,1  ,1,l), p1)
      call latlon2xyz(sph_corner(:,npx,1,l), p2)
      call vect_cross(faceEdgeNormals(:,n,l), p1, p2)
   !  print*, 'Edge: ', n
   !  print*, p1
   !  print*, p2
   !  print*, faceEdgeNormals(:,n,l)
     ! Right Edge 
      n = 2
      call latlon2xyz(sph_corner(:,npx,1  ,l), p1)
      call latlon2xyz(sph_corner(:,npx,npy,l), p2)
      call vect_cross(faceEdgeNormals(:,n,l), p1, p2)
   !  print*, 'Edge: ', n
   !  print*, p1
   !  print*, p2
   !  print*, faceEdgeNormals(:,n,l)
     ! Top Edge 
      n = 3
      call latlon2xyz(sph_corner(:,npx,npy,l), p1)
      call latlon2xyz(sph_corner(:,1  ,npy,l), p2)
      call vect_cross(faceEdgeNormals(:,n,l), p1, p2)
   !  print*, 'Edge: ', n
   !  print*, p1
   !  print*, p2
   !  print*, faceEdgeNormals(:,n,l)
     ! Left Edge 
      n = 4
      call latlon2xyz(sph_corner(:,1,npy,l), p1)
      call latlon2xyz(sph_corner(:,1,1  ,l), p2)
      call vect_cross(faceEdgeNormals(:,n,l), p1, p2)
   !  print*, 'Edge: ', n
   !  print*, p1
   !  print*, p2
   !  print*, faceEdgeNormals(:,n,l)

! Normals are based upon centers _and_ we must take care to get the sign
! right for the orientation of the tile (note swap in the xNormal case)
     do i=1,npx
          call vect_cross(xNormals(:,i,l), xyz_center(:,i,npy,l), xyz_center(:,i,1,l))
     enddo
     do j=1,npy
         call vect_cross(yNormals(:,j,l), xyz_center(:,1,j,l), xyz_center(:,npx,j,l))
     enddo

  enddo

  if (.true.) then

  do j = 1, nlat
     do i = 1, nlon

        q = xyz_latlon(:,i,j)
        cubedSphereIndex = getCubedSphereIndex(npx, npy, q, xNormals, yNormals, faceEdgeNormals)
        c2l_weight(:,i,j) = computeWeights(cubedSphereIndex, q, xyz_center, npx, npy)
        c2l_index(:,i,j) = cubedSphereIndex

     end do
  end do

  else

    do l=1,ntiles
      call latlon2xyz(sph_corner(:,2,25,l), q)
      iTile = getCubedSphereTile(q, faceEdgeNormals)
      print*, l, iTile
    enddo

  endif


contains

  function getCubedSphereIndex(npx, npy, q, xNormals, yNormals, faceEdgeNormals) result(index)
    integer, intent(in) :: npx, npy
    real(REAL8), intent(in) :: q(3)
    real(REAL8), intent(in) :: xNormals(3, npx, 6)
    real(REAL8), intent(in) :: yNormals(3, npy, 6)
    real(REAL8), intent(in) :: faceEdgeNormals(3, 4, 6)
    integer :: index(3)
    integer :: iTile

    iTile = getCubedSphereTile(q, faceEdgeNormals)
    index(1) = binCubedSphereInterval(q, xNormals(:,:,iTile), npx)
    index(2) = binCubedSphereInterval(q, yNormals(:,:,iTile), npy)
    index(3) = iTile

  end function getCubedSphereIndex

  !
  ! Use orientation with primary corners of cube to bin a given xyz point "q"
  ! into 1 of 6 faces of the cube.
  !
  integer function getCubedSphereTile(q, faceEdgeNormals) result(iTile)
    real(REAL8), intent(in) :: q(3)
    real(REAL8), intent(in) :: faceEdgeNormals(3, 4, 6)

    integer :: edge
    logical :: flag

    do iTile = 1,6
       flag = .true.
       do edge = 1,4
     !    print*, dot_product(q, faceEdgeNormals(:,edge,iTile)), edge, iTile, q, faceEdgeNormals(:,edge,iTile)
          if (dot_product(q, faceEdgeNormals(:,edge,iTile)) < 0) then
             flag = .false.
             exit
          endif
       enddo
     ! if (iTile>1) stop
       if (flag) return
    enddo 
    iTile = -1

  end function getCubedSphereTile

  !
  ! This routine, binCubedSphereInterval(), exploits the fact that
  ! each tile on the CS is subdivided by great circles.  Thus, a given
  ! point can be binned in each of the 2 horizontal directions by
  ! determing the angle the position makes with the normals for the
  ! circles along each axis. (If the circles are evenly spaced in
  ! angle, then the bin can be computed directly with modulo
  ! arithmetic on the angle.)  Actually, we only need cos(angle),
  ! which avoids expensive trig functions.
  !
  ! Note that the normals used here are for the great circles that
  ! pass through the _centers_ of the cells, _not_ the corners.
  !
  integer function binCubedSphereInterval(q, normals, n) result(index)
    integer, intent(in) :: n
    real(REAL8), intent(in) :: q(3)
    real(REAL8), intent(in) :: normals(3, n)

    integer :: i
    real(REAL8) :: cosAngle1, cosAngle2

    cosAngle1 = dot_product(q, normals(:,1))
    do i = 2, n
       cosAngle2 = dot_product(q, normals(:,i))
       if (cosAngle1*cosAngle2 < 0) then
          index = i - 1
          return
       end if
    end do
    ! If not binned yet, q is in the outer half of the last cell
    index = n
    return

  end function binCubedSphereInterval

  function computeWeights(index, q, xyz_center, npx, npy) result(weights)
    integer, intent(in) :: index(3)
    real(REAL8), intent(in) :: q(3)
    integer, intent(in) :: npx, npy
    real(REAL8), intent(in) :: xyz_center(3, npx, npy, 6)
    real(REAL8) :: weights(4)

    integer:: ic, jc, iTile
    real(REAL8) :: dist1, dist2, dist3, dist4
    real(REAL8) :: sum

    !------------------------------------------------------------!
    ! calculate shortest distance to each side of rectangle      !
    ! formed by cubed sphere cell centers                        !
    ! special corner treatment                                   !
    !------------------------------------------------------------!
    ic=index(1)
    jc=index(2)
    iTile =index(3)

	print* , ic, jc, iTile

    if (ic==npx-1 .and. jc==npy-1) then

       !------------------------------------------------------------!
       ! calculate weights for bilinear interpolation near corner   !
       !------------------------------------------------------------!
       dist1=dist2side(xyz_center(:,ic+1,jc,iTile),xyz_center(:,ic,jc+1,iTile),q(:))
       dist2=dist2side(xyz_center(:,ic+1,jc,iTile),xyz_center(:,ic,jc  ,iTile),q(:))
       dist3=dist2side(xyz_center(:,ic  ,jc,iTile),xyz_center(:,ic,jc+1,iTile),q(:))

       weights(1)=dist1      ! ic,   jc    weight
       weights(2)=dist2      ! ic,   jc+1  weight
       weights(3)=0.         ! ic+1, jc+1  weight
       weights(4)=dist3      ! ic+1, jc    weight

       sum=weights(1)+weights(2)+weights(4)
       weights(1)=weights(1)/sum
       weights(2)=weights(2)/sum
       weights(4)=weights(4)/sum

    elseif (ic==0 .and. jc==npy-1) then
       !------------------------------------------------------------!
       ! calculate weights for bilinear interpolation near corner   !
       !------------------------------------------------------------!
       dist1=dist2side(xyz_center(:,ic+1,jc+1,iTile),xyz_center(:,ic+1,jc,iTile),q(:))
       dist2=dist2side(xyz_center(:,ic+1,jc  ,iTile),xyz_center(:,ic  ,jc,iTile),q(:))
       dist3=dist2side(xyz_center(:,ic+1,jc+1,iTile),xyz_center(:,ic  ,jc,iTile),q(:))

       weights(1)=dist1      ! ic,   jc    weight
       weights(2)=0.         ! ic,   jc+1  weight
       weights(3)=dist2      ! ic+1, jc+1  weight
       weights(4)=dist3      ! ic+1, jc    weight

       sum=weights(1)+weights(3)+weights(4)
       weights(1)=weights(1)/sum
       weights(3)=weights(3)/sum
       weights(4)=weights(4)/sum

    elseif (jc==0 .and. ic==npx-1) then
       !------------------------------------------------------------!
       ! calculate weights for bilinear interpolation near corner   !
       !------------------------------------------------------------!
       dist1=dist2side(xyz_center(:,ic,jc+1,iTile),xyz_center(:,ic+1,jc+1,iTile),q(:))
       dist2=dist2side(xyz_center(:,ic,jc  ,iTile),xyz_center(:,ic+1,jc+1,iTile),q(:))
       dist3=dist2side(xyz_center(:,ic,jc  ,iTile),xyz_center(:,ic  ,jc+1,iTile),q(:))

       weights(1)=dist1      ! ic,   jc    weight
       weights(2)=dist2      ! ic,   jc+1  weight
       weights(3)=dist3      ! ic+1, jc+1  weight
       weights(4)=0.         ! ic+1, jc    weight

       sum=weights(1)+weights(2)+weights(3)
       weights(1)=weights(1)/sum
       weights(2)=weights(2)/sum
       weights(3)=weights(3)/sum

    else
       !------------------------------------------------------------!
       ! calculate weights for bilinear interpolation if no corner  !
       !------------------------------------------------------------!

       dist1=dist2side(xyz_center(:,ic  ,jc  ,iTile),xyz_center(:,ic  ,jc+1,iTile),q(:))
       dist2=dist2side(xyz_center(:,ic  ,jc+1,iTile),xyz_center(:,ic+1,jc+1,iTile),q(:))
       dist3=dist2side(xyz_center(:,ic+1,jc+1,iTile),xyz_center(:,ic+1,jc  ,iTile),q(:))
       dist4=dist2side(xyz_center(:,ic+1,jc  ,iTile),xyz_center(:,ic  ,jc  ,iTile),q(:))

       weights(1)=dist2*dist3      ! ic,   jc    weight
       weights(2)=dist3*dist4      ! ic,   jc+1  weight
       weights(3)=dist4*dist1      ! ic+1, jc+1  weight
       weights(4)=dist1*dist2      ! ic+1, jc    weight

       sum=weights(1)+weights(2)+weights(3)+weights(4)
       weights(:)=weights(:)/sum

    endif

  end function computeWeights

  end subroutine new_get_c2l_weight

  !====================================================================!
  subroutine get_c2l_weight(sph_corner, npx, npy, ntiles, xlon, ylat, nlon, nlat,  &
       c2l_index, c2l_weight, subset, elon_cubsph, elat_cubsph, elon_latlon, elat_latlon)
    !------------------------------------------------------------------!
    ! calculate weights for bilinear interpolation                     !
    ! from cubed sphere to latlon grid                                 !
    !                                                                  !
    ! input:                                                           !
    ! sph_corner      cubed sphere corner location in spherical coor   !
    ! npx, npy        number of corners per tile                       !
    ! ntiles          number of tiles                                  !
    ! xlon, ylat      latlon grid coor                                 !
    ! nlon, nlat      latlon grid dimension                            !
    !                                                                  !
    ! output:                                                          !
    ! c2l_index       cubed sphere index for latlon interpolation      !
    ! c2l_weight      weights for cubsph_to_latlon interpolation       !
    ! elon_cubsph     lon unit vector for cubed sphere center          !
    ! elat_cubsph     lat unit vector for cubed sphere center          !
    ! elon_latlon     lon unit vector for latlon grid                  !
    ! elat_latlon     lat unit vector for latlon grid                  !
    !------------------------------------------------------------------!
    use GRID_UTILS_mod, only: great_circle, dist2side, spherical_angle, &
                              xyz2latlon, latlon2xyz, unit_vect_latlon

    integer, intent(in) :: npx, npy, ntiles, nlon, nlat
    real(REAL8), dimension(2,0:npx+1, 0:npy+1, ntiles), intent(in) :: sph_corner
    real(REAL8), intent(in) :: xlon(nlon), ylat(nlat)
    logical, intent(in) :: subset

    real(REAL8), optional, dimension(3, 0:npx, 0:npy ,ntiles), intent(out) :: elon_cubsph, elat_cubsph
    real(REAL8), optional, dimension(3, nlon, nlat)          , intent(out) :: elon_latlon, elat_latlon

    real(REAL8),    dimension(4, nlon, nlat), intent(out) :: c2l_weight
    integer, dimension(3, nlon, nlat), intent(out) :: c2l_index
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    real(REAL8) :: xyz_corner(3,4)
    real(REAL8) :: sph_latlon(2), sph_center(2)
    real(REAL8) :: xyz_center(3, 0:npx  , 0:npy  , ntiles)
    real(REAL8) :: xyz_latlon(3, nlon, nlat)
    real(REAL8) :: abs_center, dcub, dlon, dlat, coslat,                       &
            distance, shortest(ntiles),                                 &
            angle_1, angle_1a, angle_1b,                                & 
            angle_2, angle_2a, angle_2b,                                & 
            angle_3, angle_3a, angle_3b,                                &
            angle_4, angle_4a, angle_4b,                                & 
            dist1, dist2, dist3, dist4, sum

    integer :: i, j, l, ic, jc, lc, icc, jcc, index(3,ntiles),          &
               i_min, i_max, j_min, j_max, iter

    logical :: found(nlon,nlat)

    !------------------------------------------------------------------!
    ! cubed sphere: cartesian coordinates of cell corners,             !
    !               cell lenghts between corners,                      !
    !               cartesian and spherical coordinates of cell centers! 
    !               calculate latlon unit vector                       !
    !------------------------------------------------------------------!
    do l=1,ntiles
       do j=0,npy
          do i=0,npx
             call latlon2xyz(sph_corner(:,i  ,j  ,l), xyz_corner(:,1))
             call latlon2xyz(sph_corner(:,i+1,j  ,l), xyz_corner(:,2))
             call latlon2xyz(sph_corner(:,i  ,j+1,l), xyz_corner(:,3))
             call latlon2xyz(sph_corner(:,i+1,j+1,l), xyz_corner(:,4))

             xyz_center(:,i,j,l)=0.25*(xyz_corner(:,1)+xyz_corner(:,2)  &
                                      +xyz_corner(:,3)+xyz_corner(:,4))
             abs_center=xyz_center(1,i,j,l)*xyz_center(1,i,j,l)         &
                       +xyz_center(2,i,j,l)*xyz_center(2,i,j,l)         &
                       +xyz_center(3,i,j,l)*xyz_center(3,i,j,l)
             if (abs_center>0.) xyz_center(:,i,j,l)=xyz_center(:,i,j,l)/sqrt(abs_center)
             if (present(elon_cubsph) .and. present(elat_cubsph)) then
                call xyz2latlon(xyz_center(:,i,j,l), sph_center(:))
                call unit_vect_latlon(sph_center(:), elon_cubsph(:,i,j,l), elat_cubsph(:,i,j,l))
             endif
          enddo
       enddo
    enddo

    !------------------------------------------------------------------!
    ! latlon: cartesian coordinates of cell centers                    !
    !         calculate latlon unit vector                             !
    !------------------------------------------------------------------!
    do j=1,nlat
       do i=1,nlon
          sph_latlon(1)=xlon(i)
          sph_latlon(2)=ylat(j)
          call latlon2xyz(sph_latlon, xyz_latlon(:,i,j))
          if (present(elon_cubsph) .and. present(elat_cubsph)) then 
             call unit_vect_latlon(sph_latlon, elon_latlon(:,i,j), elat_latlon(:,i,j))
          endif
       enddo
    enddo

    !------------------------------------------------------------------!
    ! find lower left corner on cubed sphere for given latlon location !
    !------------------------------------------------------------------!
    found(:,:)=.false.
    if (subset) then
       dlon=xlon(2)-xlon(1)
       dlat=ylat(2)-ylat(1)
    else
       dlon=(pi+pi)/real(nlon)
       dlat=pi/real(nlat-1)
    end if
    do iter=1,10
       do l=1,ntiles
          do jc=1,npy-1
             do ic=1,npx-1
                !------------------------------------------------------!
                ! guess latlon indexes for given cubed sphere cell     !
                !------------------------------------------------------!
                call xyz2latlon(xyz_center(:,ic,jc,l), sph_center(:))

                if (subset) then

                   if (sph_center(1) > pi) sph_center(1) = sph_center(1) - (pi+pi)
                   dcub=real(iter)*great_circle(xyz_center(:,ic,jc,l),xyz_center(:,ic+1,jc+1,l))

                  ! Search only if the cubed cell falls within the lat-lon subdomain
                   if ( (sph_center(1) >= xlon(1   )-dcub) .and. &
                        (sph_center(1) <= xlon(nlon)+dcub) .and. &
                        (sph_center(2) >= ylat(1   )-dcub) .and. &
                        (sph_center(2) <= ylat(nlat)+dcub) ) then
 
                      j_min=max(   1,  floor((sph_center(2)-dcub-ylat(1))/dlat)-iter+1)
                      j_max=min(nlat,ceiling((sph_center(2)+dcub-ylat(1))/dlat)+iter-1)
          
                      if (j_min==1 .or. j_max==nlat) then
                         i_min=1
                         i_max=nlon
                      else
                         i_min=max(   1,  floor((sph_center(1)-dcub-xlon(1))/dlon-iter+1))
                         i_max=min(nlon,ceiling((sph_center(1)+dcub-xlon(1))/dlon+iter-1))
                      end if

                      do j=j_min,j_max
                         do i=i_min,i_max
                            !------------------------------------------------!
                            ! for latlon cell find nearest cubed sphere cell !
                            !------------------------------------------------!
                            if (.not. found(i,j)) then
                               shortest(l)=pi+pi
                               do jcc=jc,min(npy-1,jc+1)
                                  do icc=ic,min(npx-1,ic+1)
                                     distance=great_circle(xyz_center(:,icc,jcc,l),xyz_latlon(:,i,j))
                                     if (distance<shortest(l)) then
                                        shortest(l)=distance
                                        index(1,l)=icc
                                        index(2,l)=jcc
                                        index(3,l)=l
                                     endif
                                  enddo
                               enddo
                               !------------------------------------------------!
                               ! determine lower left corner                    !
                               !------------------------------------------------!
                               call get_closest_index(index(1,l), index(2,l), index(3,l), found(i,j))
                            endif
                         enddo
                      enddo
                   end if ! search subdomain

                else ! subset else

                   dcub=real(iter)*great_circle(xyz_center(:,ic,jc,l),xyz_center(:,ic+1,jc+1,l))
  
                   j_min=max(   1,  floor((sph_center(2)-dcub+0.5*pi)/dlat)-iter+1)
                   j_max=min(nlat,ceiling((sph_center(2)+dcub+0.5*pi)/dlat)+iter-1)
                   if (j_min==1 .or. j_max==nlat) then
                      i_min=1
                      i_max=nlon
                   else
                      i_min=max(   1,  floor((sph_center(1)-dcub)/dlon-iter+1))
                      i_max=min(nlon,ceiling((sph_center(1)+dcub)/dlon+iter-1))
                   end if
                   do j=j_min,j_max
                      do i=i_min,i_max
                         !------------------------------------------------!
                         ! for latlon cell find nearest cubed sphere cell !
                         !------------------------------------------------!
                         if (.not. found(i,j)) then
                            shortest(l)=pi+pi
                            do jcc=jc,min(npy-1,jc+1)
                               do icc=ic,min(npx-1,ic+1)
                                  distance=great_circle(xyz_center(:,icc,jcc,l),xyz_latlon(:,i,j))
                                  if (distance<shortest(l)) then
                                     shortest(l)=distance
                                     index(1,l)=icc
                                     index(2,l)=jcc
                                     index(3,l)=l
                                  endif
                               enddo
                            enddo
                            !------------------------------------------------!
                            ! determine lower left corner                    !
                            !------------------------------------------------!
                            call get_closest_index(index(1,l), index(2,l), index(3,l), found(i,j))
                         endif
                      enddo 
                   enddo

                end if ! subset

             enddo
          enddo
       enddo
       if (iter>1) then
          if (all_done()) exit
       endif
    enddo
    !------------------------------------------------------------------!
    ! double check if lower left corner was found                      !
    ! calculate weights for interpolation                              !
    !------------------------------------------------------------------!
    do j=1,nlat
       do i=1,nlon
          if (.not. found(i,j)) then
             print*,"**************************************************************"
             print*,"WARNING: didn't find lower left corner for (ilon,jlat)", i,j
             print*,"will perform expensive global sweep"
             print*,"**************************************************************"
             !---------------------------------------------------------!
             ! for latlon cell find nearby cubed sphere cell           !
             !---------------------------------------------------------!
             index(:,:)=0
             shortest(:)=pi+pi
             do l=1,ntiles
                do jc=1,npx-1
                   do ic=1,npy-1
                      distance=great_circle(xyz_center(:,ic,jc,l),xyz_latlon(:,i,j))
                      if (distance<shortest(l)) then
                         shortest(l)=distance
                         index(1,l)=ic
                         index(2,l)=jc
                         index(3,l)=l
                      endif
                   enddo
                enddo
             enddo
             !---------------------------------------------------------!
             ! determine lower left corner                             !
             !---------------------------------------------------------!
             call sort_index()
             found(i,j)=.false.
             do l=1,ntiles
                if (.not. found(i,j)) then
                   call get_index(index(1,l), index(2,l), index(3,l), found(i,j))
                   if (found(i,j)) exit
                endif
             enddo
             if (.not. found(i,j)) stop "ERROR: couldn't find lower left corner"
          endif
          !------------------------------------------------------------!
          ! calculate shortest distance to each side of rectangle      !
          ! formed by cubed sphere cell centers                        !
          ! special corner treatment                                   !
          !------------------------------------------------------------!
          ic=c2l_index(1,i,j)
          jc=c2l_index(2,i,j)
          l =c2l_index(3,i,j)
          if (ic==npx-1 .and. jc==npy-1) then
             !------------------------------------------------------------!
             ! calculate weights for bilinear interpolation near corner   !
             !------------------------------------------------------------!
             dist1=dist2side(xyz_center(:,ic+1,jc,l),xyz_center(:,ic,jc+1,l),xyz_latlon(:,i,j))
             dist2=dist2side(xyz_center(:,ic+1,jc,l),xyz_center(:,ic,jc  ,l),xyz_latlon(:,i,j))
             dist3=dist2side(xyz_center(:,ic  ,jc,l),xyz_center(:,ic,jc+1,l),xyz_latlon(:,i,j))
             
             c2l_weight(1,i,j)=dist1      ! ic,   jc    weight
             c2l_weight(2,i,j)=dist2      ! ic,   jc+1  weight
             c2l_weight(3,i,j)=0.         ! ic+1, jc+1  weight
             c2l_weight(4,i,j)=dist3      ! ic+1, jc    weight
             
             sum=c2l_weight(1,i,j)+c2l_weight(2,i,j)+c2l_weight(4,i,j)
             c2l_weight(1,i,j)=c2l_weight(1,i,j)/sum
             c2l_weight(2,i,j)=c2l_weight(2,i,j)/sum
             c2l_weight(4,i,j)=c2l_weight(4,i,j)/sum

          elseif (ic==0 .and. jc==npy-1) then
             !------------------------------------------------------------!
             ! calculate weights for bilinear interpolation near corner   !
             !------------------------------------------------------------!
             dist1=dist2side(xyz_center(:,ic+1,jc+1,l),xyz_center(:,ic+1,jc,l),xyz_latlon(:,i,j))
             dist2=dist2side(xyz_center(:,ic+1,jc  ,l),xyz_center(:,ic  ,jc,l),xyz_latlon(:,i,j))
             dist3=dist2side(xyz_center(:,ic+1,jc+1,l),xyz_center(:,ic  ,jc,l),xyz_latlon(:,i,j))
             
             c2l_weight(1,i,j)=dist1      ! ic,   jc    weight
             c2l_weight(2,i,j)=0.         ! ic,   jc+1  weight
             c2l_weight(3,i,j)=dist2      ! ic+1, jc+1  weight
             c2l_weight(4,i,j)=dist3      ! ic+1, jc    weight
             
             sum=c2l_weight(1,i,j)+c2l_weight(3,i,j)+c2l_weight(4,i,j)
             c2l_weight(1,i,j)=c2l_weight(1,i,j)/sum
             c2l_weight(3,i,j)=c2l_weight(3,i,j)/sum
             c2l_weight(4,i,j)=c2l_weight(4,i,j)/sum

          elseif (jc==0 .and. ic==npx-1) then
             !------------------------------------------------------------!
             ! calculate weights for bilinear interpolation near corner   !
             !------------------------------------------------------------!
             dist1=dist2side(xyz_center(:,ic,jc+1,l),xyz_center(:,ic+1,jc+1,l),xyz_latlon(:,i,j))
             dist2=dist2side(xyz_center(:,ic,jc  ,l),xyz_center(:,ic+1,jc+1,l),xyz_latlon(:,i,j))
             dist3=dist2side(xyz_center(:,ic,jc  ,l),xyz_center(:,ic  ,jc+1,l),xyz_latlon(:,i,j))
             
             c2l_weight(1,i,j)=dist1      ! ic,   jc    weight
             c2l_weight(2,i,j)=dist2      ! ic,   jc+1  weight
             c2l_weight(3,i,j)=dist3      ! ic+1, jc+1  weight
             c2l_weight(4,i,j)=0.         ! ic+1, jc    weight
             
             sum=c2l_weight(1,i,j)+c2l_weight(2,i,j)+c2l_weight(3,i,j)
             c2l_weight(1,i,j)=c2l_weight(1,i,j)/sum
             c2l_weight(2,i,j)=c2l_weight(2,i,j)/sum
             c2l_weight(3,i,j)=c2l_weight(3,i,j)/sum

          else
             !------------------------------------------------------------!
             ! calculate weights for bilinear interpolation if no corner  !
             !------------------------------------------------------------!
             dist1=dist2side(xyz_center(:,ic  ,jc  ,l),xyz_center(:,ic  ,jc+1,l),xyz_latlon(:,i,j))
             dist2=dist2side(xyz_center(:,ic  ,jc+1,l),xyz_center(:,ic+1,jc+1,l),xyz_latlon(:,i,j))
             dist3=dist2side(xyz_center(:,ic+1,jc+1,l),xyz_center(:,ic+1,jc  ,l),xyz_latlon(:,i,j))
             dist4=dist2side(xyz_center(:,ic+1,jc  ,l),xyz_center(:,ic  ,jc  ,l),xyz_latlon(:,i,j))
             
             c2l_weight(1,i,j)=dist2*dist3      ! ic,   jc    weight
             c2l_weight(2,i,j)=dist3*dist4      ! ic,   jc+1  weight
             c2l_weight(3,i,j)=dist4*dist1      ! ic+1, jc+1  weight
             c2l_weight(4,i,j)=dist1*dist2      ! ic+1, jc    weight
             
             sum=c2l_weight(1,i,j)+c2l_weight(2,i,j)+c2l_weight(3,i,j)+c2l_weight(4,i,j)
             c2l_weight(:,i,j)=c2l_weight(:,i,j)/sum
          endif
       enddo
    enddo

#if !defined(MAPL_MODE)
    write(*,100) 
100 format (/," done calculating c2l_index and c2l_weight",/)
#endif

  contains
    !------------------------------------------------------------------!
    subroutine sort_index()
      !----------------------------------------------------------------!
      ! sort index by shortest                                         !
      !----------------------------------------------------------------!
      real(REAL8)    :: shortest_sort(ntiles)
      integer :: l, ll, lll, index_sort(3,ntiles)

      index_sort(:,:)=0
      shortest_sort(:)=pi+pi
      do l=1,ntiles
         do ll=1,ntiles
            if (shortest(l)<shortest_sort(ll)) then
               do lll=ntiles-1,ll,-1
                  index_sort(:,lll+1)=index_sort(:,lll)
                  shortest_sort(lll+1)=shortest_sort(lll)
               enddo
               index_sort(:,ll)=index(:,l)
               shortest_sort(ll)=shortest(l)
               exit
            endif
         enddo
      enddo
      shortest(:)=shortest_sort(:)
      index(:,:)=index_sort(:,:)

    end subroutine sort_index
    !------------------------------------------------------------------!
    subroutine get_closest_index(ig, jg, lg, ok)
      !----------------------------------------------------------------!
      ! determine lower left corner                                    !
      !----------------------------------------------------------------!
      integer, intent(in)  :: ig, jg, lg
      logical, intent(out) :: ok

      real(REAL8) :: angle_11, angle_11a, angle_11b,                           & 
              angle_22, angle_22a, angle_22b,                           & 
              angle_33, angle_33a, angle_33b,                           &
              angle_44, angle_44a, angle_44b


      ok=.false.
      angle_1 =spherical_angle(xyz_center(:,ig,jg,lg),xyz_center(:,ig+1,jg  ,lg),xyz_center(:,ig  ,jg+1,lg))
      angle_1a=spherical_angle(xyz_center(:,ig,jg,lg),xyz_center(:,ig+1,jg  ,lg),xyz_latlon(:,i,j))
      angle_1b=spherical_angle(xyz_center(:,ig,jg,lg),xyz_center(:,ig  ,jg+1,lg),xyz_latlon(:,i,j))
      if (max(angle_1a,angle_1b)<=angle_1) then
          if (ig+1==npx .and. jg+1==npy) then
            angle_11 =spherical_angle(xyz_center(:,ig+1,jg,lg),xyz_center(:,ig,jg+1,lg),xyz_center(:,ig,jg,lg))
            angle_11a=spherical_angle(xyz_center(:,ig+1,jg,lg),xyz_center(:,ig,jg  ,lg),xyz_latlon(:,i,j))
            angle_11b=spherical_angle(xyz_center(:,ig+1,jg,lg),xyz_center(:,ig,jg+1,lg),xyz_latlon(:,i,j))
         else
            angle_11 =spherical_angle(xyz_center(:,ig+1,jg+1,lg),xyz_center(:,ig  ,jg+1,lg),xyz_center(:,ig+1,jg,lg))
            angle_11a=spherical_angle(xyz_center(:,ig+1,jg+1,lg),xyz_center(:,ig+1,jg  ,lg),xyz_latlon(:,i,j))
            angle_11b=spherical_angle(xyz_center(:,ig+1,jg+1,lg),xyz_center(:,ig  ,jg+1,lg),xyz_latlon(:,i,j))
         endif
         if (max(angle_11a,angle_11b)<=angle_11) then
            ok=.true.
            c2l_index(1,i,j)=ig
            c2l_index(2,i,j)=jg
            c2l_index(3,i,j)=lg
         endif
      else
         angle_2 =spherical_angle(xyz_center(:,ig,jg,lg),xyz_center(:,ig,jg+1,lg),xyz_center(:,ig-1,jg,lg))
         angle_2a=angle_1b
         angle_2b=spherical_angle(xyz_center(:,ig,jg,lg),xyz_center(:,ig-1,jg,lg),xyz_latlon(:,i,j))
         if (max(angle_2a,angle_2b)<=angle_2) then
            if (ig-1==0 .and. jg+1==npy) then
               angle_22 =spherical_angle(xyz_center(:,ig,jg+1,lg),xyz_center(:,ig  ,jg,lg),xyz_center(:,ig-1,jg,lg))
               angle_22a=spherical_angle(xyz_center(:,ig,jg+1,lg),xyz_center(:,ig-1,jg,lg),xyz_latlon(:,i,j))
               angle_22b=spherical_angle(xyz_center(:,ig,jg+1,lg),xyz_center(:,ig  ,jg,lg),xyz_latlon(:,i,j))
            else
               angle_22 =spherical_angle(xyz_center(:,ig-1,jg+1,lg),xyz_center(:,ig  ,jg+1,lg),xyz_center(:,ig-1,jg,lg))
               angle_22a=spherical_angle(xyz_center(:,ig-1,jg+1,lg),xyz_center(:,ig-1,jg  ,lg),xyz_latlon(:,i,j))
               angle_22b=spherical_angle(xyz_center(:,ig-1,jg+1,lg),xyz_center(:,ig  ,jg+1,lg),xyz_latlon(:,i,j))
            endif
            if (max(angle_22a,angle_22b)<=angle_22) then
               ok=.true.
               c2l_index(1,i,j)=ig-1
               c2l_index(2,i,j)=jg
               c2l_index(3,i,j)=lg
            endif
         else
            angle_3 =spherical_angle(xyz_center(:,ig,jg,lg),xyz_center(:,ig-1,jg,lg),xyz_center(:,ig,jg-1,lg))
            angle_3a=angle_2b
            angle_3b=spherical_angle(xyz_center(:,ig,jg,lg),xyz_center(:,ig,jg-1,lg),xyz_latlon(:,i,j))
            if (max(angle_3a,angle_3b)<=angle_3 .and. ig>1 .and.jg>1) then
               angle_33 =spherical_angle(xyz_center(:,ig-1,jg-1,lg),xyz_center(:,ig  ,jg-1,lg),xyz_center(:,ig-1,jg,lg))
               angle_33a=spherical_angle(xyz_center(:,ig-1,jg-1,lg),xyz_center(:,ig-1,jg  ,lg),xyz_latlon(:,i,j))
               angle_33b=spherical_angle(xyz_center(:,ig-1,jg-1,lg),xyz_center(:,ig  ,jg-1,lg),xyz_latlon(:,i,j))
               if (max(angle_33a,angle_33b)<=angle_33) then
                  ok=.true.
                  c2l_index(1,i,j)=ig-1
                  c2l_index(2,i,j)=jg-1
                  c2l_index(3,i,j)=lg
               endif
            else
               angle_4 =spherical_angle(xyz_center(:,ig,jg,lg),xyz_center(:,ig,jg-1,lg),xyz_center(:,ig+1,jg,lg))
               angle_4a=angle_3b
               angle_4b=spherical_angle(xyz_center(:,ig,jg,lg),xyz_center(:,ig+1,jg,lg),xyz_latlon(:,i,j))
               if (max(angle_4a,angle_4b)<=angle_4) then
                  if (ig+1==npx .and. jg-1==0) then
                     angle_44 =spherical_angle(xyz_center(:,ig+1,jg,lg),xyz_center(:,ig,jg  ,lg),xyz_center(:,ig,jg-1,lg))
                     angle_44a=spherical_angle(xyz_center(:,ig+1,jg,lg),xyz_center(:,ig,jg-1,lg),xyz_latlon(:,i,j))
                     angle_44b=spherical_angle(xyz_center(:,ig+1,jg,lg),xyz_center(:,ig,jg  ,lg),xyz_latlon(:,i,j))
                  else
                     angle_44 =spherical_angle(xyz_center(:,ig+1,jg-1,lg),xyz_center(:,ig+1,jg  ,lg),xyz_center(:,ig,jg-1,lg))
                     angle_44a=spherical_angle(xyz_center(:,ig+1,jg-1,lg),xyz_center(:,ig  ,jg-1,lg),xyz_latlon(:,i,j))
                     angle_44b=spherical_angle(xyz_center(:,ig+1,jg-1,lg),xyz_center(:,ig+1,jg  ,lg),xyz_latlon(:,i,j))
                  endif
                  if (max(angle_44a,angle_44b)<=angle_44) then
                     ok=.true.
                     c2l_index(1,i,j)=ig
                     c2l_index(2,i,j)=jg-1
                     c2l_index(3,i,j)=lg
                  endif
               endif
            endif
         endif
      endif

    end subroutine get_closest_index
    !------------------------------------------------------------------!
    subroutine get_index(ig, jg, lg, ok)
      !----------------------------------------------------------------!
      ! determine lower left corner                                    !
      !----------------------------------------------------------------!
      integer, intent(in)  :: ig, jg, lg
      logical, intent(out) :: ok

      ok=.true.
      angle_1 =spherical_angle(xyz_center(:,ig,jg,lg),xyz_center(:,ig+1,jg  ,lg),xyz_center(:,ig  ,jg+1,lg))
      angle_1a=spherical_angle(xyz_center(:,ig,jg,lg),xyz_center(:,ig+1,jg  ,lg),xyz_latlon(:,i,j))
      angle_1b=spherical_angle(xyz_center(:,ig,jg,lg),xyz_center(:,ig  ,jg+1,lg),xyz_latlon(:,i,j))
      if (max(angle_1a,angle_1b)<angle_1) then
         c2l_index(1,i,j)=ig
         c2l_index(2,i,j)=jg
         c2l_index(3,i,j)=lg
      else
         angle_2 =spherical_angle(xyz_center(:,ig,jg,lg),xyz_center(:,ig,jg+1,lg),xyz_center(:,ig-1,jg,lg))
         angle_2a=angle_1b
         angle_2b=spherical_angle(xyz_center(:,ig,jg,lg),xyz_center(:,ig-1,jg,lg),xyz_latlon(:,i,j))
         if (max(angle_2a,angle_2b)<angle_2) then
            c2l_index(1,i,j)=ig-1
            c2l_index(2,i,j)=jg
            c2l_index(3,i,j)=lg
         else
            angle_3 =spherical_angle(xyz_center(:,ig,jg,lg),xyz_center(:,ig-1,jg,lg),xyz_center(:,ig,jg-1,lg))
            angle_3a=angle_2b
            angle_3b=spherical_angle(xyz_center(:,ig,jg,lg),xyz_center(:,ig,jg-1,lg),xyz_latlon(:,i,j))
            if (max(angle_3a,angle_3b)<angle_3 .and. ig>1 .and.jg>1) then
               c2l_index(1,i,j)=ig-1
               c2l_index(2,i,j)=jg-1
               c2l_index(3,i,j)=lg
            else
               angle_4 =spherical_angle(xyz_center(:,ig,jg,lg),xyz_center(:,ig,jg-1,lg),xyz_center(:,ig+1,jg,lg))
               angle_4a=angle_3b
               angle_4b=spherical_angle(xyz_center(:,ig,jg,lg),xyz_center(:,ig+1,jg,lg),xyz_latlon(:,i,j))
               if (max(angle_4a,angle_4b)<angle_4) then
                  c2l_index(1,i,j)=ig
                  c2l_index(2,i,j)=jg-1
                  c2l_index(3,i,j)=lg
               else
                  ok=.false.
               endif
            endif
         endif
      endif

    end subroutine get_index
    !------------------------------------------------------------------!
    function all_done()
      logical :: all_done

      all_done=.true.
      loop: do i=1,nlon
         do j=1,nlat
            if (.not. found(i,j)) then
               all_done=.false.
               exit loop
            endif
         enddo
      enddo loop

    end function all_done
    !------------------------------------------------------------------!
  end subroutine get_c2l_weight
  !====================================================================!
  subroutine interpolate_data(npx, npy, ntiles, data_file, data_out,     &
                              uname, vname, missing_value, fill_missing, &
                              memphis, c2l_index, c2l_weight,            &
                              elon_cubsph, elat_cubsph,                  &
                              elon_latlon, elat_latlon,                  &
                              xlon, ylat, nlon, nlat, finer_steps)
    !------------------------------------------------------------------!
    ! read cubed sphere data from file, interpolate to latlon,         !
    ! write latlon data to output file                                 !
    !------------------------------------------------------------------!
    use GHOST_CUBSPH_mod, only: A_grid, ghost_cubsph_update

#include "netcdf.inc"
    integer, intent(in) :: npx, npy, ntiles, nlon, nlat, finer_steps
    integer, intent(in) :: c2l_index(3,nlon,nlat)
    real(REAL8),    intent(in) :: c2l_weight(4,nlon,nlat), xlon(nlon), ylat(nlat)
    real(REAL8), dimension(3, 0:npx, 0:npy ,ntiles), intent(in) :: elon_cubsph, elat_cubsph
    real(REAL8), dimension(3, nlon, nlat)          , intent(in) :: elon_latlon, elat_latlon

    logical, intent(in) :: memphis, fill_missing

    character(len=120), intent(in) :: data_file, data_out,              &
                                      uname, vname, missing_value
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    real(REAL4) :: misval_r4
    real(REAL8) :: misval_r8
    real(REAL4), dimension(:), allocatable :: xlon_r4, ylat_r4
    real(REAL8), dimension(:), allocatable :: time
    real(REAL4), dimension(:,:,:), allocatable :: var_r4
    real(REAL8), dimension(:,:,:), allocatable :: var_r8
    integer*2, dimension(:,:,:), allocatable :: var_i2
    real(REAL8), dimension(:), allocatable :: xlon_crs, ylat_crs, varmisval, varscale, varoffset
    real(REAL8), dimension(:,:,:), allocatable :: var_latlon, var_latlon_crs
    real(REAL8), dimension(:,:,:,:), allocatable :: var_cubsph, ua_cubsph, va_cubsph, xyz_latlon

    real(REAL8), parameter :: todeg=180./pi

    integer :: nlon_crs, nlat_crs, ntime, icoor
    integer :: i, j, k, itile, id, iv, ia, it, nx, ny, istart, istop
    integer :: status, ndims, nvars, ngatts, ncid_out,                  &
               time_dim, xt_dim, yt_dim, lon_dim, lat_dim,              &
               time_id, xt_id, yt_id, type, start(4), count(4), n(4),   &
               dimids(4), dimid, dimlen, lon_id, lat_id,                &
               u_id, v_id, name_len, attype, attlen

    integer, dimension(:), allocatable :: ncid_in, vartype,             &
         varndims, varnatts, varid, nlev
    integer, dimension(:,:), allocatable :: vardimids

    logical :: exists
    logical, dimension(:), allocatable :: ignore, cs_static, do_c2l, misval

    character(len=120) :: time_name, filename, dimname, att_name
    character(len=120), dimension(:), allocatable :: varname

#if !defined(MAPL_MODE)

    !------------------------------------------------------------------!
    ! determine target grid resolution                                 !
    !------------------------------------------------------------------!
    if (finer_steps>0) then
       nlon_crs=nlon/2**finer_steps
       nlat_crs=(nlat-1)/2**finer_steps+1
    else
       nlon_crs=nlon
       nlat_crs=nlat
    endif
    !------------------------------------------------------------------!
    ! open netcdf files with input data                                !
    !------------------------------------------------------------------!
    allocate(ncid_in(ntiles))
    do itile=1,ntiles
       write(filename,100) trim(data_file), itile
100    format(a,'.tile',i1,'.nc')
       inquire(file=filename,exist=exists)
       if (.not. exists) then
          print 101, trim(filename)
101       format("data file ",a," doesn't exist" )
          stop
       endif

       status = nf_open(trim(filename), 0, ncid_in(itile))
       if (status /= nf_noerr) then
          print*,"nf_open: could not open file ",trim(filename)
          stop
       endif
    enddo
    write(filename,100) trim(data_file), 1
    !------------------------------------------------------------------!
    ! get basic info about dims and variables from first tile          !
    !------------------------------------------------------------------!
    status = nf_inq(ncid_in(1), ndims, nvars, ngatts, time_dim)
    status = nf_inq_dimname(ncid_in(1), time_dim, time_name)
    status = nf_inq_dimlen(ncid_in(1), time_dim, ntime)
    status = nf_inq_varid(ncid_in(1), trim(time_name), time_id)
    if (status /= nf_noerr) then
       print*,"NO time variable ", trim(time_name), " in file ",trim(filename)
       stop
    endif
    status = nf_inq_vartype(ncid_in(1), time_id, type)
    if (type /= nf_double) then
       print*," time variable ",trim(time_name), " has to be double precision"
       stop
    endif
    start(1)=1
    count(1)=ntime
    allocate(time(ntime))
    status = nf_get_vara_double(ncid_in(1), time_id, start, count, time)
    !------------------------------------------------------------------!
    ! check grid dimension and variables                               !
    !------------------------------------------------------------------!
    status = nf_inq_dimid(ncid_in(1), "grid_xt", xt_dim)
    if (status /= nf_noerr) then
       print*,"NO grid_xt dimension in file ",trim(filename)
       stop
    endif
    status = nf_inq_dimid(ncid_in(1), "grid_yt", yt_dim)
    if (status /= nf_noerr) then
       print*,"NO grid_yt dimension in file ",trim(filename)
       stop
    endif
    status = nf_inq_dimlen(ncid_in(1), xt_dim, nx)
    status = nf_inq_dimlen(ncid_in(1), yt_dim, ny)
    if (nx/=npx-1 .or. ny/=npy-1) then
       print*,"grid_xt, grid_yt don't have expected length"
       print*, nx, ny, "instead of", npx-1, npy-1
       stop
    endif
    status = nf_inq_varid(ncid_in(1), "grid_xt", xt_id)
    if (status /= nf_noerr) then
       print*,"NO grid_xt variable in file ",trim(filename)
       stop
    endif
    status = nf_inq_varid(ncid_in(1), "grid_yt", yt_id)
    if (status /= nf_noerr) then
       print*,"NO grid_xt variable in file ",trim(filename)
       stop
    endif

    !------------------------------------------------------------------!
    ! open netcdf file for output data                                 !
    !------------------------------------------------------------------!
    filename = adjustl(data_out)
    istart = index(filename,"/",  .true.)+1
    istop  = len_trim(filename)
    filename = filename(istart:istop)//".nc"
    status = nf_create(filename, nf_clobber, ncid_out)
    !------------------------------------------------------------------!
    ! define dimensions                                                !
    !------------------------------------------------------------------!
    do iv=1,ndims
       if (iv==xt_dim) then
          status = nf_def_dim(ncid_out, "lon", nlon_crs, lon_dim)
       elseif (iv==yt_dim) then
          status = nf_def_dim(ncid_out, "lat", nlat_crs, lat_dim)
       elseif (iv==time_dim) then
          status = nf_def_dim(ncid_out, trim(time_name), nf_unlimited, dimid) 
       else
          status = nf_inq_dimname(ncid_in(1), iv, dimname)
          status = nf_inq_dimlen(ncid_in(1), iv, dimlen)
          status = nf_def_dim(ncid_out, trim(dimname), dimlen, dimid)
       endif
    enddo
    !------------------------------------------------------------------!
    ! define variables and attributes                                  !
    !------------------------------------------------------------------!
    allocate(varname(nvars))
    allocate(vartype(nvars), varndims(nvars), varnatts(nvars), varid(nvars), &
             vardimids(5,nvars), nlev(nvars))
    allocate(ignore(nvars), cs_static(nvars), do_c2l(nvars))
    allocate(misval(nvars), varmisval(nvars), varscale(nvars), varoffset(nvars))
    do iv=1,nvars
       status = nf_inq_var(ncid_in(1), iv, varname(iv), vartype(iv),         &
                          varndims(iv), vardimids(1,iv), varnatts(iv))
       status = nf_inq_att(ncid_in(1), iv, trim(missing_value), attype, attlen)
       if (status == nf_noerr) then
          if (attype==nf_double) then 
             status = nf_get_att_double(ncid_in(1), iv, trim(missing_value), misval_r8)
             varmisval(iv)=misval_r8
             misval(iv) = .true.
          elseif (attype==nf_float) then
             status = nf_get_att_real(ncid_in(1), iv, trim(missing_value), misval_r4)
             varmisval(iv)=misval_r4
             misval(iv) = .true.
          else
             misval(iv) = .false.
          endif
       else
          misval(iv) = .false.
       endif

       ! check for packing
       if (vartype(iv)==nf_short) then
          ! get scale_factor
          status = nf_inq_att(ncid_in(1), iv, "scale_factor", attype, attlen)
          if (status == nf_noerr) then
             if (attype==nf_double) then 
                status = nf_get_att_double(ncid_in(1), iv, "scale_factor", misval_r8)
                varscale(iv)=misval_r8
             elseif (attype==nf_float) then
                status = nf_get_att_real(ncid_in(1), iv, "scale_factor", misval_r4)
                varscale(iv)=misval_r4
             else
                varscale(iv)=1.
             endif
          endif
          ! get add_offset
          status = nf_inq_att(ncid_in(1), iv, "add_offset", attype, attlen)
          if (status == nf_noerr) then
             if (attype==nf_double) then 
                status = nf_get_att_double(ncid_in(1), iv, "add_offset", misval_r8)
                varoffset(iv)=misval_r8
             elseif (attype==nf_float) then
                status = nf_get_att_real(ncid_in(1), iv, "add_offset", misval_r4)
                varoffset(iv)=misval_r4
             else
                varoffset(iv)=0.
             endif
          endif
       endif

       if (varndims(iv)>4) then
          ignore(iv)=.true.
          do_c2l(iv)=.false.
          cs_static(iv)=.true.
          print*, "WARNING: will ignore variable ",trim(varname(iv))
       elseif (iv==xt_id) then
          ignore(iv)=.true.
          do_c2l(iv)=.false.
          cs_static(iv)=.true.

          status = nf_def_var(ncid_out, "lon", nf_float, 1,  lon_dim,  lon_id)
          status = nf_put_att_text(ncid_out, lon_id, "long_name",     16, "T-cell longitude")
          status = nf_put_att_text(ncid_out, lon_id, "units"    ,      9, "degrees_E")
          status = nf_put_att_text(ncid_out, lon_id, "cartesian_axis", 1, "X")
       elseif (iv==yt_id) then
          ignore(iv)=.true.
          do_c2l(iv)=.false.
          cs_static(iv)=.true.

          status = nf_def_var(ncid_out, "lat", nf_float, 1,  lat_dim,  lat_id)
          status = nf_put_att_text(ncid_out, lat_id, "long_name",     15, "T-cell latitude")
          status = nf_put_att_text(ncid_out, lat_id, "units"    ,      9, "degrees_N")
          status = nf_put_att_text(ncid_out, lat_id, "cartesian_axis", 1, "Y")
       else
          ignore(iv)=.false.

          if (vardimids(varndims(iv),iv)==time_dim) then
             cs_static(iv)=.false.
          else
             cs_static(iv)=.true.
          endif

          do_c2l(iv)=.false.
          if (varndims(iv)>1) then
             if (vardimids(1,iv)==xt_dim .and. vardimids(2,iv)==yt_dim) do_c2l(iv)=.true.
             if (varndims(iv)>2 .and. vardimids(3,iv)/=time_dim) then
                status = nf_inq_dimlen(ncid_in(1), vardimids(3,iv), nlev(iv))
             else
                nlev(iv)=1
             endif
          endif

          if (.not. do_c2l(iv)) then
             status = nf_def_var(ncid_out, trim(varname(iv)), vartype(iv), &
                                 varndims(iv), vardimids(1,iv), varid(iv))
          else
             dimids(1)=lon_dim
             dimids(2)=lat_dim
             do id=3,varndims(iv)
                dimids(id)=vardimids(id,iv)
             enddo
             status = nf_def_var(ncid_out, trim(varname(iv)), vartype(iv), &
                                    varndims(iv), dimids, varid(iv))
          endif

          do ia=1,varnatts(iv)
             status = nf_inq_attname(ncid_in(1), iv, ia, att_name)
             status = nf_copy_att(ncid_in(1), iv, att_name, ncid_out, varid(iv))
          enddo
       endif
    enddo
    !------------------------------------------------------------------!
    ! copy global attributes                                           !
    !------------------------------------------------------------------!
    do ia=1,ngatts
       status = nf_inq_attname(ncid_in(1), nf_global, ia, att_name)
       if (trim(att_name)=="filename") then
          name_len=len(trim(filename))
          status = nf_put_att_text(ncid_out, nf_global, trim(att_name),  &
               name_len, trim(filename))
       else
          status = nf_copy_att(ncid_in(1), nf_global, att_name, ncid_out, nf_global)
       endif
    enddo
    status = nf_put_att_text(ncid_out, nf_global, "description", 45,     &
         "data interpolated from cubed sphere to latlon")
    status = nf_enddef(ncid_out)
    !------------------------------------------------------------------!
    ! check for vector quantities uname, vname                         !
    !------------------------------------------------------------------!
    status = nf_inq_varid(ncid_in(1), trim(uname), u_id)
    if (status /= nf_noerr) then
       u_id=0
       v_id=0
       print*,"WARNING: NO zonal flow variable", trim(uname), "found in file ",trim(filename)
    else
       status = nf_inq_varid(ncid_in(1), trim(vname), v_id)
       if (status /= nf_noerr) then
          u_id=0
          v_id=0
          print*,"WARNING: NO meridional flow variable ", trim(vname), " found in file ",trim(filename)
       else
          if (nlev(u_id)/=nlev(v_id) .or. vartype(u_id)/=vartype(v_id)) then
             print*,"WARNING: ", trim(uname), " and ", trim(vname),        &
                    "have different dimensions or types"
          elseif (cs_static(u_id) .or. cs_static(v_id)) then
             print*,"WARNING: ", trim(uname), " and/or ", trim(vname),        &
                    "are static arrays"
          elseif (varndims(u_id)/=4 .or. varndims(v_id)/=4) then
             print*,"WARNING: ", trim(uname), " and/or ", trim(vname),        &
                    "don't have 4 dimensions",  varndims(u_id), varndims(v_id)
          endif
       endif
    endif
    !------------------------------------------------------------------!
    ! write latlon grid                                                !
    !------------------------------------------------------------------!
    if (finer_steps==0) then
       allocate(xlon_r4(nlon), ylat_r4(nlat))
       xlon_r4(:)=todeg*xlon(:)
       ylat_r4(:)=todeg*ylat(:)
    else
       allocate(xlon_crs(nlon_crs), ylat_crs(nlat_crs),                 &
                xlon_r4 (nlon_crs), ylat_r4 (nlat_crs))
       call init_latlon_grid(xlon_crs, ylat_crs, nlon_crs, nlat_crs)
       xlon_r4(:)=todeg*xlon_crs(:)
       ylat_r4(:)=todeg*ylat_crs(:)
    endif
    if (memphis) then
       ylat_r4(   1)=0.5*(ylat_r4(   1)+ylat_r4(     2))
       ylat_r4(nlat)=0.5*(ylat_r4(nlat)+ylat_r4(nlat-1))
    endif
    status = nf_put_var_real(ncid_out, lon_id, xlon_r4)
    status = nf_put_var_real(ncid_out, lat_id, ylat_r4)
    !------------------------------------------------------------------!
    ! start loop over static variables                                 !
    !------------------------------------------------------------------!
    do iv=1,nvars
       if (.not.ignore(iv) .and. cs_static(iv)) then
          if (do_c2l(iv)) then
             start(1:3)=(/1,1,1/)
             count(1:3)=(/npx-1,npy-1,nlev(iv)/)
             allocate(var_cubsph(0:npx,0:npy,nlev(iv),ntiles))
             call init_corners(var_cubsph, nlev(iv))
             allocate(var_latlon(nlon,nlat,nlev(iv)))
             if (vartype(iv)==nf_double) then
                allocate(var_r8(npx-1,npy-1,nlev(iv)))
                do itile=1,ntiles
                   status = nf_get_vara_double(ncid_in(itile), iv, start, count, var_r8)
                   var_cubsph(1:npx-1,1:npy-1,1:nlev(iv),itile)=var_r8(:,:,:)
                enddo
                deallocate(var_r8)
             elseif (vartype(iv)==nf_float) then
                allocate(var_r4(npx-1,npy-1,nlev(iv)))
                do itile=1,ntiles
                   status = nf_get_vara_real(ncid_in(itile), iv, start, count, var_r4)
                   var_cubsph(1:npx-1,1:npy-1,1:nlev(iv),itile)=var_r4(:,:,:)
                enddo
                deallocate(var_r4)
             else
                print*," vartype neither nf_double nor nf_float: ",vartype(iv)
                stop
             endif
             do itile=1,ntiles
                call ghost_cubsph_update(var_cubsph, 0, npx, 0, npy, nlev(iv), 1, ntiles,  &
                                         1, nlev(iv), itile, A_grid)
                call do_c2l_interpolation(var_cubsph(:,:,:,itile), 0, npx, 0, npy, nlev(iv), itile, nlev(iv), &
                                          var_latlon, nlon, nlat, c2l_index, c2l_weight,                      &
                                          misval(iv), varmisval(iv), fill_missing)
             enddo
             if (finer_steps==0) then
                count(1:3)=(/nlon,nlat,nlev(iv)/)
                if (vartype(iv)==nf_double) then
                   allocate(var_r8(nlon,nlat,nlev(iv)))
                   var_r8(:,:,:)=var_latlon(:,:,:)
                   status = nf_put_vara_double(ncid_out, varid(iv), start, count, var_r8)
                   deallocate(var_r8)
                elseif (vartype(iv)==nf_float) then
                   allocate(var_r4(nlon,nlat,nlev(iv)))
                   var_r4(:,:,:)=var_latlon(:,:,:)
                   status = nf_put_vara_real(ncid_out, varid(iv), start, count, var_r4)
                   deallocate(var_r4)
                elseif (vartype(iv)==nf_short) then
                   allocate(var_i2(nlon,nlat,nlev(iv)))
                   var_i2(:,:,:)=(var_latlon(:,:,:)-varoffset(iv))/varscale(iv)
                   status = nf_put_vara_real(ncid_out, varid(iv), start, count, var_i2)
                   deallocate(var_i2)
                endif
             else
                allocate(var_latlon_crs(nlon_crs,nlat_crs,nlev(iv)))
                call do_latlon_coarsening(var_latlon, ylat, nlon, nlat, nlev(iv),          &
                                          var_latlon_crs, nlon_crs, nlat_crs, finer_steps, &
                                          misval(iv), varmisval(iv))
                count(1:3)=(/nlon_crs,nlat_crs,nlev(iv)/)
                if (vartype(iv)==nf_double) then
                   allocate(var_r8(nlon_crs,nlat_crs,nlev(iv)))
                   var_r8(:,:,:)=var_latlon_crs(:,:,:)
                   status = nf_put_vara_double(ncid_out, varid(iv), start, count, var_r8)
                   deallocate(var_r8, var_latlon_crs)
                elseif (vartype(iv)==nf_float) then
                   allocate(var_r4(nlon_crs,nlat_crs,nlev(iv)))
                   var_r4(:,:,:)=var_latlon_crs(:,:,:)
                   status = nf_put_vara_real(ncid_out, varid(iv), start, count, var_r4)
                   deallocate(var_r4, var_latlon_crs)
                elseif (vartype(iv)==nf_short) then
                   allocate(var_i2(nlon_crs,nlat_crs,nlev(iv)))
                   var_i2(:,:,:)=(var_latlon_crs(:,:,:)-varoffset(iv))/varscale(iv)
                   status = nf_put_vara_int2(ncid_out, varid(iv), start, count, var_i2)
                   deallocate(var_i2)
                endif
             endif
             deallocate(var_cubsph, var_latlon)
          else
             n(1:4)=1
             do id=1,varndims(iv)
                status = nf_inq_dimlen(ncid_in(1), vardimids(id,iv), n(id))
             enddo
             start(1:3)=(/1,1,1/)
             count(1:3)=(/n(1),n(2),n(3)/)
             if (vartype(iv)==nf_double) then
                allocate(var_r8(n(1),n(2),n(3)))
                status = nf_get_vara_double(ncid_in(1), iv, start, count, var_r8)
                status = nf_put_vara_double(ncid_out, varid(iv), start, count, var_r8)
                deallocate(var_r8)
             elseif (vartype(iv)==nf_float) then
                allocate(var_r4(n(1),n(2),n(3)))
                status = nf_get_vara_real(ncid_in(1), iv, start, count, var_r4)
                status = nf_put_vara_real(ncid_out, varid(iv), start, count, var_r4)
                deallocate(var_r4)
             elseif (vartype(iv)==nf_short) then
                allocate(var_i2(n(1),n(2),n(3)))
                status = nf_get_vara_int2(ncid_in(1), iv, start, count, var_i2)
                status = nf_put_vara_int2(ncid_out, varid(iv), start, count, var_i2)
                deallocate(var_i2)
             else
                print*," vartype neither nf_double nor nf_float: ",vartype(iv)
                stop
             endif
          endif
       endif
    enddo
    !------------------------------------------------------------------!
    ! start loop over time dependent arrays                            !
    !------------------------------------------------------------------!
    do it=1,ntime
       status = nf_put_var1_double(ncid_out, time_id, it, time(it))
       do iv=1,nvars
          if (iv==u_id) then
             !---------------------------------------------------------!
             ! read horizontal flow                                    !
             !---------------------------------------------------------!
             if (varndims(iv)==4) then
                start(1:4)=(/1,1,1,it/)
                count(1:4)=(/npx-1,npy-1,nlev(iv),1/)
                allocate(ua_cubsph(0:npx,0:npy,nlev(iv),ntiles), va_cubsph(0:npx,0:npy,nlev(iv),ntiles))
                call init_corners(ua_cubsph, nlev(iv))
                call init_corners(va_cubsph, nlev(iv))
                allocate(var_latlon(nlon,nlat,nlev(iv)))
                if (vartype(iv)==nf_double) then
                   allocate(var_r8(npx-1,npy-1,nlev(iv)))
                   do itile=1,ntiles
                      status = nf_get_vara_double(ncid_in(itile), u_id, start, count, var_r8)
                      ua_cubsph(1:npx-1,1:npy-1,1:nlev(iv),itile)=var_r8(:,:,:)
                      status = nf_get_vara_double(ncid_in(itile), v_id, start, count, var_r8)
                      va_cubsph(1:npx-1,1:npy-1,1:nlev(iv),itile)=var_r8(:,:,:)
                   enddo
                   deallocate(var_r8)
                elseif (vartype(iv)==nf_float) then
                   allocate(var_r4(npx-1,npy-1,nlev(iv)))
                   do itile=1,ntiles
                      status = nf_get_vara_real(ncid_in(itile), u_id, start, count, var_r4)
                      ua_cubsph(1:npx-1,1:npy-1,1:nlev(iv),itile)=var_r4(:,:,:)
                      status = nf_get_vara_real(ncid_in(itile), v_id, start, count, var_r4)
                      va_cubsph(1:npx-1,1:npy-1,1:nlev(iv),itile)=var_r4(:,:,:)
                   enddo
                   deallocate(var_r4)
                elseif (vartype(iv)==nf_short) then
                   allocate(var_i2(npx-1,npy-1,nlev(iv)))
                   do itile=1,ntiles
                      status = nf_get_vara_int2(ncid_in(itile), u_id, start, count, var_i2)
                      ua_cubsph(1:npx-1,1:npy-1,1:nlev(iv),itile)=varscale(iv)*var_i2(:,:,:)+varoffset(iv)
                      status = nf_get_vara_int2(ncid_in(itile), v_id, start, count, var_i2)
                      va_cubsph(1:npx-1,1:npy-1,1:nlev(iv),itile)=varscale(iv)*var_i2(:,:,:)+varoffset(iv)
                   enddo
                   deallocate(var_i2)
                else
                   print*," vartype neither nf_double nor nf_float of nf_short: ",vartype(iv)
                   stop
                endif
                !------------------------------------------------------!
                ! calculate and interpolate xyz flow                   !
                !------------------------------------------------------!
                allocate(var_cubsph(0:npx,0:npy,nlev(iv),ntiles), xyz_latlon(nlon,nlat,nlev(iv),3))
                call init_corners(var_cubsph, nlev(iv))
                do icoor=1,3
                   do itile=1,ntiles
                      do k=1,nlev(iv)
                         do j=1,npy-1
                            do i=1,npx-1
                               var_cubsph(i,j,k,itile)=ua_cubsph(i,j,k,itile)*elon_cubsph(icoor,i,j,itile) &
                                                      +va_cubsph(i,j,k,itile)*elat_cubsph(icoor,i,j,itile)
                            enddo
                         enddo
                      enddo
                   enddo
                   do itile=1,ntiles
                      call ghost_cubsph_update(var_cubsph, 0, npx, 0, npy, nlev(iv), 1, ntiles,  &
                                               1, nlev(iv), itile, A_grid)
                      call do_c2l_interpolation(var_cubsph(:,:,:,itile), 0, npx, 0, npy, nlev(iv), itile, nlev(iv), &
                                                xyz_latlon(:,:,:,icoor), nlon, nlat, c2l_index, c2l_weight,         &
                                                misval(iv), varmisval(iv), fill_missing)
                   enddo
                enddo
                deallocate(ua_cubsph, va_cubsph, var_cubsph)
                !------------------------------------------------------!
                ! calculate and write zonal latlon flow                !
                !------------------------------------------------------!
                do k=1,nlev(iv)
                   do j=1,nlat
                      do i=1,nlon
                         var_latlon(i,j,k)=xyz_latlon(i,j,k,1)*elon_latlon(1,i,j) &
                                          +xyz_latlon(i,j,k,2)*elon_latlon(2,i,j) &
                                          +xyz_latlon(i,j,k,3)*elon_latlon(3,i,j)
                      enddo
                   enddo
                enddo
                if (finer_steps==0) then
                   start(1:4)=(/1,1,1,it/)
                   count(1:4)=(/nlon,nlat,nlev(iv),1/)
                   if (vartype(iv)==nf_double) then
                      allocate(var_r8(nlon,nlat,nlev(iv)))
                      var_r8(:,:,:)=var_latlon(:,:,:)
                      status = nf_put_vara_double(ncid_out, varid(iv), start, count, var_r8)
                      deallocate(var_r8)
                   elseif (vartype(iv)==nf_float) then
                      allocate(var_r4(nlon,nlat,nlev(iv)))
                      var_r4(:,:,:)=var_latlon(:,:,:)
                      status = nf_put_vara_real(ncid_out, varid(iv), start, count, var_r4)
                      deallocate(var_r4)
                   elseif (vartype(iv)==nf_short) then
                      allocate(var_i2(nlon,nlat,nlev(iv)))
                      var_i2(:,:,:)=(var_latlon(:,:,:)-varoffset(iv))/varscale(iv)
                      status = nf_put_vara_int2(ncid_out, varid(iv), start, count, var_i2)
                      deallocate(var_i2)
                   endif
                else
                   allocate(var_latlon_crs(nlon_crs,nlat_crs,nlev(iv)))
                   call do_latlon_coarsening(var_latlon, ylat, nlon, nlat, nlev(iv),          &
                                             var_latlon_crs, nlon_crs, nlat_crs, finer_steps, &
                                             misval(u_id), varmisval(u_id))
                   start(1:4)=(/1,1,1,it/)
                   count(1:4)=(/nlon_crs,nlat_crs,nlev(iv),1/)
                   if (vartype(iv)==nf_double) then
                      allocate(var_r8(nlon_crs,nlat_crs,nlev(iv)))
                      var_r8(:,:,:)=var_latlon_crs(:,:,:)
                      status = nf_put_vara_double(ncid_out, varid(iv), start, count, var_r8)
                      deallocate(var_r8)
                   elseif (vartype(iv)==nf_float) then
                      allocate(var_r4(nlon_crs,nlat_crs,nlev(iv)))
                      var_r4(:,:,:)=var_latlon_crs(:,:,:)
                      status = nf_put_vara_real(ncid_out, varid(iv), start, count, var_r4)
                      deallocate(var_r4)
                   elseif (vartype(iv)==nf_short) then
                      allocate(var_i2(nlon_crs,nlat_crs,nlev(iv)))
                      var_i2(:,:,:)=(var_latlon_crs(:,:,:)-varoffset(iv))/varscale(iv)
                      status = nf_put_vara_int2(ncid_out, varid(iv), start, count, var_i2)
                      deallocate(var_i2)
                   endif
                   deallocate(var_latlon_crs)
                endif
                !------------------------------------------------------!
                ! calculate and write meridional latlon flow           !
                !------------------------------------------------------!
                do k=1,nlev(iv)
                   do j=1,nlat
                      do i=1,nlon
                         var_latlon(i,j,k)=xyz_latlon(i,j,k,1)*elat_latlon(1,i,j) &
                                          +xyz_latlon(i,j,k,2)*elat_latlon(2,i,j) &
                                          +xyz_latlon(i,j,k,3)*elat_latlon(3,i,j)
                      enddo
                   enddo
                enddo
                if (finer_steps==0) then
                   start(1:4)=(/1,1,1,it/)
                   count(1:4)=(/nlon,nlat,nlev(iv),1/)
                   if (vartype(iv)==nf_double) then
                      allocate(var_r8(nlon,nlat,nlev(iv)))
                      var_r8(:,:,:)=var_latlon(:,:,:)
                      status = nf_put_vara_double(ncid_out, varid(v_id), start, count, var_r8)
                      deallocate(var_r8)
                   elseif (vartype(iv)==nf_float) then
                      allocate(var_r4(nlon,nlat,nlev(iv)))
                      var_r4(:,:,:)=var_latlon(:,:,:)
                      status = nf_put_vara_real(ncid_out, varid(v_id), start, count, var_r4)
                      deallocate(var_r4)
                   elseif (vartype(iv)==nf_short) then
                      allocate(var_i2(nlon,nlat,nlev(iv)))
                      var_i2(:,:,:)=(var_latlon(:,:,:)-varoffset(iv))/varscale(iv)
                      status = nf_put_vara_int2(ncid_out, varid(v_id), start, count, var_i2)
                      deallocate(var_i2)
                   endif
                else
                   allocate(var_latlon_crs(nlon_crs,nlat_crs,nlev(iv)))
                   call do_latlon_coarsening(var_latlon, ylat, nlon, nlat, nlev(iv),          &
                                             var_latlon_crs, nlon_crs, nlat_crs, finer_steps, &
                                             misval(v_id), varmisval(v_id))
                   start(1:4)=(/1,1,1,it/)
                   count(1:4)=(/nlon_crs,nlat_crs,nlev(iv),1/)
                   if (vartype(iv)==nf_double) then
                      allocate(var_r8(nlon_crs,nlat_crs,nlev(iv)))
                      var_r8(:,:,:)=var_latlon_crs(:,:,:)
                      status = nf_put_vara_double(ncid_out, varid(v_id), start, count, var_r8)
                      deallocate(var_r8)
                   elseif (vartype(iv)==nf_float) then
                      allocate(var_r4(nlon_crs,nlat_crs,nlev(iv)))
                      var_r4(:,:,:)=var_latlon_crs(:,:,:)
                      status = nf_put_vara_real(ncid_out, varid(v_id), start, count, var_r4)
                      deallocate(var_r4)
                   elseif (vartype(iv)==nf_short) then
                      allocate(var_i2(nlon_crs,nlat_crs,nlev(iv)))
                      var_i2(:,:,:)=(var_latlon_crs(:,:,:)-varoffset(iv))/varscale(iv)
                      status = nf_put_vara_int2(ncid_out, varid(v_id), start, count, var_i2)
                      deallocate(var_i2)
                   endif
                   deallocate(var_latlon_crs)
                endif
                deallocate(var_latlon, xyz_latlon)
             endif
          elseif (iv/=v_id) then
             if (.not.ignore(iv) .and. .not.cs_static(iv)) then
                if (do_c2l(iv)) then
                   if (varndims(iv)==3) then
                      start(1:3)=(/1,1,it/)
                      count(1:3)=(/npx-1,npy-1,1/)
                      allocate(var_cubsph(0:npx,0:npy,nlev(iv),ntiles))
                      call init_corners(var_cubsph, nlev(iv))
                      allocate(var_latlon(nlon,nlat,nlev(iv)))
                      if (vartype(iv)==nf_double) then
                         allocate(var_r8(npx-1,npy-1,1))
                         do itile=1,ntiles
                            status = nf_get_vara_double(ncid_in(itile), iv, start, count, var_r8)
                            var_cubsph(1:npx-1,1:npy-1,1,itile)=var_r8(:,:,1)
                         enddo
                         deallocate(var_r8)
                      elseif (vartype(iv)==nf_float) then
                         allocate(var_r4(npx-1,npy-1,1))
                         do itile=1,ntiles
                            status = nf_get_vara_real(ncid_in(itile), iv, start, count, var_r4)
                            var_cubsph(1:npx-1,1:npy-1,1,itile)=var_r4(:,:,1)
                         enddo
                         deallocate(var_r4)
                      elseif (vartype(iv)==nf_short) then
                         allocate(var_i2(npx-1,npy-1,1))
                         do itile=1,ntiles
                            status = nf_get_vara_int2(ncid_in(itile), iv, start, count, var_i2)
                            var_cubsph(1:npx-1,1:npy-1,1,itile)=varscale(iv)*var_i2(:,:,1)+varoffset(iv)
                         enddo
                         deallocate(var_i2)
                      else
                         print*," vartype neither nf_double nor nf_float of nf_short: ",vartype(iv)
                         stop
                      endif
                   elseif (varndims(iv)==4) then
                      start(1:4)=(/1,1,1,it/)
                      count(1:4)=(/npx-1,npy-1,nlev(iv),1/)
                      allocate(var_cubsph(0:npx,0:npy,nlev(iv),ntiles))
                      call init_corners(var_cubsph, nlev(iv))
                      allocate(var_latlon(nlon,nlat,nlev(iv)))
                      if (vartype(iv)==nf_double) then
                         allocate(var_r8(npx-1,npy-1,nlev(iv)))
                         do itile=1,ntiles
                            status = nf_get_vara_double(ncid_in(itile), iv, start, count, var_r8)
                            var_cubsph(1:npx-1,1:npy-1,1:nlev(iv),itile)=var_r8(:,:,:)
                         enddo
                         deallocate(var_r8)
                      elseif (vartype(iv)==nf_float) then
                         allocate(var_r4(npx-1,npy-1,nlev(iv)))
                         do itile=1,ntiles
                            status = nf_get_vara_real(ncid_in(itile), iv, start, count, var_r4)
                            var_cubsph(1:npx-1,1:npy-1,1:nlev(iv),itile)=var_r4(:,:,:)
                         enddo
                         deallocate(var_r4)
                      elseif (vartype(iv)==nf_short) then
                         allocate(var_i2(npx-1,npy-1,nlev(iv)))
                         do itile=1,ntiles
                            status = nf_get_vara_int2(ncid_in(itile), iv, start, count, var_i2)
                            var_cubsph(1:npx-1,1:npy-1,1:nlev(iv),itile)=varscale(iv)*var_i2(:,:,:)+varoffset(iv)
                         enddo
                         deallocate(var_i2)
                      else
                         print*," vartype neither nf_double nor nf_float of nf_short: ",vartype(iv)
                         stop
                      endif
                   else
                      print*," unexpected number of dimensions (", varndims(iv), &
                             ") for variable ",trim(varname(iv))
                      stop
                   endif
                   
                   do itile=1,ntiles
                      call ghost_cubsph_update(var_cubsph, 0, npx, 0, npy, nlev(iv), 1, ntiles,  &
                                               1, nlev(iv), itile, A_grid)
                      call do_c2l_interpolation(var_cubsph(:,:,:,itile), 0, npx, 0, npy, nlev(iv), itile, nlev(iv), &
                                                var_latlon, nlon, nlat, c2l_index, c2l_weight,                      &
                                                misval(iv), varmisval(iv), fill_missing)
                   enddo
                
                   if (finer_steps==0) then
                      if (varndims(iv)==3) then
                         start(1:3)=(/1,1,it/)
                         count(1:3)=(/nlon,nlat,1/)
                         if (vartype(iv)==nf_double) then
                            allocate(var_r8(nlon,nlat,1))
                            var_r8(:,:,1)=var_latlon(:,:,1)
                            status = nf_put_vara_double(ncid_out, varid(iv), start, count, var_r8)
                            deallocate(var_r8)
                         elseif (vartype(iv)==nf_float) then
                            allocate(var_r4(nlon,nlat,1))
                            var_r4(:,:,1)=var_latlon(:,:,1)
                            status = nf_put_vara_real(ncid_out, varid(iv), start, count, var_r4)
                            deallocate(var_r4)
                         elseif (vartype(iv)==nf_short) then
                            allocate(var_i2(nlon,nlat,1))
                            var_i2(:,:,1)=(var_latlon(:,:,1)-varoffset(iv))/varscale(iv)
                            status = nf_put_vara_int2(ncid_out, varid(iv), start, count, var_i2)
                            deallocate(var_i2)
                         endif
                      else
                         start(1:4)=(/1,1,1,it/)
                         count(1:4)=(/nlon,nlat,nlev(iv),1/)
                         if (vartype(iv)==nf_double) then
                            allocate(var_r8(nlon,nlat,nlev(iv)))
                            var_r8(:,:,:)=var_latlon(:,:,:)
                            status = nf_put_vara_double(ncid_out, varid(iv), start, count, var_r8)
                            deallocate(var_r8)
                         elseif (vartype(iv)==nf_float) then
                            allocate(var_r4(nlon,nlat,nlev(iv)))
                            var_r4(:,:,:)=var_latlon(:,:,:)
                            status = nf_put_vara_real(ncid_out, varid(iv), start, count, var_r4)
                            deallocate(var_r4)
                         elseif (vartype(iv)==nf_short) then
                            allocate(var_i2(nlon,nlat,nlev(iv)))
                            var_i2(:,:,:)=(var_latlon(:,:,:)-varoffset(iv))/varscale(iv)
                            status = nf_put_vara_int2(ncid_out, varid(iv), start, count, var_i2)
                            deallocate(var_i2)
                         endif
                      endif
                   else
                      allocate(var_latlon_crs(nlon_crs,nlat_crs,nlev(iv)))
                      call do_latlon_coarsening(var_latlon, ylat, nlon, nlat, nlev(iv),          &
                           var_latlon_crs, nlon_crs, nlat_crs, finer_steps, &
                           misval(iv), varmisval(iv))
                      if (varndims(iv)==3) then
                         start(1:3)=(/1,1,it/)
                         count(1:3)=(/nlon_crs,nlat_crs,1/)
                         if (vartype(iv)==nf_double) then
                            allocate(var_r8(nlon_crs,nlat_crs,1))
                            var_r8(:,:,1)=var_latlon_crs(:,:,1)
                            status = nf_put_vara_double(ncid_out, varid(iv), start, count, var_r8)
                            deallocate(var_r8)
                         elseif (vartype(iv)==nf_float) then
                            allocate(var_r4(nlon_crs,nlat_crs,1))
                            var_r4(:,:,1)=var_latlon_crs(:,:,1)
                            status = nf_put_vara_real(ncid_out, varid(iv), start, count, var_r4)
                            deallocate(var_r4)
                         elseif (vartype(iv)==nf_short) then
                            allocate(var_i2(nlon_crs,nlat_crs,1))
                            var_i2(:,:,1)=(var_latlon_crs(:,:,1)-varoffset(iv))/varscale(iv)
                            status = nf_put_vara_int2(ncid_out, varid(iv), start, count, var_i2)
                            deallocate(var_i2)
                         endif
                      else
                         start(1:4)=(/1,1,1,it/)
                         count(1:4)=(/nlon_crs,nlat_crs,nlev(iv),1/)
                         if (vartype(iv)==nf_double) then
                            allocate(var_r8(nlon_crs,nlat_crs,nlev(iv)))
                            var_r8(:,:,:)=var_latlon_crs(:,:,:)
                            status = nf_put_vara_double(ncid_out, varid(iv), start, count, var_r8)
                            deallocate(var_r8)
                         elseif (vartype(iv)==nf_float) then
                            allocate(var_r4(nlon_crs,nlat_crs,nlev(iv)))
                            var_r4(:,:,:)=var_latlon_crs(:,:,:)
                            status = nf_put_vara_real(ncid_out, varid(iv), start, count, var_r4)
                            deallocate(var_r4)
                         elseif (vartype(iv)==nf_short) then
                            allocate(var_i2(nlon_crs,nlat_crs,nlev(iv)))
                            var_i2(:,:,:)=(var_latlon_crs(:,:,:)-varoffset(iv))/varscale(iv)
                            status = nf_put_vara_int2(ncid_out, varid(iv), start, count, var_i2)
                            deallocate(var_i2)
                         endif
                      endif
                      deallocate(var_latlon_crs)
                   endif
                   deallocate(var_latlon, var_cubsph)
                else
                   n(1:4)=1
                   do id=1,varndims(iv)-1
                      status = nf_inq_dimlen(ncid_in(1), vardimids(id,iv), n(id))
                   enddo
                   start(1:4)=(/1,1,1,1/)
                   count(1:4)=(/n(1),n(2),n(3),1/)
                   start(varndims(iv))=it
                   if (vartype(iv)==nf_double) then
                      allocate(var_r8(n(1),n(2),n(3)))
                      status = nf_get_vara_double(ncid_in(1), iv, start, count, var_r8)
                      status = nf_put_vara_double(ncid_out, varid(iv), start, count, var_r8)
                      deallocate(var_r8)
                   elseif (vartype(iv)==nf_float) then
                      allocate(var_r4(n(1),n(2),n(3)))
                      status = nf_get_vara_real(ncid_in(1), iv, start, count, var_r4)
                      status = nf_put_vara_real(ncid_out, varid(iv), start, count, var_r4)
                      deallocate(var_r4)
                   elseif (vartype(iv)==nf_short) then
                      allocate(var_i2(n(1),n(2),n(3)))
                      status = nf_get_vara_int2(ncid_in(1), iv, start, count, var_i2)
                      status = nf_put_vara_int2(ncid_out, varid(iv), start, count, var_i2)
                      deallocate(var_i2)
                   else
                      print*," vartype neither nf_double nor nf_float of nf_short: ",vartype(iv)
                      stop
                   endif
                endif
             endif
          endif
       enddo
    enddo
    do itile=1,ntiles
       status = nf_close(ncid_in(itile))
    enddo
    status = nf_close(ncid_out)

    deallocate(ncid_in)
    deallocate(time)

    deallocate(varname, vartype, varndims, varnatts, varid, vardimids, nlev, ignore,  &
               cs_static, do_c2l, misval, varmisval, varscale, varoffset)

    print 102, trim(data_file)//".tile?.nc",trim(filename)
102 format(/,"interpolation done for data in ",a,/,"interpolated data stored in ",a,/)

    contains
      subroutine init_corners(var, nz)
        integer, intent(in) :: nz
        real(REAL8), dimension(0:npx,0:npy,nz,ntiles), intent(inout) :: var

        var(  0,  0,:,:)=0.
        var(  0,npy,:,:)=0.
        var(npx,npy,:,:)=0.
        var(npx,  0,:,:)=0.
      end subroutine init_corners
#endif

  end subroutine interpolate_data
  !====================================================================!
  subroutine do_c2l_interpolation_r4(cubsph, n1x, nx, n1y, ny, nz, tile, nlev,   &
                                  latlon, nlon, nlat, c2l_index, c2l_weight,  &
                                  misval, varmisval, fill_missing)
    !------------------------------------------------------------------!
    ! do bilinear interpolation from cubed sphere to latlon grid       !
    ! using precalculated weights from get_c2l_weight                  !
    !                                                                  !
    ! input:                                                           !
    !                                                                  !
    ! output:                                                          !
    !                                                                  !
    !------------------------------------------------------------------!
    integer, intent(in) :: n1x, nx, n1y, ny, nz, tile, nlev, nlon, nlat

    real(REAL4), dimension(n1x:nx, n1y:ny, nz), intent(in)    :: cubsph
    real(REAL4), dimension(nlon, nlat, nz),     intent(inout) :: latlon

    real(REAL8),    dimension(4, nlon, nlat), intent(in) :: c2l_weight
    integer, dimension(3, nlon, nlat), intent(in) :: c2l_index

    real(REAL8), intent(in) :: varmisval
    logical, intent(in) :: misval, fill_missing
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    real(REAL8) :: max_weight
    integer :: i, j, k, ic, jc

    if (misval) then
       do k=1,nlev
          do i=1,nlon
             do j=1,nlat
                if (tile==c2l_index(3,i,j)) then
                   ic=c2l_index(1,i,j)
                   jc=c2l_index(2,i,j)

                   if (cubsph(ic  ,jc  ,k) == varmisval .or.               &
                       cubsph(ic  ,jc+1,k) == varmisval .or.               &
                       cubsph(ic+1,jc+1,k) == varmisval .or.               &
                       cubsph(ic+1,jc  ,k) == varmisval) then
                      if (fill_missing) then
                         max_weight=max(c2l_weight(1,i,j),c2l_weight(2,i,j), &
                                        c2l_weight(3,i,j),c2l_weight(4,i,j))
                         if (max_weight==c2l_weight(1,i,j)) latlon(i,j,k)=cubsph(ic  ,jc  ,k)
                         if (max_weight==c2l_weight(2,i,j)) latlon(i,j,k)=cubsph(ic  ,jc+1,k)
                         if (max_weight==c2l_weight(3,i,j)) latlon(i,j,k)=cubsph(ic+1,jc+1,k)
                         if (max_weight==c2l_weight(4,i,j)) latlon(i,j,k)=cubsph(ic+1,jc  ,k)
                      else
                         latlon(i,j,k)=varmisval
                      endif
                   else
                      latlon(i,j,k)=c2l_weight(1,i,j)*cubsph(ic  ,jc  ,k)  &
                                   +c2l_weight(2,i,j)*cubsph(ic  ,jc+1,k)  &
                                   +c2l_weight(3,i,j)*cubsph(ic+1,jc+1,k)  &
                                   +c2l_weight(4,i,j)*cubsph(ic+1,jc  ,k)
                   endif
                endif
             enddo
          enddo
       enddo
    else
       do k=1,nlev
          do i=1,nlon
             do j=1,nlat
                if (tile==c2l_index(3,i,j)) then
                   ic=c2l_index(1,i,j)
                   jc=c2l_index(2,i,j)

                   latlon(i,j,k)=c2l_weight(1,i,j)*cubsph(ic  ,jc  ,k)  &
                                +c2l_weight(2,i,j)*cubsph(ic  ,jc+1,k)  &
                                +c2l_weight(3,i,j)*cubsph(ic+1,jc+1,k)  &
                                +c2l_weight(4,i,j)*cubsph(ic+1,jc  ,k)
                endif
             enddo
          enddo
       enddo
    endif

  end subroutine do_c2l_interpolation_r4

  subroutine do_c2l_interpolation(cubsph, n1x, nx, n1y, ny, nz, tile, nlev,   &
                                  latlon, nlon, nlat, c2l_index, c2l_weight,  &
                                  misval, varmisval, fill_missing)
    !------------------------------------------------------------------!
    ! do bilinear interpolation from cubed sphere to latlon grid       !
    ! using precalculated weights from get_c2l_weight                  !
    !                                                                  !
    ! input:                                                           !
    !                                                                  !
    ! output:                                                          !
    !                                                                  !
    !------------------------------------------------------------------!
    integer, intent(in) :: n1x, nx, n1y, ny, nz, tile, nlev, nlon, nlat

    real(REAL8), dimension(n1x:nx, n1y:ny, nz), intent(in)    :: cubsph
    real(REAL8), dimension(nlon, nlat, nz),     intent(inout) :: latlon

    real(REAL8),    dimension(4, nlon, nlat), intent(in) :: c2l_weight
    integer, dimension(3, nlon, nlat), intent(in) :: c2l_index

    real(REAL8), intent(in) :: varmisval
    logical, intent(in) :: misval, fill_missing
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    real(REAL8) :: max_weight
    integer :: i, j, k, ic, jc

    if (misval) then
       do k=1,nlev
          do i=1,nlon
             do j=1,nlat
                if (tile==c2l_index(3,i,j)) then
                   ic=c2l_index(1,i,j)
                   jc=c2l_index(2,i,j)

                   if (cubsph(ic  ,jc  ,k) == varmisval .or.               &
                       cubsph(ic  ,jc+1,k) == varmisval .or.               &
                       cubsph(ic+1,jc+1,k) == varmisval .or.               &
                       cubsph(ic+1,jc  ,k) == varmisval) then
                      if (fill_missing) then
                         max_weight=max(c2l_weight(1,i,j),c2l_weight(2,i,j), &
                                        c2l_weight(3,i,j),c2l_weight(4,i,j))
                         if (max_weight==c2l_weight(1,i,j)) latlon(i,j,k)=cubsph(ic  ,jc  ,k)
                         if (max_weight==c2l_weight(2,i,j)) latlon(i,j,k)=cubsph(ic  ,jc+1,k)
                         if (max_weight==c2l_weight(3,i,j)) latlon(i,j,k)=cubsph(ic+1,jc+1,k)
                         if (max_weight==c2l_weight(4,i,j)) latlon(i,j,k)=cubsph(ic+1,jc  ,k)
                      else
                         latlon(i,j,k)=varmisval
                      endif
                   else
                      latlon(i,j,k)=c2l_weight(1,i,j)*cubsph(ic  ,jc  ,k)  &
                                   +c2l_weight(2,i,j)*cubsph(ic  ,jc+1,k)  &
                                   +c2l_weight(3,i,j)*cubsph(ic+1,jc+1,k)  &
                                   +c2l_weight(4,i,j)*cubsph(ic+1,jc  ,k)
                   endif
                endif
             enddo
          enddo
       enddo
    else
       do k=1,nlev
          do i=1,nlon
             do j=1,nlat
                if (tile==c2l_index(3,i,j)) then
                   ic=c2l_index(1,i,j)
                   jc=c2l_index(2,i,j)

                   latlon(i,j,k)=c2l_weight(1,i,j)*cubsph(ic  ,jc  ,k)  &
                                +c2l_weight(2,i,j)*cubsph(ic  ,jc+1,k)  &
                                +c2l_weight(3,i,j)*cubsph(ic+1,jc+1,k)  &
                                +c2l_weight(4,i,j)*cubsph(ic+1,jc  ,k)
                endif
             enddo
          enddo
       enddo
    endif

  end subroutine do_c2l_interpolation
  !====================================================================!
  subroutine do_latlon_coarsening(var_latlon, ylat, nlon, nlat, nz,     &
                                  var_latlon_crs, nlon_crs, nlat_crs,   &
                                  finer_steps, misval, varmisval)
    !------------------------------------------------------------------!
    ! calculate variable on coarser latlon grid                        !
    ! by doubling spatial resolution and preserving volume means       !
    !                                                                  !
    ! input:                                                           !
    !                                                                  !
    ! output:                                                          !
    !                                                                  !
    !------------------------------------------------------------------!
    integer, intent(in) :: nlon, nlat, nz, nlon_crs, nlat_crs, finer_steps
    real(REAL8),    intent(in) :: var_latlon(nlon,nlat,nz), ylat(nlat)
    real(REAL8), dimension(nlon_crs,nlat_crs,nz), intent(out) :: var_latlon_crs
    real(REAL8), intent(in) :: varmisval
    logical, intent(in) :: misval
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    real(REAL8), allocatable :: var_latlon_old(:,:,:), ylat_old(:),            &
                         var_latlon_new(:,:,:)
    real(REAL8) :: dlat
    integer :: nlon_old, nlat_old, nlon_new, nlat_new, steps, j

    select case (finer_steps)

    case (0)
       if (nlon_crs/=nlon .or. nlat_crs/=nlat) then
          print*," do_latlon_coarsening: grid dimensions don't match"
          stop
       endif
       var_latlon_crs(1:nlon_crs,1:nlat_crs,1:nz) = var_latlon(1:nlon_crs,1:nlat_crs,1:nz)
       
    case (1)
       call redu2x(var_latlon, ylat, nlon, nlat, var_latlon_crs, nlon_crs, nlat_crs)

    case default
       nlon_new=nlon
       nlat_new=nlat
       do steps=1,finer_steps
          nlon_old=nlon_new
          nlat_old=nlat_new
          allocate(var_latlon_old(nlon_old,nlat_old,nz),ylat_old(nlat_old))
          if (steps==1) then
             ylat_old(:)=ylat(:)
             var_latlon_old(:,:,:)=var_latlon(:,:,:)
          else
             dlat=pi/real(nlat_old-1)
             ylat_old(1)=-0.5*pi
             ylat_old(nlat_old)= 0.5*pi
             do j=2,nlat_old-1
                ylat_old(j)=ylat_old(1)+(real(j)-1.)*dlat
             enddo
             var_latlon_old(:,:,:)=var_latlon_new(:,:,:)
             deallocate(var_latlon_new)
          endif
          
          nlon_new=nlon_new/2
          nlat_new=(nlat_new-1)/2+1
          allocate(var_latlon_new(nlon_new,nlat_new,nz))
          
          call redu2x(var_latlon_old, ylat_old, nlon_old, nlat_old, var_latlon_new, nlon_new, nlat_new)
          deallocate(var_latlon_old,ylat_old)
       enddo
       var_latlon_crs(:,:,:)=var_latlon_new(:,:,:)
       deallocate(var_latlon_new)
       
    end select

  contains
    !------------------------------------------------------------------!
    subroutine redu2x(varfin, yfin, nxfin, nyfin, varcrs, nxcrs, nycrs)
      !----------------------------------------------------------------!
      ! this routine is for reducing fvccm data by a factor of 2       !
      ! volume averaging for all data except at the poles              !
      ! original developer: S.-J. Lin                                  !
      !----------------------------------------------------------------!
      integer, intent(in) :: nxfin, nyfin, nxcrs, nycrs
      real(REAL8), intent(in)  :: varfin(nxfin,nyfin,nz), yfin(nyfin)
      real(REAL8), intent(out) :: varcrs(nxcrs,nycrs,nz)
      !----------------------------------------------------------------!
      ! local variables                                                !
      !----------------------------------------------------------------!
      real(REAL8) :: cosp(nyfin), acosp(nyfin), vartmp(nxcrs,nyfin,nz)
      integer :: i, j, k, i2, j2
      !----------------------------------------------------------------!
      ! calculate cosine of latitude                                   !
      ! trick in cosp needed to maintain a constant field              !
      !----------------------------------------------------------------!
      cosp(1)=0.
      cosp(nyfin)=0.
      do j=2,nyfin-1
         cosp(j)=cos(yfin(j))
      enddo
      do j=2,nyfin-1
         acosp(j) = 1. / (cosp(j)+0.5*(cosp(j-1)+cosp(j+1)))
      enddo
      !----------------------------------------------------------------!
      ! x-sweep                                                        !
      !----------------------------------------------------------------!
      if (misval) then
         do k=1,nz
            do j=2,nyfin-1
               if (varfin(nxfin,j,k) == varmisval .or.                  &
                   varfin(1,j,k)     == varmisval .or.                  &
                   varfin(2,j,k)     == varmisval) then
                  vartmp(1,j,k) = varmisval
               else
                  vartmp(1,j,k) = 0.25*(varfin(nxfin,j,k)+2.*varfin(1,j,k)+varfin(2,j,k))
               endif
               do i=3,nxfin-1,2
                  i2 = (i+1)/2
                  if (varfin(i-1,j,k) == varmisval .or.                  &
                      varfin(i  ,j,k) == varmisval .or.                  &
                      varfin(i+1,j,k) == varmisval) then
                     vartmp(i2,j,k) = varmisval
                  else
                     vartmp(i2,j,k) = 0.25*(varfin(i-1,j,k)+2.*varfin(i,j,k)+varfin(i+1,j,k))
                  endif
               enddo
            enddo
         enddo
      else
         do k=1,nz
            do j=2,nyfin-1
               vartmp(1,j,k) = 0.25*(varfin(nxfin,j,k)+2.*varfin(1,j,k)+varfin(2,j,k))
               do i=3,nxfin-1,2
                  i2 = (i+1)/2
                  vartmp(i2,j,k) = 0.25*(varfin(i-1,j,k)+2.*varfin(i,j,k)+varfin(i+1,j,k))
               enddo
            enddo
         enddo
      endif
      !----------------------------------------------------------------!
      ! poles:                                                         !
      ! this code segment works for both the scalar and vector fields. !
      ! Winds at poles are wave-1; the follwoing is quick & dirty yet the correct way
      ! The skipping method. A more rigorous way is to                 !
      ! recompute the wave-1 components for the coarser grid.          !
      !----------------------------------------------------------------!
      do k=1,nz
         do i=1,nxcrs
            i2 = i*2-1
            varcrs(i,    1,k) = varfin(i2,    1,k)
            varcrs(i,nycrs,k) = varfin(i2,nyfin,k)
         enddo
      enddo
      !----------------------------------------------------------------!
      ! y-sweep                                                        !
      !----------------------------------------------------------------!
      if (misval) then
         do k=1,nz
            do j=2,nyfin-1
               do i=1,nxcrs
                  if (vartmp(i,j,k)/= varmisval) vartmp(i,j,k) = vartmp(i,j,k)*cosp(j)
               enddo
            enddo
         enddo
         do k=1,nz
            do j=3,nyfin-2,2
               j2 = (j+1)/2
               do i=1,nxcrs
                  if (vartmp(i,j  ,k) /= varmisval .and.                &
                      vartmp(i,j-1,k) /= varmisval .and.                &
                      vartmp(i,j+1,k) /= varmisval) then
                     varcrs(i,j2,k) = acosp(j)*(vartmp(i,j,k)+0.5*(vartmp(i,j-1,k)+vartmp(i,j+1,k)))
                  else
                     varcrs(i,j2,k) = varmisval
                  endif
               enddo
            enddo
         enddo
      else
         do k=1,nz
            do j=2,nyfin-1
               do i=1,nxcrs
                  vartmp(i,j,k) = vartmp(i,j,k)*cosp(j)
               enddo
            enddo
         enddo
         do k=1,nz
            do j=3,nyfin-2,2
               j2 = (j+1)/2
               do i=1,nxcrs
                  varcrs(i,j2,k) = acosp(j)*(vartmp(i,j,k)+0.5*(vartmp(i,j-1,k)+vartmp(i,j+1,k)))
               enddo
            enddo
         enddo
      endif

    end subroutine redu2x
    !------------------------------------------------------------------!
  end subroutine do_latlon_coarsening
  !====================================================================!
end module CUB2LATLON_mod
