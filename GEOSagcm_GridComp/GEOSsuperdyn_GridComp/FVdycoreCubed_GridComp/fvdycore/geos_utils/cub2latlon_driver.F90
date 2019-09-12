!-*- F90 -*-
program cub2latlon_driver
  !--------------------------------------------------------------------!
  ! author: Michael Herzog                                             !
  ! email:  Michael.Herzog@noaa.gov                                    !
  !                                                                    !
  ! purpose: driver for interpolation from cubed sphere to latlon      !
  !--------------------------------------------------------------------!
  use CUB2LATLON_mod, only: read_c2l_namelist,                          &
                            init_cubsph_grid, init_latlon_grid,         &
                            read_c2l_weight,  write_c2l_weight,         &
                            get_c2l_weight,   interpolate_data,         &
                            read_grid_dimensions 

  implicit none

  integer :: npx, npy, npz, ntiles, nlon, nlat, finer_steps
  logical :: read_res, write_res, memphis, fill_missing
  character(len=120) :: grid_file, data_file, data_out,                 &
                        uname, vname, missing_value

  real(REAL8), dimension(:,:,:,:), allocatable :: sph_corner
  real(REAL8), dimension(:),       allocatable :: xlon, ylat

  real(REAL8),    dimension(:,:,:), allocatable :: c2l_weight
  integer, dimension(:,:,:), allocatable :: c2l_index
  real(REAL8),    dimension(:,:,:,:), allocatable  :: elon_cubsph, elat_cubsph
  real(REAL8),    dimension(:,:,:),   allocatable  :: elon_latlon, elat_latlon

  logical :: found=.false.
  !--------------------------------------------------------------------!
  ! main setup, read in namelists                                      !
  !--------------------------------------------------------------------!
  call read_c2l_namelist(ntiles, nlon, nlat, finer_steps,               &
                         read_res, write_res, memphis,                  &
                         grid_file, data_file, data_out,                &
                         uname, vname, missing_value, fill_missing)
  !--------------------------------------------------------------------!
  ! initialize cubed sphere grid                                       !
  !--------------------------------------------------------------------!
  call read_grid_dimensions(grid_file, npx, npy)
  allocate(sph_corner(2,0:npx+1,0:npy+1,ntiles))
  call init_cubsph_grid(npx, npy, ntiles, grid_file, sph_corner)
  !--------------------------------------------------------------------!
  ! initialize latlon grid                                             !
  !--------------------------------------------------------------------!
  nlon=2**finer_steps*nlon
  nlat=2**finer_steps*(nlat-1)+1
  allocate(xlon(nlon), ylat(nlat))
  call init_latlon_grid(xlon, ylat, nlon, nlat)
  !--------------------------------------------------------------------!
  ! calculate weights for bilinear interpolation                       !
  ! from cubed sphere to latlon grid                                   !
  !--------------------------------------------------------------------!
  allocate(c2l_index(3,nlon,nlat),c2l_weight(4,nlon,nlat))
  allocate(elon_cubsph(3,0:npx,0:npy,ntiles), elon_latlon(3,nlon,nlat), &
           elat_cubsph(3,0:npx,0:npy,ntiles), elat_latlon(3,nlon,nlat))
  if (read_res) call read_c2l_weight(c2l_index, c2l_weight, nlon, nlat, npx, npy, ntiles, &
                        elon_cubsph, elat_cubsph, elon_latlon, elat_latlon, found)
  if (.not. found) then
     call get_c2l_weight(sph_corner, npx, npy, ntiles, xlon, ylat, nlon, nlat,  &
                         c2l_index, c2l_weight, elon_cubsph, elat_cubsph,       &
                         elon_latlon, elat_latlon)
     if (write_res) call write_c2l_weight(c2l_index, c2l_weight, nlon, nlat, npx, npy, ntiles, &
                                          elon_cubsph, elat_cubsph, elon_latlon, elat_latlon)
  endif
  !--------------------------------------------------------------------!
  ! do cub2latlon interpolation                                        !
  !--------------------------------------------------------------------!
  call interpolate_data(npx, npy, ntiles, data_file, data_out,              &
                        uname, vname, missing_value, fill_missing,          &
                        memphis, c2l_index, c2l_weight,                     &
                        elon_cubsph, elat_cubsph, elon_latlon, elat_latlon, &
                        xlon, ylat, nlon, nlat, finer_steps)
  !--------------------------------------------------------------------!
  ! deallocate arrays                                                  !
  !--------------------------------------------------------------------!
  deallocate(sph_corner, xlon, ylat, c2l_index, c2l_weight,             &
             elon_cubsph, elat_cubsph, elon_latlon, elat_latlon)
end program cub2latlon_driver

