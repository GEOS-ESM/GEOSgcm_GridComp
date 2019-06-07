!-*- F90 -*-
program cub2cub_driver
  !--------------------------------------------------------------------!
  ! author: Michael Herzog                                             !
  ! email:  Michael.Herzog@noaa.gov                                    !
  !                                                                    !
  ! purpose: driver for interpolation between two cubed sphere grids   !
  !          with different spatial resolution                         !
  !--------------------------------------------------------------------!
  use fv_arrays_mod,  only: REAL4, REAL8, FVPRC
  use CUB2CUB_mod,    only: read_c2c_namelist,              &
                            get_c2c_weight,                 &
                            interpolate_c2c

  use CUB2LATLON_mod, only: read_grid_dimensions,           &
                            init_cubsph_grid

  implicit none
  
  real(FVPRC), dimension(:,:,:,:), allocatable :: corner_in, corner_out, weight_c2c
  integer, dimension(:,:,:,:), allocatable :: index_c2c

  integer :: ntiles, npx_in, npy_in, npx_out, npy_out
  character(len=120) :: grid_in, grid_out, dir_in, dir_out
  !--------------------------------------------------------------------!
  ! main setup, read in namelists                                      !
  !--------------------------------------------------------------------!
  call read_c2c_namelist(ntiles, grid_in, grid_out, dir_in, dir_out)
  !--------------------------------------------------------------------!
  ! initialize cubed sphere grid: in                                   !
  !--------------------------------------------------------------------!
  call read_grid_dimensions(grid_in, npx_in, npy_in)
  allocate(corner_in(2,0:npx_in+1,0:npy_in+1,ntiles))
  call init_cubsph_grid(npx_in, npy_in, ntiles, grid_in, corner_in)
  !--------------------------------------------------------------------!
  ! initialize cubed sphere grid: out                                  !
  !--------------------------------------------------------------------!
  call read_grid_dimensions(grid_out, npx_out, npy_out)
  allocate(corner_out(2,0:npx_out+1,0:npy_out+1,ntiles))
  call init_cubsph_grid(npx_out, npy_out, ntiles, grid_out, corner_out)
  !--------------------------------------------------------------------!
  ! calculate weights and indices from bilinear interpolation          !
  ! from grid_in to grid_out                                           !
  !--------------------------------------------------------------------!
  allocate(index_c2c (3, npx_out-1, npy_out-1, ntiles))
  allocate(weight_c2c(4, npx_out-1, npy_out-1, ntiles))
  call get_c2c_weight(ntiles, npx_in, npy_in, corner_in,                &
                      npx_out, npy_out, corner_out,                     &
                      index_c2c,  weight_c2c)
  !--------------------------------------------------------------------!
  ! do cub2cub interpolation                                           !
  !--------------------------------------------------------------------!
  call interpolate_c2c(ntiles, npx_in, npy_in, npx_out, npy_out,        &
                       dir_in, dir_out, index_c2c, weight_c2c)
  !--------------------------------------------------------------------!
  ! deallocate arrays                                                  !
  !--------------------------------------------------------------------!
  deallocate(corner_in, corner_out, index_c2c, weight_c2c)

end program cub2cub_driver

