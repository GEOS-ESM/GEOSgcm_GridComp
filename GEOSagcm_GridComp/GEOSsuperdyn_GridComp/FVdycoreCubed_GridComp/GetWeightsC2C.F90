#define R8 8

subroutine GetWeightsC2C( npx_in, npy_in, npx_out, npy_out, index, weight, &
     ee1, ee2, ff1, ff2)

 use CUB2CUB_mod,     only : get_c2c_weight, get_c2c_weight_global
 use mpp_mod,         only : FATAL, mpp_error

 use ESMF
 implicit none
  integer, intent(in) :: npx_in, npy_in, npx_out, npy_out
  integer, intent(out) :: index(:,:,:,:)
  real(R8), intent(out) :: weight(:,:,:,:)
  real(R8), dimension(:,:,:), intent(out) :: ee1, ee2, ff1, ff2
! Note that the shape of ee's and ff's is different


! local vars
  integer, parameter :: ntiles=6

  integer :: npy_i, npy_o
  real(R8), pointer :: sph_in(:,:,:,:) => null()
  real(R8), pointer :: sph_out(:,:,:,:) => null()

  ! initialize cube_in and cube_out corners
  !----------------------------------------

  npy_i = npy_in/ntiles  ! should be the same as npx_in
  npy_o = npy_out/ntiles ! should be the same as npx_out
  call init_cube_corners(npx_in, npy_i, sph_in, ff1=ff1, ff2=ff2)
  call init_cube_corners(npx_out, npy_o, sph_out, ee1=ee1, ee2=ee2)

  ! calculate weights for bilinear interpolation
  ! from cubed sphere to cubed sphere
  !---------------------------------------------

  call get_c2c_weight_global(ntiles, npx_in+1, npy_i+1, sph_in,                &
                      npx_out+1, npy_o+1, sph_out,                     &
                      index,  weight)

  deallocate ( sph_in, sph_out )

  return

contains
  subroutine init_cube_corners(npx, npy, sph_corner, ee1, ee2, ff1, ff2)

    use fv_grid_utils_mod, only : gnomonic_grids, cell_center2
    use fv_grid_tools_mod, only : mirror_grid
    use GHOST_CUBSPH_mod,  only : B_grid, A_grid, ghost_cubsph_update


    integer,  intent(in   ) :: npx,  npy
    real(R8), pointer       :: sph_corner (:,:,:,:)
    real(R8), dimension(:,:,:), optional :: ee1, ee2, ff1, ff2

    ! Locals
    !-------

    integer :: npts, n, l, i, j, j1

    integer, parameter :: ntiles=6
    integer, parameter :: ndims=2


    real(R8), parameter :: PI=3.14159265358979323846
    real(R8), dimension(ndims) :: agrid
    real(R8), dimension(3) :: e1, e2, f1, f2

    ! Real*8 are needed to make fv calls.
    !-----------------------------------

    real(R8), allocatable :: grid_global(:,:,:,:)

    npts = npx + 1
    
    if (associated(sph_corner)) deallocate(sph_corner)
    allocate( sph_corner(ndims,0:npts+1,0:npts+1,ntiles))
    
    allocate( grid_global(npts,npts,ndims,ntiles))
    
    call gnomonic_grids(A_grid, npx, grid_global(:,:,1,1), grid_global(:,:,2,1))
    
    ! mirror_grid assumes that the tile=1 is centered 
    !   on equator and greenwich meridian Lon[-pi,pi]
    !------------------------------------------------
    
    call mirror_grid(grid_global, 0, npts, npts, ndims, ntiles)
    
    ! Shift the corner away from Japan.
    !  This will result in the corner 
    !  close to the east coast of China.
    !-----------------------------------
    
    grid_global(:,:,1,:) = grid_global(:,:,1,:) - PI/18.
    
    where(grid_global(:,:,1,:) < 0.) &
         grid_global(:,:,1,:) = grid_global(:,:,1,:) + 2.* PI
    
    ! Keep Equator and Greenwich exact
    !---------------------------------
    
    where(abs(grid_global(:,:,:,1)) < 1.e-10) grid_global(:,:,:,1) = 0.0
    
    
    ! Clean Up Corners
    !---------------------------------
    
    grid_global(1   ,   :,:,2)=grid_global(npts     ,:        ,:,1)
    grid_global(1   ,   :,:,3)=grid_global(npts:1:-1,npts     ,:,1)
    grid_global(:   ,npts,:,5)=grid_global(1        ,npts:1:-1,:,1)
    grid_global(:   ,npts,:,6)=grid_global(:        ,1        ,:,1)
    grid_global(:   ,   1,:,3)=grid_global(:        ,npts     ,:,2)
    grid_global(:   ,   1,:,4)=grid_global(npts     ,npts:1:-1,:,2)
    grid_global(npts,   :,:,6)=grid_global(npts:1:-1,1        ,:,2)
    grid_global(1   ,   :,:,4)=grid_global(npts     ,:        ,:,3)
    grid_global(1   ,   :,:,5)=grid_global(npts:1:-1,npts     ,:,3)
    grid_global(npts,   :,:,3)=grid_global(1        ,:        ,:,4)
    grid_global(:   ,   1,:,5)=grid_global(:        ,npts     ,:,4)
    grid_global(:   ,   1,:,6)=grid_global(npts     ,npts:1:-1,:,4)
    grid_global(1   ,   :,:,6)=grid_global(npts     ,:        ,:,5)


    sph_corner(1,1:npts,1:npts,:) = grid_global(:,:,1,:)
    sph_corner(2,1:npts,1:npts,:) = grid_global(:,:,2,:)

    ! Compute Vector Rotation arrays

    do n=1,ntiles
       do j=1,npy
          j1 = npx*(n-1) + j
          do i=1,npx
             agrid = 1.e+25
             call cell_center2(grid_global(i,  j,  1:2,n), &
                               grid_global(i+1,j,  1:2,n),   &
                               grid_global(i,  j+1,1:2,n), &
                               grid_global(i+1,j+1,1:2,n),   &
                               agrid(1:2) )
             call CreateCube2LatLonRotation( &
                  grid_global(i:i+1,j:j+1,:,n), agrid(1:2), &
                  e1(:),e2(:),f1(:),f2(:))
             if (present(ff1)) then
                ff1(i,j1,:) = f1(:)
             end if
             if (present(ff2)) then
                ff2(i,j1,:) = f2(:)
             end if
             if (present(ee1)) then
                ee1(i,j1,:) = e1(:)
             end if
             if (present(ee2)) then
                ee2(i,j1,:) = e2(:)
             end if
          enddo
       enddo
    enddo

    deallocate ( grid_global )

  ! do halo update                                                   !
  !------------------------------------------------------------------!

    do n=1,ntiles
       sph_corner(1:2,0     ,0     ,n)=0.
       sph_corner(1:2,npts+1,0     ,n)=0.
       sph_corner(1:2,0     ,npts+1,n)=0.
       sph_corner(1:2,npts+1,npts+1,n)=0.
       do L=1,2
          call ghost_cubsph_update(sph_corner(L,0:npts+1,0:npts+1,:), &
               0, npts+1, 0, npts+1,           &
               1, ntiles, 1, n, B_grid         )
       end do
    enddo

    return
  end subroutine init_cube_corners

  subroutine CreateCube2LatLonRotation(grid, center, ee1, ee2, ff1, ff2)

    use fv_grid_utils_mod, only : mid_pt_sphere
    use fv_grid_tools_mod, only : get_unit_vector

    real(R8), intent(IN)    :: grid(0:1,0:1,2), center(2)
    real(R8), intent(OUT), dimension(3) :: ee1, ee2, ff1, ff2

    real(R8), dimension(2) :: p1, p2, p3, p4
    real :: H, F


    call mid_pt_sphere(grid(0,0,:),grid(0,1,:), p1)
    call mid_pt_sphere(grid(0,0,:),grid(1,0,:), p2)
    call mid_pt_sphere(grid(1,0,:),grid(1,1,:), p3)
    call mid_pt_sphere(grid(0,1,:),grid(1,1,:), p4)

    call get_unit_vector(p3, center, p1, ee1)
    call get_unit_vector(p4, center, p2, ee2)

    H   = dot_product(ee1,ee2)
    F   = 1.0/(H**2-1.0)

    ff1 = F*(ee2*H-ee1)
    ff2 = F*(ee1*H-ee2)

    return
  end subroutine CreateCube2LatLonRotation

end subroutine GetWeightsC2C
