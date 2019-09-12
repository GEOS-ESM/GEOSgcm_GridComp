  program moist_interp_rst

  !--------------------------------------------------------------------!
  ! purpose: driver for interpolation between two cubed sphere grids   !
  !          with different spatial resolution for GEOS restarts       !
  !--------------------------------------------------------------------!
  use ESMF
  use CUB2CUB_mod,       only : get_c2c_weight, do_c2c_interpolation
  use MAPL_ConstantsMod, only : pi=> MAPL_PI
  use fv_grid_utils_mod, only : gnomonic_grids
  use fv_grid_tools_mod, only : mirror_grid
  use GHOST_CUBSPH_mod,  only : B_grid, A_grid, ghost_cubsph_update

  implicit none

  integer :: npx_in
  integer :: npy_in
  integer :: npx_out
  integer :: npy_out

  integer :: npz

  integer :: npts
  integer :: ntiles=6
  integer :: ndims=2
  real(ESMF_KIND_R8), allocatable :: xs(:,:), ys(:,:)
  real(ESMF_KIND_R8), allocatable :: grid_in(:,:,:,:)
  real(ESMF_KIND_R8), allocatable :: grid_out(:,:,:,:)
  real(ESMF_KIND_R8), allocatable :: corner_in(:,:,:,:)
  real(ESMF_KIND_R8), allocatable :: corner_out(:,:,:,:)
  real(ESMF_KIND_R8), allocatable :: weight_c2c(:,:,:,:)
  integer,            allocatable :: index_c2c(:,:,:,:)

  real(ESMF_KIND_R8), allocatable :: ak(:), bk(:)
  integer :: IUNIT=15
  integer :: OUNIT=16

  integer :: header(6)
  integer :: grid_type = 0
  integer :: i,j,l,k,j1,j2,status

  character(30) :: str_arg

  external :: getarg, iargc
  integer iargc

  if (IARGC() /= 3) then
     print*, 'ABORT: need 3 arguments input_res and output_res and levels (No vertical interp supported yet)'
     stop
  endif
 
  CALL GETARG(1, str_arg)
  read (str_arg,'(I10)') npx_in
  read (str_arg,'(I10)') npy_in
  CALL GETARG(2, str_arg)
  read (str_arg,'(I10)') npx_out
  read (str_arg,'(I10)') npy_out
  CALL GETARG(3, str_arg)
  read (str_arg,'(I10)') npz
  print*, npx_in, npx_out, npz

  !--------------------------------------------------------------------!
  ! initialize Input cubed sphere grid                                 !
  !--------------------------------------------------------------------!
  npts = npx_in+1
  allocate( xs(npts,npts) )
  allocate( ys(npts,npts) )
  allocate( grid_in(npts,npts,ndims,ntiles) )
  call gnomonic_grids(grid_type, npts-1, xs, ys)
  do j=1,npts
     do i=1,npts
        grid_in(i,j,1,1) = xs(i,j)
        grid_in(i,j,2,1) = ys(i,j)
     enddo
  enddo
  deallocate ( xs )
  deallocate ( ys )
  ! mirror_grid assumes that the tile=1 is centered on equator and greenwich meridian Lon[-pi,pi]
  call mirror_grid(grid_in, npts, npts, 2, 6)
  allocate( corner_in(ndims,0:npts+1,0:npts+1,ntiles) )
  corner_in(1,1:npts,1:npts,:) = grid_in(:,:,1,:)
  corner_in(2,1:npts,1:npts,:) = grid_in(:,:,2,:)
  !------------------------------------------------------------------!
  ! do halo update                                                   !
  !------------------------------------------------------------------!
  do l=1,ntiles
     corner_in(1:2,0     ,0     ,l)=0.
     corner_in(1:2,npts+1,0     ,l)=0.
     corner_in(1:2,0     ,npts+1,l)=0.
     corner_in(1:2,npts+1,npts+1,l)=0.
     call ghost_cubsph_update(corner_in(1,0:npts+1,0:npts+1,:), 0, npts+1, 0, npts+1, 1, &
                              1, ntiles, 1, 1, l, B_grid)
     call ghost_cubsph_update(corner_in(2,0:npts+1,0:npts+1,:), 0, npts+1, 0, npts+1, 1, &
                              1, ntiles, 1, 1, l, B_grid)
  enddo
  deallocate (grid_in)
 
  !--------------------------------------------------------------------!
  ! initialize Output cubed sphere grid                                !
  !--------------------------------------------------------------------!
  npts = npx_out+1
  allocate( xs(npts,npts) )
  allocate( ys(npts,npts) )
  allocate( grid_out(npts,npts,ndims,ntiles) )
  call gnomonic_grids(grid_type, npts-1, xs, ys)
  do j=1,npts
     do i=1,npts
        grid_out(i,j,1,1) = xs(i,j)
        grid_out(i,j,2,1) = ys(i,j)
     enddo
  enddo
  deallocate ( xs )
  deallocate ( ys )
  ! mirror_grid assumes that the tile=1 is centered on equator and greenwich meridian Lon[-pi,pi]
  call mirror_grid(grid_out, npts, npts, 2, 6)
  allocate( corner_out(ndims,0:npts+1,0:npts+1,ntiles) )
  corner_out(1,1:npts,1:npts,:) = grid_out(:,:,1,:)
  corner_out(2,1:npts,1:npts,:) = grid_out(:,:,2,:)
  !------------------------------------------------------------------!
  ! do halo update                                                   !
  !------------------------------------------------------------------!
  do l=1,ntiles
     corner_out(1:2,0     ,0     ,l)=0.
     corner_out(1:2,npts+1,0     ,l)=0.
     corner_out(1:2,0     ,npts+1,l)=0.
     corner_out(1:2,npts+1,npts+1,l)=0.
     call ghost_cubsph_update(corner_out(1,0:npts+1,0:npts+1,:), 0, npts+1, 0, npts+1, 1, &
                              1, ntiles, 1, 1, l, B_grid)
     call ghost_cubsph_update(corner_out(2,0:npts+1,0:npts+1,:), 0, npts+1, 0, npts+1, 1, &
                              1, ntiles, 1, 1, l, B_grid)
  enddo
  deallocate (grid_out)

  !--------------------------------------------------------------------!
  ! calculate weights and indices from bilinear interpolation          !
  ! from grid_in to grid_out                                           !
  !--------------------------------------------------------------------!
  allocate(index_c2c (3, npx_out, npy_out, ntiles))
  allocate(weight_c2c(4, npx_out, npy_out, ntiles))
  call get_c2c_weight(ntiles, npx_in+1, npy_in+1, corner_in,                &
                      npx_out+1, npy_out+1, corner_out,                     &
                      index_c2c,  weight_c2c)

  open(IUNIT,file='moist_internal_restart_in' ,access='sequential',form='unformatted',status='old')
  open(OUNIT,file='moist_internal_restart_out',access='sequential',form='unformatted')

! 7 3d Fields (4-bytes)
  do l=1,7
    call read_interp_write(IUNIT, OUNIT, npx_in, npy_in, npx_out, npy_out, npz, ntiles, &
                           index_c2c, weight_c2c)
  enddo
  close(IUNIT)
  close(OUNIT)

  deallocate(corner_in, corner_out, index_c2c, weight_c2c)

 contains
 
 subroutine read_interp_write(IUNIT, OUNIT, npx_in, npy_in, npx_out, npy_out, npz, ntiles, &
                              index_c2c, weight_c2c)
  integer, intent(IN) :: IUNIT, OUNIT, npx_in, npy_in, npx_out, npy_out, npz, ntiles
  integer           , intent(IN) ::  index_c2c(3, npx_out, npy_out, ntiles)
  real(ESMF_KIND_R8), intent(IN) :: weight_c2c(4, npx_out, npy_out, ntiles)

  real(ESMF_KIND_R4), allocatable :: vari(:,:)
  real(ESMF_KIND_R4), allocatable :: varo(:,:)
  real(ESMF_KIND_R8), allocatable :: var_in(:,:,:,:)
  real(ESMF_KIND_R8), allocatable :: var_out(:,:,:,:)

  integer :: i,j,l,k,j1,j2,status

  allocate( vari(npx_in,npy_in*ntiles) )
  allocate( varo(npx_out,npy_out*ntiles) )
  allocate( var_in(0:npx_in+1,0:npy_in+1,npz,ntiles) )
  allocate( var_out(npx_out,npy_out,npz,ntiles) )
! 
  do k=1,npz
    read (IUNIT, IOSTAT=status) vari
    do l=1,ntiles
       j1 = (npy_in)*(l-1) + 1
       j2 = (npy_in)*(l-1) + npx_in
       var_in(1:npx_in,1:npy_in,k,l)=vari(:,j1:j2)
    enddo

! if (k==10) then
!    print*, vari(:,1)
! endif

  enddo
  do l=1,ntiles
    call ghost_cubsph_update(var_in, 0, npx_in+1, 0, npx_in+1, npz, 1, ntiles,  &
                             1, npz, l, A_grid)
  enddo
  call do_c2c_interpolation(var_in, 0, npx_in+1, 0, npy_in+1, npz, ntiles, &
                            index_c2c, weight_c2c, npx_out, npy_out, var_out)
  do k=1,npz
    do l=1,ntiles
       j1 = (npy_out)*(l-1) + 1
       j2 = (npy_out)*(l-1) + npy_out
       varo(:,j1:j2)=var_out(:,:,k,l)
    enddo
    write(OUNIT) varo

! if (k==10) then
!    print*, '--------------'
!    print*, varo(:,1)
! endif

  enddo

  deallocate ( vari, varo, var_in, var_out )

 end subroutine read_interp_write

  end

