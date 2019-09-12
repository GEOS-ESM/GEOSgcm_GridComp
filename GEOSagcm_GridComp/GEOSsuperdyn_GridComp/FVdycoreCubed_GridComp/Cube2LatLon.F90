subroutine cube2latlon(npx, npy, nlon, nlat, data_cs, data_ll)

#define REAL8 8

 use ESMF
 use MAPL_Mod, only : MAPL_UNDEF
 use MAPL_IOMod, only : GETFILEUNIT, FREE_FILE
 use CUB2LATLON_mod,    only : init_latlon_grid, &
                           read_c2l_weight,  write_c2l_weight,         &
                           new_get_c2l_weight, get_c2l_weight, do_c2l_interpolation_r4
 use GHOST_CUBSPH_mod,   only: B_grid, A_grid, ghost_cubsph_update

 implicit none

! This EXTERNAL subroutine is in the fv directory
!  and has real*8 interfaces
  interface
     subroutine GetWeights(npx, npy, nlat, nlon, &
          index, weight, id1, id2, jdc, l2c,     &
          ee1, ee2, ff1, ff2, gg1, gg2,          &
           e1,  e2,  f1,  f2,  g1,  g2,          &
          sublons, sublats, AmNodeRoot, WriteNetcdf)
       integer    ,  intent(in   ) :: npx,  npy
       integer    ,  intent(in   ) :: nlon, nlat
       integer    ,  intent(  out) :: index(3,nlon,nlat)
       real(REAL8), intent(  out) :: weight(4,nlon,nlat)
       integer    ,  intent(  out) :: id1(npx,npy)
       integer    ,  intent(  out) :: id2(npx,npy)
       integer    ,  intent(  out) :: jdc(npx,npy)
       real(REAL8), intent(  out) :: l2c(4,npx,npy)
       real(REAL8), intent(  out) :: ee1(npx,npy,3)
       real(REAL8), intent(  out) :: ee2(npx,npy,3)
       real(REAL8), intent(  out) :: ff1(npx,npy,3)
       real(REAL8), intent(  out) :: ff2(npx,npy,3)
       real(REAL8), intent(  out) :: gg1(npx,npy,3)
       real(REAL8), intent(  out) :: gg2(npx,npy,3)
       real(REAL8), pointer       ::  e1(:,:,:)
       real(REAL8), pointer       ::  e2(:,:,:)
       real(REAL8), pointer       ::  f1(:,:,:)
       real(REAL8), pointer       ::  f2(:,:,:)
       real(REAL8), pointer       ::  g1(:,:,:)
       real(REAL8), pointer       ::  g2(:,:,:)
       real(REAL8), optional      :: sublons(:)
       real(REAL8), optional      :: sublats(:)
       logical , optional      :: AmNodeRoot
       logical , optional      :: WriteNetcdf
     end subroutine GetWeights
  end interface

 include "netcdf.inc"

 integer, intent(in) :: npx, npy, nlon, nlat
 real, dimension(npx , npy ), intent(in ) :: data_cs
 real, dimension(nlon, nlat), intent(out) :: data_ll

 integer :: ntiles=6
 integer :: npts
 integer :: ndims=2

 integer                               :: c2l_unit
 character(len=ESMF_MAXSTR)            :: c2l_fname
 logical                               :: c2l_file_exists

 integer    , pointer                          :: id1(:,:)
 integer    , pointer                          :: id2(:,:)
 integer    , pointer                          :: jdc(:,:)
 real(REAL8), pointer                          :: ee1(:,:,:)
 real(REAL8), pointer                          :: ee2(:,:,:)
 real(REAL8), pointer                          :: ff1(:,:,:)
 real(REAL8), pointer                          :: ff2(:,:,:)
 real(REAL8), pointer                          :: gg1(:,:,:)
 real(REAL8), pointer                          :: gg2(:,:,:)
 real(REAL8), pointer                          :: l2c(:,:,:)
 real(REAL8), pointer                          ::  e1(:,:,:)
 real(REAL8), pointer                          ::  e2(:,:,:)
 real(REAL8), pointer                          ::  f1(:,:,:)
 real(REAL8), pointer                          ::  f2(:,:,:)
 real(REAL8), pointer                          ::  g1(:,:,:)
 real(REAL8), pointer                          ::  g2(:,:,:)

 real(REAL8), pointer                          :: c2l_weight(:,:,:)               
 integer    , pointer                          :: c2l_index(:,:,:)

 real(REAL8) :: varmisval=MAPL_UNDEF
 real(ESMF_KIND_R4), allocatable :: var_cubsph(:,:)

 integer :: grid_type = 0

 integer :: i,j,n,l,itile, j1,j2

  npts = npx+1

   allocate( c2l_index(3,nlon,nlat) )
   allocate( c2l_weight(4,nlon,nlat) )
   allocate( id1(npx,npy) )
   allocate( id2(npx,npy) )
   allocate( jdc(npx,npy) )
   allocate( l2c(4,nlon,nlat) )
   allocate( ee1(npx,npy,3))
   allocate( ee2(npx,npy,3))
   allocate( ff1(npx,npy,3))
   allocate( ff2(npx,npy,3))
   allocate( gg1(npx,npy,3))
   allocate( gg2(npx,npy,3))

   print*, 'GetWeights:', npx, npy, nlon, nlat
   call GetWeights(npx, npy, nlat, nlon, c2l_index, c2l_weight, id1, id2, jdc, l2c, &
     ee1, ee2, ff1, ff2, gg1, gg2, e1, e2, f1, f2, g1, g2, AmNodeRoot=.true., WriteNetcdf=.true.)

#ifdef TEST
  !--------------------------------------------------------------------!
  ! perform interpolation                                              !
  !--------------------------------------------------------------------!
   print*, 'Test Interpolation'
   allocate ( var_cubsph(0:npts,0:npts) )
   do itile=1,ntiles
     j1 = (npts-1)*(itile-1) + 1
     j2 = (npts-1)*(itile-1) + npts-1
     var_cubsph(1:npts-1,1:npts-1)=data_cs(:,j1:j2)
! DSK fill the halo corners with varmisval since ghost_cubsph_update
! doesn't touch these cells
     var_cubsph(0   ,0   ) = varmisval
     var_cubsph(npts,0   ) = varmisval
     var_cubsph(0   ,npts) = varmisval
     var_cubsph(npts,npts) = varmisval
     call ghost_update(var_cubsph, data_cs, 0, npts, 0, npts, 1, ntiles,  &
                              1, itile, npx, npy, A_grid)
     call do_c2l_interpolation_r4(var_cubsph, 0, npts, 0, npts, 1, itile, 1, &
                               data_ll, nlon, nlat, c2l_index, c2l_weight,                      &
                               .true., varmisval, .true.)
   enddo
   deallocate ( var_cubsph )

! !--------------------------------------------------------------------!
! ! Lons need to follow -180:180 MAPL convention                       !
! !--------------------------------------------------------------------!
! call FlipLons(data_ll)
#endif

! NOTE: CAUTION: Code assumes interpolation is the same throughout execution !!!
! SAVE Variables to avoid costly Initialization of weights
   print*, 'Deallocate'
   deallocate( c2l_index )
   deallocate( c2l_weight )
   deallocate( id1 )
   deallocate( id2 )
   deallocate( jdc )
   deallocate( l2c )
   deallocate( ee1 )
   deallocate( ee2 )
   deallocate( ff1 )
   deallocate( ff2 )
   deallocate( gg1 )
   deallocate( gg2 )
! SAVE Variables to avoid costly Initialization of weights

   print*, 'Done'

contains

  subroutine ghost_update(var, vari, n1x, nx, n1y, ny, ng, ntiles, &
                                 ng_update, tile, npx, npy, grid_type)
    !------------------------------------------------------------------!
    ! serial ghost cell update for cubed sphere for global array var   !
    ! update ng (inner) ghost cells                                    !
    !------------------------------------------------------------------!
    integer, intent(in) :: n1x, nx, n1y, ny, ng, ntiles,            &
                                 ng_update, tile, npx, npy, grid_type

    real, dimension(npx,npy), intent(in) :: vari
    real, dimension(n1x:nx,n1y:ny), intent(inout) :: var
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    integer :: i, j, i_in, j_in, i_out, j_out, ig,                   &
               n1x_comp, nx_comp, n1y_comp, ny_comp

    n1x_comp=n1x+ng
    nx_comp=nx-ng
    n1y_comp=n1y+ng
    ny_comp=ny-ng
    if (grid_type==B_grid) then
       select case (tile)

       case (1)
             do i=n1x_comp,nx_comp
                do ig=1,ng_update
                   i_in =i
                   j_in =n1y_comp-ig
                   i_out=i
                   j_out=ny_comp-ig + 5*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)

                   i_in =i
                   j_in =ny_comp+ig
                   i_out=n1x_comp+ig
                   j_out=ny_comp-i+1 + 2*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)
                enddo
             enddo
             do j=n1y_comp,ny_comp
                do ig=1,ng_update
                   i_in =n1x_comp-ig
                   j_in =j
                   i_out=nx_comp-j+1
                   j_out=ny_comp-ig + 4*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)

                   i_in =nx_comp+ig
                   j_in =j
                   i_out=n1x_comp+ig
                   j_out=j + 1*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)
                enddo
             enddo

       case (2)
             do i=n1x_comp,nx_comp
                do ig=1,ng_update
                   i_in =i
                   j_in =n1y_comp-ig
                   i_out=nx_comp-ig
                   j_out=ny_comp-i+1 + 5*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)

                   i_in =i
                   j_in =ny_comp+ig
                   i_out=i
                   j_out=n1y_comp+ig + 2*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)
                enddo
             enddo
             do j=n1y_comp,ny_comp
                do ig=1,ng_update
                   i_in =n1x_comp-ig
                   j_in =j
                   i_out=nx_comp-ig
                   j_out=j
                   var(i_in,j_in)=vari(i_out,j_out)

                   i_in =nx_comp+ig
                   j_in =j
                   i_out=nx_comp-j+1
                   j_out=n1y_comp+ig + 3*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)
                enddo
             enddo

       case (3)
             do i=n1x_comp,nx_comp
                do ig=1,ng_update
                   i_in =i
                   j_in =n1y_comp-ig
                   i_out=i
                   j_out=ny_comp-ig + 1*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)

                   i_in =i
                   j_in =ny_comp+ig
                   i_out=n1x_comp+ig
                   j_out=ny_comp-i+1 + 4*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)
                enddo
             enddo
             do j=n1y_comp,ny_comp
                do ig=1,ng_update
                   i_in =n1x_comp-ig
                   j_in =j
                   i_out=nx_comp-j+1
                   j_out=ny_comp-ig
                   var(i_in,j_in)=vari(i_out,j_out)

                   i_in =nx_comp+ig
                   j_in =j
                   i_out=n1x_comp+ig
                   j_out=j + 3*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)
                enddo
             enddo

       case (4)
             do i=n1x_comp,nx_comp
                do ig=1,ng_update
                   i_in =i
                   j_in =n1y_comp-ig
                   i_out=nx_comp-ig
                   j_out=ny_comp-i+1 + 1*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)

                   i_in =i
                   j_in =ny_comp+ig
                   i_out=i
                   j_out=n1y_comp+ig + 4*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)
                enddo
             enddo
             do j=n1y_comp,ny_comp
                do ig=1,ng_update
                   i_in =n1x_comp-ig
                   j_in =j
                   i_out=nx_comp-ig
                   j_out=j + 2*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)

                   i_in =nx_comp+ig
                   j_in =j
                   i_out=nx_comp-j+1
                   j_out=n1y_comp+ig + 5*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)
                enddo
             enddo

       case (5)
             do i=n1x_comp,nx_comp
                do ig=1,ng_update
                   i_in =i
                   j_in =n1y_comp-ig
                   i_out=i
                   j_out=ny_comp-ig + 3*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)

                   i_in =i
                   j_in =ny_comp+ig
                   i_out=n1x_comp+ig
                   j_out=ny_comp-i+1
                   var(i_in,j_in)=vari(i_out,j_out)
                enddo
             enddo
             do j=n1y_comp,ny_comp
                do ig=1,ng_update
                   i_in =n1x_comp-ig
                   j_in =j
                   i_out=nx_comp-j+1
                   j_out=ny_comp-ig + 2*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)

                   i_in =nx_comp+ig
                   j_in =j
                   i_out=n1x_comp+ig
                   j_out=j + 5*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)
                enddo
             enddo

       case (6)
             do i=n1x_comp,nx_comp
                do ig=1,ng_update
                   i_in =i
                   j_in =n1y_comp-ig
                   i_out=nx_comp-ig
                   j_out=ny_comp-i+1 + 3*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)

                   i_in =i
                   j_in =ny_comp+ig
                   i_out=i
                   j_out=n1y_comp+ig
                   var(i_in,j_in)=vari(i_out,j_out)
                enddo
             enddo
             do j=n1y_comp,ny_comp
                do ig=1,ng_update
                   i_in =n1x_comp-ig
                   j_in =j
                   i_out=nx_comp-ig
                   j_out=j + 4*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)

                   i_in =nx_comp+ig
                   j_in =j
                   i_out=nx_comp-j+1
                   j_out=n1y_comp+ig + 1*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)
                enddo
             enddo
       end select

    elseif (grid_type==A_grid) then

       select case (tile)

       case (1)
             do i=n1x_comp,nx_comp
                do ig=1,ng_update
                   i_in =i
                   j_in =n1y_comp-ig
                   i_out=i
                   j_out=ny_comp-ig+1 + 5*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)

                   i_in =i
                   j_in =ny_comp+ig
                   i_out=n1x_comp+ig-1
                   j_out=ny_comp-i+1 + 2*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)
                enddo
             enddo
             do j=n1y_comp,ny_comp
                do ig=1,ng_update
                   i_in =n1x_comp-ig
                   j_in =j
                   i_out=nx_comp-j+1
                   j_out=ny_comp-ig+1 + 4*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)

                   i_in =nx_comp+ig
                   j_in =j
                   i_out=n1x_comp+ig-1
                   j_out=j + 1*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)
                enddo
             enddo

       case (2)
             do i=n1x_comp,nx_comp
                do ig=1,ng_update
                   i_in =i
                   j_in =n1y_comp-ig
                   i_out=nx_comp-ig+1
                   j_out=ny_comp-i+1 + 5*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)

                   i_in =i
                   j_in =ny_comp+ig
                   i_out=i
                   j_out=n1y_comp+ig-1 + 2*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)
                enddo
             enddo
             do j=n1y_comp,ny_comp
                do ig=1,ng_update
                   i_in =n1x_comp-ig
                   j_in =j
                   i_out=nx_comp-ig+1
                   j_out=j
                   var(i_in,j_in)=vari(i_out,j_out)

                   i_in =nx_comp+ig
                   j_in =j
                   i_out=nx_comp-j+1
                   j_out=n1y_comp+ig-1 + 3*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)
                enddo
             enddo

       case (3)
             do i=n1x_comp,nx_comp
                do ig=1,ng_update
                   i_in =i
                   j_in =n1y_comp-ig
                   i_out=i
                   j_out=ny_comp-ig+1 + 1*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)

                   i_in =i
                   j_in =ny_comp+ig
                   i_out=n1x_comp+ig-1
                   j_out=ny_comp-i+1 + 4*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)
                enddo
             enddo
             do j=n1y_comp,ny_comp
                do ig=1,ng_update
                   i_in =n1x_comp-ig
                   j_in =j
                   i_out=nx_comp-j+1
                   j_out=ny_comp-ig+1
                   var(i_in,j_in)=vari(i_out,j_out)

                   i_in =nx_comp+ig
                   j_in =j
                   i_out=n1x_comp+ig-1
                   j_out=j + 3*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)
                enddo
             enddo

       case (4)
             do i=n1x_comp,nx_comp
                do ig=1,ng_update
                   i_in =i
                   j_in =n1y_comp-ig
                   i_out=nx_comp-ig+1
                   j_out=ny_comp-i+1 + 1*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)

                   i_in =i
                   j_in =ny_comp+ig
                   i_out=i
                   j_out=n1y_comp+ig-1 + 4*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)
                enddo
             enddo
             do j=n1y_comp,ny_comp
                do ig=1,ng_update
                   i_in =n1x_comp-ig
                   j_in =j
                   i_out=nx_comp-ig+1
                   j_out=j + 2*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)

                   i_in =nx_comp+ig
                   j_in =j
                   i_out=nx_comp-j+1
                   j_out=n1y_comp+ig-1 + 5*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)
                enddo
             enddo

       case (5)
             do i=n1x_comp,nx_comp
                do ig=1,ng_update
                   i_in =i
                   j_in =n1y_comp-ig
                   i_out=i
                   j_out=ny_comp-ig+1 + 3*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)

                   i_in =i
                   j_in =ny_comp+ig
                   i_out=n1x_comp+ig-1
                   j_out=ny_comp-i+1
                   var(i_in,j_in)=vari(i_out,j_out)
                enddo
             enddo
             do j=n1y_comp,ny_comp
                do ig=1,ng_update
                   i_in =n1x_comp-ig
                   j_in =j
                   i_out=nx_comp-j+1
                   j_out=ny_comp-ig+1 + 2*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)

                   i_in =nx_comp+ig
                   j_in =j
                   i_out=n1x_comp+ig-1
                   j_out=j + 5*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)
                enddo
             enddo

       case (6)
             do i=n1x_comp,nx_comp
                do ig=1,ng_update
                   i_in =i
                   j_in =n1y_comp-ig
                   i_out=nx_comp-ig+1
                   j_out=ny_comp-i+1 + 3*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)

                   i_in =i
                   j_in =ny_comp+ig
                   i_out=i
                   j_out=n1y_comp+ig-1
                   var(i_in,j_in)=vari(i_out,j_out)
                enddo
             enddo
             do j=n1y_comp,ny_comp
                do ig=1,ng_update
                   i_in =n1x_comp-ig
                   j_in =j
                   i_out=nx_comp-ig+1
                   j_out=j + 4*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)

                   i_in =nx_comp+ig
                   j_in =j
                   i_out=nx_comp-j+1
                   j_out=n1y_comp+ig-1 + 1*(npy/6)
                   var(i_in,j_in)=vari(i_out,j_out)
                enddo
             enddo
       end select
    else
       stop 'grid_type incorrect'
    end if

  end subroutine ghost_update
  !====================================================================!

!-------------------------------------------------------------------------
   subroutine FlipLons(q)
!-------------------------------------------------------------------------
     implicit none
     real,dimension(:,:),intent(inout) :: q
     integer :: im, j
     real,dimension(size(q,1)/2) :: d
     im=size(q,1)
     do j=1,size(q,2)
       d(     1:im/2  ) = q(     1:im/2,j)
       q(     1:im/2,j) = q(im/2+1:im  ,j)
       q(im/2+1:im  ,j) = d(     1:im/2  )
     end do
   end subroutine FlipLons

end subroutine cube2latlon

