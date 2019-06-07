
module CubeHaloMod
  use mpp_domains_mod,   only: domain2d
  use mpp_mod,           only: mpp_error, FATAL, NOTE
  implicit none

#define MAX_HALOTYPES 2
  type(domain2D) :: domainList(MAX_HALOTYPES)

  public :: mpp_domain_decomp
  public :: domainList

  contains


!ALT: Next subroutine is a straight copy from the following file:
!     fvdycore/tools/external_ic.F90
      subroutine mpp_domain_decomp(domain,npx,npy,nregions,ng,grid_type, &
!@                                   is,ie,js,je,isd,ied,jsd,jed,tile, &
                                   npes, npes_x, npes_y)
         use mpp_domains_mod, only: mpp_domains_init, MPP_DOMAIN_TIME, &
                              mpp_define_mosaic, mpp_get_compute_domain, mpp_get_data_domain, &
                              mpp_domains_set_stack_size, mpp_define_layout
         use mpp_mod,         only : mpp_pe
         type(domain2D), intent(OUT) :: domain
         integer, intent(IN)  :: npx,npy,nregions,ng,grid_type
!@         integer, intent(OUT) :: is,ie,js,je,isd,ied,jsd,jed,tile
         integer, intent(IN)  :: npes, npes_x, npes_y

         integer :: layout(2)
         integer, allocatable :: pe_start(:), pe_end(:)

         integer :: num_contact, ntiles, npes_per_tile
         integer, allocatable, dimension(:)       :: npes_tile, tile1, tile2
         integer, allocatable, dimension(:)       :: istart1, iend1, jstart1, jend1
         integer, allocatable, dimension(:)       :: istart2, iend2, jstart2, jend2
         integer, allocatable, dimension(:,:)     :: layout2D, global_indices

         character*80 :: evalue
         integer :: ios,nx,ny,n,num_alloc
         character(len=32) :: type = "unknown"
         logical :: is_symmetry
         logical :: debug=.false.

         nx = npx-1
         ny = npy-1

!@         call print_memuse_stats('external_ic:mpp_domain_decomp: top')

         call mpp_domains_init(MPP_DOMAIN_TIME)

!@         call print_memuse_stats('external_ic:mpp_domain_decomp: mpp_domains_init')

       ! call mpp_domains_set_stack_size(10000)
       ! call mpp_domains_set_stack_size(900000)
       ! call mpp_domains_set_stack_size(1500000)
         call mpp_domains_set_stack_size(3000000)

         select case(nregions)
         case ( 1 )  ! Lat-Lon "cyclic"

            select case (grid_type)
            case (3)   ! Lat-Lon "cyclic"
               type="Lat-Lon: cyclic"
               ntiles = 4
               num_contact = 8
               if( mod(npes,ntiles) .NE. 0 ) then
                  call mpp_error(NOTE,'TEST_MPP_DOMAINS: for Cyclic mosaic, npes should be multiple of ntiles. ' // &
                                       'No test is done for Cyclic mosaic. ' )
                  return
               end if
               npes_per_tile = npes/ntiles
               call mpp_define_layout( (/1,npx-1,1,npy-1/), npes_per_tile, layout )
               layout = (/1,npes_per_tile/) ! force decomp only in Lat-Direction
            case (4)   ! Cartesian, double periodic
               type="Cartesian: double periodic"
               ntiles = 1
               num_contact = 2
               npes_per_tile = npes/ntiles
               call mpp_define_layout( (/1,npx-1,1,npy-1/), npes_per_tile, layout )
            case (5)   ! latlon patch
               type="Lat-Lon: patch"
               ntiles = 1
               num_contact = 0
               npes_per_tile = npes/ntiles
               call mpp_define_layout( (/1,npx-1,1,npy-1/), npes_per_tile, layout )
            case (6)   ! latlon strip
               type="Lat-Lon: strip"
               ntiles = 1
               num_contact = 1
               npes_per_tile = npes/ntiles
               call mpp_define_layout( (/1,npx-1,1,npy-1/), npes_per_tile, layout )
            case (7)   ! Cartesian, channel
               type="Cartesian: channel"
               ntiles = 1
               num_contact = 1
               npes_per_tile = npes/ntiles
               call mpp_define_layout( (/1,npx-1,1,npy-1/), npes_per_tile, layout )
            end select

         case ( 6 )  ! Cubed-Sphere
            type="Cubic: cubed-sphere"
            ntiles = 6
            num_contact = 12
            !--- cubic grid always have six tiles, so npes should be multiple of 6
            if( mod(npes,ntiles) .NE. 0 .OR. npx-1 .NE. npy-1) then
               call mpp_error(NOTE,'mpp_domain_decomp: for Cubic_grid mosaic, npes should be multiple of ntiles(6) ' // &
                                   'and npx-1 should equal npy-1, mpp_domain_decomp is NOT done for Cubic-grid mosaic. ' )
               return
            end if
            npes_per_tile = npes/ntiles
            call  mpp_define_layout( (/1,npx-1,1,npy-1/), npes_per_tile, layout )

!@            if ( npes_x == 0 ) then
!@               npes_x = layout(1)
!@            endif
!@            if ( npes_y == 0 ) then
!@               npes_y = layout(2)
!@            endif

            if ( (npx/npes_x < ng) .or. (npy/npes_y < ng) ) then
               write(*,310) npes_x, npes_y, npx/npes_x, npy/npes_y
 310           format('Invalid layout, NPES_X:',i4.4,'NPES_Y:',i4.4,'ncells_X:',i4.4,'ncells_Y:',i4.4)
               call mpp_error(FATAL, 'mpp_domain_decomp: bad decomp')
            endif

            layout = (/npes_x,npes_y/)
         case default
            call mpp_error(FATAL, 'mpp_domain_decomp: no such test: '//type)
         end select

!@         call print_memuse_stats('external_ic:mpp_domain_decomp: mpp_define_layout')

         allocate(layout2D(2,ntiles), global_indices(4,ntiles), npes_tile(ntiles) )
         allocate(pe_start(ntiles),pe_end(ntiles))
         npes_tile = npes_per_tile
         do n = 1, ntiles
            global_indices(:,n) = (/1,npx-1,1,npy-1/)
            layout2D(:,n)         = layout
            pe_start(n) = (n-1)*layout(1)*layout(2)
            pe_end(n)   = pe_start(n) + layout(1)*layout(2) -1
         end do
         num_alloc=max(1,num_contact)
         allocate(tile1(num_alloc), tile2(num_alloc) )
         allocate(istart1(num_alloc), iend1(num_alloc), jstart1(num_alloc), jend1(num_alloc) )
         allocate(istart2(num_alloc), iend2(num_alloc), jstart2(num_alloc), jend2(num_alloc) )

         is_symmetry = .true.

!@         call print_memuse_stats('external_ic:mpp_domain_decomp: allocates 1')

         select case(nregions)
         case ( 1 )

            select case (grid_type)
            case (3)   ! Lat-Lon "cyclic"
               !--- Contact line 1, between tile 1 (EAST) and tile 2 (WEST)
               tile1(1) = 1; tile2(1) = 2
               istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
               istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
               !--- Contact line 2, between tile 1 (SOUTH) and tile 3 (NORTH)  --- cyclic
               tile1(2) = 1; tile2(2) = 3
               istart1(2) = 1;  iend1(2) = nx; jstart1(2) = 1;   jend1(2) = 1
               istart2(2) = 1;  iend2(2) = nx; jstart2(2) = ny;  jend2(2) = ny
               !--- Contact line 3, between tile 1 (WEST) and tile 2 (EAST) --- cyclic
               tile1(3) = 1; tile2(3) = 2
               istart1(3) = 1;  iend1(3) = 1;  jstart1(3) = 1;  jend1(3) = ny
               istart2(3) = nx; iend2(3) = nx; jstart2(3) = 1;  jend2(3) = ny
               !--- Contact line 4, between tile 1 (NORTH) and tile 3 (SOUTH)
               tile1(4) = 1; tile2(4) = 3
               istart1(4) = 1;  iend1(4) = nx; jstart1(4) = ny;  jend1(4) = ny
               istart2(4) = 1;  iend2(4) = nx; jstart2(4) = 1;   jend2(4) = 1
               !--- Contact line 5, between tile 2 (SOUTH) and tile 4 (NORTH) --- cyclic
               tile1(5) = 2; tile2(5) = 4
               istart1(5) = 1;  iend1(5) = nx; jstart1(5) = 1;  jend1(5) = 1
               istart2(5) = 1;  iend2(5) = nx; jstart2(5) = ny; jend2(5) = ny
               !--- Contact line 6, between tile 2 (NORTH) and tile 4 (SOUTH)
               tile1(6) = 2; tile2(6) = 4
               istart1(6) = 1;  iend1(6) = nx; jstart1(6) = ny;  jend1(6) = ny
               istart2(6) = 1;  iend2(6) = nx; jstart2(6) = 1;   jend2(6) = 1
               !--- Contact line 7, between tile 3 (EAST) and tile 4 (WEST)
               tile1(7) = 3; tile2(7) = 4
               istart1(7) = nx; iend1(7) = nx; jstart1(7) = 1;  jend1(7) = ny
               istart2(7) = 1;  iend2(7) = 1;  jstart2(7) = 1;  jend2(7) = ny
               !--- Contact line 8, between tile 3 (WEST) and tile 4 (EAST) --- cyclic
               tile1(8) = 3; tile2(8) = 4
               istart1(8) = 1;  iend1(8) = 1;  jstart1(8) = 1;  jend1(8) = ny
               istart2(8) = nx; iend2(8) = nx; jstart2(8) = 1;  jend2(8) = ny
               is_symmetry = .false.
            case (4)   ! Cartesian, double periodic
               !--- Contact line 1, between tile 1 (EAST) and tile 1 (WEST)
               tile1(1) = 1; tile2(1) = 1
               istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
               istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
               !--- Contact line 2, between tile 1 (SOUTH) and tile 1 (NORTH)  --- cyclic
               tile1(2) = 1; tile2(2) = 1
               istart1(2) = 1;  iend1(2) = nx; jstart1(2) = 1;   jend1(2) = 1
               istart2(2) = 1;  iend2(2) = nx; jstart2(2) = ny;  jend2(2) = ny
            case (5)   ! latlon patch

            case (6)   !latlon strip
               !--- Contact line 1, between tile 1 (EAST) and tile 1 (WEST)
               tile1(1) = 1; tile2(1) = 1
               istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
               istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
            case (7)   ! Cartesian, channel
               !--- Contact line 1, between tile 1 (EAST) and tile 1 (WEST)
               tile1(1) = 1; tile2(1) = 1
               istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
               istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
            end select

         case ( 6 )  ! Cubed-Sphere
            !--- Contact line 1, between tile 1 (EAST) and tile 2 (WEST)
            tile1(1) = 1; tile2(1) = 2
            istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
            istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
            !--- Contact line 2, between tile 1 (NORTH) and tile 3 (WEST)
            tile1(2) = 1; tile2(2) = 3
            istart1(2) = 1;  iend1(2) = nx; jstart1(2) = ny; jend1(2) = ny
            istart2(2) = 1;  iend2(2) = 1;  jstart2(2) = ny; jend2(2) = 1
            !--- Contact line 3, between tile 1 (WEST) and tile 5 (NORTH)
            tile1(3) = 1; tile2(3) = 5
            istart1(3) = 1;  iend1(3) = 1;  jstart1(3) = 1;  jend1(3) = ny
            istart2(3) = nx; iend2(3) = 1;  jstart2(3) = ny; jend2(3) = ny
            !--- Contact line 4, between tile 1 (SOUTH) and tile 6 (NORTH)
            tile1(4) = 1; tile2(4) = 6
            istart1(4) = 1;  iend1(4) = nx; jstart1(4) = 1;  jend1(4) = 1
            istart2(4) = 1;  iend2(4) = nx; jstart2(4) = ny; jend2(4) = ny
            !--- Contact line 5, between tile 2 (NORTH) and tile 3 (SOUTH)
            tile1(5) = 2; tile2(5) = 3
            istart1(5) = 1;  iend1(5) = nx; jstart1(5) = ny; jend1(5) = ny
            istart2(5) = 1;  iend2(5) = nx; jstart2(5) = 1;  jend2(5) = 1
            !--- Contact line 6, between tile 2 (EAST) and tile 4 (SOUTH)
            tile1(6) = 2; tile2(6) = 4
            istart1(6) = nx; iend1(6) = nx; jstart1(6) = 1;  jend1(6) = ny
            istart2(6) = nx; iend2(6) = 1;  jstart2(6) = 1;  jend2(6) = 1
            !--- Contact line 7, between tile 2 (SOUTH) and tile 6 (EAST)
            tile1(7) = 2; tile2(7) = 6
            istart1(7) = 1;  iend1(7) = nx; jstart1(7) = 1;  jend1(7) = 1
            istart2(7) = nx; iend2(7) = nx; jstart2(7) = ny; jend2(7) = 1
            !--- Contact line 8, between tile 3 (EAST) and tile 4 (WEST)
            tile1(8) = 3; tile2(8) = 4
            istart1(8) = nx; iend1(8) = nx; jstart1(8) = 1;  jend1(8) = ny
            istart2(8) = 1;  iend2(8) = 1;  jstart2(8) = 1;  jend2(8) = ny
            !--- Contact line 9, between tile 3 (NORTH) and tile 5 (WEST)
            tile1(9) = 3; tile2(9) = 5
            istart1(9) = 1;  iend1(9) = nx; jstart1(9) = ny; jend1(9) = ny
            istart2(9) = 1;  iend2(9) = 1;  jstart2(9) = ny; jend2(9) = 1
            !--- Contact line 10, between tile 4 (NORTH) and tile 5 (SOUTH)
            tile1(10) = 4; tile2(10) = 5
            istart1(10) = 1;  iend1(10) = nx; jstart1(10) = ny; jend1(10) = ny
            istart2(10) = 1;  iend2(10) = nx; jstart2(10) = 1;  jend2(10) = 1
            !--- Contact line 11, between tile 4 (EAST) and tile 6 (SOUTH)
            tile1(11) = 4; tile2(11) = 6
            istart1(11) = nx; iend1(11) = nx; jstart1(11) = 1;  jend1(11) = ny
            istart2(11) = nx; iend2(11) = 1;  jstart2(11) = 1;  jend2(11) = 1
            !--- Contact line 12, between tile 5 (EAST) and tile 6 (WEST)
            tile1(12) = 5; tile2(12) = 6
            istart1(12) = nx; iend1(12) = nx; jstart1(12) = 1;  jend1(12) = ny
            istart2(12) = 1;  iend2(12) = 1;  jstart2(12) = 1;  jend2(12) = ny
         end select

       call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, num_contact, tile1, tile2, &
                              istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
                              pe_start=pe_start, pe_end=pe_end, symmetry=is_symmetry,              &
                              shalo = ng, nhalo = ng, whalo = ng, ehalo = ng, name = type)
!@       call print_memuse_stats('external_ic:mpp_domain_decomp: mpp_define_mosaic')

       deallocate(pe_start,pe_end)

!@        !--- find the tile number
!@         tile = mpp_pe()/npes_per_tile+1
!@         call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
!@         call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )

!@         call print_memuse_stats('external_ic:mpp_domain_decomp: mpp_get domains')

      end subroutine mpp_domain_decomp

end module CubeHaloMod

subroutine CubeHalo(domainIdx, input)
  use cubehalomod
  use mpp_domains_mod, only: mpp_update_domains
  implicit none
  integer :: domainIdx
  real*4  :: input(:, :)

  type(domain2d) :: mydomain

  mydomain = domainList(domainIdx)
  call mpp_update_domains(input, mydomain)
  return
end subroutine CubeHalo

subroutine CubeHaloInit(comm, IM_world, npes, NX, NY, domainIdx)
  use cubehalomod
  use fms_mod,         only: fms_init
  implicit none
  integer :: comm
  integer :: im_world
  integer :: npes
  integer :: nx
  integer :: ny
  integer :: domainIdx
  
  type(domain2d) :: mydomain
  integer        :: npx, npy
  integer        :: ng, nregions, grid_type

  nregions = 6
  npx = IM_WORLD+1
  npy = npx
  ng = 1        ! halowidth
  grid_type = 0 ! cubed-sphere

  call fms_init(comm)
  call mpp_domain_decomp(mydomain,npx,npy,nregions,ng,grid_type, &
       npes,nx,ny/nregions)
  domainList(domainIdx) = mydomain
end subroutine CubeHaloInit
