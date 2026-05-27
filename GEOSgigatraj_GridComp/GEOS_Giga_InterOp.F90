! This module define the interface bewteen GEOS and gigatraj
! The functions are defined in gigatraj

module GEOS_Giga_InterOpMod
   use, intrinsic :: iso_c_binding, only : c_double, c_int, c_ptr, c_null_char, c_associated
   use, intrinsic :: iso_c_binding, only : c_loc, c_null_ptr
   use mpi
   implicit none
   private

   public :: initMetGEOSDistributedLatLonData
   public :: initMetGEOSDistributedCubedData
   public :: updateFields
   public :: RK4_advance
   public :: setData
   public :: getData
   public :: getData2d
 
   public :: test_Field3D
   public :: test_dataflow
   public :: test_metData

   interface

     function initMetGEOSDistributedCubedData(comm, ijToRank, Ig, lev, i1, i2, j1, j2, nzs, lons_ptr, lats_ptr, eta_ptr, ctime_ptr) result (metdata_ptr) bind(C, name="initGigaGridDistributedCubedData")
       import :: c_int, c_ptr
       implicit none
       integer(c_int), intent(in), value :: comm, Ig, lev, i1,i2,j1,j2, nzs
       type(c_ptr), intent(in), value    :: ijToRank, lons_ptr, lats_ptr, eta_ptr, ctime_ptr
       type(c_ptr) :: metdata_ptr
     end function

     function initMetGEOSDistributedLatLonData(comm, ijToRank, Ig, Jg,lev, nlon_local, nlat_local, nzs, lons_ptr, lats_ptr, eta_ptr, ctime_ptr) result (metdata_ptr) bind(C, name="initGigaGridDistributedLatLonData")
       import :: c_int, c_ptr
       implicit none
       integer(c_int), intent(in), value :: comm, Ig, Jg, lev, nlon_local, nlat_local, nzs
       type(c_ptr), intent(in), value    :: ijToRank, lons_ptr, lats_ptr, eta_ptr, ctime_ptr
       type(c_ptr) :: metdata_ptr
     end function

     subroutine updateFields( metSrc_ptr, ctime_ptr, u_ptr, v_ptr, w_ptr, p_ptr) bind(C, name="updateFields")
       import :: c_ptr
       implicit none
       type(c_ptr), intent(in), value    :: metSrc_ptr, ctime_ptr, u_ptr, v_ptr, w_ptr, p_ptr
     end subroutine

     subroutine RK4_advance(metsrc_ptr, ctime_ptr, dt, n, lons_ptr, lats_ptr, levs_ptr) bind( C, name='RK4_advance')
       import :: c_ptr, c_int, c_double
       type(c_ptr), intent(in), value    :: metsrc_ptr
       real(c_double), intent(in), value :: dt
       integer(c_int), intent(in), value :: n
       type(c_ptr), intent(in), value    :: ctime_ptr, lons_ptr, lats_ptr, levs_ptr
     end subroutine 

     subroutine test_Field3d(obj_ptr) bind(C, name="test_Field3D")
       import :: c_ptr
       implicit none
       type(c_ptr), intent(in), value :: obj_ptr
     end subroutine

     subroutine test_metData(obj_ptr, time, n, lons_ptr, lats_ptr, levs_ptr, u_ptr, v_ptr, w_ptr) bind(C, name="test_metData")
       import :: c_ptr,c_int, c_double
       type(c_ptr), intent(in), value :: obj_ptr
       real(c_double), intent(in), value :: time
       integer(c_int), intent(in), value :: n
       type(c_ptr), intent(in), value :: lons_ptr, lats_ptr, levs_ptr, u_ptr, v_ptr, w_ptr
     end subroutine
   
     subroutine setData ( metSrc_ptr, ctime, quantity_ptr, data_ptr) bind(C, name="setData")
       import :: c_ptr
       type(c_ptr), intent(in), value :: metSrc_ptr, ctime, quantity_ptr, data_ptr
     end subroutine setData

     subroutine getData ( metSrc_ptr, ctime, quantity_ptr, n, lons_ptr, lats_ptr, levs_ptr, values_ptr) bind(C, name="getData")
       import :: c_ptr, c_int
       integer(c_int), intent(in), value :: n
       type(c_ptr), intent(in), value :: metSrc_ptr, ctime, quantity_ptr, lons_ptr, lats_ptr, levs_ptr, values_ptr
     end subroutine getData

     subroutine getData2d ( metSrc_ptr, ctime, quantity_ptr, n, lons_ptr, lats_ptr, values_ptr) bind(C, name="getData2d")
       import :: c_ptr, c_int
       integer(c_int), intent(in), value :: n
       type(c_ptr), intent(in), value :: metSrc_ptr, ctime, quantity_ptr, lons_ptr, lats_ptr, values_ptr
     end subroutine getData2d
   end interface

contains

  subroutine test_dataflow(num_parcels, lons, lats, zs, CellToRank, DIMS, comm)
    integer :: num_parcels, comm, DIMS(3)
    real, dimension(:), intent(in) :: lons, lats,zs
    integer, dimension(:,:), intent(in) :: CellToRank

    integer :: i, npes, ierror, rank, my_rank
    real :: dlon, dlat
    real, allocatable :: lons_positive(:)

    real, allocatable :: lons_send(:), lats_send(:), zs_send(:)
    real, allocatable :: lons_recv(:), lats_recv(:), zs_recv(:)
    real, allocatable :: U_recv(:), U_send(:)
    real, allocatable :: U(:), V(:), W(:), pos(:)

    integer, allocatable :: counts_send(:),counts_recv(:), II(:), JJ(:), ranks(:)
    integer, allocatable :: disp_send(:), disp_recv(:), tmp_position(:)
      
    dlon = 360.0 / DIMS(1)
    dlat = 180.0 / DIMS(2)

    lons_positive = lons
    where (lons_positive < 0) lons_positive=lons_positive + 360.0
    II = min( max(ceiling (lons_positive/dlon),1), DIMS(1))
    JJ = min( max(ceiling ((lats + 90.0)/dlat),1), DIMS(2))

    call MPI_Comm_size(comm, npes, ierror)
    call MPI_Comm_rank(comm, my_rank, ierror)
    
    allocate(ranks(num_parcels))
    allocate(counts_send(npes))
    allocate(counts_recv(npes))
    allocate(disp_send(npes))
    allocate(disp_recv(npes))

    do i = 1, num_parcels
       ranks(i) = CellToRank(II(i), JJ(i))
    enddo

!-- -------------------
!step 4) Pack the location data and send them to where the metData sit
!-- -------------------

    do rank = 0, npes-1
       counts_send(rank+1) = count(ranks == rank)
    enddo

    call MPI_AllToALL(counts_send, 1, MPI_INTEGER, counts_recv, 1, MPI_INTEGER, comm, ierror)

    disp_send = 0
    do rank = 1, npes-1
       disp_send(rank+1) = disp_send(rank)+ counts_send(rank)
    enddo
    disp_recv = 0
    do rank = 1, npes-1
       disp_recv(rank+1) = disp_recv(rank)+ counts_recv(rank)
    enddo


    ! re-arranged lats lons, and ids
    tmp_position = disp_send
    allocate(lons_send(num_parcels))
    allocate(lons_recv(sum(counts_recv)))
    allocate(lats_send(num_parcels))
    allocate(lats_recv(sum(counts_recv)))
    allocate(zs_send(num_parcels))
    allocate(zs_recv(sum(counts_recv)))

    allocate(pos(num_parcels))
    do i = 1, num_parcels
       rank   = ranks(i)
       pos(i) = tmp_position(rank+1) +1
       lons_send(pos(i)) = lons(i)
       lats_send(pos(i)) = lats(i)
       zs_send(pos(i))   = zs(i)
       tmp_position(rank+1) = tmp_position(rank+1) + 1
    enddo

    call MPI_AllToALLv(lons_send, counts_send, disp_send, MPI_REAL, lons_recv, counts_recv, disp_recv, MPI_REAL, comm, ierror)
    call MPI_AllToALLv(lats_send, counts_send, disp_send, MPI_REAL, lats_recv, counts_recv, disp_recv, MPI_REAL, comm, ierror)
    call MPI_AllToALLv(zs_send,   counts_send, disp_send, MPI_REAL, zs_recv, counts_recv, disp_recv, MPI_REAL, comm, ierror)
!-- -------------------
!step 5) Interpolate the data ( horiontally and vertically) and send back where they are from
!-- -------------------
     allocate(U_recv(sum(counts_recv)), source = my_rank*1.0)
     allocate(U_send(num_parcels), source = -1.0)
     !
     ! Horizontal and vertical interpolator here
     !
     call MPI_AllToALLv(U_recv, counts_recv, disp_recv, MPI_REAL, U_send, counts_send, disp_send, MPI_REAL, comm, ierror)

!---------------------
!step 6) Rearrange data ( not necessary if ids was rearranged ins step 4)
!---------------------

     allocate(U(num_parcels))
     allocate(V(num_parcels))
     allocate(W(num_parcels))
     U(:) = U_send(pos(:))

  end subroutine

end module GEOS_Giga_InterOpMod
