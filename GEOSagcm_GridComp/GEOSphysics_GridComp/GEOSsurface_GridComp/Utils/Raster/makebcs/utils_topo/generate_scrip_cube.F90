!#define _VERIFY(A) if(A/=0) write(*,*)__LINE__;call MPI_ABort(MPI_COMM_WORLD,error_code,status)
#define _VERIFY(A) if(A/=0) call local_abort(A,__LINE__)

    program ESMF_GenerateCSGridDescription

! ESMF Framework module
    use ESMF
    use mpi
    use netcdf
    use, intrinsic :: iso_fortran_env, only: REAL64,REAL32

    implicit none
    
! Local variables to be used in the Regrid method calls.
! The code creating and filling these variables is not included in the
! example documentation because those interfaces are not specific to
! Regrid.
    type(ESMF_Grid) ::  dstgrid
    integer :: rc

!EOC

    real(REAL64), parameter           :: PI = 3.14159265358979323846
    integer                           :: npets, localPet
    integer                           :: i, j, k
    type(ESMF_VM)                     :: vm
    integer                           :: IM_World, JM_World, scrip_size
    integer, parameter                :: grid_type = 0
    integer, parameter                :: KM_WORLD=1
    integer, parameter                :: NX=1
    integer, parameter                :: NY=6
    integer, parameter                :: ntiles=6
    integer, parameter                :: ndims=2
    integer                           :: N
    integer                           :: info
    integer                           :: start(2), cnt(2), hull_num, hull(4)
    integer                           :: UNIT
    real(ESMF_KIND_R8), allocatable   :: grid_global(:,:,:,:)
    real(ESMF_KIND_R8), allocatable   :: SCRIP_CenterLat(:), SCRIP_CenterLon(:)
    real(ESMF_KIND_R8), allocatable   :: SCRIP_CornerLat(:,:), SCRIP_CornerLon(:,:)
    real(ESMF_KIND_R8), allocatable   :: SCRIP_Area(:)
    real(ESMF_KIND_R8), allocatable   :: SCRIP_rrfac(:)
    real(ESMF_KIND_R8)                :: node_xy(2,4), node_xy_tmp(2,4), lon_e, lon_w
    integer                           :: cornerdim, gridsize, griddim, rankdim, mask
    integer                           :: cornerlon, cornerlat, centerlon, centerlat,cellarea, cellrrfac
    integer, allocatable              :: IMS(:,:), JMS(:,:), sendData(:), GlobalCounts(:), recvCounts(:), recvOffsets(:)
    integer, allocatable              :: grid_imask(:)
    character(len=ESMF_MAXSTR)        :: gridname, FMT, title
    integer                           :: NPX, NPY
    integer                           :: myTile
    integer                           :: npts, tmp, mpiC
    integer                           :: IG, JG, rrfac_max
    logical                           :: do_schmidt
    logical, allocatable              :: fallback_mask(:)
    integer                           :: global_start
    integer                           :: start_mask(1), cnt_mask(1)
    integer                           :: varid_mask_fallback  ! NetCDF variable ID
    integer, allocatable              :: mask_fallback(:)     ! Array to hold data
    real(ESMF_KIND_R8)                :: p1(2),p2(2),p3(2),p4(2)
    real(ESMF_KIND_R8)                :: local_max_length, max_length, local_min_length, min_length
    real(ESMF_KIND_R4)                :: target_lon, target_lat, stretch_factor
    type(ESMF_HConfig)                 :: CF
    integer                           :: status,error_code
    real(ESMF_KIND_R8), pointer :: center_lons(:,:),center_lats(:,:),corner_lons(:,:),corner_lats(:,:)
    type(ESMF_CubedSphereTransform_Args) :: transformArgument
    character(len=ESMF_MAXPATHLEN) :: output_scrip, output_geos
    type(ESMF_Field) :: field
    type(ESMF_FieldBundle) :: field_bundle
    type(ESMF_Clock) :: clock
    type(ESMF_Time) :: start_time
    type(ESMF_TimeInterval) :: time_interval
    real(ESMF_KIND_R4), pointer  :: ptr2d(:,:)
    integer :: failed_cells
    real(ESMF_KIND_R8) :: min_area, max_area


    call ESMF_Initialize(logKindFlag=ESMF_LOGKIND_NONE,rc=status)
    _VERIFY(status)
    call ESMF_CalendarSetDefault(ESMF_CALKIND_GREGORIAN,rc=status)
    _VERIFY(status)

    call ESMF_VMGetGlobal(vm, rc=status)
    _VERIFY(status)

! Get number of PETs we are running with
! --------------------------------------
    call ESMF_VMGet(vm, localPet=localPet, petCount=npets, mpiCommunicator=mpiC, rc=status)
    _VERIFY(status)
    if (npets /= 6) call local_abort(1,__LINE__)

    cf = ESMF_HConfigCreate(filename='GenScrip.yaml',rc=status)
    _VERIFY(STATUS)

    im_world = ESMF_HConfigAsI4(cf,keyString='CUBE_DIM',rc=status)
    _VERIFY(STATUS)
    JM_WORLD = 6 * IM_WORLD
    gridname = "cube_grid"

    output_scrip = ESMF_HConfigAsString(cf,keyString='output_scrip',rc=status)
    _VERIFY(STATUS)
    output_geos = ESMF_HConfigAsString(cf,keyString='output_geos',rc=status)
    _VERIFY(STATUS)

    do_schmidt=.false.
    if (ESMF_HConfigIsDefined(cf,keyString='DO_SCHMIDT')) then
       do_schmidt = ESMF_HConfigAsLogical(cf,keystring='DO_SCHMIDT',rc=status)
       _VERIFY(status)
    end if
    if (do_schmidt) then
       target_lon = ESMF_HConfigAsR4(cf,keyString='TARGET_LON',rc=status)
       _VERIFY(status)
       target_lat = ESMF_HConfigAsR4(cf,keyString='TARGET_LAT',rc=status)
       _VERIFY(status)
       stretch_factor = ESMF_HConfigAsR4(cf,keyString='STRETCH_FACTOR',rc=status)
       _VERIFY(status)
    end if

    allocate(ims(1,6),jms(1,6))
    do i=1,ntiles
       ims(1,i)=im_world
       jms(1,i)=im_world
    enddo

    if (do_schmidt) then
       transformArgument%stretch_factor = stretch_factor
       transformArgument%target_lon = target_lon * pi/180.0d0
       transformArgument%target_lat = target_lat * pi/180.0d0 
       dstgrid = ESMF_GridCreateCubedSphere(im_world,countsPerDEDim1PTile=ims, &
            countsPerDEDim2PTile=jms, name='bobo',  &
            staggerLocList=[ESMF_STAGGERLOC_CENTER,ESMF_STAGGERLOC_CORNER], &
            transformArgs=transformArgument, &
            coordSys=ESMF_COORDSYS_SPH_RAD,rc=status)
            _VERIFY(status)
    else
       dstgrid = ESMF_GridCreateCubedSphere(im_world,countsPerDEDim1PTile=ims, &
            countsPerDEDim2PTile=jms, name='bobo', &
            staggerLocList=[ESMF_STAGGERLOC_CENTER,ESMF_STAGGERLOC_CORNER], &
            coordSys=ESMF_COORDSYS_SPH_RAD, rc=status)
            _VERIFY(status)
    end if
    call create_gmao_file(dstgrid,im_world,output_geos,rc=status)
    _VERIFY(status)
    
    call ESMF_GridGetCoord(dstgrid, coordDim=1, localDE=0, &
        staggerloc=ESMF_STAGGERLOC_CENTER, &
        farrayPtr=center_lons, rc=status)
    _VERIFY(status)
    call ESMF_GridGetCoord(dstgrid, coordDim=2, localDE=0, &
        staggerloc=ESMF_STAGGERLOC_CENTER, &
        farrayPtr=center_lats, rc=status)
    _VERIFY(status)
    call ESMF_GridGetCoord(dstgrid, coordDim=1, localDE=0, &
        staggerloc=ESMF_STAGGERLOC_CORNER, &
        farrayPtr=corner_lons, rc=status)
    _VERIFY(status)
    call ESMF_GridGetCoord(dstgrid, coordDim=2, localDE=0, &
        staggerloc=ESMF_STAGGERLOC_CORNER, &
        farrayPtr=corner_lats, rc=status)
    _VERIFY(status)
    tmp = im_world*im_world
    allocate(SCRIP_CenterLon(tmp),stat=status)
    _VERIFY(status)
    allocate(SCRIP_CenterLat(tmp),stat=status)
    _VERIFY(status)
    allocate(SCRIP_CornerLat(4,tmp),stat=status)
    _VERIFY(status)
    allocate(SCRIP_CornerLon(4,tmp),stat=status)
    _VERIFY(status)
    allocate(SCRIP_Area(tmp),stat=status)
    _VERIFY(status)
    allocate(SCRIP_rrfac(tmp),stat=status)
    _VERIFY(status)
    allocate(fallback_mask(tmp), stat=status)
    _VERIFY(status)

    local_max_length = 0.0d0
    local_min_length = 1.e15
    fallback_mask = .false.
    failed_cells = 0
    min_area = 1.d30
    max_area = -1.d30
    n = 1

    do jg=1,im_world
       do ig=1,im_world
          i=ig
          j=jg
          mytile=localPet
          SCRIP_CenterLon(n) = modulo(center_lons(i,j)*(180._8/pi), 360.0_8)
          SCRIP_CenterLat(n) = center_lats(i,j)*(180._8/pi)
    
          node_xy(1,1:4) = [corner_lons(i  ,j),corner_lons(i+1,j),corner_lons(i+1,j+1),corner_lons(i,j+1)]
          node_xy(2,1:4) = [corner_lats(i  ,j),corner_lats(i+1,j),corner_lats(i+1,j+1),corner_lats(i,j+1)]
          node_xy_tmp = node_xy
    
          ! Correct for the periodic boundary at 0/360
          lon_w = min( corner_lons(i  ,j),corner_lons(i+1,j),corner_lons(i+1,j+1),corner_lons(i,j+1) )
          lon_e = max( corner_lons(i  ,j),corner_lons(i+1,j),corner_lons(i+1,j+1),corner_lons(i,j+1) )
          if ( abs(lon_e - lon_w) > 1.5_8*pi .and. (SCRIP_CenterLon(n) < pi) ) then
             where(node_xy(1,:) > pi) node_xy_tmp(1,:) = node_xy(1,:) - 2._8*pi
          elseif ( abs(lon_e - lon_w) > 1.5_8*pi .and. (SCRIP_CenterLon(n) > pi) ) then
             where(node_xy(1,:) < pi) node_xy_tmp(1,:) = node_xy(1,:) + 2._8*pi
          endif
    
          !———————————————————————————————————————————————————————————————
          !   Only stretched grids need the convex-hull fix
          !———————————————————————————————————————————————————————————————
          if (do_schmidt) then
             call points_hull_2d(4, node_xy_tmp, hull_num, hull)
    
             if (any(hull == 0)) then
                write(*,100) 'Zero Hull ', corner_lons(i  ,j), corner_lons(i+1,j), &
                                        corner_lons(i+1,j+1), corner_lons(i,j+1)
                write(*,100) 'Zero Hull ', node_xy_tmp(1,:)
             endif
    
             do k = 1, 4
                SCRIP_CornerLon(k,n) = modulo(node_xy_tmp(1,hull(k)) * (180._8/pi), 360.0_8)
                SCRIP_CornerLat(k,n) = node_xy_tmp(2,hull(k)) * (180._8/pi)
             end do
          else
             ! Uniform PE faces keep the original ordering
             do k = 1, 4
                SCRIP_CornerLon(k,n) = node_xy(1,k) * (180._8/pi)
                SCRIP_CornerLat(k,n) = node_xy(2,k) * (180._8/pi)
             end do
          endif
    
          p1 = [corner_lons(i,j),corner_lats(i,j)]
          p2 = [corner_lons(i,j+1),corner_lats(i,j+1)]
          p3 = [corner_lons(i+1,j),corner_lats(i+1,j)]
          p4 = [corner_lons(i+1,j+1),corner_lats(i+1,j+1)]
    
          SCRIP_Area(n) = get_area_spherical_polygon(p1, p2, p3, p4)
    
          !------------------------------------------------------------
          !--- ONLY ADDITION REQUIRED FOR MASK DETECTION STARTS HERE --
          !------------------------------------------------------------
          if (SCRIP_Area(n) <= 0.d0) then
             failed_cells = failed_cells + 1
             fallback_mask(n) = .true.
             ! Print first 10 failures, then every 10-millionth cell afterward
             if (failed_cells <= 10 .or. mod(n,10000000) == 0) then
                write(*,*) 'Negative/zero area at cell ', n, &
                           ' Area:', SCRIP_Area(n), &
                           ' Lon:', SCRIP_CenterLon(n), &
                           ' Lat:', SCRIP_CenterLat(n)
             endif
          endif          
    
          ! Update global min/max area statistics
          min_area = min(min_area, SCRIP_Area(n))
          max_area = max(max_area, SCRIP_Area(n))
          !------------------------------------------------------------
          !--- ONLY ADDITION REQUIRED FOR MASK DETECTION ENDS HERE ----
          !------------------------------------------------------------
    
          call get_grid_length(p1, p2, p3, SCRIP_rrfac(n), local_max_length, local_min_length)
    
          n=n+1
       enddo
    enddo
    
    call MPI_AllReduce(local_max_length, max_length, 1, MPI_DOUBLE, MPI_MAX, mpiC,status)
    call MPI_AllReduce(local_min_length, min_length, 1, MPI_DOUBLE, MPI_MIN, mpiC,status)
    SCRIP_rrfac = max_length/SCRIP_rrfac
    rrfac_max = int(ceiling(max_length/min_length))
    
    !------------------------------------------------------------
    !--- ADDITIONAL REQUIRED DIAGNOSTIC (END OF LOOP) ----------
    !------------------------------------------------------------
    write(*,*) 'SUMMARY: INVALID CELLS:', failed_cells
    write(*,*) 'SUMMARY: AREA RANGE:', min_area, max_area
    !------------------------------------------------------------

 100     format(a,4f20.15)
 101     format(a,f20.15)
 102     format(2f20.15)
 103     format(a)

    deallocate( IMS )
    deallocate( JMS )

    scrip_size = IM_World*JM_World
    call MPI_Info_create(info, status)
    _VERIFY(status)
    call MPI_Info_set(info, "cb_buffer_size", "1048576", status)
    _VERIFY(status)

    status = nf90_create(trim(output_scrip), IOR(NF90_MPIIO,IOR(NF90_CLOBBER,NF90_NETCDF4)), unit, comm=mpiC, info=info)
    _VERIFY(status)

    FMT = '(A,' // 'A,' //'A)'
    write(title,trim(FMT)) 'GMAO ',trim(gridname),' Grid'
    status = nf90_put_att(UNIT, NF90_GLOBAL, 'title',trim(title))
    _VERIFY(status)
    status = nf90_put_att(UNIT, NF90_GLOBAL, 'GridDescriptionFormat', 'SCRIP')
    _VERIFY(status)
    if (do_schmidt) then
       status = nf90_put_att(UNIT, NF90_GLOBAL, 'rrfac_max', rrfac_max)
       _VERIFY(status)
    endif

    status = NF90_DEF_DIM(UNIT, 'grid_size'  , SCRIP_size, gridsize)
    _VERIFY(status)
    status = NF90_DEF_DIM(UNIT, 'grid_corners'  , 4, cornerdim)
    _VERIFY(status)
!! Peggy Li suggested setting grid_rank=1 and grid_dims=1
!! so that ESMF will treat the grid as unstructured
! ------------------------------------------------------
    status = NF90_DEF_DIM(UNIT, 'grid_rank'  , 1, rankdim)
    _VERIFY(status)

!! Grid dimensions
!! ---------------
    status = nf90_def_var(UNIT, "grid_dims", NF90_INT, [rankdim], griddim)
    _VERIFY(status)

!! Grid mask
!! ---------
    status = nf90_def_var(UNIT, "grid_imask", NF90_INT, [gridsize], mask)
    _VERIFY(status)
    status = nf90_put_att(UNIT, mask, "units"    , "unitless")
    _VERIFY(status)

!! Inavalid cells stretched grid mask
!! ----------------------------------
   status = nf90_def_var(UNIT, "grid_fallback_mask", NF90_INT, [gridsize], varid_mask_fallback)
   _VERIFY(status)
   status = nf90_put_att(UNIT, varid_mask_fallback, "description", "1 = fallback area used, 0 = valid")
   _VERIFY(status)
!! cell center Longitude variable
!! ------------------------------
    status = nf90_def_var(UNIT, "grid_center_lon", NF90_DOUBLE, [gridsize], centerlon)
    _VERIFY(status)
    status = nf90_put_att(UNIT, centerlon, "units"    , "degrees")
    _VERIFY(status)

!! cell center Latitude variable
!! -----------------------------
    status = nf90_def_var(UNIT, "grid_center_lat", NF90_DOUBLE,  [gridsize], centerlat)
    _VERIFY(status)
    status = nf90_put_att(UNIT, centerlat, "units"    , "degrees")
    _VERIFY(status)

!! cell corner Longitude variable
!! ------------------------------
    status = nf90_def_var(UNIT, "grid_corner_lon", NF90_DOUBLE, [cornerdim,gridsize], cornerlon)
    _VERIFY(status)
    status = nf90_put_att(UNIT, cornerlon, "units"    , "degrees")
    _VERIFY(status)

!! cell corner Latitude variable
!! -----------------------------
    status = nf90_def_var(UNIT, "grid_corner_lat", NF90_DOUBLE, [cornerdim,gridsize], cornerlat)
    _VERIFY(status)
    status = nf90_put_att(UNIT, cornerlat, "units"    ,"degrees")
    _VERIFY(status)
    status = nf90_def_var(UNIT, "grid_area", NF90_DOUBLE, [gridsize], cellarea)
    _VERIFY(status)
    status = nf90_put_att(UNIT, cellarea, "units"    ,"radians^2")
    _VERIFY(status)
    if (do_schmidt) then
       status = nf90_def_var(UNIT, "rrfac", NF90_DOUBLE, [gridsize], cellrrfac)
       _VERIFY(status)
       status = nf90_put_att(UNIT, cellrrfac, "units"    ,"unitless")
       _VERIFY(status)
    endif

    status = nf90_enddef(UNIT)
    _VERIFY(status)

    rc = NF90_PUT_VAR(UNIT, griddim, [1])

    allocate (sendData(1),GlobalCounts(npets), recvCounts(npets), recvOffsets(npets), stat=status)
    _VERIFY(status)
    sendData = tmp
    recvCounts=1
    recvOffsets=0
    do i=2, npets
      recvOffsets(i) = recvOffsets(i-1) + recvCounts(i-1)
    enddo
    call ESMF_VMGatherV(vm,sendData=sendData,sendCount=1,recvData=GlobalCounts,recvCounts=recvCounts,recvOffsets=recvOffsets,rootPet=0,rc=status)
    _VERIFY(status)
    call ESMF_VMBroadcast(vm,bcstData=GlobalCounts,count=npets,rootPet=0, rc=status)
    _VERIFY(status)

    start=1
    do i=1,localPet
       start(1) = start(1)+GlobalCounts(i)
    enddo

    cnt(1) = tmp; cnt(2)=1
    status = NF90_PUT_VAR(UNIT, centerlon, SCRIP_CenterLon, start, cnt)
    _VERIFY(status)

    status = NF90_PUT_VAR(UNIT, centerlat, SCRIP_CenterLat, start, cnt)
    _VERIFY(status)

    status = NF90_PUT_VAR(UNIT, cellarea, SCRIP_Area,start, cnt)
    _VERIFY(status)

    if (do_schmidt) then
       status = NF90_PUT_VAR(UNIT, cellrrfac, SCRIP_rrfac,start, cnt)
       _VERIFY(status)
    endif
!! Grid mask
!! ---------
    allocate(grid_imask(tmp), stat=status)
    _VERIFY(status)
    grid_imask = 1
    status = NF90_PUT_VAR(UNIT, mask, grid_imask, start, cnt)
    _VERIFY(status)
    deallocate(grid_imask)

    global_start = start(1)    
    start(2)=start(1)
    start(1)=1
    cnt(1)=4
    cnt(2)=tmp

    ! Fallback mask write (1D version, safe)
    allocate(mask_fallback(tmp), stat=status)
    _VERIFY(status)
    mask_fallback = 0
    where (fallback_mask)
       mask_fallback = 1
    end where

    ! Use local start/count for 1D fallback mask
    start_mask(1) = global_start
    cnt_mask(1)   = tmp

    print *, 'DEBUG PET', localPet, 'start_mask=', start_mask(1), 'cnt_mask=', cnt_mask(1) 

    if (cnt_mask(1) > 0 .and. size(mask_fallback) == cnt_mask(1)) then
       status = nf90_put_var(UNIT, varid_mask_fallback, mask_fallback, start_mask, cnt_mask)
       _VERIFY(status)
    else
       print *, 'WARNING: PET', localPet, 'skipping fallback mask write (size mismatch or empty)'
    endif

    status = NF90_PUT_VAR(UNIT, cornerlat, SCRIP_CornerLat, start, cnt)
    _VERIFY(STATUS)

    status = NF90_PUT_VAR(UNIT, cornerlon, SCRIP_CornerLon, start, cnt)
    _VERIFY(STATUS)

    status = NF90_CLOSE(UNIT)
    _VERIFY(STATUS)
    call ESMF_VMBarrier(vm, rc=status)

    deallocate(SCRIP_CenterLat)
    deallocate(SCRIP_CenterLon)
    deallocate(SCRIP_CornerLat)
    deallocate(SCRIP_CornerLon)
    deallocate(SCRIP_Area)
    deallocate(SCRIP_rrfac)
    deallocate(sendData)
    deallocate(GlobalCounts)
    deallocate(recvCounts)
    deallocate(recvOffsets)
    deallocate(fallback_mask)

    call ESMF_Finalize ( rc=status)
    _VERIFY(status)

    contains

subroutine angle_rad_2d ( p1, p2, p3, res )
!*****************************************************************************80
!
!! ANGLE_RAD_2D returns the angle swept out between two rays in 2D.
!
!  Discussion:
!
!    Except for the zero angle case, it should be true that
!
!      ANGLE_RAD_2D ( P1, P2, P3 ) + ANGLE_RAD_2D ( P3, P2, P1 ) = 2 * PI
!
!        P1
!        /
!       /
!      /
!     /
!    P2--------->P3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( REAL64 ) P1(2), P2(2), P3(2), define the rays
!    P1 - P2 and P3 - P2 which define the angle.
!
!    Output, real ( REAL64 ) ANGLE_RAD_2D, the angle swept out by the rays,
!    in radians.  0 <= ANGLE_RAD_2D < 2 * PI.  If either ray has zero
!    length, then ANGLE_RAD_2D is set to 0.
  implicit none

  integer, parameter :: dim_num = 2

  real    (REAL64), parameter :: pi = 3.141592653589793D+00
  real    (REAL64) p(dim_num)
  real    (REAL64) p1(dim_num)
  real    (REAL64) p2(dim_num)
  real    (REAL64) p3(dim_num)
  real    (REAL64) res

  p(1) = ( p3(1) - p2(1) ) * ( p1(1) - p2(1) ) &
       + ( p3(2) - p2(2) ) * ( p1(2) - p2(2) )


  p(2) = ( p3(1) - p2(1) ) * ( p1(2) - p2(2) ) &
       - ( p3(2) - p2(2) ) * ( p1(1) - p2(1) )

  if ( p(1) == 0.0D+00 .and. p(2) == 0.0D+00 ) then
    res = 0.0D+00
    return
  end if

  res = atan2 ( p(2), p(1) )

  if ( res < 0.0D+00 ) then
    res = res + 2.0D+00 * pi
  end if

  return
end

subroutine points_hull_2d ( node_num, node_xy, hull_num, hull )

!*****************************************************************************80
!
!! POINTS_HULL_2D computes the convex hull of 2D points.
!
!  Discussion:
!
!    The work involved is N*log(H), where N is the number of points, and H is
!    the number of points that are on the hull.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, real ( REAL64 ) NODE_XY(2,NODE_NUM), the coordinates of the nodes.
!
!    Output, integer  HULL_NUM, the number of nodes that lie on
!    the convex hull.
!
!    Output, integer  HULL(NODE_NUM).  Entries 1 through HULL_NUM
!    contain the indices of the nodes that form the convex hull, in order.
!
  implicit none

  integer node_num
  real    (REAL64) angle
  real    (REAL64) angle_max
  real    (REAL64) di
  real    (REAL64) dr
  integer  first
  integer  hull(node_num)
  integer  hull_num
  integer  i
  real    (REAL64) node_xy(2,node_num)
  real    (REAL64) p_xy(2)
  integer  q
  real    (REAL64) q_xy(2)
  integer  r
  real    (REAL64) r_xy(2)

  if ( node_num < 1 ) then
    hull_num = 0
    return
  end if
!
!  If NODE_NUM = 1, the hull is the point.
!
  if ( node_num == 1 ) then
    hull_num = 1
    hull(1) = 1
    return
  end if
!
!  If NODE_NUM = 2, then the convex hull is either the two distinct points,
!  or possibly a single (repeated) point.
!
  if ( node_num == 2 ) then

    if ( node_xy(1,1) /= node_xy(1,2) .or. node_xy(2,1) /= node_xy(2,2) ) then
      hull_num = 2
      hull(1) = 1
      hull(2) = 2
    else
      hull_num = 1
      hull(1) = 1
    end if

    return

  end if
!
!  Find the leftmost point and call it "Q".
!  In case of ties, take the bottom-most.
!
  q = 1
  do i = 2, node_num
    if ( node_xy(1,i) < node_xy(1,q) .or. &
       ( node_xy(1,i) == node_xy(1,q) .and. node_xy(2,i) < node_xy(2,q) ) ) then
      q = i
    end if
  end do

  q_xy(1:2) = node_xy(1:2,q)
!
!  Remember the starting point, so we know when to stop!
!
  first = q
  hull_num = 1
  hull(1) = q
!
!  For the first point, make a dummy previous point, 1 unit south,
!  and call it "P".
!
  p_xy(1) = q_xy(1)
  p_xy(2) = q_xy(2) - 1.0D+00
!
!  Now, having old point P, and current point Q, find the new point R
!  so the angle PQR is maximal.
!
!  Watch out for the possibility that the two nodes are identical.
!
  do

    r = 0
    angle_max = 0.0D+00

    do i = 1, node_num

      if ( i /= q .and. &
           ( node_xy(1,i) /= q_xy(1) .or. node_xy(2,i) /= q_xy(2) ) ) then

        call angle_rad_2d(p_xy, q_xy, node_xy(1:2,i),angle)

        if ( r == 0 .or. angle_max < angle ) then

          r = i
          r_xy(1:2) = node_xy(1:2,r)
          angle_max = angle
!
!  In case of ties, choose the nearer point.
!
        else if ( r /= 0 .and. angle == angle_max ) then

          di = ( node_xy(1,i) - q_xy(1) )**2 + ( node_xy(2,i) - q_xy(2) )**2
          dr = ( r_xy(1)      - q_xy(1) )**2 + ( r_xy(2)      - q_xy(2) )**2

          if ( di < dr ) then
            r = i
            r_xy(1:2) = node_xy(1:2,r)
            angle_max = angle
          end if

        end if

      end if

    end do
!
!  We are done when we have returned to the first point on the convex hull.
!
    if ( r == first ) then
      exit
    end if

    hull_num = hull_num + 1
    if ( node_num < hull_num ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POINTS_HULL_2D - Fatal error!'
      write ( *, '(a)' ) '  The algorithm has failed.'
      stop
    end if
!
!  Add point R to convex hull.
!
    hull(hull_num) = r
!
!  Set P := Q, Q := R, and prepare to search for next point R.
!
    q = r

    p_xy(1:2) = q_xy(1:2)
    q_xy(1:2) = r_xy(1:2)

  end do

  return
end subroutine

subroutine create_gmao_file(grid,im_world,filename,rc)
   type(ESMF_Grid), intent(in) :: grid
   integer, intent(in) :: im_world
   character(len=*), intent(in) :: filename
   integer, intent(out), optional :: rc

   integer :: status
   real(ESMF_KIND_R8), pointer :: coords(:,:)

   integer :: ncid,info,lat_id,lon_id,clat_id,clon_id,nf_id,x_id,y_id,rank
   integer :: xp1_id,yp1_id
   real(ESMF_KIND_R8), allocatable :: temp_var(:,:)


   call MPI_Info_create(info, status)
   _VERIFY(status)
   call MPI_Info_set(info, "cb_buffer_size", "1048576", status)
   _VERIFY(status)

   status = nf90_create(filename,NF90_NETCDF4,ncid,comm=MPI_COMM_WORLD,info=info)
   _VERIFY(status)
   
   status = nf90_def_dim(ncid,"nf",6,nf_id)
   _VERIFY(status)
   status = nf90_def_dim(ncid,"Xdim",im_world,x_id)
   _VERIFY(status)
   status = nf90_def_dim(ncid,"Ydim",im_world,y_id)
   _VERIFY(status)
   status = nf90_def_dim(ncid,"XCdim",im_world+1,xp1_id)
   _VERIFY(status)
   status = nf90_def_dim(ncid,"YCdim",im_world+1,yp1_id)
   _VERIFY(status)

   status = nf90_def_var(ncid,"lons",NF90_DOUBLE,[x_id,y_id,nf_id],lon_id)
   _VERIFY(status)
   status = nf90_def_var(ncid,"lats",NF90_DOUBLE,[x_id,y_id,nf_id],lat_id)
   _VERIFY(status)
   status = nf90_def_var(ncid,"corner_lons",NF90_DOUBLE,[xp1_id,yp1_id,nf_id],clon_id)
   _VERIFY(status)
   status = nf90_def_var(ncid,"corner_lats",NF90_DOUBLE,[xp1_id,yp1_id,nf_id],clat_id)
   _VERIFY(status)

   call MPI_COMM_RANK(MPI_COMM_WORLD,rank,status)
   ! centers
   call ESMF_GridGetCoord(grid, coordDim=1, localDE=0, &
       staggerloc=ESMF_STAGGERLOC_CENTER, &
       farrayPtr=coords, rc=status)
   _VERIFY(status)
   allocate(temp_var(im_world,im_world))
   temp_var = coords*180.d0/pi
   status = NF90_put_var(ncid,lon_id,temp_var,start=[1,1,rank+1],count=[im_world,im_world,1])
   _VERIFY(status)
   call ESMF_GridGetCoord(grid, coordDim=2, localDE=0, &
       staggerloc=ESMF_STAGGERLOC_CENTER, &
       farrayPtr=coords, rc=status)
   _VERIFY(status)
   temp_var = coords*180.d0/pi
   status = NF90_put_var(ncid,lat_id,temp_var,start=[1,1,rank+1],count=[im_world,im_world,1])
   _VERIFY(status)
   deallocate(temp_var)
   ! corners
   call ESMF_GridGetCoord(grid, coordDim=1, localDE=0, &
       staggerloc=ESMF_STAGGERLOC_CORNER, &
       farrayPtr=coords, rc=status)
   _VERIFY(status)
   allocate(temp_var(im_world+1,im_world+1))
   temp_var = coords*180.d0/pi
   status = NF90_put_var(ncid,clon_id,temp_var,start=[1,1,rank+1],count=[im_world+1,im_world+1,1])
   _VERIFY(status)
   call ESMF_GridGetCoord(grid, coordDim=2, localDE=0, &
       staggerloc=ESMF_STAGGERLOC_CORNER, &
       farrayPtr=coords, rc=status)
   _VERIFY(status)
   temp_var = coords*180.d0/pi
   status = NF90_put_var(ncid,clat_id,temp_var,start=[1,1,rank+1],count=[im_world+1,im_world+1,1])
   _VERIFY(status)
   deallocate(temp_var)
   if (present(rc)) rc=0

   end subroutine   


!------------------------------------------------------------
! Girard's Formula computes the area of a polygon on a sphere by 
! calculating the sum of the polygon's corner angles (interior angles) 
! and subtracting the expected sum of angles for a planar polygon. 
! On a sphere, polygons have slightly larger angle sums compared 
! to their planar counterparts, and this difference known as 
! spherical excess directly gives the polygon's area.
! stretched grid has very small areas and need a safeguard. 
!------------------------------------------------------------
 function get_area_spherical_polygon(p1,p4,p2,p3) result(area)
    real(REAL64), parameter           :: PI = 3.14159265358979323846
    real(real64) :: area
    real(real64), intent(in) :: p1(2),p2(2),p3(2),p4(2)

    real(real64) :: e1(3),e2(3),e3(3)
    real(real64) :: ang1,ang2,ang3,ang4


    e1 = convert_to_cart(p1)
    e2 = convert_to_cart(p2)
    e3 = convert_to_cart(p4)
    ang1 = spherical_angles(e1, e2, e3)

    e1 = convert_to_cart(p2)
    e2 = convert_to_cart(p3)
    e3 = convert_to_cart(p1)
    ang2 = spherical_angles(e1, e2, e3)

    e1 = convert_to_cart(p3)
    e2 = convert_to_cart(p4)
    e3 = convert_to_cart(p2)
    ang3 = spherical_angles(e1, e2, e3)

    e1 = convert_to_cart(p4)
    e2 = convert_to_cart(p3)
    e3 = convert_to_cart(p1)
    ang4 = spherical_angles(e1, e2, e3)

    area = ang1 + ang2 + ang3 + ang4 - 2.0d0*PI

    ! stretched grid safeguard
   if (area <= 0.d0 .and. area > -1e-10) then
      area = abs(area) ! minor numerical fix
   endif

 end function get_area_spherical_polygon

 subroutine get_grid_length(p1,p2,p3, local,max_length, min_length)
   real(real64), intent(in) :: p1(2),p2(2),p3(2)
   real(REAL64), intent(out) :: local
   real(REAL64), intent(inout) :: max_length
   real(REAL64), intent(inout) :: min_length
   real(REAL64) :: dx, dy

   dx = great_circle_dist(p1,p2)
   dy = great_circle_dist(p2,p3)
   local = 0.5d0*(dx+dy)
   max_length = max(local, max_length)
   min_length = min(local, min_length)
 end subroutine

 function convert_to_cart(v) result(xyz)
    real(real64), intent(in) :: v(2)
    real(real64) :: xyz(3)

    xyz(1)=cos(v(2))*cos(v(1))
    xyz(2)=cos(v(2))*sin(v(1))
    xyz(3)=sin(v(2))

 end function convert_to_cart

function spherical_angles(p1,p2,p3) result(spherical_angle)
   real(real64) :: spherical_angle
   real(real64), intent(in) :: p1(3),p2(3),p3(3)

   real (real64):: e1(3), e2(3), e3(3)
   real (real64):: px, py, pz
   real (real64):: qx, qy, qz
   real (real64):: angle, ddd
   integer n
   real(REAL64), parameter           :: PI = 3.14159265358979323846

   do n=1,3
      e1(n) = p1(n)
      e2(n) = p2(n)
      e3(n) = p3(n)
   enddo

   !-------------------------------------------------------------------
   ! Page 41, Silverman's book on Vector Algebra; spherical trigonmetry
   !-------------------------------------------------------------------
   ! Vector P:
   px = e1(2)*e2(3) - e1(3)*e2(2)
   py = e1(3)*e2(1) - e1(1)*e2(3)
   pz = e1(1)*e2(2) - e1(2)*e2(1)
   ! Vector Q:
   qx = e1(2)*e3(3) - e1(3)*e3(2)
   qy = e1(3)*e3(1) - e1(1)*e3(3)
   qz = e1(1)*e3(2) - e1(2)*e3(1)

   ddd = (px*px+py*py+pz*pz)*(qx*qx+qy*qy+qz*qz)

   if ( ddd <= 0.0d0 ) then
      angle = 0.d0
   else
      ddd = (px*qx+py*qy+pz*qz) / sqrt(ddd)
      ! Added numerical safeguard for acos domain
      if (ddd >  1.0d0) ddd =  1.0d0
      if (ddd < -1.0d0) ddd = -1.0d0
      if ( abs(ddd)>1.d0) then
         angle = 0.5d0 * PI
         if (ddd < 0.d0) then
            angle = PI
         else
            angle = 0.d0
         end if
      else
         angle = acos( ddd )
      endif
   endif

   spherical_angle = angle
end function

subroutine local_abort(rc,line_number)
   integer, intent(in) :: rc
   integer, intent(in) :: line_number
   integer :: status, error_code
   write(*,*)"aborted on with stat ",line_number,rc
   call MPI_Abort(MPI_COMM_WORLD,error_code,status)
end subroutine

real(REAL64) function great_circle_dist( q1, q2, radius )
     real(REAL64), intent(IN)           :: q1(2), q2(2)
     real(REAL64), intent(IN), optional :: radius

     real (REAL64):: p1(2), p2(2)
     real (REAL64):: beta
     integer n

     do n=1,2
        p1(n) = q1(n)
        p2(n) = q2(n)
     enddo

     beta = asin( sqrt( sin((p1(2)-p2(2))/2.)**2 + cos(p1(2))*cos(p2(2))*   &
                        sin((p1(1)-p2(1))/2.)**2 ) ) * 2.

     if ( present(radius) ) then
          great_circle_dist = radius * beta
     else
          great_circle_dist = beta   ! Returns the angle
     endif

end function great_circle_dist

    end program ESMF_GenerateCSGridDescription
