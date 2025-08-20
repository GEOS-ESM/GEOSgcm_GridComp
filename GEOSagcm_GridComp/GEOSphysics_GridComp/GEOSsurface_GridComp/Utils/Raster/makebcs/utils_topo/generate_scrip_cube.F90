!#define _VERIFY(A) if(A/=0) write(*,*)__LINE__;call MPI_ABort(MPI_COMM_WORLD,error_code,status)
#define _VERIFY(A) if(A/=0) call local_abort(A,__LINE__)

    program ESMF_GenerateCSGridDescription

! ESMF Framework module
    use ESMF
    use mpi
    use netcdf
    use, intrinsic :: iso_fortran_env, only: REAL64,REAL32

    implicit none
    
! The code creating and filling these variables is not included in the
! example documentation because those interfaces are not specific to
! Regrid.
    type(ESMF_Grid) ::  dstgrid
    type(ESMF_Grid) ::  tempgrid
    type(ESMF_Grid) ::  unigrid
    type(ESMF_CubedSphereTransform_Args) :: transformArgument

!EOC

    real(REAL64), parameter           :: PI = 3.14159265358979323846
    integer                           :: npets, localPet
    integer                           :: i, j, k
    integer :: rc
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
    real(ESMF_KIND_R8), allocatable   :: SCRIP_CenterLat(:), SCRIP_CenterLon(:)
    real(ESMF_KIND_R8), allocatable   :: SCRIP_CornerLat(:,:), SCRIP_CornerLon(:,:)
    real(ESMF_KIND_R8), allocatable   :: SCRIP_Area(:)
    real(ESMF_KIND_R8), allocatable   :: SCRIP_rrfac(:)
    real(ESMF_KIND_R8)                :: node_xy(2,4), node_xy_tmp(2,4), lon_e, lon_w
    integer                           :: gridsize, griddim, rankdim, mask
    integer                           :: cornerlon, cornerlat, centerlon, centerlat,cellarea, cellrrfac
    integer, allocatable              :: IMS(:,:), JMS(:,:), sendData(:), GlobalCounts(:), recvCounts(:), recvOffsets(:)
    integer, allocatable              :: grid_imask(:)
    character(len=ESMF_MAXSTR)        :: gridname, FMT, title
    integer                           :: myTile
    integer                           :: tmp, mpiC
    integer                           :: IG, JG, rrfac_max
    logical                           :: do_schmidt
    logical, allocatable              :: fallback_mask(:)
    integer                           :: start_mask(1), cnt_mask(1)
    integer                           :: varid_mask_fallback  
    integer, allocatable              :: mask_fallback(:)     
    real(ESMF_KIND_R8)                :: p1(2),p2(2),p3(2),p4(2)
    real(ESMF_KIND_R8)                :: local_max_length, max_length, local_min_length, min_length
    real(ESMF_KIND_R4)                :: target_lon, target_lat, stretch_factor
    type(ESMF_HConfig)                :: CF
    integer                           :: status
    character(len=ESMF_MAXPATHLEN)    :: output_scrip, output_geos
    integer                           :: failed_cells
    real(ESMF_KIND_R8)                :: max_area
    real(ESMF_KIND_R8), allocatable   :: local_max_length_all(:), local_min_length_all(:)
    integer, allocatable              :: localDEList(:)
    integer                           :: cell, num_cells
    integer                           :: localDECount
    integer                           :: mpi_err
    integer                           :: de
    integer                           :: global_chosen_pet, chosen_pet, chosen_de
    integer                           :: chosen_de_stretch, chosen_i_stretch, chosen_j_stretch, chosen_pet_stretch
    real(ESMF_KIND_R8), pointer       :: tmp_center_lons(:,:) => null()
    real(ESMF_KIND_R8), pointer       :: tmp_center_lats(:,:) => null()
    real(ESMF_KIND_R8), pointer       :: tmp_corner_lons(:,:) => null()
    real(ESMF_KIND_R8), pointer       :: tmp_corner_lats(:,:) => null()
    real(ESMF_KIND_R8), pointer       :: uni_corner_lons(:,:) => null()
    real(ESMF_KIND_R8), pointer       :: uni_corner_lats(:,:) => null()
    real(ESMF_KIND_R8)                :: tiny_dlon, tiny_dlat
    real(ESMF_KIND_R8)                :: clon, clat, area_signed
    real(ESMF_KIND_R8)                :: global_max_length, global_min_length, min_allowed_length
    real(ESMF_KIND_R8)                :: max_rrfac
    real(ESMF_KIND_R8)                :: midpoint_lon, midpoint_lat
    real(ESMF_KIND_R8)                :: dist, min_dist_local, min_dist_global, lon_diff, lat_diff
    real(ESMF_KIND_R8)                :: chosen_center_lon, chosen_center_lat
    real(ESMF_KIND_R8)                :: global_chosen_lon, global_chosen_lat
    real(ESMF_KIND_R8), parameter     :: min_length_threshold = 1.0d-12
    integer           , parameter     :: grid_corners = 4
    logical                           :: panel_found 
    integer                           :: num_rrfac_max
    integer                           :: j_offset, jc_offset
    integer                           :: num_local_cells, n_start
    integer                           :: kk
    integer                           :: best_idx
    integer                           :: num_cells_global
    integer                           :: n_end
    integer                           :: i_local, j_local
    integer                           :: num_rrfac_max_local
    integer                           :: local_idx_max_rrfac, global_idx_max_rrfac
    integer                           :: valid_count_local
    integer                           :: n_midpoint
    integer                           :: cornerdimID
    integer                           :: local_idx_in_slice
    integer                           :: clamp_count_local, clamp_count_global
    integer                           :: fallback_count_local, fallback_count_global
    integer                           :: owner_pet
    integer                           :: start_file, cnt_file, mem_start, mem_end
    integer                           :: j_offset_chosen
    integer                           :: pet
    real(ESMF_KIND_R8)                :: local_max_rrfac, global_max_rrfac
    real(ESMF_KIND_R8)                :: global_lon_max_rrfac, global_lat_max_rrfac
    real(ESMF_KIND_R8)                :: eps
    logical                           :: found_degenerate
    logical                           :: bad_corner
    logical                           :: has_midpoint
    real(ESMF_KIND_R8)                :: epsilon
    integer                           :: l
    real(ESMF_KIND_R8)                :: global_max_area, global_min_area,ratio
    real(ESMF_KIND_R8), allocatable   :: my_corner_lat(:,:), my_corner_lon(:,:)
    real(ESMF_KIND_R8), allocatable   :: A_uniform(:)
    real(ESMF_KIND_R8)                :: lon_tmp(4), lat_tmp(4)
    real(ESMF_KIND_R8)                :: max_rrfac_allowed, rrfac_tmp
    real(ESMF_KIND_R8)                :: midpoint_lon_deg, midpoint_lat_deg
    real(ESMF_KIND_R8)                :: max_area_local, min_area_local
    real(ESMF_KIND_R8), dimension(2)  :: local_val_idx, global_val_idx    
    real(ESMF_KIND_R8)                :: local_min_rrfac, global_min_rrfac_print, global_max_rrfac_print
    real(ESMF_KIND_R8)                :: owner_lon, owner_lat
    real(ESMF_KIND_R8)                :: u1(2), u2(2), u3(2), u4(2)
    real(ESMF_KIND_R8)                :: swap_p(2)
    real(ESMF_KIND_R8)                :: deg2rad, lonc, latc, half_power_radius_deg, theta0
    real(ESMF_KIND_R8)                :: desired_peak, dlon, lat1, lat2, dlam, cc, dtheta, w
    real(ESMF_KIND_R8)                :: best_dist
    real(ESMF_KIND_R8)                :: pair_local(2), pair_global(2)
    real(ESMF_KIND_R8)                :: tmp_len
    real(ESMF_KIND_R8)                :: dummy_max, dummy_min
    real(ESMF_KIND_R8)                :: max_len_local, min_len_local
    character(len=3)                  :: orient


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
       if (target_lon < 0.0) target_lon = target_lon + 360.0
    end if

    allocate(ims(1,6),jms(1,6))
    do i=1,ntiles
       ims(1,i)=im_world
       jms(1,i)=im_world
    enddo

    !======================== GRID CREATION & SETUP ========================
    
    ! Build grid(s) based on do_schmidt
    if (.not. do_schmidt) then
      !------------------------------ REGULAR GRID (legacy) ------------------------------
      dstgrid = ESMF_GridCreateCubedSphere(im_world, countsPerDEDim1PTile=ims, &
           countsPerDEDim2PTile=jms, name='bobo',                           &
           staggerLocList=[ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER], &
           coordSys=ESMF_COORDSYS_SPH_RAD, rc=status); _VERIFY(status)
    
      call ESMF_GridGet(dstgrid, localDECount=localDECount, rc=status); _VERIFY(status)
    
    else
      !------------------------------ STRETCHED GRID (Schmidt) ---------------------------
    
      ! 1) Pick the midpoint (closest cell center to requested target) on an
      !    unstretched grid in degrees, so the search is simple.
      tempgrid = ESMF_GridCreateCubedSphere(im_world, countsPerDEDim1PTile=ims, &
           countsPerDEDim2PTile=jms, name='temp',                               &
           staggerLocList=[ESMF_STAGGERLOC_CENTER],                             &
           coordSys=ESMF_COORDSYS_SPH_DEG, rc=status); _VERIFY(status)
    
      call ESMF_GridGet(tempgrid, localDECount=localDECount, rc=status); _VERIFY(status)
    
      allocate(localDEList(localDECount))
      do de = 0, localDECount - 1
         localDEList(de+1) = de
      enddo
    
      min_dist_local = huge(1.0d0)
      do de = 1, localDECount
         call ESMF_GridGetCoord(tempgrid, coordDim=1, localDE=localDEList(de), &
              staggerloc=ESMF_STAGGERLOC_CENTER, farrayPtr=tmp_center_lons, rc=status); _VERIFY(status)
         call ESMF_GridGetCoord(tempgrid, coordDim=2, localDE=localDEList(de), &
              staggerloc=ESMF_STAGGERLOC_CENTER, farrayPtr=tmp_center_lats, rc=status); _VERIFY(status)
    
         do j = 1, size(tmp_center_lats,2)
           do i = 1, size(tmp_center_lons,1)
             lon_diff = modulo(tmp_center_lons(i,j) - target_lon + 180.0d0, 360.0d0) - 180.0d0
             lat_diff = tmp_center_lats(i,j) - target_lat
             dist = sqrt(lon_diff**2 + lat_diff**2)
    
             if (dist < min_dist_local) then
               min_dist_local    = dist
               chosen_center_lon = tmp_center_lons(i,j)
               chosen_center_lat = tmp_center_lats(i,j)
               chosen_pet        = localPet
               chosen_de         = localDEList(de)
             endif
           enddo
         enddo
      enddo
    
      ! 2) Find global min distance and broadcast the midpoint in degrees
      call MPI_Allreduce(min_dist_local, min_dist_global, 1, MPI_DOUBLE_PRECISION, MPI_MIN, mpiC, mpi_err)
    
      if (abs(min_dist_local - min_dist_global) < 1.0d-10) then
         global_chosen_pet = localPet
      else
         global_chosen_pet = -1
      endif
      call MPI_Allreduce(MPI_IN_PLACE, global_chosen_pet, 1, MPI_INTEGER, MPI_MAX, mpiC, mpi_err)
    
      global_chosen_lon = chosen_center_lon
      global_chosen_lat = chosen_center_lat
      call MPI_Bcast(global_chosen_lon, 1, MPI_DOUBLE_PRECISION, global_chosen_pet, mpiC, mpi_err)
      call MPI_Bcast(global_chosen_lat, 1, MPI_DOUBLE_PRECISION, global_chosen_pet, mpiC, mpi_err)
    
      midpoint_lon = global_chosen_lon
      midpoint_lat = global_chosen_lat
    
      ! 3) Create the actual stretched grid at the chosen midpoint (radians)
      transformArgument%stretch_factor = stretch_factor
      transformArgument%target_lon     = midpoint_lon * pi/180.0d0
      transformArgument%target_lat     = midpoint_lat * pi/180.0d0
    
      dstgrid = ESMF_GridCreateCubedSphere(im_world, countsPerDEDim1PTile=ims,     &
           countsPerDEDim2PTile=jms, name='bobo',                                   &
           staggerLocList=[ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER],         &
           transformArgs=transformArgument,                                         &
           coordSys=ESMF_COORDSYS_SPH_RAD, rc=status); _VERIFY(status)
    
      call ESMF_GridGet(dstgrid, localDECount=localDECount, rc=status); _VERIFY(status)
    
      ! 4) uniform “reference” corners used for diagnostics/comparisons
      unigrid = ESMF_GridCreateCubedSphere(im_world, countsPerDEDim1PTile=ims, &
           countsPerDEDim2PTile=jms, name='uniform',                           &
           staggerLocList=[ESMF_STAGGERLOC_CORNER],                            &
           coordSys=ESMF_COORDSYS_SPH_RAD, rc=status); _VERIFY(status)
    
      call ESMF_GridGetCoord(unigrid, coordDim=1, localDE=0, &
           staggerloc=ESMF_STAGGERLOC_CORNER, farrayPtr=uni_corner_lons, rc=status); _VERIFY(status)
      call ESMF_GridGetCoord(unigrid, coordDim=2, localDE=0, &
           staggerloc=ESMF_STAGGERLOC_CORNER, farrayPtr=uni_corner_lats, rc=status); _VERIFY(status)
    
      ! 5) Find the exact (PET,DE,i,j) of the midpoint on the *stretched* grid
      min_dist_local = huge(1.0d0)
      do de = 0, localDECount - 1
         call ESMF_GridGetCoord(dstgrid, coordDim=1, localDE=de, &
              staggerloc=ESMF_STAGGERLOC_CENTER, farrayPtr=tmp_center_lons, rc=status); _VERIFY(status)
         call ESMF_GridGetCoord(dstgrid, coordDim=2, localDE=de, &
              staggerloc=ESMF_STAGGERLOC_CENTER, farrayPtr=tmp_center_lats, rc=status); _VERIFY(status)
    
         do j = 1, size(tmp_center_lats,2)
           do i = 1, size(tmp_center_lons,1)
             lon_diff = modulo(tmp_center_lons(i,j)*180.d0/pi - midpoint_lon + 180.0d0, 360.0d0) - 180.0d0
             lat_diff = tmp_center_lats(i,j)*180.d0/pi - midpoint_lat
             dist = sqrt(lon_diff**2 + lat_diff**2)
             if (dist < min_dist_local) then
               min_dist_local    = dist
               chosen_de_stretch = de
               chosen_i_stretch  = i
               chosen_j_stretch  = j
             endif
           enddo
         enddo
      enddo
    
      call MPI_Allreduce(min_dist_local, min_dist_global, 1, MPI_DOUBLE_PRECISION, MPI_MIN, mpiC, mpi_err)
      if (abs(min_dist_local - min_dist_global) < 1.0d-10) then
         chosen_pet_stretch = localPet
      else
         chosen_pet_stretch = -1
      endif
      call MPI_Allreduce(MPI_IN_PLACE, chosen_pet_stretch, 1, MPI_INTEGER, MPI_MAX, mpiC, mpi_err)
    
      call MPI_Bcast(chosen_pet_stretch, 1, MPI_INTEGER, chosen_pet_stretch, mpiC, mpi_err)
      call MPI_Bcast(chosen_de_stretch,  1, MPI_INTEGER, chosen_pet_stretch, mpiC, mpi_err)
      call MPI_Bcast(chosen_i_stretch,   1, MPI_INTEGER, chosen_pet_stretch, mpiC, mpi_err)
      call MPI_Bcast(chosen_j_stretch,   1, MPI_INTEGER, chosen_pet_stretch, mpiC, mpi_err)
    
      midpoint_lon_deg = -999.d0
      midpoint_lat_deg = -999.d0
      if (localPet == chosen_pet_stretch) then
         call ESMF_GridGetCoord(dstgrid, coordDim=1, localDE=chosen_de_stretch, &
              staggerloc=ESMF_STAGGERLOC_CENTER, farrayPtr=tmp_center_lons, rc=status); _VERIFY(status)
         call ESMF_GridGetCoord(dstgrid, coordDim=2, localDE=chosen_de_stretch, &
              staggerloc=ESMF_STAGGERLOC_CENTER, farrayPtr=tmp_center_lats, rc=status); _VERIFY(status)
         midpoint_lon_deg = modulo(tmp_center_lons(chosen_i_stretch,chosen_j_stretch)*180.d0/pi, 360.d0)
         midpoint_lat_deg =        tmp_center_lats(chosen_i_stretch,chosen_j_stretch)*180.d0/pi
      endif
      call MPI_Bcast(midpoint_lon_deg, 1, MPI_DOUBLE_PRECISION, chosen_pet_stretch, mpiC, mpi_err)
      call MPI_Bcast(midpoint_lat_deg, 1, MPI_DOUBLE_PRECISION, chosen_pet_stretch, mpiC, mpi_err)
    
    end if
    !-------------------------- end split regular / stretched --------------------------
    
    ! 6) Pointers to THIS PET’s face on dstgrid (used by both paths later)
    call ESMF_GridGetCoord(dstgrid, coordDim=1, localDE=0, &
         staggerloc=ESMF_STAGGERLOC_CENTER, farrayPtr=tmp_center_lons, rc=status); _VERIFY(status)
    call ESMF_GridGetCoord(dstgrid, coordDim=2, localDE=0, &
         staggerloc=ESMF_STAGGERLOC_CENTER, farrayPtr=tmp_center_lats, rc=status); _VERIFY(status)
    
    call ESMF_GridGetCoord(dstgrid, coordDim=1, localDE=0, &
         staggerloc=ESMF_STAGGERLOC_CORNER, farrayPtr=tmp_corner_lons, rc=status); _VERIFY(status)
    call ESMF_GridGetCoord(dstgrid, coordDim=2, localDE=0, &
         staggerloc=ESMF_STAGGERLOC_CORNER, farrayPtr=tmp_corner_lats, rc=status); _VERIFY(status)
    
    ! 7) Indexing helpers (REGULAR keeps legacy local indexing; STRETCH uses global)
    jm_world         = im_world * npets
    num_local_cells  = im_world * im_world
    num_cells_global = im_world * jm_world
    
    if (.not. do_schmidt) then
      n_start = 1
      n_end   = num_local_cells
    else
      n_start = localPet * num_local_cells + 1
      n_end   = n_start + num_local_cells - 1
    endif
    
    !====================== END GRID CREATION & SETUP ======================
    ! --- Optional: quick sanity print from locals (in degrees) ---
    if (localPet == 0) then
      write(*,*) "[LOCAL FACE] center(1,1) lon/lat (deg):", &
           modulo(tmp_center_lons(1,1)*180.d0/pi,360.d0), tmp_center_lats(1,1)*180.d0/pi
      write(*,*) "[LOCAL FACE] corner(1,1) lon/lat (deg):", &
           modulo(tmp_corner_lons(1,1)*180.d0/pi,360.d0), tmp_corner_lats(1,1)*180.d0/pi
    end if
    
    ! wrap local tile corners to [0,2π) – applies to both modes
    do j = 1, size(tmp_corner_lons,2)
      do i = 1, size(tmp_corner_lons,1)
        if (tmp_corner_lons(i,j) >= 2.d0*pi) tmp_corner_lons(i,j) = tmp_corner_lons(i,j) - 2.d0*pi
        if (tmp_corner_lons(i,j) <  0.d0)    tmp_corner_lons(i,j) = tmp_corner_lons(i,j) + 2.d0*pi
      end do
    end do
    write(*,*) "[POST PERIODICITY] PET:", localPet, &
         " corner_lons(1,1):", tmp_corner_lons(1,1)*180.d0/pi, &
         " corner_lats(1,1):", tmp_corner_lats(1,1)*180.d0/pi
    
    ! Extra midpoint diagnostics ONLY for stretched mode
    if (do_schmidt) then
      j_offset_chosen = chosen_pet_stretch * im_world   
      midpoint_lon_deg = -999.d0
      midpoint_lat_deg = -999.d0
      if (localPet == chosen_pet_stretch) then
        midpoint_lon_deg = modulo(tmp_center_lons(chosen_i_stretch, chosen_j_stretch)*180.d0/pi, 360.d0)
        midpoint_lat_deg =        tmp_center_lats(chosen_i_stretch, chosen_j_stretch)*180.d0/pi
      endif
      call MPI_Bcast(midpoint_lon_deg, 1, MPI_DOUBLE_PRECISION, chosen_pet_stretch, MPI_COMM_WORLD, mpi_err)
      call MPI_Bcast(midpoint_lat_deg, 1, MPI_DOUBLE_PRECISION, chosen_pet_stretch, MPI_COMM_WORLD, mpi_err)
      if (localPet == chosen_pet_stretch) write(*,*) "[MIDPOINT CONFIRMED] PET:", localPet
    end if
    
    ! ---- Allocate SCRIP arrays (size depends on mode) ----
    if (.not. do_schmidt) then
      tmp = im_world * im_world          ! local face (legacy behavior)
    else
      tmp = num_cells_global             ! global (stretched workflow)
    endif
    
    allocate(SCRIP_CenterLon(tmp),stat=status); _VERIFY(status)
    allocate(SCRIP_CenterLat(tmp),stat=status); _VERIFY(status)
    allocate(SCRIP_CornerLat(4,tmp),stat=status); _VERIFY(status)
    allocate(SCRIP_CornerLon(4,tmp),stat=status); _VERIFY(status)
    allocate(SCRIP_Area(tmp),stat=status); _VERIFY(status)
    allocate(SCRIP_rrfac(tmp),stat=status); _VERIFY(status)
    SCRIP_rrfac = 0.0d0
    allocate(fallback_mask(tmp),stat=status); _VERIFY(status)
    
    if (do_schmidt) then
      allocate(A_uniform(tmp), stat=status); _VERIFY(status)
    end if
    
    num_cells = tmp
    allocate(local_max_length_all(num_cells))
    allocate(local_min_length_all(num_cells))
    local_max_length_all(:) = 0.0d0
    local_min_length_all(:) = 1.0d15
    fallback_mask = .false.
    failed_cells  = 0
    
    mytile   = localPet + 1
    j_offset = (mytile-1)*im_world
    jc_offset= (mytile-1)*im_world
    
    write(*,'("[CHECK PET SLICE] PET:",I2," j_offset:",I8," jc_offset:",I8)') localPet,j_offset,jc_offset
    
    if (localPet == 0) then
      write(*,*) "PET 0 local corner (1,1) lon/lat:", &
                 tmp_corner_lons(1,1)*180.d0/pi, tmp_corner_lats(1,1)*180.d0/pi
    end if
    if (localPet == npets-1) then
      write(*,*) "PET last local corner (end,end) lon/lat:", &
                 tmp_corner_lons(im_world+1,im_world+1)*180.d0/pi, &
                 tmp_corner_lats(im_world+1,im_world+1)*180.d0/pi
    end if
    
    num_local_cells = im_world*im_world
    
    ! *** CRITICAL: indexing differs by mode ***
    if (.not. do_schmidt) then
      n_start = 1
      n_end   = num_local_cells
    else
      n_start = num_local_cells * localPet + 1
      n_end   = n_start + num_local_cells - 1
    endif
    n = n_start
    
    write(*,*) "[PRE MAIN LOOP] PET:", localPet, " j_offset:", j_offset, " jc_offset:", jc_offset, &
     " corner_lons(1,1):", tmp_corner_lons(1,1)*180.d0/pi, &
     " corner_lats(1,1):", tmp_corner_lats(1,1)*180.d0/pi
    
    do j = 1, im_world
      do i = 1, im_world
        ! centers (local tile → degrees)
        SCRIP_CenterLon(n) = modulo(tmp_center_lons(i,j)*(180.d0/pi), 360.d0)
        SCRIP_CenterLat(n) =        tmp_center_lats(i,j)*(180.d0/pi)
    
        ! corners (local tile → radians here, you wrap before area or convert when writing)
        node_xy(1,:) = [ tmp_corner_lons(i  ,j  ), tmp_corner_lons(i+1,j  ), &
                         tmp_corner_lons(i+1,j+1), tmp_corner_lons(i  ,j+1) ]
        node_xy(2,:) = [ tmp_corner_lats(i  ,j  ), tmp_corner_lats(i+1,j  ), &
                         tmp_corner_lats(i+1,j+1), tmp_corner_lats(i  ,j+1) ]
    
        node_xy_tmp = node_xy
    
        ! Correct for longitude periodicity
        lon_w = minval(node_xy(1,:))
        lon_e = maxval(node_xy(1,:))
        if (abs(lon_e - lon_w) > 1.5_8*pi) then
          if (tmp_center_lons(i,j) < pi) then
            where (node_xy(1,:) > pi) node_xy_tmp(1,:) = node_xy(1,:) - 2._8*pi
          else
            where (node_xy(1,:) < pi) node_xy_tmp(1,:) = node_xy(1,:) + 2._8*pi
          endif
        endif
    
        node_xy = node_xy_tmp
    
        ! === Uniform-grid reference area: STRETCH ONLY ===
        if (do_schmidt) then
          u1 = [uni_corner_lons(i  ,j  ), uni_corner_lats(i  ,j  )]
          u2 = [uni_corner_lons(i+1,j  ), uni_corner_lats(i+1,j  )]
          u3 = [uni_corner_lons(i+1,j+1), uni_corner_lats(i+1,j+1)]
          u4 = [uni_corner_lons(i  ,j+1), uni_corner_lats(i  ,j+1)]
          A_uniform(n) = get_area_spherical_polygon(u1,u2,u3,u4)
        end if
    
        if (n == n_start) then
          write(*,*) "[FIRST CELL CORNERS] PET:", localPet, " n:", n, &
                     " corners_lon(deg):", node_xy(1,:)*180/pi, &
                     " corners_lat(deg):", node_xy(2,:)*180/pi
        endif

         !------------------------------------------------------------------------------------
         !  STRETCHED-ONLY geometry repair (WHY & ORDER)
         !
         !  WHY:
         !    Schmidt stretching can distort panel quads near the target and panel edges:
         !    - corners can wrap across 0/360, 
         !    - quads can become non-convex or nearly collinear,
         !    - numerical noise can yield NaNs/Inf or negative/tiny areas,
         !    - corner ordering can be inconsistent.
         !    SCRIP expects a simple quad with CCW ordering and positive area.
         !
         !  WHAT & ORDER (progressively stronger interventions):
         !    1) Degeneracy check: if any pair of corners nearly coincides, build a tiny
         !       CCW square around the cell center (fallback) and flag fallback_mask.
         !    2) Clamp corners: wrap lon to [0,2π), clamp lon away from exact 0/2π, clamp
         !       lat to [-π/2, π/2] to avoid out-of-range/NaN geometry.
         !    3) Convex hull: reorder the (possibly scrambled) corners into a proper quad.
         !       If hull fails, jitter a tiny ε box around the center and retry.
         !       If it still fails, build the explicit tiny CCW square fallback.
         !    4) Finite check: if any corner is NaN/Inf, rebuild the tiny CCW square fallback.
         !    5) Consistent ordering: call safe_reorder_hull → (p1,p2,p3,p4).
         !    6) Area guard: compute signed spherical area; if NaN/Inf/|area| tiny/huge,
         !       rebuild tiny CCW square and recompute; mark fallback_mask.
         !    7) CCW enforcement: if signed area < 0, swap p2↔p4 to make it CCW so area>0.
         !    8) Write corners/area: store in SCRIP arrays (degrees; lon wrapped to 0..360).
         !
         !  NOTE:
         !    - All of the above is under do_schmidt=.true.; the regular grid path preserves
         !      legacy behavior, with only CCW enforcement to guarantee positive area.
         !    - fallback_mask marks cells where we had to invent a tiny polygon; later logic
         !      (e.g., rrfac) can treat these conservatively.
         !------------------------------------------------------------------------------------         
         if (do_schmidt) then
           ! 1) Detect truly degenerate quads and do tiny‐polygon fallback
           found_degenerate = .false.
           epsilon = 5.0d-6        ! ≈ 1″; stops “everything looks degenerate”
           do k = 1,4
             do l = k+1,4
               if (abs(node_xy_tmp(1,k)-node_xy_tmp(1,l)) < epsilon  .and.  &
                   abs(node_xy_tmp(2,k)-node_xy_tmp(2,l)) < epsilon) then
                 found_degenerate = .true.
               end if
             end do
           end do
         
           if (found_degenerate) then
             ! build a tiny fallback polygon around the center
             tiny_dlon = 1.0d-2 * pi/180._8
             tiny_dlat = 1.0d-2 * pi/180._8
             clon = tmp_center_lons(i,j)   ! radians
             clat = tmp_center_lats(i,j)   ! radians
             p1 = [clon - tiny_dlon, clat - tiny_dlat]
             p2 = [clon + tiny_dlon, clat - tiny_dlat]
             p3 = [clon + tiny_dlon, clat + tiny_dlat]
             p4 = [clon - tiny_dlon, clat + tiny_dlat]
         
             ! wrap each longitude of the tiny fallback polygon
             p1(1) = modulo(p1(1), 2.0d0*pi)
             p2(1) = modulo(p2(1), 2.0d0*pi)
             p3(1) = modulo(p3(1), 2.0d0*pi)
             p4(1) = modulo(p4(1), 2.0d0*pi)
         
             ! compute area
             area_signed = get_signed_area_spherical_polygon(p1,p2,p3,p4)
             if (area_signed <= 0.d0 .or. abs(area_signed) < 1.0d-12) then
               area_signed = 1.0d-12
             end if
             SCRIP_Area(n) = abs(area_signed)
         
             ! write fallback into SCRIP arrays
             SCRIP_CornerLon(:,n) = modulo([p1(1),p2(1),p3(1),p4(1)]*(180._8/pi),360.0_8)
             SCRIP_CornerLat(:,n) =        [p1(2),p2(2),p3(2),p4(2)]*(180._8/pi)
             fallback_mask(n)     = .true.
             failed_cells         = failed_cells + 1
         
             n = n + 1
             cycle
           end if
         
           ! 2) Clamp and correct the “real” node_xy_tmp before hull
           do k = 1,4
             node_xy_tmp(1,k) = modulo(node_xy_tmp(1,k), 2.0d0*pi)
         
             ! longitude clamp
             if (node_xy_tmp(1,k) <  1.0d-8)             node_xy_tmp(1,k) =  1.0d-8
             if (node_xy_tmp(1,k) >  2.0d0*pi - 1.0d-8)  node_xy_tmp(1,k) =  2.0d0*pi - 1.0d-8
         
             ! latitude clamp
             if (node_xy_tmp(2,k) < -pi/2.0d0) node_xy_tmp(2,k) = -pi/2.0d0
             if (node_xy_tmp(2,k) >  pi/2.0d0) node_xy_tmp(2,k) =  pi/2.0d0
           end do
         
           ! 3) Attempt the real hull
           call points_hull_2d(4, node_xy_tmp, hull_num, hull)
         
           ! If hull fails, try jitter fallback near cell center
           if (hull_num /= 4) then
             eps = 1.0d-5*pi/180.0d0
             node_xy_tmp(:,1) = [tmp_center_lons(i,j)-eps, tmp_center_lats(i,j)-eps]
             node_xy_tmp(:,2) = [tmp_center_lons(i,j)+eps, tmp_center_lats(i,j)-eps]
             node_xy_tmp(:,3) = [tmp_center_lons(i,j)+eps, tmp_center_lats(i,j)+eps]
             node_xy_tmp(:,4) = [tmp_center_lons(i,j)-eps, tmp_center_lats(i,j)+eps]
             call points_hull_2d(4, node_xy_tmp, hull_num, hull)
           end if
         
           ! If still fails, fallback to explicit tiny cell at real center
           if (hull_num /= 4) then
             tiny_dlon = 1.0d-2 * pi/180.0d0
             tiny_dlat = 1.0d-2 * pi/180.0d0
             clon = tmp_center_lons(i,j)
             clat = tmp_center_lats(i,j)
             node_xy_tmp(:,1) = [modulo(clon-tiny_dlon,2*pi), clat-tiny_dlat]
             node_xy_tmp(:,2) = [modulo(clon+tiny_dlon,2*pi), clat-tiny_dlat]
             node_xy_tmp(:,3) = [modulo(clon+tiny_dlon,2*pi), clat+tiny_dlat]
             node_xy_tmp(:,4) = [modulo(clon-tiny_dlon,2*pi), clat+tiny_dlat]
             hull_num = 4
             hull = [1,2,3,4]
           end if
         
           ! Ensure corners are finite; if not, rebuild tiny polygon
           bad_corner = .false.
           do k=1,4
             if (node_xy_tmp(1,k) /= node_xy_tmp(1,k) .or. abs(node_xy_tmp(1,k)) > 1.0d10) bad_corner = .true.
             if (node_xy_tmp(2,k) /= node_xy_tmp(2,k) .or. abs(node_xy_tmp(2,k)) > 1.0d10) bad_corner = .true.
           end do
           if (bad_corner) then
             clon = tmp_center_lons(i,j)
             clat = tmp_center_lats(i,j)
             tiny_dlon = 1.0d-4 * pi/180._8
             tiny_dlat = 1.0d-4 * pi/180._8
             node_xy_tmp(:,1) = [modulo(clon-tiny_dlon,2*pi), clat-tiny_dlat]
             node_xy_tmp(:,2) = [modulo(clon+tiny_dlon,2*pi), clat-tiny_dlat]
             node_xy_tmp(:,3) = [modulo(clon+tiny_dlon,2*pi), clat+tiny_dlat]
             node_xy_tmp(:,4) = [modulo(clon-tiny_dlon,2*pi), clat+tiny_dlat]
             hull_num = 4
             hull = [1,2,3,4]
             fallback_mask(n) = .true.
           endif
         
           ! 4) Reorder hull consistently -> p1,p2,p3,p4
           call safe_reorder_hull(node_xy_tmp, hull, p1, p2, p3, p4, n, i, j)
         
           area_signed = get_signed_area_spherical_polygon(p1, p2, p3, p4)
         
           ! If NaN/inf/tiny/huge area, rebuild a tiny CCW square and recompute
           if (area_signed /= area_signed .or. abs(area_signed) < 1.0d-12 .or. abs(area_signed) > 1.0d10) then
             clon = tmp_center_lons(i,j)
             clat = tmp_center_lats(i,j)
             tiny_dlon = 1.0d-4 * pi/180._8
             tiny_dlat = 1.0d-4 * pi/180._8
             p1 = [modulo(clon-tiny_dlon,2*pi), clat-tiny_dlat]
             p2 = [modulo(clon+tiny_dlon,2*pi), clat-tiny_dlat]
             p3 = [modulo(clon+tiny_dlon,2*pi), clat+tiny_dlat]
             p4 = [modulo(clon-tiny_dlon,2*pi), clat+tiny_dlat]
             area_signed = get_signed_area_spherical_polygon(p1, p2, p3, p4)
             fallback_mask(n) = .true.
           end if
         
           ! Enforce CCW orientation for SCRIP (positive signed area)
           if (area_signed < 0.0d0) then
             swap_p = p2; p2 = p4; p4 = swap_p
             area_signed = -area_signed
           endif
         
           ! Write CCW corners and POSITIVE area
           SCRIP_CornerLon(:,n) = modulo([p1(1),p2(1),p3(1),p4(1)]*(180._8/pi), 360.0_8)
           SCRIP_CornerLat(:,n) =        [p1(2),p2(2),p3(2),p4(2)]*(180._8/pi)
           SCRIP_Area(n)        = area_signed
         
         else
          ! ----- Regular grid path: enforce CLOCKWISE corners -----
          ! node_xy is (lon,lat) in radians with canonical perimeter order:
          !   1:(i,j), 2:(i+1,j), 3:(i+1,j+1), 4:(i,j+1)
          
          p1 = node_xy(:,1)
          p2 = node_xy(:,2)
          p3 = node_xy(:,3)
          p4 = node_xy(:,4)

          area_signed = get_signed_area_spherical_polygon(p1, p2, p3, p4)  ! >0 CCW, <0 CW
          
          if (area_signed > 0.d0) then
            swap_p = p2;  p2 = p4;  p4 = swap_p   ! make CW
          end if
          
          SCRIP_CornerLon(:,n) = modulo([p1(1),p2(1),p3(1),p4(1)]*(180._8/pi), 360.0_8)
          SCRIP_CornerLat(:,n) =        [p1(2),p2(2),p3(2),p4(2)]*(180._8/pi)
          
          SCRIP_Area(n) = sph_tri_area_rad(p1,p2,p3) + sph_tri_area_rad(p1,p3,p4)  ! steradians
          if (SCRIP_Area(n) <= 0.d0) SCRIP_Area(n) = 1.d-12          
         end if           
          !==============================================================
          ! Per-cell metrics & rrfac inputs
          ! - Regular: rrfac(n) holds the per-cell length (legacy); we
          !            convert to final rrfac after the loop.
          ! - Stretched: rrfac is computed later (Gaussian profile);
          !              here we only collect a robust per-cell length.
          !==============================================================
          if (.not. do_schmidt) then
             !----- ORIGINAL REGULAR GRID LOGIC (legacy) -----
             dummy_max = 0.0d0
             dummy_min = huge(1.0d0)
             call get_grid_length(p1, p2, p3, tmp_len, dummy_max, dummy_min)

             ! Legacy behavior: store per-cell length in SCRIP_rrfac(n)
             SCRIP_rrfac(n) = tmp_len

             ! Also track the same per-cell length into our arrays so the
             ! post-loop reduction works uniformly for both modes
             local_max_length_all(n) = tmp_len
             local_min_length_all(n) = tmp_len

          else
             !----- STRETCHED GRID PATH -----
             ! Keep the area guard
             if (SCRIP_Area(n) /= SCRIP_Area(n) .or. SCRIP_Area(n) <= 0.0d0) then
                write(*,*) '[CRITICAL FIX] Non-positive or NaN area at cell:', n, ' original area:', SCRIP_Area(n)
                SCRIP_Area(n)   = 1.0d-12
                fallback_mask(n) = .true.
             endif

             ! Use the same geometric length metric as regular grids
             dummy_max = 0.0d0
             dummy_min = huge(1.0d0)
             call get_grid_length(p1, p2, p3, tmp_len, dummy_max, dummy_min)

             ! Store per-cell length for later global reductions
             local_max_length_all(n) = tmp_len
             local_min_length_all(n) = tmp_len
          endif

          ! Minimal diagnostics (both modes)
          if (local_min_length_all(n) <= 1.0d-12) then
             write(*,*) "CRITICAL DEBUG cell:", n, i, j
             write(*,*) "local_min_length very small or zero:", local_min_length_all(n)
             write(*,*) "Polygon lengths (km):", &
                        great_circle_dist(p1,p2,6371.d0), &
                        great_circle_dist(p2,p3,6371.d0), &
                        great_circle_dist(p3,p4,6371.d0), &
                        great_circle_dist(p4,p1,6371.d0)
             write(*,*) "Polygon coords (deg):", &
                        "p1", p1*180/pi, "p2", p2*180/pi, &
                        "p3", p3*180/pi, "p4", p4*180/pi
          endif

          n = n + 1
        end do  ! closes i loop
      end do    ! closes j loop

      write(*,*) 'Finished per-cell geometry/length pass'
      call MPI_Barrier(mpiC, mpi_err)

      !---------------- Global min/max of per-cell length ----------------
      local_max_length = maxval(local_max_length_all)
      local_min_length = minval(local_min_length_all)

      call MPI_Allreduce(local_max_length, global_max_length, 1, MPI_DOUBLE_PRECISION, MPI_MAX, mpiC, mpi_err)
      call MPI_Allreduce(local_min_length, global_min_length, 1, MPI_DOUBLE_PRECISION, MPI_MIN, mpiC, mpi_err)

      if (global_min_length <= 0.0d0 .or. global_min_length /= global_min_length) then
         global_min_length = max(1.0d-12, global_max_length * 1.0d-4)
      end if
      min_allowed_length = max(1.0d-12, global_max_length * 1.0d-3)


      if (do_schmidt) then
        !--- Guard areas and compute global area bounds (diagnostics only)
        do cell = n_start, n_end
          if (SCRIP_Area(cell) < 1.0d-12 .or. SCRIP_Area(cell) /= SCRIP_Area(cell)) then
            SCRIP_Area(cell)    = 1.0d-12
            fallback_mask(cell) = .true.
          end if
        end do
      
        valid_count_local = count(.not. fallback_mask(n_start:n_end))
        if (valid_count_local > 0) then
          max_area_local = maxval(SCRIP_Area(n_start:n_end), mask=.not. fallback_mask(n_start:n_end))
          min_area_local = minval(SCRIP_Area(n_start:n_end), mask=.not. fallback_mask(n_start:n_end))
        else
          max_area_local = 0.0d0
          min_area_local = huge(1.0d0)
        end if
        call MPI_Allreduce(max_area_local, global_max_area, 1, MPI_DOUBLE_PRECISION, MPI_MAX, mpiC, mpi_err)
        call MPI_Allreduce(min_area_local, global_min_area, 1, MPI_DOUBLE_PRECISION, MPI_MIN, mpiC, mpi_err)
        if (global_max_area <= 0.0d0 .or. global_max_area /= global_max_area) global_max_area = 1.0d-12
        if (global_min_area <= 0.0d0 .or. global_min_area /= global_min_area) global_min_area = 1.0d-12
      
        !--- Build the Gaussian radial rrfac about (target_lon, target_lat)
        max_rrfac_allowed    = 100.0d0
        clamp_count_local    = 0
        fallback_count_local = count(fallback_mask(n_start:n_end))
      
        deg2rad = acos(-1.0d0)/180.0d0
        lonc    = target_lon; if (lonc < 0.d0) lonc = lonc + 360.d0
        latc    = target_lat
      
        ! Use global length contrast computed earlier to set the peak
        valid_count_local = count(.not. fallback_mask(n_start:n_end))
        if (valid_count_local > 0) then
          max_len_local = maxval(local_max_length_all(n_start:n_end), mask=.not. fallback_mask(n_start:n_end))
          min_len_local = minval(local_min_length_all(n_start:n_end), mask=.not. fallback_mask(n_start:n_end))
        else
          max_len_local = 0.0d0
          min_len_local = huge(1.0d0)
        end if
        call MPI_Allreduce(max_len_local, global_max_length, 1, MPI_DOUBLE_PRECISION, MPI_MAX, mpiC, mpi_err)
        call MPI_Allreduce(min_len_local, global_min_length, 1, MPI_DOUBLE_PRECISION, MPI_MIN, mpiC, mpi_err)
        if (global_max_length <= 0.d0 .or. global_max_length /= global_max_length) global_max_length = 1.0d0
        if (global_min_length <= 0.d0 .or. global_min_length /= global_min_length) global_min_length = global_max_length
      
          !==============================================================
          ! Width knob: choose a *continental* footprint.
          ! For sf=2.5 this gives half-power ~25° (covers most of CONUS).
          ! If you want even broader, bump 40 -> 45 or 50. 
          ! (FWHM radius = half_power_radius_deg by our definition.)          
          ! Width knob, narrower if stretch_factor is larger
          ! USING LENGTH: Peak magnitude from *edge-length* ratio (matches legacy behavior), capped
          ! --- Radial “bullseye” centered on (target_lon, target_lat) ---
          ! Continental footprint: ~25° half-power for sf≈2.5; scales gently with stretch         
          !==============================================================
        half_power_radius_deg = 40.0d0 / sqrt(max(1.0d0, dble(stretch_factor)))
        theta0 = (half_power_radius_deg*deg2rad) / sqrt(log(2.d0))
      
        ! Peak magnitude from *edge-length* ratio (legacy-aligned), capped
        ratio        = global_max_length / max(global_min_length, 1.0d-12)
        if (ratio < 1.d0) ratio = 1.d0
        desired_peak = min(max_rrfac_allowed, ratio)
        if (localPet == 0) then
          write(*,*) 'RRFAC length ratio:', ratio, ' desired_peak:', desired_peak, &
                     ' half_power_radius_deg:', half_power_radius_deg
        end if
      
        do cell = n_start, n_end
          if (fallback_mask(cell)) then
            SCRIP_rrfac(cell) = 1.d0
          else
            dlon = modulo(SCRIP_CenterLon(cell) - lonc + 180.d0, 360.d0) - 180.d0
            lat1 = SCRIP_CenterLat(cell) * deg2rad
            lat2 = latc                 * deg2rad
            dlam = dlon * deg2rad
            cc   = sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(dlam)
            cc   = max(-1.d0, min(1.d0, cc))
            dtheta = acos(cc)
      
            ! Gaussian radial profile: 1 far away, ~desired_peak at center
            w = exp( - (dtheta/theta0)**2 )
            SCRIP_rrfac(cell) = 1.d0 + (desired_peak - 1.d0) * w
      
            if (SCRIP_rrfac(cell) > max_rrfac_allowed) then
              SCRIP_rrfac(cell) = max_rrfac_allowed
              clamp_count_local = clamp_count_local + 1
            end if
          end if
        end do
      
        ! --- Peak/diagnostics block (global max, apex, counts) ---
        if (any(.not. fallback_mask(n_start:n_end))) then
          local_max_rrfac = maxval(SCRIP_rrfac(n_start:n_end), mask=.not. fallback_mask(n_start:n_end))
        else
          local_max_rrfac = -huge(1.0d0)
        end if
        call MPI_Allreduce(local_max_rrfac, global_max_rrfac, 1, MPI_DOUBLE_PRECISION, MPI_MAX, mpiC, mpi_err)
      
        ! Choose unique apex: the near-max cell closest to (lonc,latc)
        best_dist = huge(1.0d0); best_idx = -1
        lat2 = latc * deg2rad
        do cell = n_start, n_end
          if (.not. fallback_mask(cell)) then
            if (SCRIP_rrfac(cell) >= (1.0d0 - 1.0d-6)*global_max_rrfac) then
              dlon = modulo(SCRIP_CenterLon(cell) - lonc + 180.d0, 360.d0) - 180.d0
              lat1 = SCRIP_CenterLat(cell) * deg2rad
              dlam = dlon * deg2rad
              cc   = sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(dlam)
              cc   = max(-1.d0, min(1.d0, cc))
              dtheta = acos(cc)
              if (dtheta < best_dist) then
                best_dist = dtheta
                best_idx  = cell
              end if
            end if
          end if
        end do
      
        pair_local(1) = best_dist
        pair_local(2) = real(max(0,best_idx), 8)
        call MPI_Allreduce(pair_local, pair_global, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, mpiC, mpi_err)
        global_idx_max_rrfac = int(pair_global(2))
      
        ! Nudge chosen cell to be strictly largest, then recompute the peak
        if (global_idx_max_rrfac >= n_start .and. global_idx_max_rrfac <= n_end) then
          SCRIP_rrfac(global_idx_max_rrfac) = max(SCRIP_rrfac(global_idx_max_rrfac), global_max_rrfac*(1.0d0+1.0d-12))
        end if
        if (any(.not. fallback_mask(n_start:n_end))) then
          local_max_rrfac = maxval(SCRIP_rrfac(n_start:n_end), mask=.not. fallback_mask(n_start:n_end))
        else
          local_max_rrfac = -huge(1.0d0)
        end if
        call MPI_Allreduce(local_max_rrfac, global_max_rrfac, 1, MPI_DOUBLE_PRECISION, MPI_MAX, mpiC, mpi_err)
      
        ! rrfac_max for downstream
        rrfac_max = int(ceiling(min(max_rrfac_allowed, global_max_rrfac)))
        write(*,*) 'Computed rrfac_max:', rrfac_max
      
        ! Location of the max (owner PET + broadcast of lon/lat)
        if (global_idx_max_rrfac < 1 .or. global_idx_max_rrfac > size(SCRIP_CenterLon)) then
          owner_pet = 0; owner_lon = -999.d0; owner_lat = -999.d0
        else
          owner_pet = (global_idx_max_rrfac - 1) / (im_world*im_world)
          owner_lon = -999.d0; owner_lat = -999.d0
          if (localPet == owner_pet) then
            owner_lon = SCRIP_CenterLon(global_idx_max_rrfac)
            owner_lat = SCRIP_CenterLat(global_idx_max_rrfac)
          end if
        end if
        call MPI_Bcast(owner_lon, 1, MPI_DOUBLE_PRECISION, owner_pet, mpiC, mpi_err)
        call MPI_Bcast(owner_lat, 1, MPI_DOUBLE_PRECISION, owner_pet, mpiC, mpi_err)
        global_lon_max_rrfac = owner_lon
        global_lat_max_rrfac = owner_lat
      
        ! Diagnostics (min/max/clamp/fallback counts)
        if (any(.not. fallback_mask(n_start:n_end))) then
          local_min_rrfac = minval(SCRIP_rrfac(n_start:n_end), mask=.not. fallback_mask(n_start:n_end))
          local_max_rrfac = maxval(SCRIP_rrfac(n_start:n_end), mask=.not. fallback_mask(n_start:n_end))
        else
          local_min_rrfac = huge(1.0d0)
          local_max_rrfac = -huge(1.0d0)
        end if
        call MPI_Allreduce(local_min_rrfac,  global_min_rrfac_print, 1, MPI_DOUBLE_PRECISION, MPI_MIN, mpiC, mpi_err)
        call MPI_Allreduce(local_max_rrfac,  global_max_rrfac_print, 1, MPI_DOUBLE_PRECISION, MPI_MAX, mpiC, mpi_err)
        call MPI_Allreduce(clamp_count_local, clamp_count_global, 1, MPI_INTEGER, MPI_SUM, mpiC, mpi_err)
        call MPI_Allreduce(fallback_count_local, fallback_count_global, 1, MPI_INTEGER, MPI_SUM, mpiC, mpi_err)
      
        if (localPet == 0) then
          write(*,*) "Final RRFAC diagnostics:"
          write(*,*) "  Min rrfac:", global_min_rrfac_print
          write(*,*) "  Max rrfac:", global_max_rrfac_print
          write(*,*) "  rrfac_max:", rrfac_max
          write(*,*) "  [GLOBAL] clamped cells:", clamp_count_global, "  geometry fallbacks:", fallback_count_global
          write(*,*) "  MPI GLOBAL MAX RRFAC RESULTS:"
          write(*,*) "  global_max_rrfac:", global_max_rrfac
          write(*,*) "  global_idx_max_rrfac:", global_idx_max_rrfac
          write(*,*) "  owner PET:", owner_pet
          write(*,*) "  global_lon_max_rrfac:", global_lon_max_rrfac
          write(*,*) "  global_lat_max_rrfac:", global_lat_max_rrfac
        end if
      
        num_rrfac_max_local = count( (.not. fallback_mask(n_start:n_end)) .and. &
             (SCRIP_rrfac(n_start:n_end) >= (1.0d0 - 1.0d-5) * global_max_rrfac) )
        call MPI_Reduce(num_rrfac_max_local, num_rrfac_max, 1, MPI_INTEGER, MPI_SUM, 0, mpiC, mpi_err)
        if (localPet == 0) then
          write(*,*) ">>> Number of cells near max rrfac (0.001%):", num_rrfac_max
          write(*,*) ">>> Intended center: Lon =", target_lon, "Lat =", target_lat
          write(*,*) ">>> Actual max rrfac at Lon =", global_lon_max_rrfac, "Lat =", global_lat_max_rrfac
        end if
      
      else
        !---------------- Regular (non-stretch) final rrfac ----------------
        ! Legacy rule: rrfac = global_max_length / per_cell_length
        SCRIP_rrfac(n_start:n_end) = global_max_length / max(SCRIP_rrfac(n_start:n_end), 1.0d-12)
      end if
      !======================== END FINAL RRFAC =========================

 100     format(a,4f20.15)
 101     format(a,f20.15)
 102     format(2f20.15)
 103     format(a)
    deallocate( IMS )
    deallocate( JMS )

    scrip_size = IM_World*JM_World
    call MPI_Info_create(info, status);                                   _VERIFY(status)
    call MPI_Info_set(info, "cb_buffer_size", "1048576", status);         _VERIFY(status)

    if (len_trim(output_scrip) == 0) then
       write(*,*) 'ERROR: output_scrip is blank!'
       call ESMF_Finalize(rc=status)
       stop 1
    endif

    status = nf90_create(trim(output_scrip), IOR(NF90_MPIIO,IOR(NF90_CLOBBER,NF90_NETCDF4)), unit, comm=mpiC, info=info)
    _VERIFY(status)

    FMT = '(A,' // 'A,' //'A)'
    write(title,trim(FMT)) 'GMAO ',trim(gridname),' Grid'
    status = nf90_put_att(UNIT, NF90_GLOBAL, 'title',trim(title));        _VERIFY(status)
    status = nf90_put_att(UNIT, NF90_GLOBAL, 'GridDescriptionFormat','SCRIP'); _VERIFY(status)
    if (do_schmidt) then
       status = nf90_put_att(UNIT, NF90_GLOBAL, 'rrfac_max', rrfac_max);  _VERIFY(status)
    endif

    status = NF90_DEF_DIM(UNIT, 'grid_size'   , scrip_size,  gridsize);   _VERIFY(status)
    status = NF90_DEF_DIM(UNIT, 'grid_corners', grid_corners, cornerdimID); _VERIFY(status)
    status = NF90_DEF_DIM(UNIT, 'grid_rank'   , 1,          rankdim);     _VERIFY(status)

    ! grid_dims (unstructured marker)
    status = nf90_def_var(UNIT, "grid_dims", NF90_INT, [rankdim], griddim); _VERIFY(status)

    ! grid_imask
    status = nf90_def_var(UNIT, "grid_imask", NF90_INT, [gridsize], mask);   _VERIFY(status)
    status = nf90_put_att(UNIT, mask, "units", "unitless");                  _VERIFY(status)

    ! stretched/geometry fallback mask (1D)
    status = nf90_def_var(UNIT, "grid_fallback_mask", NF90_INT, [gridsize], varid_mask_fallback); _VERIFY(status)
    status = nf90_put_att(UNIT, varid_mask_fallback, "description", "1 = fallback area used, 0 = valid"); _VERIFY(status)
    status = nf90_put_att(UNIT, varid_mask_fallback, "units", "unitless"); _VERIFY(status)

    ! centers
    status = nf90_def_var(UNIT, "grid_center_lon", NF90_DOUBLE, [gridsize], centerlon); _VERIFY(status)
    status = nf90_put_att(UNIT, centerlon, "units", "degrees");                          _VERIFY(status)
    status = nf90_def_var(UNIT, "grid_center_lat", NF90_DOUBLE, [gridsize], centerlat); _VERIFY(status)
    status = nf90_put_att(UNIT, centerlat, "units", "degrees");                          _VERIFY(status)

    ! corners
    status = nf90_def_var(UNIT, "grid_corner_lon", NF90_DOUBLE, [cornerdimID,gridsize], cornerlon); _VERIFY(status)
    status = nf90_put_att(UNIT, cornerlon, "units", "degrees");                                      _VERIFY(status)
    status = nf90_def_var(UNIT, "grid_corner_lat", NF90_DOUBLE, [cornerdimID,gridsize], cornerlat); _VERIFY(status)
    status = nf90_put_att(UNIT, cornerlat, "units", "degrees");                                      _VERIFY(status)

    ! area
    status = nf90_def_var(UNIT, "grid_area", NF90_DOUBLE, [gridsize], cellarea); _VERIFY(status)
    status = nf90_put_att(UNIT, cellarea, "units", "radians^2");                 _VERIFY(status)

    ! rrfac only for stretched
    if (do_schmidt) then
       status = nf90_def_var(UNIT, "rrfac", NF90_DOUBLE, [gridsize], cellrrfac); _VERIFY(status)
       status = nf90_put_att(UNIT, cellrrfac, "units", "unitless");              _VERIFY(status)
    endif

    status = nf90_enddef(UNIT); _VERIFY(status)

    ! grid_dims value
    rc = NF90_PUT_VAR(UNIT, griddim, (/ scrip_size /))

    allocate (sendData(1), GlobalCounts(npets), recvCounts(npets), recvOffsets(npets), stat=status); _VERIFY(status)
    sendData    = tmp
    recvCounts  = 1
    recvOffsets = 0
    do i=2, npets
      recvOffsets(i) = recvOffsets(i-1) + recvCounts(i-1)
    end do
    call ESMF_VMGatherV(vm, sendData=sendData, sendCount=1, recvData=GlobalCounts, recvCounts=recvCounts, recvOffsets=recvOffsets, rootPet=0, rc=status); _VERIFY(status)
    call ESMF_VMBroadcast(vm, bcstData=GlobalCounts, count=npets, rootPet=0, rc=status); _VERIFY(status)
    
    ! ------------------- File vs memory spans -------------------
    ! FILE: global slice for this PET (must differ by PET)
    if (do_schmidt) then
      start_file = n_start                     ! already global index in stretched mode
    else
      start_file = localPet * num_local_cells + 1   
    end if
    cnt_file   = num_local_cells
    
    ! MEMORY: regular uses local-sized arrays; stretched uses global-sized
    if (do_schmidt) then
      mem_start = n_start
    else
      mem_start = 1
    end if
    mem_end = mem_start + cnt_file - 1    
    ! -----------------------------------------------------------
    
    ! 1D writes: centers / area / (rrfac only if stretched) / imask
    start(1) = start_file
    cnt(1)   = cnt_file
    
    status = NF90_PUT_VAR(UNIT, centerlon, SCRIP_CenterLon(mem_start:mem_end), start, cnt); _VERIFY(status)
    status = NF90_PUT_VAR(UNIT, centerlat, SCRIP_CenterLat(mem_start:mem_end), start, cnt); _VERIFY(status)
    status = NF90_PUT_VAR(UNIT, cellarea , SCRIP_Area     (mem_start:mem_end), start, cnt); _VERIFY(status)
    
    if (do_schmidt) then
      status = NF90_PUT_VAR(UNIT, cellrrfac, SCRIP_rrfac(mem_start:mem_end), start, cnt); _VERIFY(status)
    end if
    
    allocate(grid_imask(cnt_file), stat=status); _VERIFY(status)
    grid_imask = 1
    status = NF90_PUT_VAR(UNIT, mask, grid_imask, start, cnt); _VERIFY(status)
    deallocate(grid_imask)
    
    ! Fallback mask (1D)
    allocate(mask_fallback(cnt_file), stat=status); _VERIFY(status)
    mask_fallback = 0
    where (fallback_mask(mem_start:mem_end)) mask_fallback = 1
    start(1) = start_file
    cnt(1)   = cnt_file
    status = NF90_PUT_VAR(UNIT, varid_mask_fallback, mask_fallback, start, cnt); _VERIFY(status)
    deallocate(mask_fallback)
    
    ! 2D corners: (corner=1..4, cells = this PET's global slice)
    start(1) = 1
    start(2) = start_file
    cnt(1)   = grid_corners
    cnt(2)   = cnt_file
    
    allocate(my_corner_lat(grid_corners, cnt_file), stat=status); _VERIFY(status)
    allocate(my_corner_lon(grid_corners, cnt_file), stat=status); _VERIFY(status)
    
    my_corner_lat = SCRIP_CornerLat(:, mem_start:mem_end)
    my_corner_lon = SCRIP_CornerLon(:, mem_start:mem_end)
    
    status = NF90_PUT_VAR(UNIT, cornerlat, my_corner_lat, start, cnt); _VERIFY(status)
    status = NF90_PUT_VAR(UNIT, cornerlon, my_corner_lon, start, cnt); _VERIFY(status)
    
    deallocate(my_corner_lat, my_corner_lon)
    
    if (localPet == 0) write(*,*) 'WRITE slices  file:[', start_file, start_file+cnt_file-1, ']  mem:[', mem_start, mem_end, ']'      
    
    call MPI_Barrier(mpiC, mpi_err); _VERIFY(mpi_err)
    status = NF90_CLOSE(UNIT);       _VERIFY(status)
    call MPI_Barrier(mpiC, mpi_err); _VERIFY(mpi_err)
    call ESMF_VMBarrier(vm, rc=status); _VERIFY(status)
    call ESMF_Finalize ( rc=status); _VERIFY(status)
    
    if (allocated(SCRIP_CenterLat)) deallocate(SCRIP_CenterLat)
    if (allocated(SCRIP_CenterLon)) deallocate(SCRIP_CenterLon)
    if (allocated(SCRIP_CornerLat)) deallocate(SCRIP_CornerLat)
    if (allocated(SCRIP_CornerLon)) deallocate(SCRIP_CornerLon)
    if (allocated(SCRIP_Area))      deallocate(SCRIP_Area)
    if (allocated(SCRIP_rrfac))     deallocate(SCRIP_rrfac)
    if (allocated(sendData))        deallocate(sendData)
    if (allocated(GlobalCounts))    deallocate(GlobalCounts)
    deallocate(recvCounts)
    deallocate(recvOffsets)
    deallocate(fallback_mask)

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
!   Jun,2025 Robustness enhancements added for numerical stability and polar cases.
!
  implicit none

  integer, intent(in) :: node_num
  real(REAL64) angle
  real(REAL64) angle_max
  real(REAL64) di
  real(REAL64) dr
  real(REAL64) eps, lon_min, lon_max
  integer, intent(out) :: hull_num
  integer, intent(out) :: hull(node_num)
  integer :: i, q, r, first
  real(REAL64) node_xy(2,node_num)
  real(REAL64) p_xy(2)
  real(REAL64) q_xy(2)
  real(REAL64) r_xy(2)
  real(REAL64) centroid(2)
  logical :: polar_cell
  integer :: node_num_unique
  real(REAL64), allocatable :: unique_xy(:,:)
  real(REAL64), parameter :: pi = 3.141592653589793d0

  if ( node_num < 1 ) then
     hull_num = 0
     return
  elseif ( node_num <= 2 ) then
     hull_num = node_num
     hull(1:node_num) = (/ (i, i=1,node_num) /)
     return
  endif

  ! Step 1: Remove duplicate points  (numerical degeneracy fix)

  eps = 1.0d-12
  allocate(unique_xy(2, node_num))

  call remove_duplicate_points(node_xy, node_num, eps, unique_xy, node_num_unique)

  ! Update the original array:
  node_xy(:,1:node_num_unique) = unique_xy(:,1:node_num_unique)

  deallocate(unique_xy)

  ! Step 2: Handle longitude periodicity  

  do i = 1, node_num_unique
     node_xy(1,i) = modulo(node_xy(1,i), 2.0d0*pi)
  enddo

  lon_min = minval(node_xy(1,1:node_num_unique))
  lon_max = maxval(node_xy(1,1:node_num_unique))

  if (lon_max - lon_min > pi) then
     do i = 1, node_num_unique
        if (node_xy(1,i) < lon_min + pi) then
           node_xy(1,i) = node_xy(1,i) + 2.0d0*pi
        endif
     enddo
  endif

  ! Step 3: Detect and handle polar cells

    polar_cell = all(abs(node_xy(2,1:node_num_unique))*180.d0/pi >= 89.0d0)
    if (polar_cell) then
       call handle_polar_cell(node_num_unique, node_xy(:,1:node_num_unique), hull_num, hull)
       return
    endif

  ! Step 4: Original convex hull algorithm begins (unchanged original logic)

!  Find the leftmost point and call it "Q".
!  In case of ties, take the bottom-most.
!
  q = 1
  do i = 2, node_num_unique
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

    do i = 1, node_num_unique

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
    if ( node_num_unique < hull_num ) then
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

  pure function central_angle(v1, v2) result(ang)
    real(ESMF_KIND_R8), intent(in) :: v1(3), v2(3)
    real(ESMF_KIND_R8) :: cp(3), dotp, ang
    cp   = [ v1(2)*v2(3)-v1(3)*v2(2),  v1(3)*v2(1)-v1(1)*v2(3),  v1(1)*v2(2)-v1(2)*v2(1) ]
    dotp = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)
    ang  = atan2( sqrt(cp(1)**2 + cp(2)**2 + cp(3)**2), dotp )
  end function central_angle

   real(ESMF_KIND_R8) pure function sph_tri_area_rad(pA, pB, pC) result(area)
     implicit none
     real(ESMF_KIND_R8), intent(in) :: pA(2), pB(2), pC(2)
     real(ESMF_KIND_R8) :: a, b, c, s, t
     real(ESMF_KIND_R8) :: A3(3), B3(3), C3(3)
   
     A3 = [cos(pA(2))*cos(pA(1)), cos(pA(2))*sin(pA(1)), sin(pA(2))]
     B3 = [cos(pB(2))*cos(pB(1)), cos(pB(2))*sin(pB(1)), sin(pB(2))]
     C3 = [cos(pC(2))*cos(pC(1)), cos(pC(2))*sin(pC(1)), sin(pC(2))]
   
     a = central_angle(B3, C3)
     b = central_angle(C3, A3)
     c = central_angle(A3, B3)
   
     s = 0.5d0*(a+b+c)
     t = tan(0.5d0*s)*tan(0.5d0*(s-a))*tan(0.5d0*(s-b))*tan(0.5d0*(s-c))
   
     if (t <= 0.d0) then
        area = 0.d0
     else
        area = 4.d0*atan( sqrt(t) )
     end if
   end function sph_tri_area_rad

  function get_area_spherical_polygon(p1,p2,p3,p4) result(area)
    ! p1..p4 are (lon,lat) in radians, ordered around the cell
    real(ESMF_KIND_R8), intent(in) :: p1(2), p2(2), p3(2), p4(2)
    real(ESMF_KIND_R8) :: area
    area = sph_tri_area_rad(p1,p2,p3) + sph_tri_area_rad(p1,p3,p4)   ! steradians
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
   real (real64):: angle, ddd, threshold
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

    ddd = (px*px + py*py + pz*pz)*(qx*qx + qy*qy + qz*qz)
    
    threshold = 1.d-24
    if (ddd <= threshold) then
        angle = 0.d0
    else
        ddd   = (px*qx + py*qy + pz*qz) / sqrt(ddd)
        ddd   = min(1.d0, max(-1.d0, ddd))
        angle = acos(ddd)
    endif

   spherical_angle = angle
end function

subroutine local_abort(rc,line_number)
  integer, intent(in) :: rc
  integer, intent(in) :: line_number
  integer :: status
  write(*,*) 'Aborting at line', line_number, 'rc=', rc
  call MPI_Abort(MPI_COMM_WORLD, rc, status)
end subroutine local_abort


real(REAL64) function great_circle_dist( q1, q2, radius )
     real(REAL64), intent(IN)           :: q1(2), q2(2)
     real(REAL64), intent(IN), optional :: radius

     real (REAL64):: p1(2), p2(2)
     real (REAL64):: beta
     real(REAL64) :: dlon, dlat
     real(REAL64), parameter           :: pi = 3.141592653589793d0
     integer n

     do n=1,2
        p1(n) = q1(n)
        p2(n) = q2(n)
     enddo

        dlon = modulo( (p1(1)-p2(1)) + pi, 2.0d0*pi ) - pi
        dlat = p1(2) - p2(2)
        beta = sin(0.5d0*dlat)**2 + cos(p1(2))*cos(p2(2))*sin(0.5d0*dlon)**2
        beta = max(0.d0, min(1.d0, beta))
        beta = 2.d0*atan2( sqrt(beta), sqrt(max(0.d0, 1.d0 - beta)) )

     if ( present(radius) ) then
          great_circle_dist = radius * beta
     else
          great_circle_dist = beta   ! Returns the angle
     endif

end function great_circle_dist

subroutine remove_duplicate_points(points, num_points, eps, unique_points, num_unique)
   implicit none
   integer, intent(in) :: num_points
   real(REAL64), intent(in) :: points(2,num_points), eps
   real(REAL64), intent(out) :: unique_points(2,num_points)
   integer, intent(out) :: num_unique
   integer :: i, j
   logical :: duplicate

   num_unique = 0
   do i = 1, num_points
      duplicate = .false.
      do j = 1, num_unique
         if (sum(abs(points(:,i) - unique_points(:,j))) < eps) then
            duplicate = .true.
            exit
         endif
      enddo
      if (.not. duplicate) then
         num_unique = num_unique + 1
         unique_points(:,num_unique) = points(:,i)
      endif
   enddo
end subroutine

subroutine safe_reorder_hull(node_xy_tmp, hull, p1, p2, p3, p4, n, i, j)
  implicit none
  real(8), intent(in)  :: node_xy_tmp(2,4)
  integer, intent(in)  :: hull(4)
  real(8), intent(out) :: p1(2), p2(2), p3(2), p4(2)
  integer, intent(in)  :: n, i, j
  integer :: idx
  real(8) :: tmp_nodes(2,4)

  do idx = 1,4
     if (hull(idx) < 1 .or. hull(idx) > 4) then
        write(*,*) "CRITICAL INDEX ERROR at idx=", idx, &
                   " hull(idx)=", hull(idx), " cell (n,i,j)=", n,i,j
        stop "Exiting due to invalid hull index"
     endif
     tmp_nodes(:,idx) = node_xy_tmp(:,hull(idx))
  enddo

  call reorder_hull_quad(tmp_nodes, p1, p2, p3, p4)

end subroutine safe_reorder_hull

subroutine handle_polar_cell(node_num, points, hull_num, hull)
   implicit none
   integer, intent(in) :: node_num
   real(REAL64), intent(in) :: points(2,node_num)
   integer, intent(out) :: hull_num, hull(node_num)
   real(REAL64) :: angles(node_num), centroid(2)
   integer :: i, idx(node_num)

   centroid = [sum(points(1,:)), sum(points(2,:))] / node_num

   do i = 1, node_num
      angles(i) = atan2(points(2,i)-centroid(2), points(1,i)-centroid(1))
      idx(i) = i
   enddo

   call sort_by_angles(node_num, angles, idx)

   hull_num = node_num
   hull(1:node_num) = idx(1:node_num)
end subroutine

subroutine sort_by_angles(n, angles, idx)
   implicit none
   integer, intent(in) :: n
   real(REAL64), intent(inout) :: angles(n)
   integer, intent(inout) :: idx(n)
   integer :: i, j, tmp_idx
   real(REAL64) :: tmp_angle
   do i = 1, n-1
      do j = i+1, n
         if (angles(j) < angles(i)) then
            tmp_angle = angles(i); angles(i) = angles(j); angles(j) = tmp_angle
            tmp_idx = idx(i); idx(i) = idx(j); idx(j) = tmp_idx
         endif
      enddo
   enddo
end subroutine

subroutine reorder_hull_quad(node_xy, p1, p2, p3, p4)
  implicit none
  real(8), intent(in) :: node_xy(2,4)
  real(8), intent(out) :: p1(2), p2(2), p3(2), p4(2)
  real(8) :: centroid(2), angles(4)
  integer :: i, order(4), temp_order
  real(8) :: temp_angle, temp_x(4), temp_y(4)
  logical :: swapped
  real(8) :: xyz1(3), xyz2(3), xyz3(3), xyz4(3), normal(3), xyz_centroid(3), orientation, tmp_p(2)
  logical, parameter :: verbose = .false.

  ! Compute centroid
  centroid = [sum(node_xy(1,:))/4.0_8, sum(node_xy(2,:))/4.0_8]

  ! Compute angles from centroid to each corner
  do i = 1,4
      angles(i) = atan2(node_xy(2,i)-centroid(2), node_xy(1,i)-centroid(1))
  enddo
  ! Immediately at start of routine:
   if (verbose) then
      write(*,*) "INSIDE reorder_hull_quad (INPUT):"
       do i = 1,4
         write(*,*) " Corner ", i, ": Lon,Lat(deg):", node_xy(1,i)*180.0d0/pi, node_xy(2,i)*180.0d0/pi
       end do
   end if

  ! Initialize corner order
  order = [1, 2, 3, 4]

  ! Robust bubble sort angles (ascending order)
  do
      swapped = .false.
      do i = 1, 3
          if (angles(i) > angles(i+1)) then
              ! Swap angles
              temp_angle = angles(i)
              angles(i) = angles(i+1)
              angles(i+1) = temp_angle

              ! Swap order indices
              temp_order = order(i)
              order(i) = order(i+1)
              order(i+1) = temp_order

              swapped = .true.
          endif
      enddo
      if (.not. swapped) exit  ! Sorted correctly
  enddo

  ! Print explicitly sorted angles and corresponding corners
  temp_x = node_xy(1,:)
  temp_y = node_xy(2,:)

  ! Assign sorted corners explicitly
  p1 = [temp_x(order(1)), temp_y(order(1))]
  p2 = [temp_x(order(2)), temp_y(order(2))]
  p3 = [temp_x(order(3)), temp_y(order(3))]
  p4 = [temp_x(order(4)), temp_y(order(4))]

   xyz1 = lonlat_to_xyz(p1(1), p1(2))
   xyz2 = lonlat_to_xyz(p2(1), p2(2))
   xyz3 = lonlat_to_xyz(p3(1), p3(2))
   xyz4 = lonlat_to_xyz(p4(1), p4(2))

   normal = cross(xyz1, xyz2) + cross(xyz2, xyz3) + cross(xyz3, xyz4) + cross(xyz4, xyz1)
   xyz_centroid = (xyz1 + xyz2 + xyz3 + xyz4) / 4.0_8
   orientation = dot_product(normal, xyz_centroid)

   if (orientation > 0.0_8) then
       tmp_p = p2
       p2 = p4
       p4 = tmp_p
   endif
   ! Just before exiting:
   if (verbose) then
   write(*,*) "INSIDE reorder_hull_quad (OUTPUT):"
   write(*,*) " p1(deg):", p1*180.0d0/pi
   write(*,*) " p2(deg):", p2*180.0d0/pi
   write(*,*) " p3(deg):", p3*180.0d0/pi
   write(*,*) " p4(deg):", p4*180.0d0/pi
   end if

end subroutine reorder_hull_quad

! Calculate signed spherical polygon area using Girard's formula:
function get_signed_area_spherical_polygon(p1, p2, p3, p4) result(area_signed)
  implicit none
  real(8), intent(in) :: p1(2), p2(2), p3(2), p4(2)
  real(8) :: area_signed
  real(8) :: angles(4), excess
  real(8), dimension(3) :: xyz1, xyz2, xyz3, xyz4

  xyz1 = lonlat_to_xyz(p1(1), p1(2))
  xyz2 = lonlat_to_xyz(p2(1), p2(2))
  xyz3 = lonlat_to_xyz(p3(1), p3(2))
  xyz4 = lonlat_to_xyz(p4(1), p4(2))

  angles(1) = vertex_angle(xyz4, xyz1, xyz2)
  angles(2) = vertex_angle(xyz1, xyz2, xyz3)
  angles(3) = vertex_angle(xyz2, xyz3, xyz4)
  angles(4) = vertex_angle(xyz3, xyz4, xyz1)

  excess = sum(angles) - 2.0_8*pi
  area_signed = excess  ! negative for CW, positive for CCW
end function

! Helper function to calculate vertex angle:
function vertex_angle(xyzA, xyzB, xyzC) result(angle)
  implicit none
  real(8), intent(in) :: xyzA(3), xyzB(3), xyzC(3)
  real(8) :: angle
  real(8), dimension(3) :: AB, CB, cross_prod
  real(8) :: norm_cross, dot_prod

  AB = xyzA - xyzB
  CB = xyzC - xyzB

  cross_prod = cross(AB, CB)
  norm_cross = sqrt(dot_product(cross_prod, cross_prod))
  dot_prod = dot_product(AB, CB)

  angle = atan2(norm_cross, dot_prod)
end function

! Convert spherical coordinates (lon,lat) to Cartesian coordinates
function lonlat_to_xyz(lon, lat) result(xyz)
  real(kind=8), intent(in) :: lon, lat
  real(kind=8) :: xyz(3)

  xyz(1) = cos(lat) * cos(lon)
  xyz(2) = cos(lat) * sin(lon)
  xyz(3) = sin(lat)
end function

function cross(a, b) result(c)
  real(kind=8), intent(in) :: a(3), b(3)
  real(kind=8) :: c(3)

  c(1) = a(2)*b(3) - a(3)*b(2)
  c(2) = a(3)*b(1) - a(1)*b(3)
  c(3) = a(1)*b(2) - a(2)*b(1)
end function cross


    end program ESMF_GenerateCSGridDescription
