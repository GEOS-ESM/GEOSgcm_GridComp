program main
!Main purpose: Reads SMAP L4 runoff data (2016–2023) from a NetCDF file and computes the climatological mean discharge for each catchment.

  use omp_lib
  use river_read
  use constant, only : nlat=>nlat09, nlon=>nlon09, nc, nupmax

  implicit none

  ! Define variables:
  real, allocatable    :: runoff(:,:), qrunf(:), temp(:,:), qri(:), qin(:)
  integer, allocatable :: nts(:), downid(:), upstream(:,:)
  integer              :: i, j, nmax, did

  character(len=900)   :: file_path !"input/SMAPL4_OL7000_runoff_mean_2016_2023.nc"   

  if (command_argument_count() /= 1) then
    print *, "no <file_path> found"
    stop
  endif
  call get_command_argument(1, file_path)

  ! Allocate arrays for runoff (grid), catchment runoff, and a temporary grid array:
  allocate(runoff(nlon, nlat), qrunf(nc), temp(nlon, nlat))
  
  ! Read the "mean_runoff_flux" variable from the NetCDF file:
  call read_ncfile_real2d(trim(file_path), "mean_runoff_flux", runoff, nlon, nlat)
  ! Replace missing values (-9999) with 0:
  where(runoff == -9999.) runoff = 0.
  
  ! Flip the grid vertically (reverse the latitude order) and assign back to runoff:
  temp = runoff(:, nlat:1:-1)
  runoff = temp
  
  ! Convert runoff from [mm/s] to [mm/d]
  runoff = runoff * 86400.
  
  ! Map runoff from the M09 grid to catchments using the function M09_to_cat.
  ! The result is in kg/s; then convert to m^3/s by dividing by 1.e3.
  qrunf = M09_to_cat(runoff, nlon, nlat, nc)  ! kg/s
  qrunf = qrunf / 1.e3                        ! m^3/s

  ! Write catchment runoff (qrunf):
  open(88, file="output/Pfaf_qstr.txt")
  do i = 1, nc
    write(88, *) qrunf(i)
  end do

  ! Allocate arrays for "steps to sink" (nts), downstream id (downid) and aggregated runoff (qri):
  allocate(nts(nc), downid(nc), qri(nc))
  ! Read the number of steps to sink for each catchment from file:
  open(77, file="output/Pfaf_tosink.txt")
  read(77, *) nts
  ! Read the downstream connectivity (immediate downstream catchment id) from file:
  open(77, file="output/downstream_1D_new_noadj.txt")  
  read(77, *) downid

  ! Get the maximum number of steps among all catchments:
  nmax = maxval(nts)
  ! Initialize qri with the catchment runoff values:
  qri = qrunf
  ! Aggregate runoff upstream: For each catchment with a given number of steps j,
  ! add its runoff to its downstream catchment.
  do j = nmax, 1, -1
    do i = 1, nc  
      if (nts(i) == j) then
        did = downid(i)
        qri(did) = qri(did) + qri(i)
      endif
    end do
  end do

  ! Write the aggregated runoff (qri) to file "Pfaf_qri.txt":
  open(88, file="output/Pfaf_qri.txt")
  do i = 1, nc
    write(88, *) qri(i)
  end do  

  ! Allocate arrays for upstream connectivity and inlet discharge (qin):
  allocate(upstream(nupmax, nc), qin(nc))
  ! Read upstream connectivity information from file "upstream_1D.txt":
  open(77, file="output/upstream_1D.txt")
  read(77, *) upstream
  ! Initialize qin to -9999:
  qin = -9999.
  ! For catchments that have upstream connectivity (upstream(1,:) /= -1),
  ! set qin as the difference between outlet discharge (qri) and runoff (qrunf);
  ! for catchments with no upstream (upstream(1,:) == -1), set qin to half of direct runoff.
  where(upstream(1,:) /= -1) qin = qri - qrunf
  where(upstream(1,:) == -1) qin = qrunf / 2.
  
  ! Write the inlet discharge (qin):
  open(88, file="output/Pfaf_qin.txt")
  do i = 1, nc
    write(88, *) qin(i)
  end do   

contains
  !------------------------------------------------------------------------------
  ! Function: M09_to_cat
  ! Purpose : Maps runoff data from the M09 grid resolution to catchments using
  !           sub-area information. It aggregates runoff from sub-areas weighted by
  !           their area fractions.
  !
  ! Input:
  !   runoff - Runoff array of size (nlon, nlat) [in mm/d]
  !   nlon   - Number of longitude grid cells.
  !   nlat   - Number of latitude grid cells.
  !   ncat   - Number of catchments.
  !
  ! Output:
  !   Qrunf  - Runoff mapped to catchments (in kg/s, then converted to m^3/s).
  !------------------------------------------------------------------------------
  function M09_to_cat(runoff, nlon, nlat, ncat) result(Qrunf)
    integer, intent(in) :: nlon, nlat, ncat    ! Grid dimensions and number of catchments
    real, intent(in) :: runoff(nlon, nlat)       ! Input runoff array at grid resolution
    real :: Qrunf(ncat)                          ! Output catchment runoff array

    real, parameter :: small = 1.e-12         

    ! Define sub-area parameters (same as in the M09 dataset)
    integer, parameter :: nmax = 458             ! Maximum number of sub-areas per catchment
    integer, parameter :: nc = 291284            ! Total number of catchments

    ! Declare allocatable arrays to hold sub-area data:
    real, allocatable, dimension(:,:) :: subarea, frac  ! subarea: area of each sub-area, frac: fraction of total
    integer, allocatable, dimension(:,:) :: subx, suby    ! Coordinates of sub-areas in the grid
    real, allocatable, dimension(:) :: tot, runfC, fracA    ! tot: total catchment area; runfC: aggregated runoff; fracA: fraction sum
    integer, allocatable, dimension(:) :: nsub             ! nsub: number of sub-areas per catchment

    integer :: i, j, sx, sy  ! Loop counters and sub-area grid coordinates

    ! Allocate arrays for sub-area information and total area:
    allocate(nsub(nc), subarea(nmax, nc), subx(nmax, nc), suby(nmax, nc), tot(nc))

    ! Read sub-area data from text files:
    open(77, file="output/Pfaf_nsub_M09.txt"); read(77, *) nsub
    open(77, file="output/Pfaf_asub_M09.txt"); read(77, *) subarea
    open(77, file="output/Pfaf_xsub_M09.txt"); read(77, *) subx
    open(77, file="output/Pfaf_ysub_M09.txt"); read(77, *) suby
    open(77, file="output/Pfaf_area.txt"); read(77, *) tot

    ! Allocate fraction array (fraction of sub-area relative to total catchment area)
    allocate(frac(nmax, nc))

    ! Compute the fraction for each sub-area:
    do i = 1, nc
      frac(:, i) = subarea(:, i) / tot(i)
    end do

    ! Allocate arrays to accumulate runoff and fraction sums per catchment:
    allocate(runfC(nc), fracA(nc))
    runfC = 0.  ! Initialize aggregated runoff for each catchment to zero
    fracA = 0.  ! Initialize fraction accumulation to zero

    !$OMP PARALLEL default(shared) private(i,j,sx,sy)
    !$OMP DO
    ! Loop over all catchments and their sub-areas:
    do i = 1, nc
      do j = 1, nsub(i)
        sy = suby(j, i)  ! Get y-coordinate of the sub-area
        sx = subx(j, i)  ! Get x-coordinate of the sub-area
        ! Only consider valid sub-areas (non-zero fraction and valid runoff values)
        if (frac(j, i) > 0. .and. runoff(sx, sy) < 1.e14 .and. runoff(sx, sy) >= 0.) then
          runfC(i) = runfC(i) + frac(j, i) * runoff(sx, sy)
          fracA(i) = fracA(i) + frac(j, i)
        endif
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! Convert aggregated runoff to kg/s by multiplying by total catchment area (in m²)
    ! and dividing by the number of seconds per day (86400):
    Qrunf = runfC * (tot * 1.e6) / 86400.
    
    ! Deallocate allocated arrays to free memory:
    deallocate(subarea, subx, suby, tot, frac, &
               runfC, fracA, nsub)

  end function M09_to_cat
  !------------------------------------------------------------------------------
  
end program main
