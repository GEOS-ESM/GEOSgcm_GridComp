module k_module
!module for K parameter calculations.

  use river_read
  use constant, only: nga

  implicit none
  private
  public :: read_usgs_data, process_usgs_data, find_nearest_coords, get_station_inf, regression
  public :: filter_station, cal_Kmodel, get_valide_stations_gageii

contains
!------------------------------------------------------------
  subroutine read_usgs_data(file_vel, file_dis, nl, data)
    !------------------------------------------------------------
    ! Subroutine: read_usgs_data
    ! Purpose   : Reads USGS velocity and discharge data from text files
    !             and stores the data in a 2D array.
    !
    ! Input:
    !   nl   - Total number of records (lines) to read.
    !
    ! Output:
    !   data - 2D array (nl x 2) where column 1 contains velocity and
    !          column 2 contains discharge.
    !------------------------------------------------------------
    character(len=900),intent(in)  :: file_vel, file_dis
    integer, intent(in)            :: nl
    real, allocatable, intent(out) :: data(:,:)
    
    character(len=100)             :: var(2)
    character(len=900)             :: filename
    character(len=100)             :: line
    character(len=100)             :: x(100)

    integer                        :: i, j, l, io, k
    integer, allocatable           :: nv(:)

    ! Define the variable names for the two data files
    var = (/ "velocity", "discharge" /)

    ! Allocate the data array with nl rows and 2 columns
    allocate(data(nl, 2))

    ! Loop over both "velocity" and "discharge" files
    do l = 1, 2
      ! Construct the file name 
      if(l==1)filename = file_vel
      if(l==2)filename = file_dis
      open(unit=77, file=trim(filename), status='old')

      ! Allocate a temporary array to count the number of valid tokens per line
      allocate(nv(nl))
      
      ! Read each line from the file
      do k = 1, nl
        read(77, '(A)', iostat=io) line
        if (io /= 0) then
          print *, "Error reading line ", k, " from file: ", trim(filename)
          exit
        endif
        
        ! Read tokens from the line and store them in array x
        do i = 1, 100
          read(line, *, iostat=io) (x(j), j = 1, i)
          if (io == -1) then
            exit
          endif
        end do
        
        ! Record the number of valid tokens read from this line
        nv(k) = i - 1
        
        ! If valid data is found, extract the first value for the data array;
        ! otherwise assign a missing value (-9999)
        if (nv(k) >= 1) then
          read(x(1), *, iostat=io) data(k, l)
        else
          data(k, l) = -9999
        end if
      end do

      ! Deallocate the temporary token count array
      deallocate(nv)

      ! Close the current file
      close(77)
    end do
  end subroutine read_usgs_data

!------------------------------------------------------------
  subroutine process_usgs_data(file_usid, nl, ns, data, nv, nt, vel, dis)
    !------------------------------------------------------------
    ! Subroutine: process_usgs_data
    ! Purpose   : Processes the raw USGS data by reading unique station IDs,
    !             counting valid records per station.
    !
    ! Input:
    !   nl   - Total number of records.
    !   data - 2D array with raw velocity and discharge data.
    !
    ! Output:
    !   ns   - Number of unique stations.
    !   nv   - Array with the count of valid records per station.
    !   nt   - Total number of valid records.
    !   vel  - Array of velocity values.
    !   dis  - Array of discharge values.
    !------------------------------------------------------------
    character(len=900),intent(in)     :: file_usid    
    integer, intent(in)               :: nl
    integer, intent(out)              :: ns
    real, intent(inout), allocatable  :: data(:,:)
    real, allocatable, intent(out)    :: vel(:), dis(:)
    integer, allocatable, intent(out) :: nv(:)
    integer, intent(out)              :: nt

    character(len=20), allocatable    :: id(:)
    integer, allocatable              :: nu(:)
    character(len=20), allocatable    :: idu(:)
    integer                           :: i, k, ii

    ! Allocate array to hold station IDs for each record
    allocate(id(nl))

    ! Read station IDs from file "input/USGSID.txt"
    open(unit=11, file=trim(file_usid), status="old")
    read(11, *) id
    close(11)

    ! Initialize station count and count unique IDs
    k = 1
    do i = 2, nl
      if (.not.(trim(id(i)) == trim(id(i-1)))) then
        k = k + 1
      end if
    end do
    ns = k
    allocate(nu(ns), nv(ns))
    allocate(idu(ns))
    
    nu(1) = 1
    idu(1) = id(1)
    k = 1
    do i = 2, nl
      if (trim(id(i)) == trim(id(i-1))) then
        nu(k) = nu(k) + 1
      else
        k = k + 1
        nu(k) = 1
        idu(k) = id(i)
      end if
    end do

    ! Write unique station IDs to files (with and without commas)
    open(unit=13, file="temp/id_for_site.txt")
    do i = 1, ns
      write(13, '(A)') trim(idu(i)) // ","
    end do
    close(13)
    open(unit=13, file="temp/id_for_site_nocomma.txt")
    do i = 1, ns
      write(13, '(A)') trim(idu(i))
    end do
    close(13)    

    ! Read record
    nv = 0
    nv(1) = 1  
    k = 1
    ii = 0
    k = k - 1
    do i = 2, nl
      if (id(i) == id(i - 1)) then
        k = k + 1
        if (data(i,1) <= 0.0) then
          k = k - 1
        else if (data(i,2) <= 0.0) then
          k = k - 1
        end if
      else
        nv(ii + 1) = k  
        k = 1
        ii = ii + 1
        if (data(i,1) <= 0.0) then
          k = k - 1
        else if (data(i,2) <= 0.0) then
          k = k - 1
        end if
      end if
    end do
    nv(ii + 1) = k  
    nt = sum(nv)
    allocate(vel(nt), dis(nt))
    k = 0
    do i = 1, nl
      if (data(i,1) > 0.0 .and. data(i,2) > 0.0) then
        k = k + 1
        vel(k) = data(i,1)
        dis(k) = data(i,2)
      endif
    enddo

    ! Deallocate temporary arrays
    deallocate(id)
    deallocate(nu)
    deallocate(idu)
    deallocate(data)

  end subroutine process_usgs_data  

!------------------------------------------------------------
  subroutine find_nearest_coords(file_lats, file_lons, file_lat1m, file_lon1m, ns, nlat, nlon, lats, lons, lati, loni)
    !------------------------------------------------------------
    ! Subroutine: find_nearest_coords
    ! Purpose   : For each station, finds the nearest grid point in a 1-minute
    !             resolution grid and returns the corresponding indices.
    !
    ! Input:
    !   ns    - Number of stations.
    !   nlat  - Number of latitude grid points in the high-resolution grid.
    !   nlon  - Number of longitude grid points in the high-resolution grid.
    !
    ! In/Out:
    !   lats, lons - Arrays to store the station latitude and longitude values.
    !
    ! Output:
    !   lati, loni - Arrays of indices corresponding to the nearest grid points.
    !------------------------------------------------------------
    character(len=900),intent(in)     :: file_lats, file_lons, file_lat1m, file_lon1m    
    integer, intent(in)               :: ns, nlat, nlon
    real, allocatable, intent(inout)  :: lats(:), lons(:)
    integer, allocatable, intent(out) :: lati(:), loni(:)
    real, allocatable                 :: lat1m(:), lon1m(:)
    real                              :: min_dist_lat, min_dist_lon, dist
    integer                           :: i, j, idx_min_lat, idx_min_lon

    ! Allocate output arrays for grid indices for each station
    allocate(lati(ns), loni(ns))
    ! Allocate arrays for station coordinates
    allocate(lats(ns), lons(ns))
    ! Allocate arrays for 1-minute grid coordinates
    allocate(lat1m(nlat), lon1m(nlon))

    ! Read station latitudes from file "input/lat_for_site_200.txt"
    open(unit=10, file=trim(file_lats), status='old')
    do i = 1, ns
       read(10, *) lats(i)
    end do
    close(10)

    ! Read station longitudes from file "input/lon_for_site_200.txt"
    open(unit=11, file=trim(file_lons), status='old')
    do i = 1, ns
       read(11, *) lons(i)
    end do
    close(11)

    ! Read high-resolution latitude grid from file "input/lat_1m.txt"
    open(unit=12, file=trim(file_lat1m), status='old')
    do i = 1, nlat
       read(12, *) lat1m(i)
    end do
    close(12)

    ! Read high-resolution longitude grid from file "input/lon_1m.txt"
    open(unit=13, file=trim(file_lon1m), status='old')
    do i = 1, nlon
       read(13, *) lon1m(i)
    end do
    close(13)

    ! For each station, determine the nearest latitude and longitude indices
    do i = 1, ns
       min_dist_lat = 1.0e20
       min_dist_lon = 1.0e20
       idx_min_lat = -1
       idx_min_lon = -1

       ! Find nearest latitude index
       do j = 1, nlat
          dist = abs(lats(i) - lat1m(j))
          if (dist < min_dist_lat) then
             min_dist_lat = dist
             idx_min_lat = j
          end if
       end do
       lati(i) = idx_min_lat

       ! Find nearest longitude index
       do j = 1, nlon
          dist = abs(lons(i) - lon1m(j))
          if (dist < min_dist_lon) then
             min_dist_lon = dist
             idx_min_lon = j
          end if
       end do
       loni(i) = idx_min_lon
    end do

    ! Deallocate high-resolution grid arrays
    deallocate(lat1m)
    deallocate(lon1m)
  end subroutine find_nearest_coords
!------------------------------------------------------------
  subroutine get_station_inf(file_pfafmap, ns, nc, nlat, nlon, lati, loni, catid, Qclmt, slp, KImodel_all, exp_slp, exp_clmt, fac_str)
    !------------------------------------------------------------
    ! Subroutine: get_station_inf
    ! Purpose   : Retrieves station catchment information from a NetCDF file,
    !             assigns climate runoff (Qclmt) and slope values for each station,
    !             and computes modeled K values for all catchments.
    !
    ! Input:
    !   ns         - Number of stations.
    !   nc         - Total number of catchments.
    !   nlat, nlon - Dimensions of the grid.
    !   lati, loni - Grid indices for each station.
    !   exp_slp, exp_clmt - Exponents for slope and climatology discharge.
    !   fac_str    - Scaling factor for stream.
    !
    ! Output:
    !   catid      - Array of catchment IDs for each station.
    !   Qclmt      - Array of climatology discharge values for stations.
    !   slp        - Array of slope values for stations.
    !   KImodel_all - Array of modeled K values for all catchments.
    !------------------------------------------------------------
    character(len=900),intent(in)     :: file_pfafmap    
    integer, intent(in)               :: ns, nc, nlat, nlon
    integer, intent(in)               :: lati(nlon), loni(nlon)
    integer, allocatable, intent(out) :: catid(:)
    real, allocatable, intent(out)    :: Qclmt(:), slp(:)
    real, allocatable, intent(out)    :: KImodel_all(:)
    real, intent(in)                  :: exp_slp, exp_clmt, fac_str

    integer, allocatable              :: catchind(:,:)
    real, allocatable, dimension(:)   :: Qclmt_all, slp_all, Kstr_all, Qstr_all
    integer                           :: i

    ! Allocate arrays for the catchment index and station outputs
    allocate(catchind(nlon, nlat), catid(ns))
    allocate(Qclmt_all(nc), slp_all(nc))
    allocate(Qclmt(ns), slp(ns))
    allocate(KImodel_all(nc), Kstr_all(nc), Qstr_all(nc))    
    
    ! Read catchment index data from the NetCDF file "input/SRTM_PfafData.nc"
    call read_ncfile_int2d(trim(file_pfafmap), "CatchIndex", catchind, nlon, nlat)

    ! For each station, assign the catchment ID based on its grid location
    do i = 1, ns
      catid(i) = catchind(loni(i), lati(i))
    end do

    ! Write station catchment IDs to a temporary file
    open(88, file="temp/catid_for_site_200.txt")
    do i = 1, ns
      write(88, *) catid(i)
    end do    
    close(88)

    ! Read climate runoff data from file "output/Pfaf_qri.txt"
    open(77, file="output/Pfaf_qri.txt")
    read(77, *) Qclmt_all
    where(Qclmt_all < 1.e-8) Qclmt_all = 1.e-8
    ! Read slope data from file "temp/Pfaf_slope.txt"
    open(77, file="temp/Pfaf_slope.txt")
    read(77, *) slp_all           
    ! Read clmt discharge data from file "output/Pfaf_qstr.txt"
    open(77, file="output/Pfaf_qstr.txt")
    read(77, *) Qstr_all    
    where(Qstr_all < 1.e-8) Qstr_all = 1.e-8

    ! For each station, assign Qclmt and slope using the catchment ID
    do i = 1, ns
      if (catid(i) /= -9999) then
        Qclmt(i) = Qclmt_all(catid(i))
        slp(i) = slp_all(catid(i))
      else
        Qclmt(i) = -9999
        slp(i) = -9999
      endif
    end do 

    ! Calculate modeled K values for all catchments
    KImodel_all = (Qclmt_all**(exp_clmt)) * (slp_all**(exp_slp))

    ! Calculate stream K values using the scaling factor
    Kstr_all = fac_str * (Qstr_all**(exp_clmt)) * (slp_all**(exp_slp))

    ! Write stream K values to an output file
    open(88, file="output/Pfaf_Kstr_PR_fac1_0p35_0p45_0p2_n0p2.txt")
    do i = 1, nc
      write(88, *) Kstr_all(i)
    end do    
    close(88)

    ! Deallocate temporary arrays
    deallocate(catchind, Qclmt_all, slp_all, Kstr_all, Qstr_all)

  end subroutine get_station_inf
!------------------------------------------------------------
  subroutine get_valide_stations_gageii(file_gage_id, file_gage_acar, ns, nc, catid_sta, flag_thres)
    !------------------------------------------------------------
    ! Subroutine: get_valide_stations_gageii
    ! Purpose   : Compares station drainage area with GAGE-II dataset and applies an
    !             area ratio threshold to determine valid stations.
    !
    ! Input:
    !   ns        - Number of stations.
    !   nc        - Total number of catchments.
    !   catid_sta - Array of catchment IDs for stations.
    !
    ! Output:
    !   flag_thres - Array indicating valid stations (1 for valid, 0 otherwise).
    !------------------------------------------------------------
    character(len=900),intent(in)     :: file_gage_id, file_gage_acar    
    integer, intent(in)               :: ns, nc
    integer, intent(in)               :: catid_sta(ns)
    integer, allocatable, intent(out) :: flag_thres(:)

    real                              :: thr_sel = 0.3  ! Threshold selection factor

    real, dimension(:), allocatable   :: acar_pfaf
    integer                           :: i, j, k, cid
    character(len=20)                 :: id_gages(nga)
    character(len=20)                 :: id_sta(ns)
    integer                           :: flag_gageii(ns)
    real                              :: acar_gages(nga)
    real                              :: acar_gages_sta(ns), acar_sta(ns)
    character(len=20)                 :: line
    integer                           :: ios

    allocate(flag_thres(ns))
  
    ! Initialize station area ratios with a missing value
    acar_sta = -9999.0
    k = 0
  
    ! Read GAGE-II station IDs from "input/id_gagesii.txt"
    open(unit=10, file=trim(file_gage_id), status="old", action="read")
    do j = 1, nga
      read(10, '(A)', iostat=ios) id_gages(j)
      if (ios /= 0) then
        print *, "Error reading id_gagesii.txt"
        stop
      end if
    end do
    close(10)
  
    ! Read station IDs for the sites from "temp/id_for_site_nocomma.txt"
    open(unit=11, file="temp/id_for_site_nocomma.txt", status="old", action="read")
    do i = 1, ns
      read(11, '(A)', iostat=ios) id_sta(i)
      if (ios /= 0) then
        print *, "Error reading id_for_site_nocomma.txt"
        stop
      end if
    end do
    close(11)
  
    ! Read area ratios for GAGE-II stations from "input/acar_gagesii.txt"
    open(unit=12, file=trim(file_gage_acar), status="old", action="read")
    do j = 1, nga
      read(12, *, iostat=ios) acar_gages(j)
      if (ios /= 0) then
        print *, "Error reading acar_gagesii.txt"
        stop
      end if
    end do
    close(12)
  
    ! Initialize the GAGE-II flag array to zero (no match)
    flag_gageii = 0
    ! Compare station IDs with GAGE-II IDs and mark matches
    do i = 1, ns
      do j = 1, nga
        if (trim(id_gages(j)) == trim(id_sta(i))) then
          acar_gages_sta(i) = acar_gages(j)
          flag_gageii(i) = 1
          k = k + 1
          exit  ! Exit loop after a match is found
        end if
      end do
    end do

    print *, "Number of matches:", sum(flag_gageii)

    allocate(acar_pfaf(nc))
    open(77, file="temp/Pfaf_acar.txt")
    read(77, *) acar_pfaf
    close(77)

    ! For each station, assign the area ratio based on its catchment ID
    do i = 1, ns
      if (catid_sta(i) /= -9999) then
        cid = catid_sta(i)
        acar_sta(i) = acar_pfaf(cid)
      else
        acar_sta(i) = -9999.
      end if
    end do  

    ! Apply threshold criteria to flag valid stations
    flag_thres = 0
    do i = 1, ns
      if (flag_gageii(i) == 1 .and. catid_sta(i) /= -9999) then
        if (acar_sta(i) .ge. (1. - thr_sel) * acar_gages_sta(i) .and. &
            acar_sta(i) .le. (1. + thr_sel) * acar_gages_sta(i)) then
          flag_thres(i) = 1
        endif
      endif
    end do

    print *, "Number of valid:", sum(flag_thres)

    deallocate(acar_pfaf)
  end subroutine get_valide_stations_gageii
!------------------------------------------------------------
  subroutine regression(nt, vel_ori, dis_ori, nv, ns, Qclmt, slp, KKobs, KImodel, exp_slp, exp_clmt, mm, MU)
    !------------------------------------------------------------
    ! Subroutine: regression
    ! Purpose   : For each station with sufficient valid records, performs a
    !             regression between discharge and velocity to obtain a calibration
    !             factor, and then computes the observed K value (KKobs) for that station.
    !
    ! Input:
    !   nt      - Total number of valid records.
    !   ns      - Number of stations.
    !   nv      - Array containing the count of valid records per station.
    !   vel_ori - Original velocity data (in ft/s, will be converted).
    !   dis_ori - Original discharge data (in ft^3/s, will be converted).
    !   Qclmt   - Climatology discharge data for each station.
    !   slp     - Slope data for each station.
    !   exp_slp, exp_clmt - Exponents for slope and climatology discharge.
    !   mm, MU  - Model parameters.
    !
    ! Output:
    !   KKobs   - Array of observed K values for each station.
    !   KImodel - Array of modeled K values (init guess) for each station.
    !------------------------------------------------------------
    integer, intent(in)              :: nt, ns
    real, intent(inout), allocatable :: vel_ori(:), dis_ori(:)
    integer, intent(in)              :: nv(ns)
    real, intent(inout), allocatable :: Qclmt(:), slp(:)
    real, intent(out), allocatable   :: KKobs(:), KImodel(:)
    real, intent(in)                 :: exp_slp, exp_clmt, mm, MU

    real, allocatable, dimension(:)  :: x, y, yest
    integer                          :: thres = 100
    integer :: i, j
    real :: k(ns), cdtm(ns), med
    integer :: acc(ns)
    real, allocatable :: vel(:), dis(:)

    ! Convert velocity from ft/s to m/s and discharge from ft^3/s to m^3/s
    allocate(vel(nt), dis(nt))
    vel = vel_ori * 0.3048
    dis = dis_ori * 0.0283168

    ! Calculate cumulative counts to index into the valid records for each station
    acc(1) = nv(1)
    do i = 2, ns
      acc(i) = acc(i - 1) + nv(i)
    end do

    ! For each station with enough valid records, perform regression
    do i = 1, ns
      if (nv(i) >= thres) then
        allocate(x(nv(i)), y(nv(i)), yest(nv(i))) 
        x = dis(acc(i) - nv(i) + 1 : acc(i))**mm
        y = vel(acc(i) - nv(i) + 1 : acc(i))
        k(i) = sum(x * y) / sum(x * x)
        yest = k(i) * x
        cdtm(i) = cal_cdtm(y, yest)
        deallocate(x, y, yest)
      else
        k(i) = -9999.
        cdtm(i) = -9999.
      endif
    end do
    med = median(cdtm)

    ! Invalidate calibration factors for stations with low determination coefficient
    where(cdtm < 0.5) k = -9999.

    allocate(KKobs(ns))
    do i = 1, ns
      if (k(i) /= -9999. .and. Qclmt(i) /= -9999.) then
        KKobs(i) = k(i) / (Qclmt(i)**(MU - mm))
      else
        KKobs(i) = -9999.
      endif
    end do
    
    ! Calculate modeled K values (init guess) using the provided exponents
    allocate(KImodel(ns))
    KImodel = (Qclmt**(exp_clmt)) * (slp**(exp_slp))

    deallocate(vel, dis)
  end subroutine regression
!------------------------------------------------------------
  subroutine filter_station(nc, ns, np, lats_full, lons_full, Qclmt_full, slp_full, catid_full, KKobs_full, KImodel_full, Qclmt, slp, catid, KKobs, KImodel, flag_gageii)
    !------------------------------------------------------------
    ! Subroutine: filter_station
    ! Purpose   : Filters out stations that do not meet several criteria:
    !             valid catchment ID, valid K values, minimum slope threshold,
    !             and a positive GAGE-II flag. It then outputs the filtered data.
    !
    ! Input:
    !   nc                - Total number of catchments.
    !   ns                - Number of stations.
    !   lats_full, lons_full - Full arrays of station latitudes and longitudes.
    !   Qclmt_full, slp_full - Full climatology discharge and slope data for stations.
    !   catid_full        - Full catchment ID array for stations.
    !   KKobs_full, KImodel_full - Full observed and modeled K values (initial guess).
    !   flag_gageii       - GAGE-II validation flags.
    !
    ! Output:
    !   np      - Number of stations that passed the filter.
    !   Qclmt, slp, KKobs, KImodel - Filtered arrays for clmt discharge, slope, observed and modeled K (init guess).
    !   catid   - Filtered catchment IDs for the valid stations.
    !------------------------------------------------------------
    integer, intent(in)                 :: ns, nc
    integer, intent(out)                :: np
    real, intent(inout), allocatable    :: lats_full(:), lons_full(:), Qclmt_full(:), slp_full(:), KKobs_full(:), KImodel_full(:)
    real, intent(out), allocatable      :: Qclmt(:), slp(:), KKobs(:), KImodel(:)
    integer, intent(inout), allocatable :: catid_full(:)
    integer, intent(out), allocatable   :: catid(:)
    integer, intent(inout), allocatable :: flag_gageii(:)

    integer, allocatable                :: flag_slp(:)
    real, allocatable                   :: lats(:), lons(:)
    integer                             :: i, k
    integer, allocatable                :: flag_7065(:)

    ! Allocate and read slope flag data from file "temp/Pfaf_slope_flag.txt"
    allocate(flag_slp(nc))
    open(77, file="temp/Pfaf_slope_flag.txt")
    read(77, *) flag_slp 

    allocate(flag_7065(ns))
    flag_7065 = 0

    k = 0
    ! Count stations that meet all filtering criteria
    do i = 1, ns
      if (catid_full(i) .ne. -9999 .and. KKobs_full(i) /= -9999. .and. &
          slp_full(i) > 1.e-5 .and. flag_slp(catid_full(i)) == 1 .and. &
          flag_gageii(i) == 1) then
        k = k + 1
      endif
    end do
    np = k
    print *, "number of valid stations: ", np

    ! Allocate filtered output arrays
    allocate(Qclmt(np), slp(np), catid(np), KKobs(np), KImodel(np))
    allocate(lats(np), lons(np))
    k = 0
    do i = 1, ns
      if (catid_full(i) .ne. -9999 .and. KKobs_full(i) /= -9999. .and. &
          slp_full(i) > 1.e-5 .and. flag_slp(catid_full(i)) == 1 .and. &
          flag_gageii(i) == 1) then
        k = k + 1
        Qclmt(k) = Qclmt_full(i)
        slp(k) = slp_full(i)
        KKobs(k) = KKobs_full(i)
        KImodel(k) = KImodel_full(i)
        catid(k) = catid_full(i)
        lats(k) = lats_full(i)
        lons(k) = lons_full(i)
        flag_7065(i) = 1
      endif
    end do    

    ! Deallocate temporary full arrays that are no longer needed
    deallocate(Qclmt_full, slp_full, KKobs_full, KImodel_full, flag_slp, flag_gageii, lats, lons)

  end subroutine filter_station
!------------------------------------------------------------
  subroutine cal_Kmodel(ns, np, nc, MU, exp_slp, exp_clmt, Qclmt, slp, KKobs, KImodel, KImodel_all, catid, catid_full, ccr, rms)
    !------------------------------------------------------------
    ! Subroutine: cal_Kmodel
    ! Purpose   : Calibrates the model by adjusting catchment K values with a scaling
    !             factor computed from the percentiles of observed and modeled K values.
    !             It then computes the correlation coefficient (ccr) and RMS error.
    !
    ! Input/Output:
    !   ns        - Number of stations.
    !   np        - Number of valid stations.
    !   nc        - Total number of catchments.
    !   MU, exp_slp, exp_clmt - Model parameters.
    !   Qclmt, slp, KKobs, KImodel - Arrays for station data.
    !   KImodel_all - Modeled K values for all catchments.
    !   catid, catid_full - Filtered and full catchment ID arrays.
    !
    ! Output:
    !   ccr       - Correlation coefficient between observed and calibrated K.
    !   rms       - RMS error between observed and calibrated K.
    !------------------------------------------------------------
    integer, intent(in)                 :: ns, np, nc
    real, intent(in)                    :: MU, exp_slp, exp_clmt
    real, intent(inout), allocatable    :: Qclmt(:), slp(:), KKobs(:), KImodel(:)
    real, intent(inout), allocatable    :: KImodel_all(:)
    integer, intent(inout), allocatable :: catid(:), catid_full(:)
    real, intent(inout)                 :: ccr, rms
    
    real, allocatable                   :: KKobs_sort(:), KImodel_sort(:), KKmodel_full(:)
    real, allocatable, dimension(:)     :: dis, sca, Kv, KKmodel
    integer, allocatable, dimension(:)  :: gear

    character(len=50)                   :: MU_s, exp_slp_s, exp_clmt_s

    integer                             :: bulk, i, lev
    real                                :: Kper(11), KMper(11), rat(11), dis_full(11)

    ! Format model parameters into strings for output naming purposes
    write(MU_s, '(f4.2)') MU
    write(exp_slp_s, '(f4.2)') exp_slp
    if (exp_clmt >= 0.) then
      write(exp_clmt_s, '(f4.2)') exp_clmt
    else
      write(exp_clmt_s, '(f4.2)') -1.*exp_clmt
      exp_clmt_s = "n" // trim(exp_clmt_s)
    endif

    ! Allocate arrays for sorted K values
    allocate(KKobs_sort(np), KImodel_sort(np))
    call sort(np, KKobs, KKobs_sort)
    call sort(np, KImodel, KImodel_sort)

    ! Compute percentile thresholds by dividing sorted arrays into 10 equal parts
    bulk = np / 10
    Kper(1) = KKobs_sort(1)
    KMper(1) = KImodel_sort(1)
    do i = 2, 10
      Kper(i) = KKobs_sort(bulk * (i - 1))
      KMper(i) = KImodel_sort(bulk * (i - 1))
    end do
    Kper(11) = KKobs_sort(np)
    KMper(11) = KImodel_sort(np)
    rat = Kper / KMper

    ! Allocate arrays for scaling calculations over all catchments
    allocate(gear(nc), dis(nc), sca(nc), Kv(nc))
    
    ! Initialize gear to default (12) and compute distance to percentile thresholds
    gear = 12
    dis = -9999.
    do i = 1, nc
      do lev = 1, 11
        if (KImodel_all(i) <= KMper(lev)) then
          gear(i) = lev
          dis(i) = KMper(lev) - KImodel_all(i)
          exit
        endif
      end do
    end do    

    ! Calculate differences between consecutive percentile thresholds
    dis_full(1) = KMper(1)
    do i = 2, 11
      dis_full(i) = KMper(i) - KMper(i - 1)
    end do

    ! Compute scaling factors for each catchment based on its percentile position
    do i = 1, nc
      if (gear(i) == 1) then
        sca(i) = rat(1)
      elseif (gear(i) == 12) then
        sca(i) = rat(11)
      else
        sca(i) = ( rat(gear(i)-1) * dis(i) + rat(gear(i)) * (dis_full(gear(i)) - dis(i)) ) / dis_full(gear(i))
      endif 
      Kv(i) = KImodel_all(i) * sca(i)
    end do    

    ! Write scaled K values for each catchment to an output file
    open(88, file="output/Pfaf_Kv_PR_0p35_0p45_0p2_n0p2.txt")
    do i = 1, nc
      write(88, *) Kv(i)
    end do
    close(88)

    ! For each station, assign the corresponding scaled K value from its catchment
    allocate(KKmodel_full(ns))
    do i = 1, ns
      if (catid_full(i) /= -9999) then
        KKmodel_full(i) = Kv(catid_full(i))
      else
        KKmodel_full(i) = -9999.
      endif
    end do

    ! For filtered stations, extract the modeled K values
    allocate(KKmodel(np))
    do i = 1, np
      KKmodel(i) = Kv(catid(i))
    end do

    ! Compute correlation coefficient and RMS error between observed and modeled K values
    ccr = cal_ccr(KKobs, KKmodel)
    rms = cal_rms(KKobs, KKmodel, np)

    ! Deallocate temporary arrays and full data arrays
    deallocate(KKobs_sort, KImodel_sort)
    deallocate(KImodel_all, gear, dis, sca, Kv)
    deallocate(Qclmt, slp, KKobs, KImodel, catid, KKmodel, catid_full, KKmodel_full)
  end subroutine cal_Kmodel

  subroutine sort(np, data, data_sort)
    !------------------------------------------------------------
    ! Subroutine: sort
    ! Purpose   : Sorts an array of real numbers in ascending order using bubble sort.
    !
    ! Input:
    !   np        - Number of elements in the array.
    !   data      - Input array to be sorted.
    !
    ! Output:
    !   data_sort - Sorted array.
    !------------------------------------------------------------
    integer, intent(in) :: np           ! Size of the array
    real, intent(in)    :: data(np)          ! Input array
    real, intent(out)   :: data_sort(np)    ! Output sorted array
    integer             :: i, j
    real                :: temp

    ! Copy the input array to the output array
    data_sort = data

    ! Bubble sort algorithm
    do i = 1, np - 1
      do j = 1, np - i
        if (data_sort(j) > data_sort(j + 1)) then
          temp = data_sort(j)
          data_sort(j) = data_sort(j + 1)
          data_sort(j + 1) = temp
        end if
      end do
    end do
  end subroutine sort

  function cal_ccr(y, yest) result(ccr)
    !------------------------------------------------------------
    ! Function: cal_ccr
    ! Purpose : Calculates the correlation coefficient between observed and
    !           estimated arrays.
    !
    ! Input:
    !   y     - Observed data array.
    !   yest  - Estimated (modeled) data array.
    !
    ! Output:
    !   ccr   - Correlation coefficient.
    !------------------------------------------------------------
    real, intent(in) :: y(:)
    real, intent(in) :: yest(:)
    real             :: ccr
    real             :: mean_y, mean_yest
    real             :: sum_y, sum_yest
    real             :: sum_num, sum_den_y, sum_den_yest
    integer          :: n
    integer          :: i

    n = size(y)
    if (n /= size(yest)) then
      print *, "Error: Arrays must have the same length"
      ccr = 0.0
      return
    endif

    ! Compute means
    sum_y = sum(y)
    sum_yest = sum(yest)
    mean_y = sum_y / n
    mean_yest = sum_yest / n

    ! Compute numerator and denominators for the correlation coefficient
    sum_num = 0.0
    sum_den_y = 0.0
    sum_den_yest = 0.0
    do i = 1, n
      sum_num = sum_num + (y(i) - mean_y) * (yest(i) - mean_yest)
      sum_den_y = sum_den_y + (y(i) - mean_y)**2
      sum_den_yest = sum_den_yest + (yest(i) - mean_yest)**2
    end do

    if (sum_den_y == 0.0 .or. sum_den_yest == 0.0) then
      print *, "Error: Zero variance in input arrays"
      ccr = 0.0
    else
      ccr = sum_num / sqrt(sum_den_y * sum_den_yest)
    end if

  end function cal_ccr  

  function cal_rms(k_obs, k_model, n) result(rms)
    !------------------------------------------------------------
    ! Function: cal_rms
    ! Purpose : Calculates the relative root mean square error between observed
    !           and modeled K values.
    !
    ! Input:
    !   k_obs   - Observed K values array.
    !   k_model - Modeled K values array.
    !   n       - Number of elements.
    !
    ! Output:
    !   rms     - Relative RMS error.
    !------------------------------------------------------------
    implicit none
    integer, intent(in) :: n
    real, intent(in)    :: k_obs(n), k_model(n)
    real                :: rms
    real                :: sum_sq_diff
    integer             :: i

    sum_sq_diff = 0.0

    do i = 1, n
      sum_sq_diff = sum_sq_diff + ((k_model(i) - k_obs(i)) / k_obs(i))**2
    end do

    rms = sqrt(sum_sq_diff / n)
  end function cal_rms 

  function cal_cdtm(y, yest) result(dtmc)
    !------------------------------------------------------------
    ! Function: cal_cdtm
    ! Purpose : Computes the coefficient of determination (R^2) between observed
    !           and estimated data.
    !
    ! Input:
    !   y     - Observed data array.
    !   yest  - Estimated data array.
    !
    ! Output:
    !   dtmc  - Coefficient of determination (R^2).
    !------------------------------------------------------------
    real, intent(in) :: y(:)
    real, intent(in) :: yest(:)
    real             :: dtmc
    real             :: ss_tot, ss_res
    real             :: mean_y
    integer          :: n, i

    n = size(y)
    if (n /= size(yest)) then
      print *, "Error: Arrays must have the same length"
      dtmc = 0.0
      return
    endif

    mean_y = sum(y) / n
    ss_tot = sum((y - mean_y)**2)
    ss_res = sum((y - yest)**2)

    if (ss_tot == 0.0) then
      print *, "Error: Zero total sum of squares"
      dtmc = 0.0
    else
      dtmc = 1.0 - (ss_res / ss_tot)
    endif

  end function cal_cdtm

  function median(data) result(med)
    !------------------------------------------------------------
    ! Function: median
    ! Purpose : Computes the median of an array, ignoring values equal to -9999.0.
    !
    ! Input:
    !   data  - Array of real numbers.
    !
    ! Output:
    !   med   - Median value.
    !------------------------------------------------------------
    implicit none
    real, intent(in) :: data(:)  
    real             :: med                  
    real             :: sorted_data(size(data))  
    integer          :: n_valid           
    integer          :: i                 
    
    n_valid = 0
    do i = 1, size(data)
        if (data(i) /= -9999.0) then
            n_valid = n_valid + 1
            sorted_data(n_valid) = data(i)
        end if
    end do

    if (n_valid == 0) then
        med = -9999.0  
        return
    end if

    call sort2(sorted_data(1:n_valid))
    
    if (mod(n_valid, 2) == 0) then
        med = (sorted_data(n_valid/2) + sorted_data(n_valid/2 + 1)) / 2.0
    else
        med = sorted_data((n_valid + 1) / 2)
    end if

  end function median

  subroutine sort2(arr)
    !------------------------------------------------------------
    ! Subroutine: sort2
    ! Purpose   : Sorts an array of real numbers in ascending order using
    !             insertion sort.
    !
    ! Input/Output:
    !   arr - Array to be sorted.
    !------------------------------------------------------------
    implicit none
    real, intent(inout) :: arr(:)
    integer             :: i, j
    real                :: temp
    
    do i = 2, size(arr)
        temp = arr(i)
        j = i - 1
        do while (j >= 1 .and. arr(j) > temp)
            arr(j + 1) = arr(j)
            j = j - 1
        end do
        arr(j + 1) = temp
    end do
  end subroutine sort2

!------------------------------------------------------------
end module k_module