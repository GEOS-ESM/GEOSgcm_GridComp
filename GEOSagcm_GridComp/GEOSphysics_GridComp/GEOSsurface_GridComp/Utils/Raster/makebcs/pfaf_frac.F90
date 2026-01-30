#include "MAPL_ErrLog.h"

module pfaf_fracMod

  ! Main purpose: Assigns a catchment‐tile index from Pfaf catchment definition files to each tile for global cylindrical EASE grid;
  !               writes information into EASEv[x]_tile2pfaf.nc4

  use MAPL
  use LogRectRasterizeMod, only: nc=>SRTM_maxcat

  implicit none
  private

  public :: get_Pfaf_frac

contains 
  
  subroutine get_Pfaf_frac(file_out, BCS_PATH, GridName)
    
    character(*), intent(in) :: file_out, BCS_PATH    
    character(*), intent(in) :: GridName
    
    integer, parameter :: nlon=21600         ! Number of longitude grid cells in raster grid
    integer, parameter :: nlat=10800         ! Number of latitude  grid cells in raster grid

    ! Variable declarations:
    integer              :: id, xi, yi, i,j, flag, nmax,nmax2,ntot,it,ns_tot
    integer, allocatable :: nsub(:)                            ! Array to store the number of sub-areas for each catchment
    integer, allocatable :: xsub(:,:), ysub(:,:), isub(:,:)
    ! xsub and ysub: 2D arrays to store mapped x and y indices for sub-catchments (not using subi_global in this code)
    real,    allocatable :: asub(:,:)                          ! 2D array to store aggregated area for each sub-catchment
    
    real*8,  allocatable :: lon(  :), lat(  :)                 ! Arrays to hold longitude and latitude values from the NetCDF file
    real*8,  allocatable :: lons( :), lats( :)                 ! Arrays to hold longitude and latitude values from the NetCDF file
    integer, allocatable :: loni( :), lati( :)
    integer, allocatable :: loni2(:), lati2(:)
    ! loni and lati: Arrays holding mapping indices from 1-minute resolution data files
    integer, allocatable :: catchind(:,:)                      ! 2D array holding catchment indices for each grid cell
    real,    allocatable :: cellarea(:,:)                      ! 2D array containing the area of each grid cell
    
    ! Declare allocatable arrays for catchment ID, parent catchment ID, and boundary coordinates
    integer, allocatable, dimension(:)   :: tid, catid, lati_tile, loni_tile
    real*8,  allocatable, dimension(:)   :: lon_left, lon_right, lat_bottom, lat_top, latc,lonc
    integer, allocatable, dimension(:,:) :: map_tile
    
    real,    allocatable, dimension(:)   :: area_cat, pfaf_frac(:)
    real,    allocatable, dimension(:,:) :: frac,frac_tile
    integer, allocatable, dimension(:)   :: nsub_tile
    integer, allocatable                 :: csub(:,:),flag2(:,:), tile_id(:), pfaf_index(:)

    type(NetCDF4_FileFormatter) :: formatter
    type (FileMetadata)         :: meta
    type(Variable)              :: v
    integer                     :: nc_ease, nr_ease
    real                        :: tmp_lat, tmp_lon
    
    ! Define file path for input routing data:
    character(len=256)   :: pfafData_file    
    character(len=256)   :: cellarea_file

    pfafData_file = trim(BCS_PATH)//"/land/topo/v1/SRTM-TopoData/SRTM_PfafData.nc"
    cellarea_file = trim(BCS_PATH)//"/route/routing_model/v1/cellarea.nc"
    call MAPL_ease_extent( trim(GridName), nc_ease, nr_ease)   
 
    ! Allocate arrays with the specified dimensions:
    allocate(catchind(nlon, nlat), cellarea(nlon, nlat))
    allocate(lon(nlon), lat(nlat))
    allocate(loni(nlon), lati(nlat),loni2(nlon), lati2(nlat))
    allocate(lons(nc_ease),lats(nr_ease))
    allocate(nsub(nc))
    
    call formatter%open(trim(pfafData_file), PFIO_READ)
    call formatter%get_var("latitude", lat)
    call formatter%get_var("longitude", lon)
    call formatter%get_var("CatchIndex", catchind)
    call formatter%close()
    
    do i = 1, nc_ease
       call MAPL_ease_inverse( trim(GridName), real(i-1), 0.0, tmp_lat, tmp_lon)
       lons(i) = tmp_lon
    enddo
    do j = 1, nr_ease
       call MAPL_ease_inverse( trim(GridName), 0.0, real(j-1), tmp_lat, tmp_lon)
       lats(j) = tmp_lat
    enddo
    call nearest_index_vector(lat, lats, lati)
    call nearest_index_vector(lon, lons, loni)
    
    call formatter%open(trim(cellarea_file), PFIO_READ)
    call formatter%get_var("data", cellarea)
    call formatter%close()
    cellarea = cellarea / 1.e6  ! Convert cell area (e.g., from m^2 to km^2)
    ! Initialize aggregation arrays to zero:
    allocate(xsub(9999, nc), ysub(9999, nc))
    nsub = 0
    ! Loop over all raster grid cells to aggregate cell areas by catchment and sub-area:
    do xi = 1, nlon
      do yi = 1, nlat
        if (catchind(xi, yi) >= 1) then
          ! The raster grid cell belongs to a catchment
          id = catchind(xi, yi)  ! Get the catchment id for the current cell
          flag = 0              ! Reset flag to indicate whether a matching sub-area is found      
          ! If the catchment already has one or more sub-areas, check for a matching sub-area:
          if (nsub(id) >= 1) then
            do i = 1, nsub(id)
              if (loni(xi) == xsub(i, id) .and. lati(yi) == ysub(i, id)) then
                flag = 1
                exit  ! Exit the inner loop since a matching sub-area has been found
              endif
            end do
          endif      
          ! If no matching sub-area was found, create a new sub-area:
          if (flag == 0) then
            nsub(id) = nsub(id) + 1
            xsub(nsub(id), id) = loni(xi)
            ysub(nsub(id), id) = lati(yi)        
          endif
        endif
      end do
    end do
    nmax = maxval(nsub)
    !print *,nmax
    deallocate(xsub,ysub)
    allocate(xsub(nmax, nc), ysub(nmax, nc), asub(nmax, nc))
    ! Initialize aggregation arrays to zero:
    nsub = 0;xsub = 0;ysub = 0;asub = 0.
    ! Loop over all raster grid cells to aggregate cell areas by catchment and sub-area:
    do xi = 1, nlon
      do yi = 1, nlat
        if (catchind(xi, yi) >= 1) then
          ! The raster grid cell belongs to a catchment
          id = catchind(xi, yi)  ! Get the catchment id for the current cell
          flag = 0              ! Reset flag to indicate whether a matching sub-area is found      
          ! If the catchment already has one or more sub-areas, check for a matching sub-area:
          if (nsub(id) >= 1) then
            do i = 1, nsub(id)
              if (loni(xi) == xsub(i, id) .and. lati(yi) == ysub(i, id)) then
                flag = 1
                ! If a match is found, accumulate the cell area into the existing sub-area:
                asub(i, id) = asub(i, id) + cellarea(xi, yi)
                exit  ! Exit the inner loop since a matching sub-area has been found
              endif
            end do
          endif      
          ! If no matching sub-area was found, create a new sub-area:
          if (flag == 0) then
            nsub(id) = nsub(id) + 1
            xsub(nsub(id), id) = loni(xi)
            ysub(nsub(id), id) = lati(yi)
            asub(nsub(id), id) = cellarea(xi, yi)
          endif
        endif
      end do
    end do
    
    ! Open the catchment definition file for the EASE grid and read the total number of tiles (header)
    open(77, file="clsm/catchment.def");read(77, *) ntot
    ! Allocate arrays with size nt
    allocate(tid(ntot), catid(ntot), lon_left(ntot), lon_right(ntot), lat_bottom(ntot), lat_top(ntot),latc(ntot),lonc(ntot))
    allocate(lati_tile(ntot),loni_tile(ntot),map_tile(nc_ease,nr_ease))
    ! Loop over each catchment and read: id, catchment id, left/right longitudes, bottom/top latitudes
    do i = 1, ntot
      read(77, *) tid(i), catid(i), lon_left(i), lon_right(i), lat_bottom(i), lat_top(i)
    end do
    latc = (lat_bottom + lat_top) / 2.
    lonc = (lon_left + lon_right) / 2.
    call nearest_index_vector(latc,lats,lati_tile)
    call nearest_index_vector(lonc,lons,loni_tile)
    map_tile=-9999
    do i =1, ntot
      map_tile(loni_tile(i),lati_tile(i)) = i
    enddo
    
    allocate(isub(nmax, nc))
    isub=0
    ! Loop over each catchment
    do i = 1, nc
      ! Loop over each potential sub-area within the current catchment
      do j = 1, nmax
        xi = xsub(j, i)  ! Retrieve the x-coordinate for the sub-area
        yi = ysub(j, i)  ! Retrieve the y-coordinate for the sub-area
        if (xi /= 0) then  ! Check if a valid sub-area exists (non-zero x-coordinate)
          isub(j, i) = map_tile(xi, yi)  ! Map the sub-area coordinates to a global tile index using map_tile
        endif
      enddo
    enddo
    
    allocate(area_cat(nc),frac(nmax,nc),nsub_tile(ntot),flag2(nmax,nc))
    where(isub==-9999) asub=0.
    area_cat=0.
    do i=1,nc
      area_cat(i)=sum(asub(:,i))
    enddo
    nsub_tile=0
    frac=0.
    do i=1,nc
      do j=1,nmax
        if(isub(j,i)>0)then
          it=isub(j,i)
          nsub_tile(it)=nsub_tile(it)+1
          frac(j,i)=asub(j,i)/area_cat(i)
        endif
        if(isub(j,i)==0)exit
      enddo
    enddo
    nmax2=maxval(nsub_tile)
    allocate(csub(nmax2,ntot),frac_tile(nmax2,ntot))
    csub=0
    nsub_tile=0
    frac_tile=0.
    do i=1,nc
      do j=1,nmax
        if(isub(j,i)>0)then
          it=isub(j,i)
          nsub_tile(it)=nsub_tile(it)+1
          csub(nsub_tile(it),it)=i
          frac_tile(nsub_tile(it),it)=frac(j,i)
        endif
        if(isub(j,i)==0)exit
      enddo
    enddo
    ns_tot=sum(nsub_tile)
    allocate(tile_id(ns_tot), pfaf_index(ns_tot), pfaf_frac(ns_tot))
    it = 0
    do i=1,ntot
      do j=1,nmax2
        if(csub(j,i)>0)then
          it = it+1
          tile_id(it) = i
          pfaf_index(it) = csub(j,i)
          pfaf_frac(it)  = frac_tile(j,i)
        endif
        if(csub(j,i)==0)exit
      enddo
    enddo

    call meta%add_dimension('tile', ns_tot)

    v = Variable(type=PFIO_INT32, dimensions='tile')
    call v%add_attribute('units', '1')
    call v%add_attribute('long_name', 'tile_id')
    call meta%add_variable('tile_id', v)

    v = Variable(type=pFio_INT32, dimensions='tile')
    call v%add_attribute('units', '1')
    call v%add_attribute('long_name', 'pfaf_index')
    call meta%add_variable('pfaf_index', v)

    v = Variable(type=pFio_REAL32, dimensions='tile')
    call v%add_attribute('units', '1')
    call v%add_attribute('long_name', 'area_fraction_of_catchment')
    call meta%add_variable('pfaf_frac', v)

    call formatter%create(file_out, mode=PFIO_CLOBBER)
    call formatter%write(meta)
    call formatter%put_var('tile_id', tile_id)
    call formatter%put_var('pfaf_index', pfaf_index)
    call formatter%put_var('pfaf_frac',  pfaf_frac)
    call formatter%close()

  end subroutine get_Pfaf_frac

  subroutine nearest_index_vector(candidates, targets, idx)
    ! For each targets(k), find argmin_j |candidates(j) - targets(k)|
    real*8, intent(in)  :: targets(:)
    real*8, intent(in)  :: candidates(:)
    integer, intent(out) :: idx(size(candidates))   ! 1-based indices
    integer :: k, j, nT, nC, best_j
    real*8 :: best_d, d

    nT = size(targets)
    nC = size(candidates)

    do k = 1, nC
       best_d = huge(1.0d0)
       best_j = 1
       do j = 1, nT
          d = abs(candidates(k) - targets(j))
          if (d < best_d) then
             best_d = d
             best_j = j
          end if
       end do
       idx(k) = best_j    ! already 1-based
    end do
  end subroutine nearest_index_vector

end module pfaf_fracMod
