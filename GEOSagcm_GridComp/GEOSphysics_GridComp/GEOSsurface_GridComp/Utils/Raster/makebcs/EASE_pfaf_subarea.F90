#include "MAPL_ErrLog.h"

module EASE_pfaf_subareaMod

  ! Main purpose: Assigns a catchment‐tile index from Pfaf catchment definition files to each tile for global cylindrical EASE grid;
  !               writes information into EASEv[x]_tile2pfaf.nc4

  use MAPL
  use LogRectRasterizeMod, only: nPfaf=>SRTM_maxcat

  implicit none
  private

  public :: EASE_get_Pfaf_subarea
  public :: EASE_find_subs

contains 
  
  subroutine EASE_get_Pfaf_subarea(file_out, BCS_PATH, GridName, gfile)
    
    character(*), intent(in) :: file_out, BCS_PATH    
    character(*), intent(in) :: GridName
    character(*), intent(in) :: gfile                 ! name of (EASE) til file (minus file name extension)
    
    integer, parameter :: nlon=21600         ! Number of longitude grid cells in raster grid
    integer, parameter :: nlat=10800         ! Number of latitude  grid cells in raster grid

    integer, parameter :: nsub_max_init=9999
    
    ! Variable declarations:
    integer              :: id, xi, yi, i,j, flag, nmax,nmax2,ntot,it,ns_tot
    integer, allocatable :: nsub(:)                            ! Array to store the number of sub-areas for each catchment
    integer, allocatable :: xsub( :,:), ysub( :,:), isub( :,:)
    integer, allocatable :: xsub0(:,:), ysub0(:,:)
    ! xsub and ysub: 2D arrays to store mapped x and y indices for sub-catchments (not using subi_global in this code)
    real,    allocatable :: asub( :,:), asub0(:,:)             ! 2D array to store aggregated area for each sub-catchment
    
    real*8,  allocatable :: lon(  :), lat(  :)                 ! Arrays to hold longitude and latitude values
    integer, allocatable :: loni( :), lati( :)
    ! loni and lati: Arrays holding mapping indices from 1-minute resolution data files
    integer, allocatable :: catchind(:,:)                      ! 2D array holding catchment indices for each grid cell
    
    ! Declare allocatable arrays for catchment ID, parent catchment ID, and boundary coordinates
    integer, allocatable, dimension(:)   :: lati_tile, loni_tile
    real*8,  allocatable, dimension(:)   :: latc,lonc
    integer, allocatable, dimension(:,:) :: map_tile
    
    real,    allocatable, dimension(:)   :: area_cat, subtile_area(:)
    real,    allocatable, dimension(:,:) :: frac,frac_tile
    integer, allocatable, dimension(:)   :: nsub_tile
    integer, allocatable                 :: csub(:,:),flag2(:,:), tile_id(:), pfaf_index(:)

    type(NetCDF4_FileFormatter) :: formatter
    type (FileMetadata)         :: meta
    type(Variable)              :: v
    integer                     :: nc_ease, nr_ease
    real                        :: tmp_lati, tmp_loni
    integer                     :: ntile, npfaf0, nx, ny, type, pfaf_ind 
    
    ! Define file path for input routing data:
    character(len=256)   :: pfafData_file    

    pfafData_file = trim(BCS_PATH)//"/land/topo/v1/SRTM-TopoData/SRTM_PfafData.nc"
    call MAPL_ease_extent( trim(GridName), nc_ease, nr_ease)   
 
    ! Allocate arrays with the specified dimensions:
    allocate(catchind(nlon, nlat))
    allocate(lon(nlon), lat(nlat))
    allocate(loni(nlon), lati(nlat))
    allocate(nsub(nPfaf))
    
    call formatter%open(trim(pfafData_file), PFIO_READ)
    call formatter%get_var("latitude",   lat)
    call formatter%get_var("longitude",  lon)
    call formatter%get_var("CatchIndex", catchind)
    call formatter%close()

    ! for each raster grid cell longitude, get longitude index of overlying EASE grid cell
    do i = 1, nlon
       call MAPL_ease_convert( trim(GridName), 0.0, real(lon(i)), tmp_loni, tmp_lati)
       loni(i) = nint(tmp_loni) + 1
    enddo

    ! for each raster grid cell latitude,  get latitude  index of overlying EASE grid cell
    do i = 1, nlat
       call MAPL_ease_convert( trim(GridName), real(lat(i)), 0.0, tmp_loni, tmp_lati)
       lati(i) = nint(tmp_lati) + 1
    enddo
    
    ! allocate aggregation arrays
    allocate(xsub0(nsub_max_init, nPfaf), ysub0(nsub_max_init, nPfaf), asub0(nsub_max_init, nPfaf))

    ! collect aggregation info
    call EASE_Find_subs(catchind, loni, lati, lat, nsub, xsub0, ysub0, asub0)

    ! reduce array sizes as much as possible
    
    nmax = maxval(nsub)
    !print *,nmax
    allocate(xsub(nmax, nPfaf), ysub(nmax, nPfaf), asub(nmax, nPfaf))
    xsub=xsub0(1:nmax,:)
    ysub=ysub0(1:nmax,:)
    asub=asub0(1:nmax,:)
    deallocate(xsub0,ysub0,asub0)
    
    ! Open the catchment definition file for the EASE grid and header (ntot = total number of *land* tiles)
    ! NOTE: This approach works only when the routed runoff is from land tiles. 
    !       When we add the routing of runoff from landice tiles, this will have to be revisited.
    open(77, file="clsm/catchment.def");read(77, *) ntot

    ! get center lat/lon of EASE tiles from *.til file
    allocate(latc(ntot),lonc(ntot),lati_tile(ntot),loni_tile(ntot))
    open(10,file="til/"//trim(gfile)//".til", form="formatted", status="old", action='read')
    read(10,*) ntile, npfaf0, nx, ny
    do i=1,4
      read(10,*) 
    enddo
    do i=1,ntot                                                                ! read *land* tile center coords --> assumes land is first in til file!
      read(10,*) type, pfaf_ind, lonc(i), latc(i), loni_tile(i), lati_tile(i)  ! til file format here is specific to EASE grid!
    end do  
    close(10) 

    lati_tile = lati_tile + 1
    loni_tile = loni_tile + 1
    ! record into map_tile
    allocate(map_tile(nc_ease,nr_ease))
    map_tile=-9999
    do i =1, ntot
      map_tile(loni_tile(i),lati_tile(i)) = i
    enddo
    
    allocate(isub(nmax, nPfaf))
    isub=0
    ! Loop over each catchment
    do i = 1, nPfaf
      ! Loop over each potential sub-area within the current catchment
      do j = 1, nmax
        xi = xsub(j, i)                  ! Retrieve the x-coordinate for the sub-area
        yi = ysub(j, i)                  ! Retrieve the y-coordinate for the sub-area
        if (xi /= 0) then                ! Check if a valid sub-area exists (non-zero x-coordinate)
          isub(j, i) = map_tile(xi, yi)  ! Map the sub-area coordinates to a global tile index using map_tile
        endif
      enddo
    enddo

    allocate(area_cat(nPfaf),frac(nmax,nPfaf),nsub_tile(ntot),flag2(nmax,nPfaf))
    where(isub==-9999) asub=0.
    area_cat=0.
    do i=1,nPfaf
      area_cat(i)=sum(asub(:,i))
    enddo
    nsub_tile=0
    frac=0.
    do i=1,nPfaf
      do j=1,nmax
        if(isub(j,i)>0)then
          it=isub(j,i)
          nsub_tile(it)=nsub_tile(it)+1
          frac(j,i)=asub(j,i)/area_cat(i)
        endif
        if(isub(j,i)==-9999) nsub(i) = nsub(i) - 1
        if(isub(j,i)==0)exit
      enddo
    enddo

    nmax2=maxval(nsub_tile)
    allocate(csub(nmax2,ntot),frac_tile(nmax2,ntot))
    csub=0
    nsub_tile=0
    frac_tile=0.
    do i=1,nPfaf
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
    allocate(tile_id(ns_tot), pfaf_index(ns_tot), subtile_area(ns_tot))
    it = 0
    do i=1,ntot
      do j=1,nmax2
        if(csub(j,i)>0)then
          it = it+1
          tile_id(it) = i
          pfaf_index(it) = csub(j,i)
          subtile_area(it)  = frac_tile(j,i)*area_cat(pfaf_index(it))*1.e6
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
    call v%add_attribute('units', 'm+2')
    call v%add_attribute('long_name', 'area_of_subtile')
    call meta%add_variable('subtile_area', v)

    call formatter%create(file_out, mode=PFIO_CLOBBER)
    call formatter%write(meta)
    call formatter%put_var('tile_id', tile_id)
    call formatter%put_var('pfaf_index', pfaf_index)
    call formatter%put_var('subtile_area',  subtile_area)
    call formatter%close()

  end subroutine EASE_get_Pfaf_subarea

  ! ----------------------------------------------------------------------------------
  
  subroutine EASE_Find_subs(catchind, loni, lati, lat, nsub, xsub0, ysub0, asub0)
    
    integer,           intent(in)  :: catchind(:,:)
    integer,           intent(in)  :: loni(:)          ! lon index of overlying EASE grid cell;                 size = nlon
    integer,           intent(in)  :: lati(:)          ! lat index of overlying EASE grid cell;                 size = nlat
    real(kind=REAL64), intent(in)  :: lat(:)           ! lat of raster grid cell;                               size = nlat
    integer,           intent(out) :: nsub(:)          ! # raster grid cells that contribute to Pfaf catchment; size = nPfaf
    integer,           intent(out) :: xsub0(:,:)       ! lon index of overlying EASE grid cell;                 size = nsub_max_init x nPfaf
    integer,           intent(out) :: ysub0(:,:)       ! lat index of overlying EASE grid cell;                 size = nsub_max_init x nPfaf
    real,              intent(out) :: asub0(:,:)       ! area of raster grid cell;                              size = nsub_max_init x nPfaf

    integer            :: xi, yi, flag, id, i, nlon, nlat
    real(kind=REAL64)  :: cellarea, delta

    nlon  = size(loni)
    nlat  = size(lati)
    nsub  = 0 
    xsub0 = 0
    ysub0 = 0
    asub0 = 0.

    delta = 2*MAPL_PI_R8/nlon * MAPL_PI_R8/nlat                    ! pre-compute term for raster grid cell area
    
    ! Loop through all raster grid cells to aggregate cell areas by catchment and sub-area:
    do yi = 1, nlat
      cellarea = cos(lat(yi) * MAPL_DEGREES_TO_RADIANS_R8)*delta   ! area of raster grid cell for latitude yi [radians^2]
      cellarea = cellarea*MAPL_RADIUS/1.e3*MAPL_RADIUS/1.e3        ! convert to km^2
      do xi = 1, nlon
        if (catchind(xi, yi) >= 1) then
          ! The raster grid cell belongs to a catchment
          id = catchind(xi, yi)  ! Get the catchment id for the current cell
          flag = 0               ! Reset flag to indicate whether a matching sub-area is found      
          ! If the catchment already has one or more sub-areas, check for a matching sub-area:
          if (nsub(id) >= 1) then
            do i = 1, nsub(id)
              if (loni(xi) == xsub0(i, id) .and. lati(yi) == ysub0(i, id)) then
                flag = 1
                ! If a match is found, accumulate the cell area into the existing sub-area:
                asub0(i, id) = asub0(i, id) + cellarea
                exit  ! Exit the inner loop since a matching sub-area has been found
              endif
            end do
          endif
          ! If no matching sub-area was found, create a new sub-area:
          if (flag == 0) then
            nsub(           id) = nsub(id) + 1
            xsub0(nsub(id), id) = loni(xi)
            ysub0(nsub(id), id) = lati(yi)
            asub0(nsub(id), id) = cellarea
          endif
        endif
      end do
    end do
  end subroutine

end module EASE_pfaf_subareaMod

