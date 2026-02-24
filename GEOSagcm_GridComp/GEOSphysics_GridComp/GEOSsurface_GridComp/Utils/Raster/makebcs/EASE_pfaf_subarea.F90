#include "MAPL_ErrLog.h"

module EASE_pfaf_subareaMod

  ! Main purpose: Calculates info needed to map from EASE tile space to Pfafstetter catchment space;
  !               writes information into EASEv[x]_tile2pfaf.nc4.
  !
  ! Background:   The standard EASE tile space is created such that there is at most one tile per EASE grid cell.
  !               Each EASE tile is nominally assigned to the Pfafstetter catchment that covers the largest portion
  !               of the tile. Information about other Pfafstetter catchments that intersect with the EASE grid cell
  !               is lost in the process.
  !               For routing, we need to know about all of the Pfafstetter catchments that intersect with the tile,
  !               including the area (or fraction) of the intersect.
  !               Given an EASE grid resolution, this program calculates the full mapping information from the 1-arcmin
  !               raster file on which the Pfafstetter catchments are defined.
  !               Note that this is essentially equivalent to creating an EASE tile space *without* the restriction of 
  !               at most one tile per EASE grid cell. 

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
    integer              :: id, xi, yi, i,j, flag, nmax,nmax2,ntotland,it,ns_tot
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
    integer                     :: ntile, npfaf0, nx, ny, type, pfaf_ind, nlon30s, nlat30s 
    
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
    ! Open the catchment definition file for the EASE grid and header (ntotland = total number of *land* tiles)
    ! NOTE: This approach works only when the routed runoff is from land tiles. 
    !       When we add the routing of runoff from landice tiles, this will have to be revisited.
    open(77, file="clsm/catchment.def");read(77, *) ntotland
    nlon30s = 2*nlon; nlat30s = 2*nlat
    ! collect aggregation info
    call EASE_Find_subs(catchind, loni, lati, lat, ntotland, nlon30s, nlat30s, "rst/"//trim(gfile)//".rst", &
                        nsub, xsub0, ysub0, asub0)

    ! reduce array sizes as much as possible
    
    nmax = maxval(nsub)
    !print *,nmax
    allocate(xsub(nmax, nPfaf), ysub(nmax, nPfaf), asub(nmax, nPfaf))
    xsub=xsub0(1:nmax,:)
    ysub=ysub0(1:nmax,:)
    asub=asub0(1:nmax,:)
    deallocate(xsub0,ysub0,asub0)   
    ! get center lat/lon of EASE tiles from *.til file
    allocate(latc(ntotland),lonc(ntotland),lati_tile(ntotland),loni_tile(ntotland))

    open(10,file="til/"//trim(gfile)//".til", form="formatted", status="old", action='read')
    read(10,*) ntile, npfaf0, nx, ny
    do i=1,4
      read(10,*) 
    enddo
    do i=1,ntotland                                                                ! read *land* tile center coords --> assumes land is first in til file!
      read(10,*) type, pfaf_ind, lonc(i), latc(i), loni_tile(i), lati_tile(i)  ! til file format here is specific to EASE grid!
    end do  
    close(10) 

    lati_tile = lati_tile + 1
    loni_tile = loni_tile + 1
    ! record into map_tile
    allocate(map_tile(nc_ease,nr_ease))
    map_tile=-9999
    do i =1, ntotland
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

    allocate(area_cat(nPfaf),frac(nmax,nPfaf),nsub_tile(ntotland),flag2(nmax,nPfaf))
    where(isub==-9999) asub=0.
    area_cat=0.
    ! Calculate the area for each catchment
    do i=1,nPfaf
      area_cat(i)=sum(asub(:,i))
    enddo
    ! Calculate the fraction of each sub-area to the entire catchment
    nsub_tile=0
    frac=0.
    do i=1,nPfaf
      do j=1,nmax
        if(isub(j,i)>0)then
          it=isub(j,i)
          nsub_tile(it)=nsub_tile(it)+1   ! count on the number of sub-area within a tile
          frac(j,i)=asub(j,i)/area_cat(i) ! fraction of each sub-area to the entire catchment, stored in (sub-area, catchment) array
        endif
        if(isub(j,i)==-9999) nsub(i) = nsub(i) - 1
        if(isub(j,i)==0)exit
      enddo
    enddo

    nmax2=maxval(nsub_tile)
    allocate(csub(nmax2,ntotland),frac_tile(nmax2,ntotland))
    csub=0      
    nsub_tile=0  
    frac_tile=0. 
    do i=1,nPfaf
      do j=1,nmax
        if(isub(j,i)>0)then
          it=isub(j,i)
          nsub_tile(it)=nsub_tile(it)+1          ! count on the number of sub-area within a tile
          csub(nsub_tile(it),it)=i               ! to store the catchment ids (Pfaf ids) in a tile
          frac_tile(nsub_tile(it),it)=frac(j,i)  ! frac_tile is the fraction of each sub-area to the entire catchment, stored in (sub-area, tile) array
        endif
        if(isub(j,i)==0)exit
      enddo
    enddo 
    ns_tot=sum(nsub_tile)                        ! the number of sub-area
    allocate(tile_id(ns_tot), pfaf_index(ns_tot), subtile_area(ns_tot))
    it = 0
    do i=1,ntotland
      do j=1,nmax2
        if(csub(j,i)>0)then
          it = it+1
          tile_id(it) = i                                                  ! tile_id of each sub-area  
          pfaf_index(it) = csub(j,i)                                       ! pfaf_id of each sub-area  
          subtile_area(it)  = frac_tile(j,i)*area_cat(pfaf_index(it))*1.e6 ! the area of each sub-area  
        endif
        if(csub(j,i)==0)exit
      enddo
    enddo 

    ! The length of the vectors in EASEv[x]_tile2pfaf.nc4 is ns_tot, which is the
    ! total number of EASE land tiles that we would have if the EASE tile space was
    ! generated without the restriction of at most one tile per EASE grid cell.
    
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
  
  subroutine EASE_Find_subs(catchind, loni, lati, lat, nland, nlon30s, nlat30s, rstfile, nsub, xsub0, ysub0, asub0)
    
    integer,           intent(in)  :: catchind(:,:)
    integer,           intent(in)  :: loni(:)          ! lon index of overlying EASE grid cell;                 size = nlon
    integer,           intent(in)  :: lati(:)          ! lat index of overlying EASE grid cell;                 size = nlat
    real(kind=REAL64), intent(in)  :: lat(:)           ! lat of raster grid cell;                               size = nlat
    integer,           intent(in)  :: nland            ! number of land tiles
    integer,           intent(in)  :: nlon30s          ! number of lon coordinate in 30-s grid
    integer,           intent(in)  :: nlat30s          ! number of lat coordinate in 30-s grid
    character(*),      intent(in)  :: rstfile          ! path of *.rst file
    integer,           intent(out) :: nsub(:)          ! # raster grid cells that contribute to Pfaf catchment; size = nPfaf
    integer,           intent(out) :: xsub0(:,:)       ! lon index of overlying EASE grid cell;                 size = nsub_max_init x nPfaf
    integer,           intent(out) :: ysub0(:,:)       ! lat index of overlying EASE grid cell;                 size = nsub_max_init x nPfaf
    real,              intent(out) :: asub0(:,:)       ! area of raster grid cell;                              size = nsub_max_init x nPfaf

    integer :: xi, yi, flag, id, i, j, nlon, nlat
    real(kind=REAL64)             :: cellarea, delta, area30s(2,2)
    real(kind=REAL64),allocatable :: lat30s(:), cellarea30s(:)
    integer,allocatable           :: rst(:,:)

    nlon  = size(loni)
    nlat  = size(lati)
    nsub  = 0 
    xsub0 = 0
    ysub0 = 0
    asub0 = 0.

    allocate(rst(nlon30s,nlat30s),lat30s(nlat30s),cellarea30s(nlat30s)) 
    delta=180./nlat30s
    ! Calculate the lat of 30-s grid.
    do j=1,nlat  
      if(lat(1)<0.d0)then
        lat30s(2*j-1)=lat(j)-delta/2.d0
        lat30s(2*j)=lat(j)+delta/2.d0
      else
        lat30s(2*j-1)=lat(j)+delta/2.d0
        lat30s(2*j)=lat(j)-delta/2.d0
      endif      
    enddo
    delta = 2*MAPL_PI_R8/nlon30s * MAPL_PI_R8/nlat30s
    open(20,file=trim(rstfile),form="unformatted",status="old")
    ! Read rst file and calculate the area of 30-s gridcell.
    do j=1,nlat30s
      read(20) rst(:,j)
      cellarea30s(j) = cos(lat30s(j) * MAPL_DEGREES_TO_RADIANS_R8)*delta * MAPL_RADIUS/1.e3*MAPL_RADIUS/1.e3 !km^2
    enddo
    where (rst<1.or.rst>nland) rst=0 !For no-land tiles in the rst file, set rst to 0.

    delta = 2*MAPL_PI_R8/nlon * MAPL_PI_R8/nlat                    ! pre-compute term for raster grid cell area
    ! Loop through all raster grid cells to aggregate cell areas by catchment and sub-area:
    do yi = 1, nlat
      cellarea = cos(lat(yi) * MAPL_DEGREES_TO_RADIANS_R8)*delta   ! area of raster grid cell for latitude yi [radians^2]
      cellarea = cellarea*MAPL_RADIUS/1.e3*MAPL_RADIUS/1.e3        ! convert to km^2
      do xi = 1, nlon
        if (catchind(xi, yi) >= 1) then
          area30s(:,1)=cellarea30s(2*yi-1) !area30s(2,2) is the area of 30s cells corresponding to the (xi,yi) of 1-min cell.
          area30s(:,2)=cellarea30s(2*yi)
          !set land area to 0 in the area30s(2,2), so the sum(area30s(2,2)) is the no-land area in the (xi,yi) of 1-min cell.
          where(rst(2*xi-1:2*xi,2*yi-1:2*yi)/=0) area30s=0.d0          
          ! The raster grid cell belongs to a catchment
          id = catchind(xi, yi)  ! Get the catchment id for the current cell
          flag = 0               ! Reset flag to indicate whether a matching sub-area is found      
          ! If the catchment already has one or more sub-areas, check for a matching sub-area:
          if (nsub(id) >= 1) then
            do i = 1, nsub(id)
              if (loni(xi) == xsub0(i, id) .and. lati(yi) == ysub0(i, id)) then
                flag = 1
                ! If a match is found, accumulate the cell area into the existing sub-area:
                ! subtract the no-land area sum(area30s) here, so only land area is summed to sub-area
                asub0(i, id) = asub0(i, id) + max(0.,cellarea - sum(area30s)) 
                exit  ! Exit the inner loop since a matching sub-area has been found
              endif
            end do
          endif
          ! If no matching sub-area was found, create a new sub-area:
          if (flag == 0) then
            nsub(           id) = nsub(id) + 1
            xsub0(nsub(id), id) = loni(xi)
            ysub0(nsub(id), id) = lati(yi)
            ! subtract the no-land area sum(area30s) here, so only land area is summed to sub-area
            asub0(nsub(id), id) = max(0.,cellarea - sum(area30s))
          endif
        endif
      end do
    end do

    deallocate(rst,lat30s,cellarea30s)

  end subroutine EASE_Find_subs

end module EASE_pfaf_subareaMod

