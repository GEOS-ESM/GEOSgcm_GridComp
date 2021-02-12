
module easeV1_conv
  
  ! ==========================================================================
  !
  ! easeV1_conv.F90 - FORTRAN routines for conversion of azimuthal 
  !                   equal area and equal area cylindrical grid coordinates
  ! 
  ! 30-Jan-1992 H.Maybee
  ! 20-Mar-1992 Ken Knowles  303-492-0644  knowles@kryos.colorado.edu
  ! 16-Dec-1993 MJ Brodzik   303-492-8263  brodzik@jokull.colorado.edu
  !              Copied from nsmconv.f, changed resolutions from 
  !              40-20-10 km to 25-12.5 km
  ! 21-Dec-1993 MJ Brodzik   303-492-8263  brodzik@jokull.colorado.edu
  !              Fixed sign of Southern latitudes in ease_inverse.
  ! 12-Sep-1994 David Hoogstrate 303-492-4116 hoogstra@jokull.colorado.edu
  ! 	       Changed grid cell size. Changed "c","f" to "l","h"
  ! 25-Oct-1994 David Hoogstrate 303-492-4116 hoogstra@jokull.colorado.edu
  ! 	       Changed row size from 587 to 586 for Mercator projection
  ! 11-May-2011 reichle: Changed "smap" to "easeV1".  
  !                      Added SSM/I and AMSR-E "M25" grid.
  !                      So far ONLY for cylindrical grids.
  !                      Converted from *.f to *.F90 module
  ! 
  ! $Log$
  ! Revision 1.1.4.2  2012/04/18 18:04:57  smahanam
  ! updated on 4-18-2012 to process global data on native grids and time steps
  !
  ! Revision 1.1.2.1  2012-01-27 20:22:10  smahanam
  ! Added Richard's equation solver, replaced SMAP f77 files with F90 module
  !
  ! Revision 1.1  2011-05-11 21:58:46  rreichle
  !
  ! Adding utilities to map between EASE grids and lat/lon coordinates.
  !
  ! Revision 1.3  1994/11/01 23:40:43  brodzik
  ! Replaced all references to 'ease' with 'smap'
  ! Replaced all references to 'smap' with 'easeV1' -- reichle
  ! 
  ! ==========================================================================

  implicit none
  
  private

  public :: easeV1_convert
  public :: easeV1_inverse


  ! ***NEVER*** change these constants to GEOS-5 MAPL constants!!!!
    
  ! radius of the earth (km), authalic sphere based on International datum 
  
  real*8, parameter :: RE_km = 6371.228
  
  ! scale factor for standard paralles at +/-30.00 degrees
  
  real*8, parameter :: COS_PHI1 = .866025403
  
  real*8, parameter :: PI = 3.141592653589793
  

contains  
  
  ! *******************************************************************
  
  subroutine easeV1_convert (grid, lat, lon, r, s)
    
    ! convert geographic coordinates (spherical earth) to 
    ! azimuthal equal area or equal area cylindrical grid coordinates
    ! 
    ! status = easeV1_convert (grid, lat, lon, r, s)
    ! 
    ! input : grid - projection name '[M][xx]'
    !            where xx = approximate resolution [km]
    !               ie xx = "01", "03", "09", "36"       (SMAP)
    !               or xx = "12", "25"                   (SSM/I, AMSR-E)
    ! 	    lat, lon = geo. coords. (decimal degrees)
    ! 
    ! output: r, s - column, row coordinates
    ! 
    ! result: status = 0 indicates normal successful completion
    ! 		-1 indicates error status (point not on grid)
    ! 
    ! --------------------------------------------------------------------------
        
    character*(*), intent(in)  :: grid
    real,          intent(in)  :: lat, lon
    real,          intent(out) :: r, s

    ! local variables
    
    integer :: cols, rows, scale
    real*8  :: Rg, phi, lam, rho, CELL_km, r0, s0
    
    ! ---------------------------------------------------------------------
    
    call easeV1_get_params( grid, CELL_km, cols, rows, r0, s0, Rg )
    
    phi = lat*PI/180.   ! convert from degree to radians
    lam = lon*PI/180.   ! convert from degree to radians
    
    if (grid(1:1).eq.'N') then
       rho = 2 * Rg * sin(PI/4. - phi/2.)
       r = r0 + rho * sin(lam)
       s = s0 + rho * cos(lam)
       
    else if (grid(1:1).eq.'S') then
       rho = 2 * Rg * cos(PI/4. - phi/2.)
       r = r0 + rho * sin(lam)
       s = s0 - rho * cos(lam)
       
    else if (grid(1:1).eq.'M') then
       r = r0 + Rg * lam * COS_PHI1
       s = s0 - Rg * sin(phi) / COS_PHI1
       
    endif
        
  end subroutine easeV1_convert
  
  ! *******************************************************************
  
  subroutine easeV1_inverse (grid, r, s, lat, lon)
    
    ! convert azimuthal equal area or equal area cylindrical 
    ! grid coordinates to geographic coordinates (spherical earth)
    ! 
    ! status = easeV1_inverse (grid, r, s, lat, lon)
    ! 
    ! input : grid - projection name '[M][xx]'
    !            where xx = approximate resolution [km]
    !               ie xx = "01", "03", "09", "36"       (SMAP)
    !               or xx = "12", "25"                   (SSM/I, AMSR-E)
    ! 	    r, s - column, row coordinates
    ! 
    ! output: lat, lon = geo. coords. (decimal degrees)
    ! 
    ! result: status = 0 indicates normal successful completion
    ! 		-1 indicates error status (point not on grid)
    ! 
    ! --------------------------------------------------------------------------

    character*(*), intent(in)  :: grid
    real,          intent(in)  :: r, s
    real,          intent(out) :: lat, lon

    ! local variables
    
    integer :: cols, rows, scale
    real*8    :: Rg, phi, lam, rho, CELL_km, r0, s0
    real*8    :: gamma, beta, epsilon, x, y, c
    real*8    :: sinphi1, cosphi1

    ! ---------------------------------------------------------------------
    
    call easeV1_get_params( grid, CELL_km, cols, rows, r0, s0, Rg )
        
    x = r - r0
    y = -(s - s0)
    
    if ((grid(1:1).eq.'N').or.(grid(1:1).eq.'S')) then 
       rho = sqrt(x*x + y*y)
       if (rho.eq.0.0) then
          if (grid(1:1).eq.'N') lat = 90.0 
          if (grid(1:1).eq.'S') lat = -90.0 
          lon = 0.0
       else
          if (grid(1:1).eq.'N') then
             sinphi1 = sin(PI/2.)
             cosphi1 = cos(PI/2.)
             if (y.eq.0.) then
                if (r.le.r0) lam = -PI/2.
                if (r.gt.r0) lam = PI/2.
             else
                lam = atan2(x,-y)
             endif
          else if (grid(1:1).eq.'S') then
             sinphi1 = sin(-PI/2.)
             cosphi1 = cos(-PI/2.)
             if (y.eq.0.) then
                if (r.le.r0) lam = -PI/2.
                if (r.gt.r0) lam = PI/2.
             else
                lam = atan2(x,y)
             endif
          endif
          gamma = rho/(2 * Rg)
          if (abs(gamma) .gt. 1.) return
          c = 2 * asin(gamma)
          beta = cos(c) * sinphi1 + y * sin(c) * (cosphi1/rho)
          if (abs(beta).gt.1.) return
          phi = asin(beta)
          lat = phi*180./PI   ! convert from radians to degree
          lon = lam*180./PI   ! convert from radians to degree
       endif
       
    else if (grid(1:1).eq.'M') then
       
       ! 	  allow .5 cell tolerance in arcsin function
       ! 	  so that grid coordinates which are less than .5 cells
       ! 	  above 90.00N or below 90.00S are given a lat of 90.00
       
       epsilon = 1 + 0.5/Rg
       beta = y*COS_PHI1/Rg
       if (abs(beta).gt.epsilon) return
       if (beta.le.-1.) then
          phi = -PI/2.
       else if (beta.ge.1.) then
          phi = PI/2.
       else
          phi = asin(beta)
       endif
       lam = x/COS_PHI1/Rg
       lat = phi*180./PI   ! convert from radians to degree
       lon = lam*180./PI   ! convert from radians to degree
    endif
        
  end subroutine easeV1_inverse
  
  ! *******************************************************************

  subroutine easeV1_get_params( grid, CELL_km, cols, rows, r0, s0, Rg )
    
    implicit none
    
    character*(*), intent(in)  :: grid
    real*8,        intent(out) :: CELL_km, r0, s0, Rg
    integer,       intent(out) :: cols, rows
    
    ! --------------------------------------------------------
    !
    ! r0,s0 are defined such that cells at all scales have 
    ! coincident center points
    ! 
    !c        r0 = (cols-1)/2. * scale
    !c        s0 = (rows-1)/2. * scale
    !
    ! --------------------------------------------------------
    
    if ((grid(1:1).eq.'N').or.(grid(1:1).eq.'S')) then

       print *,'Polar projections not implemented yet'
       stop

    else if (grid(1:1).eq.'M') then

       if      (grid .eq. 'M36') then ! SMAP 36 km grid
          CELL_km = 36.00040279063    ! nominal cell size in kilometers
          cols = 963
          rows = 408
          r0 = 481.0
          s0 = 203.5
       
       else if (grid .eq. 'M25') then ! SSM/I, AMSR-E 25 km grid
          CELL_km = 25.067525         ! nominal cell size in kilometers
          cols = 1383
          rows = 586
          r0 = 691.0
          s0 = 292.5
       
       else if (grid .eq. 'M09') then ! SMAP  9 km grid
          CELL_km = 9.00010069766     ! nominal cell size in kilometers
          cols = 3852
          rows = 1632
          r0 = 1925.5
          s0 = 815.5
       
       else if (grid .eq. 'M03') then ! SMAP  3 km grid
          CELL_km = 3.00003356589     ! nominal cell size in kilometers
          cols = 11556
          rows = 4896
          r0 = 5777.5
          s0 = 2447.5
       
       else if (grid .eq. 'M01') then ! SMAP  1 km grid
          CELL_km = 1.00001118863     ! nominal cell size in kilometers
          cols = 34668
          rows = 14688
          r0 = 17333.5
          s0 = 7343.5
       
       else
       
          print *,'easeV1_convert: unknown resolution: ',grid
          stop
       
       endif

    else
       
       print *, 'easeV1_convert: unknown projection: ', grid
       stop
       
    endif
        
    Rg = RE_km/CELL_km
    
  end subroutine easeV1_get_params
  
  ! *******************************************************************
  
end module easeV1_conv

! =============================== EOF =================================

