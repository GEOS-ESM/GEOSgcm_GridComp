
module CubedSphere_GridMod

#define r8 kind=8

  implicit none

  real(r8) :: pi = 3.14159265358979323846 

  public Get_CubedSphere_Grid
  private

#ifdef EIGHT_BYTE
 integer, parameter:: f_p = 8!selected_real_kind(15)   ! same as 12 on Altix
#else
! Higher precisions for grid geometrical factors:
 integer, parameter:: f_p = 8!selected_real_kind(20)
#endif

contains
!!!!!!!!!!!!!!!%%%%%%%%%%%%%%%%%%%%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine Get_CubedSphere_Grid(npx, npy, xlocs, ylocs, grid_type, shift_west)
 integer, intent(in)   :: npx, npy
 real(r8), intent(out) :: xlocs(npx,npy), ylocs(npx,npy)
 integer, intent(in)   :: grid_type
 logical, optional, intent(in) :: shift_west

! ErrLog variables
!-----------------

 integer                      :: STATUS
 integer                      :: RC

! Local variables
!-----------------
  integer                       :: isg, ieg
  integer                       :: jsg, jeg
  integer                       :: npts
  integer                       :: ntiles=6
  integer                       :: ndims=2
  integer                       :: I, J, N
  integer                       :: IG, JG
  integer                       :: myTile

  real(r8), allocatable :: xs(:,:)
  real(r8), allocatable :: ys(:,:)
  real(r8), allocatable :: grid_global(:,:,:,:)
  real(r8), pointer     :: lons(:,:)
  real(r8), pointer     :: lats(:,:)
  real(r8)              :: alocs(2)

  logical doShiftWest

 isg = 1
 ieg = npx
 jsg = 1
 jeg = npy

 npts = (npy/ntiles)
 if (npts /= npx) then
    print*, 'Error npts /= npx', npts, npx
    stop
 endif

 allocate( xs(npts,npts) )
 allocate( ys(npts,npts) )
 allocate( grid_global(npts,npts,ndims,ntiles) )

  if(grid_type==0) call gnomonic_ed(  npts-1, xs, ys)
  if(grid_type==1) call gnomonic_dist(npts-1, xs, ys)
  if(grid_type==2) call gnomonic_angl(npts-1, xs, ys)
  call symm_ed(npts-1, xs, ys)
  do j=1,npts
     do i=1,npts
        grid_global(i,j,1,1) = xs(i,j) - pi
        grid_global(i,j,2,1) = ys(i,j)
     enddo
  enddo
! mirror_grid assumes that the tile=1 is centered on equator and greenwich meridian Lon[-pi,pi] 
  call mirror_grid(grid_global, 0, npts, npts, 2, 6)
! Cell Vertices
  doShiftWest = .false.
  if (present(shift_west)) doShiftWest = shift_west

  if (doShiftWest) then 
    do myTile=1,ntiles
    do jg=jsg,jeg/ntiles
       do ig=isg,ieg
          i=ig
          j=jg

! This will result in a corner close to east coast of China avoiding mountains in Japan
          grid_global(i,j,1,myTile) = grid_global(i,j,1,myTile) - pi/18.
          if ( grid_global(i,j,1,myTile) < 0. )              &
               grid_global(i,j,1,myTile) = grid_global(i,j,1,myTile) + 2.*pi
          if (ABS(grid_global(i,j,1,myTile)) < 1.e-10) grid_global(i,j,1,myTile) = 0.0
          if (ABS(grid_global(i,j,2,myTile)) < 1.e-10) grid_global(i,j,2,myTile) = 0.0

          alocs(1) = grid_global(i,j,1,myTile)
          alocs(2) = grid_global(i,j,2,myTile)
          i=ig-isg+1
          j=jg+(myTile-1)*npts
          xlocs(i,j) = alocs(1)*180.0/pi
          ylocs(i,j) = alocs(2)*180.0/pi
       enddo
    enddo
    enddo
  else
    do myTile=1,ntiles
    do jg=jsg,jeg/ntiles
       do ig=isg,ieg
          i=ig
          j=jg
          if ( grid_global(i,j,1,myTile) < 0. )              &
               grid_global(i,j,1,myTile) = grid_global(i,j,1,myTile) + 2.*pi
          if (ABS(grid_global(i,j,1,myTile)) < 1.e-10) grid_global(i,j,1,myTile) = 0.0
          if (ABS(grid_global(i,j,2,myTile)) < 1.e-10) grid_global(i,j,2,myTile) = 0.0
          alocs(1) = grid_global(i,j,1,myTile)
          alocs(2) = grid_global(i,j,2,myTile)
          i=ig-isg+1
          j=jg+(myTile-1)*npts
          xlocs(i,j) = alocs(1)*180.0/pi
          ylocs(i,j) = alocs(2)*180.0/pi
       enddo
    enddo
    enddo
  endif

  deallocate( xs )
  deallocate( ys )
  deallocate( grid_global )

 end subroutine Get_CubedSphere_Grid


 subroutine gnomonic_ed(im, lamda, theta)
!-----------------------------------------------------
! Equal distance along the 4 edges of the cubed sphere
!-----------------------------------------------------
! Properties: 
!            * defined by intersections of great circles
!            * max(dx,dy; global) / min(dx,dy; global) = sqrt(2) = 1.4142
!            * Max(aspect ratio) = 1.06089
!            * the N-S coordinate curves are const longitude on the 4 faces with equator 
! For C2000: (dx_min, dx_max) = (3.921, 5.545)    in km unit
! This is the grid of choice for global cloud resolving
 
 integer, intent(in):: im
 real(r8), intent(out):: lamda(im+1,im+1)
 real(r8), intent(out):: theta(im+1,im+1)
       
! Local:
 real(r8) pp(3,im+1,im+1)
! real(f_p):: pi, rsq3, alpha, delx, dely
 real(r8):: rsq3, alpha, delx, dely
 integer i, j, k
 
  rsq3 = 1./sqrt(3.) 
 alpha = asin( rsq3 )
 
! Ranges:
! lamda = [0.75*pi, 1.25*pi]
! theta = [-alpha, alpha]
 
    dely = 2.*alpha / real(im)
    
! Define East-West edges:
 do j=1,im+1
    lamda(1,   j) = 0.75*pi                  ! West edge
    lamda(im+1,j) = 1.25*pi                  ! East edge
    theta(1,   j) = -alpha + dely*real(j-1)  ! West edge
    theta(im+1,j) = theta(1,j)               ! East edge
 enddo

! Get North-South edges by symmetry:

 do i=2,im
    call mirror_latlon(lamda(1,1), theta(1,1), lamda(im+1,im+1), theta(im+1,im+1), &
                       lamda(1,i), theta(1,i), lamda(i,1),       theta(i,      1) )
    lamda(i,im+1) =  lamda(i,1)
    theta(i,im+1) = -theta(i,1)
 enddo

! Set 4 corners:
    call latlon2xyz2(lamda(1    ,  1), theta(1,      1), pp(1,   1,   1))
    call latlon2xyz2(lamda(im+1,   1), theta(im+1,   1), pp(1,im+1,   1))
    call latlon2xyz2(lamda(1,   im+1), theta(1,   im+1), pp(1,   1,im+1))
    call latlon2xyz2(lamda(im+1,im+1), theta(im+1,im+1), pp(1,im+1,im+1))

! Map edges on the sphere back to cube:
! Intersections at x=-rsq3

 i=1
 do j=2,im
    call latlon2xyz2(lamda(i,j), theta(i,j), pp(1,i,j))
    pp(2,i,j) = -pp(2,i,j)*rsq3/pp(1,i,j)
    pp(3,i,j) = -pp(3,i,j)*rsq3/pp(1,i,j)
 enddo

 j=1
 do i=2,im
    call latlon2xyz2(lamda(i,j), theta(i,j), pp(1,i,1))
    pp(2,i,1) = -pp(2,i,1)*rsq3/pp(1,i,1)
    pp(3,i,1) = -pp(3,i,1)*rsq3/pp(1,i,1)
 enddo

 do j=1,im+1
    do i=1,im+1
       pp(1,i,j) = -rsq3
    enddo
 enddo

 do j=2,im+1
    do i=2,im+1
! Copy y-z face of the cube along j=1
       pp(2,i,j) = pp(2,i,1)
! Copy along i=1
       pp(3,i,j) = pp(3,1,j)
    enddo
 enddo

 call cart_to_latlon( (im+1)*(im+1), pp, lamda, theta)

 end subroutine gnomonic_ed

subroutine gnomonic_angl(im, lamda, theta)
! This is the commonly known equi-angular grid
 integer im
 real(r8) lamda(im+1,im+1)
 real(r8) theta(im+1,im+1)
 real(r8) p(3,im+1,im+1)
! Local
 real(r8) rsq3, xf, y0, z0, y, x, z, ds
 real(r8) dy, dz
 integer j,k
 real(r8) dp

 dp = 0.5*pi/real(im)

 rsq3 = 1./sqrt(3.)
 do k=1,im+1
    do j=1,im+1
       p(1,j,k) =-rsq3               ! constant
       p(2,j,k) =-rsq3*tan(-0.25*pi+(j-1)*dp)
       p(3,j,k) = rsq3*tan(-0.25*pi+(k-1)*dp)
    enddo
 enddo

 call cart_to_latlon( (im+1)*(im+1), p, lamda, theta)

 end subroutine gnomonic_angl

 subroutine gnomonic_dist(im, lamda, theta)
! This is the commonly known equi-distance grid
 integer im
 real(r8) lamda(im+1,im+1)
 real(r8) theta(im+1,im+1)
 real(r8) p(3,im+1,im+1)
! Local
 real(r8) rsq3, xf, y0, z0, y, x, z, ds
 real(r8) dy, dz
 integer j,k

! Face-2

 rsq3 = 1./sqrt(3.)
 xf = -rsq3
 y0 =  rsq3;  dy = -2.*rsq3/im
 z0 = -rsq3;  dz =  2.*rsq3/im

 do k=1,im+1
    do j=1,im+1
       p(1,j,k) = xf
       p(2,j,k) = y0 + (j-1)*dy
       p(3,j,k) = z0 + (k-1)*dz
    enddo
 enddo
 call cart_to_latlon( (im+1)*(im+1), p, lamda, theta)

 end subroutine gnomonic_dist

      subroutine mirror_grid(grid_global,ng,npx,npy,ndims,nregions)
         integer, intent(IN)    :: ng,npx,npy,ndims,nregions
         real(r8)   , intent(INOUT) :: grid_global(1-ng:npx  +ng,1-ng:npy  +ng,ndims,1:nregions)
         integer :: i,j,n,n1,n2,nreg
         real(r8) :: x1,y1,z1, x2,y2,z2, ang
!
!    Mirror Across the 0-longitude
!
         nreg = 1
         do j=1,ceiling(npy/2.)
            do i=1,ceiling(npx/2.)

            x1 = 0.25 * (ABS(grid_global(i        ,j        ,1,nreg)) + &
                         ABS(grid_global(npx-(i-1),j        ,1,nreg)) + &
                         ABS(grid_global(i        ,npy-(j-1),1,nreg)) + &
                         ABS(grid_global(npx-(i-1),npy-(j-1),1,nreg)))
            grid_global(i        ,j        ,1,nreg) = SIGN(x1,grid_global(i        ,j        ,1,nreg))
            grid_global(npx-(i-1),j        ,1,nreg) = SIGN(x1,grid_global(npx-(i-1),j        ,1,nreg))
            grid_global(i        ,npy-(j-1),1,nreg) = SIGN(x1,grid_global(i        ,npy-(j-1),1,nreg))
            grid_global(npx-(i-1),npy-(j-1),1,nreg) = SIGN(x1,grid_global(npx-(i-1),npy-(j-1),1,nreg))

            y1 = 0.25 * (ABS(grid_global(i        ,j        ,2,nreg)) + &
                         ABS(grid_global(npx-(i-1),j        ,2,nreg)) + &
                         ABS(grid_global(i        ,npy-(j-1),2,nreg)) + &
                         ABS(grid_global(npx-(i-1),npy-(j-1),2,nreg)))
            grid_global(i        ,j        ,2,nreg) = SIGN(y1,grid_global(i        ,j        ,2,nreg))
            grid_global(npx-(i-1),j        ,2,nreg) = SIGN(y1,grid_global(npx-(i-1),j        ,2,nreg))
            grid_global(i        ,npy-(j-1),2,nreg) = SIGN(y1,grid_global(i        ,npy-(j-1),2,nreg))
            grid_global(npx-(i-1),npy-(j-1),2,nreg) = SIGN(y1,grid_global(npx-(i-1),npy-(j-1),2,nreg))

           ! force dateline/greenwich-meridion consitency
            if (mod(npx,2) /= 0) then
              if ( (i==1+(npx-1)/2.0) ) then
                 grid_global(i,j        ,1,nreg) = 0.0
                 grid_global(i,npy-(j-1),1,nreg) = 0.0
              endif
            endif

            enddo
         enddo

         do nreg=2,nregions
           do j=1,npy
             do i=1,npx

               x1 = grid_global(i,j,1,1)
               y1 = grid_global(i,j,2,1)
               z1 = 1.0

               if (nreg == 2) then
                  ang = -90.
                  call rot_3d( 3, x1, y1, z1, ang, x2, y2, z2, 1, 1)  ! rotate about the z-axis
               elseif (nreg == 3) then
                  ang = -90.
                  call rot_3d( 3, x1, y1, z1, ang, x2, y2, z2, 1, 1)  ! rotate about the z-axis
                  ang = 90.
                  call rot_3d( 1, x2, y2, z2, ang, x1, y1, z1, 1, 1)  ! rotate about the x-axis
                  x2=x1
                  y2=y1
                  z2=z1

           ! force North Pole and dateline/greenwich-meridion consitency
                  if (mod(npx,2) /= 0) then
                     if ( (i==1+(npx-1)/2.0) .and. (i==j) ) then
                        x2 = 0.0
                        y2 = pi/2.0
                     endif
                     if ( (j==1+(npy-1)/2.0) .and. (i < 1+(npx-1)/2.0) ) then
                        x2 = 0.0
                     endif
                     if ( (j==1+(npy-1)/2.0) .and. (i > 1+(npx-1)/2.0) ) then
                        x2 = pi
                     endif
                  endif

               elseif (nreg == 4) then
                  ang = -180.
                  call rot_3d( 3, x1, y1, z1, ang, x2, y2, z2, 1, 1)  ! rotate about the z-axis
                  ang = 90.
                  call rot_3d( 1, x2, y2, z2, ang, x1, y1, z1, 1, 1)  ! rotate about the x-axis
                  x2=x1
                  y2=y1
                  z2=z1

               ! force dateline/greenwich-meridion consitency
                  if (mod(npx,2) /= 0) then
                    if ( (j==1+(npy-1)/2.0) ) then
                       x2 = pi
                    endif
                  endif

               elseif (nreg == 5) then
                  ang = 90.
                  call rot_3d( 3, x1, y1, z1, ang, x2, y2, z2, 1, 1)  ! rotate about the z-axis
                  ang = 90.
                  call rot_3d( 2, x2, y2, z2, ang, x1, y1, z1, 1, 1)  ! rotate about the y-axis
                  x2=x1
                  y2=y1
                  z2=z1
               elseif (nreg == 6) then
                  ang = 90.
                  call rot_3d( 2, x1, y1, z1, ang, x2, y2, z2, 1, 1)  ! rotate about the y-axis
                  ang = 0.
                  call rot_3d( 3, x2, y2, z2, ang, x1, y1, z1, 1, 1)  ! rotate about the z-axis
                  x2=x1
                  y2=y1
                  z2=z1

           ! force South Pole and dateline/greenwich-meridion consitency
                  if (mod(npx,2) /= 0) then
                     if ( (i==1+(npx-1)/2.0) .and. (i==j) ) then
                        x2 = 0.0
                        y2 = -pi/2.0
                     endif
                     if ( (i==1+(npx-1)/2.0) .and. (j > 1+(npy-1)/2.0) ) then
                        x2 = 0.0
                     endif
                     if ( (i==1+(npx-1)/2.0) .and. (j < 1+(npy-1)/2.0) ) then
                        x2 = pi
                     endif
                  endif

               endif

               grid_global(i,j,1,nreg) = x2
               grid_global(i,j,2,nreg) = y2

              enddo
            enddo
          enddo

  end subroutine mirror_grid

      subroutine rot_3d(axis, x1in, y1in, z1in, angle, x2out, y2out, z2out, degrees, convert)

         integer, intent(IN) :: axis         ! axis of rotation 1=x, 2=y, 3=z
         real(r8) , intent(IN)    :: x1in, y1in, z1in
         real(r8) , intent(INOUT) :: angle        ! angle to rotate in radians
         real(r8) , intent(OUT)   :: x2out, y2out, z2out
         integer, intent(IN), optional :: degrees ! if present convert angle 
                                                  ! from degrees to radians
         integer, intent(IN), optional :: convert ! if present convert input point
                                                  ! from spherical to cartesian, rotate, 
                                                  ! and convert back
                     
         real(r8)  :: c, s  
         real(r8)  :: x1,y1,z1, x2,y2,z2
                  
         if ( present(convert) ) then
           call spherical_to_cartesian(x1in, y1in, z1in, x1, y1, z1)
         else
           x1=x1in
           y1=y1in
           z1=z1in
         endif
            
         if ( present(degrees) ) then
            angle = angle*pi/180.d0
         endif

         c = COS(angle)
         s = SIN(angle)

         SELECT CASE(axis)
                  
            CASE(1)
               x2 =  x1
               y2 =  c*y1 + s*z1
               z2 = -s*y1 + c*z1
            CASE(2)
               x2 = c*x1 - s*z1
               y2 = y1
               z2 = s*x1 + c*z1
            CASE(3)
               x2 =  c*x1 + s*y1
               y2 = -s*x1 + c*y1
               z2 = z1
            CASE DEFAULT
              write(*,*) "Invalid axis: must be 1 for X, 2 for Y, 3 for Z."
 
         END SELECT
 
         if ( present(convert) ) then
           call cartesian_to_spherical(x2, y2, z2, x2out, y2out, z2out)
         else
           x2out=x2
           y2out=y2
           z2out=z2
         endif

      end subroutine rot_3d

      subroutine cartesian_to_spherical(x, y, z, lon, lat, r)
      real(r8) , intent(IN)  :: x, y, z
      real(r8) , intent(OUT) :: lon, lat, r

      r = SQRT(x*x + y*y + z*z)
      if ( (abs(x) + abs(y)) < 1.E-10 ) then       ! poles:
           lon = 0.
      else
           lon = ATAN2(y,x)    ! range: [-pi,pi]
      endif

#ifdef RIGHT_HAND
      lat = asin(z/r)
#else
      lat = ACOS(z/r) - pi/2.
#endif

      end subroutine cartesian_to_spherical

      subroutine spherical_to_cartesian(lon, lat, r, x, y, z)

         real(r8) , intent(IN)  :: lon, lat, r
         real(r8) , intent(OUT) :: x, y, z

         x = r * COS(lon) * cos(lat)
         y = r * SIN(lon) * cos(lat)

#ifdef RIGHT_HAND
         z =  r * SIN(lat)
#else
         z = -r * sin(lat)
#endif

      end subroutine spherical_to_cartesian

 subroutine symm_ed(im, lamda, theta)
! Make grid symmetrical to i=im/2+1
 integer im
 real(r8) lamda(im+1,im+1)
 real(r8) theta(im+1,im+1)
 integer i,j,ip,jp
 real(r8) avg

 do j=2,im+1
    do i=2,im
       lamda(i,j) = lamda(i,1)
    enddo
 enddo

 do j=1,im+1
    do i=1,im/2
       ip = im + 2 - i
       avg = 0.5*(lamda(i,j)-lamda(ip,j))
       lamda(i, j) = avg + pi
       lamda(ip,j) = pi - avg
       avg = 0.5*(theta(i,j)+theta(ip,j))
       theta(i, j) = avg
       theta(ip,j) = avg
    enddo
 enddo

! Make grid symmetrical to j=im/2+1
 do j=1,im/2
       jp = im + 2 - j
    do i=2,im
       avg = 0.5*(lamda(i,j)+lamda(i,jp))
       lamda(i, j) = avg
       lamda(i,jp) = avg
       avg = 0.5*(theta(i,j)-theta(i,jp))
       theta(i, j) =  avg
       theta(i,jp) = -avg
    enddo
 enddo

 end subroutine symm_ed

 subroutine cart_to_latlon(np, q, xs, ys)
! vector version of cart_to_latlon1
  integer, intent(in):: np
  real(r8), intent(inout):: q(3,np)
  real(r8), intent(inout):: xs(np), ys(np)
! local
  real(r8), parameter:: esl=1.e-10
  real (f_p):: p(3)
  real (f_p):: dist, lat, lon
  integer i,k
 
  do i=1,np
     do k=1,3
        p(k) = q(k,i)
     enddo
     dist = sqrt(p(1)**2 + p(2)**2 + p(3)**2)
     do k=1,3
        p(k) = p(k) / dist
     enddo
 
     if ( (abs(p(1))+abs(p(2)))  < esl ) then
          lon = 0.
     else
          lon = atan2( p(2), p(1) )   ! range [-pi,pi]
     endif

     if ( lon < 0.) lon = 2.*pi + lon
     lat = asin(p(3))
      
     xs(i) = lon
     ys(i) = lat
! q Normalized:
     do k=1,3
        q(k,i) = p(k)
     enddo
  enddo
 
 end  subroutine cart_to_latlon

 subroutine mirror_latlon(lon1, lat1, lon2, lat2, lon0, lat0, lon, lat)
!   
! Given the "mirror" as defined by (lon1, lat1), (lon2, lat2), and center 
! of the sphere, compute the mirror image of (lon0, lat0) as  (lon, lat)
 
 real(r8), intent(in):: lon1, lat1, lon2, lat2, lon0, lat0
 real(r8), intent(inout):: lon(1), lat(1)
!   
 real(r8) p0(3), p1(3), p2(3), nb(3), pp(3)
 real(r8) pdot
 integer k
    
 call latlon2xyz2(lon0, lat0, p0)
 call latlon2xyz2(lon1, lat1, p1)
 call latlon2xyz2(lon2, lat2, p2)
 call vect_cross(nb, p1, p2)
    pdot = sqrt(nb(1)**2+nb(2)**2+nb(3)**2)
 do k=1,3
    nb(k) = nb(k) / pdot
 enddo
 
 pdot = p0(1)*nb(1) + p0(2)*nb(2) + p0(3)*nb(3)
 do k=1,3
    pp(k) = p0(k) - 2.*pdot*nb(k)
 enddo
 
 call cart_to_latlon(1, pp, lon, lat)
 
 end subroutine  mirror_latlon

 subroutine latlon2xyz2(lon, lat, p3)
 real(r8), intent(in):: lon, lat
 real(r8), intent(out):: p3(3)
 real(r8) e(2)

    e(1) = lon;    e(2) = lat
    call latlon2xyz(e, p3)
 
 end subroutine latlon2xyz2
 
    
 subroutine latlon2xyz(p, e)
!   
! Routine to map (lon, lat) to (x,y,z)
!
 real(r8), intent(in) :: p(2)
 real(r8), intent(out):: e(3)
    
 integer n
 real (f_p):: q(2)
 real (f_p):: e1, e2, e3
 
    do n=1,2
       q(n) = p(n)
    enddo

    e1 = cos(q(2)) * cos(q(1))
    e2 = cos(q(2)) * sin(q(1))
    e3 = sin(q(2))
!-----------------------------------
! Truncate to the desired precision:
!-----------------------------------
    e(1) = e1
    e(2) = e2
    e(3) = e3
 
 end subroutine latlon2xyz

 subroutine cell_center2(q1, q2, q3, q4, e2)
      real(r8) , intent(in ) :: q1(2), q2(2), q3(2), q4(2)
      real(r8) , intent(OUT) :: e2(2)
! Local
      real(r8) p1(3), p2(3), p3(3), p4(3)
      real(r8) ec(3)
      real(r8) dd
      integer k

      call latlon2xyz(q1, p1)
      call latlon2xyz(q2, p2)
      call latlon2xyz(q3, p3)
      call latlon2xyz(q4, p4)

      do k=1,3
         ec(k) = p1(k) + p2(k) + p3(k) + p4(k)
      enddo
      dd = sqrt( ec(1)**2 + ec(2)**2 + ec(3)**2 )

      do k=1,3
         ec(k) = ec(k) / dd
      enddo

      call cart_to_latlon(1, ec, e2(1), e2(2))

 end subroutine cell_center2

 subroutine vect_cross(e, p1, p2)
 real(r8), intent(in) :: p1(3), p2(3)
 real(r8), intent(out):: e(3)
!       
! Perform cross products of 3D vectors: e = P1 X P2
!    
      e(1) = p1(2)*p2(3) - p1(3)*p2(2)
      e(2) = p1(3)*p2(1) - p1(1)*p2(3) 
      e(3) = p1(1)*p2(2) - p1(2)*p2(1)
        
 end subroutine vect_cross

end module CubedSphere_GridMod
