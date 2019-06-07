
!  $Id$


#ifdef MESH
#define LRRasterize2Mem  LRRasterize2Mem0
#define LRRasterize2File  LRRasterize2File0
#define POLYSIZE 4
#else
#define POLYSIZE size(xv,3)
#endif

!BOP
! !IROUTINE: LRRasterize -- Rasterizes a 2-D array of polygons.

  subroutine LRRasterize2File(GridName, xv,yv,nc,nr,xmn,xmx,ymn,ymx,     &
                              SurfaceType, Verb, Here,jseg               )

! This routine rasterizes a grid defined by a 2-D array of polygons.
! The raster value assigned each pixel is either sign((I*pushleft + J),sgn),
! where I and J are the inner and outer indeces in the 2-D array of 
! polygons of the polygon containing the pixel, or an index into a 
! table describing the polygon. The former option is kept for backward
! compatibility, and is the default; the latter is the prefered option,
! and it is invoked by specifying TileFile, the name of the file where
! the table is to be stored. The sgn argument has no effect is the table
! option is chosen, while SurfaceType and GridName are used only with 
! the table option.
!
! A value is assigned only if the pixel is RASTERUNDEF. As a result,
! if polygons overlap, the pixel is assigned to the first polygon in
! Fortran index order, i.e, (1,1), (2,1), etc.
!
! The polygons are specified with the arguments xv and yv, which 
! contain the two coordinates of each of the vertices.
! The xv and yv are modified on exit. The third  dimension of xv and yv
! contains the vertices. Vertices should be in counterclockwise order
! as viewed from the top or from "outside" the sphere. 
! A polygon array with varying number of vertices can be accomodated
! by repeating the last vertex as necessary.
!   
! The x and y refer to the coordinates of the
! the two dimensions of the raster field, whose bounds can be specified.
! with the optional arguments xmn, xmx, ymn, ymx.
!
! By default, polygon vertices are connected by straight lines assuming 
! cartesian coordinates. If sphere is specified as .true. the coordinates
! are assumed to be longitude (X) and latitude (Y) in degrees and the
! vertices are connected with great circles.
! 
! The x coordinate is assumed to be cyclic.
!
! A 2-D array of polygons can be rasterized in sections by specifying 
! the beginning ideces of xv and yv with the optional arguments imin and jmin.
! The defaults are (1,1). 

    character*(*),   intent(INOUT)  :: GridName  ! Raster file name
#ifdef MESH
    REAL_,  intent(INOUT)  :: xv(:,:) ! X coordinates of vertices
    REAL_,  intent(INOUT)  :: yv(:,:) ! Y coordinates of vertices
#else
    REAL_,  intent(INOUT)  :: xv(:,:,:) ! X coordinates of vertices
    REAL_,  intent(INOUT)  :: yv(:,:,:) ! Y coordinates of vertices 
#endif
    integer, optional, intent(IN) :: nc,nr ! Raster field sizes
    REAL_,   optional, intent(IN) :: xmn   ! LL x of LL raster cell (-180)
    REAL_,   optional, intent(IN) :: ymn   ! LL y of LL raster cell ( -90)
    REAL_,   optional, intent(IN) :: xmx   ! UR x of UR raster cell ( 180)
    REAL_,   optional, intent(IN) :: ymx   ! UR y of UR raster cell (  90)
    logical, optional, intent(IN) :: verb  ! Verbose
    logical, optional, intent(IN) :: here  ! write here
    integer, optional, intent(IN) :: SurfaceType
    integer, optional, intent(IN) :: jseg

    character*(128)       :: TileFile
    character*(128)       :: TilDir='', RstDir=''
    character*(128)       :: RasterFile
    integer, allocatable  :: Raster(:,:)
    integer               :: idx, j, unc,  unr
    integer               :: STATUS,i,n
    logical               :: DoZip, failed

! Begin
!------

! Process optiona raster sizes. Defaults are 8640 x 4320
!-------------------------------------------------------

    if(present(nc)) then
       unc = nc
    else
       unc = NX
    endif

    if(present(nr)) then
       unr = nr
    else
       unr = NY
    endif

    if(present(Here)) then
       if(.not.Here) then
          TilDir = 'til/'
          RstDir = 'rst/'
       end if
    end if

! Allocate interger*4 raster array and initialize to RASTERUNDEF
!---------------------------------------------------------------

    allocate(Raster(Unc,Unr),stat=STATUS)
    ASSERT_(STATUS==0)

    Raster = RASTERUNDEF

! Check if the gridname ends in .gz, which
! is the signal for creating gzipped files
!-----------------------------------------

    idx    = index(GridName,'.',back=.true.)
    DoZip  = idx>0 .and. GridName(idx+1:idx+2)=="gz"

! Create the file names
!----------------------

    if(DoZip) then
       GridName   = GridName(1:idx-1)
       TileFile   = trim(Tildir)//trim(GridName)//'.til.gz'
       RasterFile = trim(Rstdir)//trim(GridName)//'.rst.gz'
    else
       TileFile   = trim(Tildir)//trim(GridName)//'.til'
       RasterFile = trim(Rstdir)//trim(GridName)//'.rst'
    end if

! Do the Rasterization and write the TileFile (.til)
!---------------------------------------------------

    call LRRasterize2Mem (Raster, xv, yv, TileFile,   &
            xmn, xmx, ymn, ymx,         &
            SurfaceType, verb, jseg           )

! Make sure it worked
!--------------------

    failed=.false.
    do j=1,unr
       n = count(raster(:,j)==RASTERUNDEF)
       if(n>0) then
          print *, ' At j=',j,' there are ',n,' undefined pixels'
          failed=.true.
       end if
    enddo

    if(.not.failed) print *, 'Rasterization completed successfully'

! Write the Raster File (.rst)
!-----------------------------

    call WriteRaster(RasterFile,Raster,DoZip)

    deallocate(Raster)

    return
  end subroutine LRRasterize2File

!===================================================================


!BOP

! !IROUTINE: LRRasterize -- Rasterizes a 2-D array of polygons.

  subroutine LRRasterize2Mem(Raster,xv,yv,Tilefile,    &
                             xmn,xmx,ymn,ymx,          &
                             SurfaceType,verb,jseg     )

! This routine rasterizes a grid defined by a 2-D array of polygons.
! The raster value assigned is an index to the table.
! A value is assigned only if the pixel is RASTERUNDEF. As a result,
! if polygons overlap the pixel is assigned to the first one in
! Fortran index order.
!
! The polygons are specified with the arguments xv and yv, which 
! contain the two coordinates of each of the vertices.
! The xv and yv are modified on exit.
! 
! The x and y are the Cartesian coordinates of the
! the two dimensions of the raster field, whose bounds can be specified.
! The x coordinate is assumed
! to be cyclic. 

    integer,       intent(INOUT) :: Raster(:,:)  ! Raster field to be filled
#ifdef MESH
    REAL_,         intent(INOUT) :: xv(:,:  )    ! X coordinates of vertices
    REAL_,         intent(INOUT) :: yv(:,:  )    ! Y coordinates of vertices 
#else
    REAL_,         intent(INOUT) :: xv(:,:,:)    ! X coordinates of vertices
    REAL_,         intent(INOUT) :: yv(:,:,:)    ! Y coordinates of vertices 
#endif
    character*(*), intent(IN   ) :: TileFile

    REAL_,   optional, intent(IN) :: xmn   ! LL x of LL raster cell (-180)
    REAL_,   optional, intent(IN) :: ymn   ! LL y of LL raster cell ( -90)
    REAL_,   optional, intent(IN) :: xmx   ! UR x of UR raster cell ( 180)
    REAL_,   optional, intent(IN) :: ymx   ! UR y of UR raster cell (  90)
    logical, optional, intent(IN) :: verb  ! Verbose
    integer, optional, intent(IN) :: SurfaceType
    integer, optional, intent(IN) :: jseg

    !  X abd Y bounds of each polygon

    REAL_ ::  xmin, xmax
    REAL_ ::  ymin, ymax
    REAL_ ::  minx, miny
    REAL_ ::  maxx, maxy

    !  x and y coordinates of the Raster grid

    REAL_, dimension(size(Raster,1))        :: xs, xcs, xss
    REAL_, dimension(size(Raster,2))        :: ys, ycs, yss, da

    integer :: IM, JM, NV                   ! Shape of input grid
    REAL_   :: dx, dy, dxi, dyi             ! Grid spacing of raster grid
    integer :: xsize, ysize                 ! Dimensions of Raster grid

    integer :: i, j, jn, n, ib, jb, fill, uType, js, k
    REAL_   :: range, d2r, r2d, ddx, grid_ymin, grid_ymax, xc, yc, Area, xx
    logical :: DoZip, uVerb
    integer :: idx, ct
    integer :: count0,count1,count_rate
    integer :: maxtiles

    character*(128)  :: GridName, TilFile

    integer, allocatable :: iTable(:,:)
    REAL_,   allocatable :: rTable(:,:)

    integer :: useg, unit, fq


    integer, dimension(POLYSIZE) &
            :: nxt
    REAL_,   dimension(POLYSIZE) &
            :: xvc, yvc, xvs, yvs, xrd, yrd, x3, y3, z3, &
               dx3, dy3, dz3, x31, x32, y31, y32, z31, z32, &
               dx4, dy4, x4, y4,   xu,  yu

! Process optionals

    if(present(jseg)) then
       useg = jseg
    else
       useg = 1
    endif

    if(present(SurfaceType)) then
       uType = SurfaceType
    else
       uType = -1
    endif

    if(present(Verb       )) then
       uVerb = Verb
    else
       uVerb = .false.
    endif

! Domain to rasterize

    if(present(xmn )) then
       xmin = xmn
    else
       xmin = -180.
    end if

    if(present(xmx )) then
       xmax = xmx
    else
       xmax =  180.
    end if

    if(present(ymn )) then
       ymin = ymn
    else
       ymin = -90.
    end if

    if(present(ymx )) then
       ymax = ymx
    else
       ymax =  90.
    end if

!  Input grid sizes
    
    ib = 0
    jb = 0

    print *, 'In Rasterize...'

! Polygon sizes

    NV     = POLYSIZE
    ASSERT_(NV > 1)

#ifdef MESH
    IM     = size(xv,1)-1
    JM     = size(xv,2)-1*useg

    ASSERT_(IM == size(yv,1)-1)
    ASSERT_(JM == size(yv,2)-1*useg)
#else
    IM     = size(xv,1)
    JM     = size(xv,2)

    ASSERT_(IM == size(yv,1))
    ASSERT_(JM == size(yv,2))
    ASSERT_(NV == size(yv,3))
#endif


    do i=1,NV
       nxt(i) = mod(i,NV) + 1
    enddo

    !  Raster sizes

    xsize = size(Raster,1)
    ysize = size(Raster,2)

    !  Raster grid sizes

    range = (xmax-xmin)
    dx    = (xmax-xmin)/xsize
    dy    = (ymax-ymin)/ysize
    dxi   = 1.0/dx
    dyi   = 1.0/dy

    d2r    = (2._8*RASTER_PI)/range
    r2d    = range/(2._8*RASTER_PI)

    !  Report

    if(uVerb) then
       print *, 'Input grid extents are within Raster limits:'
       print *, '  Raster limits:'
       print *, '   Min and Max X: ', xmin, xmax
       print *, '   Min and Max Y: ', ymin, ymax
    end if

    ! Compute ys and xs at centers of raster cell

    xs(1) = xmin + dx*0.5
    do i = 2,xsize
       xs(i) = xs(1) + dx*(i-1)
    enddo

    ys(1) = ymin + dy*0.5
    do j = 1,ysize
       ys(j) = ys(1) + dy*(j-1)
!       if(ys(j)<grid_ymin) Raster(:,j) = 0
!       if(ys(j)>grid_ymax) Raster(:,j) = 0
       if(ys(j)<ymin) Raster(:,j) = 0
       if(ys(j)>ymax) Raster(:,j) = 0
    enddo

    xss = xs*d2r
    yss = ys*d2r
    xcs = cos(xss)
    ycs = cos(yss)
    xss = sin(xss)
    yss = sin(yss)
    da  = (sin(d2r*(ys+0.5*dy)) -   &
           sin(d2r*(ys-0.5*dy)) )*(dx*d2r)

! Process tilefile name
!----------------------

    TilFile = TileFile

    idx = index(TilFile,'.',back=.true.)
    DoZip = idx>0 .and. TilFile(idx+1:idx+2)=="gz"
    if(DoZip) TilFile = TilFile(1:idx-1)

    idx = index(TilFile,'.til',back=.true.)
    GridName = TilFile(1:idx-1)

! Polygon counter; also value with which raster is filled
!  for pixels with that polygon.
!--------------------------------------------------------

    fill = 1
    fq = 1
    ASSERT_(mod(jm,useg)==0)	
    jn = jm/useg	

     if (useg > 1) then
       call OpenTiling(Unit,TilFile,Gridname,im,jm,im*jm,  &
                     xsize,ysize,Dozip,uVerb                  )
       maxtiles = 1
     else	
       maxtiles = im*jm
     endif

! Allocate Table for a single grid
!---------------------------------

    allocate(iTable(0:3,maxtiles),stat=k)

    if(k/=0) then 
       print *, 'allocation error ', k, ' at line ',__LINE__
       call exit(1)
    endif		

    allocate(rTable(1:4,maxtiles),stat=k)

    if(k/=0) then 
       print *, 'allocation error ', k, ' at line ',__LINE__
       call exit(1)
    endif		

    if(uVerb) then
       call system_clock(count0,count_rate)
       write (6, '(A)', advance='NO') 'Begin looping over polygons'
    end if

    Polygons: do k=0,useg-1 
	do j = 1, jn
          js = (jn+1)*k + j
    if(uVerb .and. mod(j,50)==0) write (6, '(A)', advance='NO') '.'
       do i = 1, im


#ifdef MESH
          xu(1) = xv(i  ,js  )
          xu(2) = xv(i+1,js  )
          xu(3) = xv(i+1,js+1)
          xu(4) = xv(i  ,js+1)
          yu(1) = yv(i  ,js  )
          yu(2) = yv(i+1,js  )
          yu(3) = yv(i+1,js+1)
          yu(4) = yv(i  ,js+1)
#else
          xu = xv(i,j,:)
          yu = yv(i,j,:)
#endif

! Fill-in Raster values for each polygon
!---------------------------------------

          ct   = 0
          Area = 0.
          xc   = 0.
          yc   = 0.
          xx   = 0.

          do n=1,nv 
             if(xu(n)-xu(nxt(n)) > 90.001 ) xu(n) = xu(n)-range
             if(xu(n)-xu(nxt(n)) <-90.001 ) xu(n) = xu(n)+range
          enddo

          miny  = minval(yu)
          maxy  = maxval(yu)

          if(miny<-86.) miny = -90.
          if(maxy> 86.) maxy =  90.

          minx  = minval(xu)
          maxx  = maxval(xu)

          call FillPoly(Zero)

          if    (maxx>=xmax) then
             call FillPoly( range)
          elseif(minx<=xmin) then
             call FillPoly(-range)
          end if

! Skip empty polygons
!--------------------

          if(ct==0) then
             if(uVerb) then
                print *, 'Empty Poly:', i, js, j, k
                print *, '           ', xu
                print *, '           ', yu
             end if
             cycle
          end if

! Compute lat and lon
!--------------------

          yc = yc/area
          xc = atan2(xx,xc)*r2d

! Update tabled for .til file
!---------------------------- 


          iTable(0 ,fq) = uType
          iTable(1 ,fq) = ct
          iTable(2 ,fq) = i
          iTable(3 ,fq) = jn*k + j

          rTable(1 ,fq) = xc
          rTable(2 ,fq) = yc
          rTable(3 ,fq) = area
          rTable(4 ,fq) = area

! Next polygon
!-------------

          if(useg>1) then
            call writeline(TilFile,Unit,iTable(:,1), rTable(:,1),fill,Dozip,uVerb              )
            fill = fill + 1
          else
            fill = fill + 1
            fq = fill
          endif

       enddo
       enddo

    enddo Polygons

    if(uVerb) then
       call system_clock(count1)
       print *, 'Done. Time = ', (count1-count0)/float(count_rate)
       count0 = count1
    end if

    if(useg>1) then
      call CloseTiling(TilFile, Unit, im*jm,  DoZip, uVerb)
    else
      fill = fill - 1
      call Writetiling(TilFile,(/Gridname/),(/im/),(/jm/),(/fill/),  &
                     xsize,ysize,iTable(:,1:fill),rTable(:,1:fill),&
                     Dozip,uVerb                                   )
    endif

    if(uVerb) then
       call system_clock(count1)
       print *, 'Wrote Tiling. Time = ', (count1-count0)/float(count_rate)
       count0 = count1
    end if
    return

  contains

  subroutine FillPoly(sh)

    REAL_, intent(IN) :: sh

    logical :: IsIn
    integer :: i1, i2, jj1
    integer :: ii, jj, n1, n2, jx
    integer, save :: j1, j2
    integer :: HALO=10
    REAL_   :: x0, y0


    if    (abs(miny+90._8) < 1.e-10_8) then
       i1 = 1
       i2 = xsize
    elseif(abs(maxy-90._8) < 1.e-10_8) then
       i1 = 1 
       i2 = xsize
    else
       i1 = max(nint((minx-(xmin+sh))*dxi),1    )
       i2 = min(nint((maxx-(xmin+sh))*dxi),xsize)
    end if

    if(sh==Zero) then

       xrd = xu*d2r
       yrd = yu*d2r
       xvc = cos(xrd)
       xvs = sin(xrd)
       yvc = cos(yrd)
       yvs = sin(yrd)
       z3  = yvs
       x3  = yvc*xvc
       y3  = yvc*xvs

       do n=1,NV
          n1 = n
          n2 = nxt(n)
          dx3(n1) = x3(n2)-x3 (n1)
          dy3(n1) = y3(n2)-y3 (n1)
          dz3(n1) = z3(n2)-z3 (n1)
          x31(n1) = x3(n1)*dz3(n1)
          x32(n1) = x3(n1)*dy3(n1)
          y31(n1) = y3(n1)*dx3(n1)
          y32(n1) = y3(n1)*dz3(n1)
          z31(n1) = z3(n1)*dy3(n1)
          z32(n1) = z3(n1)*dx3(n1)

          z31(n1) = (z31(n1)-y32(n1))
          x31(n1) = (x31(n1)-z32(n1))
          y31(n1) = (y31(n1)-x32(n1))

          dx3(n1) = x3(n1)*z31(n1) + &
                    y3(n1)*x31(n1) + &
                    z3(n1)*y31(n1)

          halo = max(2,nint(10.*(1.-minval(yvc))))
          j1 = max(nint((miny- ymin    )*dyi)-halo,1    )
          j2 = min(nint((maxy- ymin    )*dyi)+halo,ysize)
       end do
    end if

    LonLoop: do ii=i1,i2

       jj1 = -1
       jj  = j1

       LatLoop: do

          EmptyPixel: if(Raster(ii,jj) == RASTERUNDEF) then

             x0 = ycs(jj)*xcs(ii)
             y0 = ycs(jj)*xss(ii)
             
             do n=1,NV
                Isin = (x0*z31(n) + y0*x31(n) + yss(jj)*y31(n) - dx3(n)) <= 1.e-12
                if(.not.IsIn) exit
             end do

             if(IsIn) then
                if(jj1<0) then
                   jj1 = jj
                   jj  = j2
                   cycle
                else
                   do jx=jj1,jj
                      Raster(ii,jx) = fill
                      ct   = ct   + 1
                      xc   = xc   + xcs(ii)*da(jx)
                      xx   = xx   + xss(ii)*da(jx)
                      Area = Area +         da(jx)
                      yc   = yc   + ys (jx)*da(jx)
                   enddo
                   exit
                endif
             end if

          end if EmptyPixel

          if(jj1<0) then
             if(jj==j2) exit
             jj = jj + 1
          else
             jj = jj - 1
          endif

       enddo LatLoop
    enddo LonLoop

  end subroutine FillPoly


end subroutine LRRasterize2Mem


#ifdef MESH
#undef LRRasterize2Mem
#undef LRRasterize2File
#endif


#undef POLYSIZE

