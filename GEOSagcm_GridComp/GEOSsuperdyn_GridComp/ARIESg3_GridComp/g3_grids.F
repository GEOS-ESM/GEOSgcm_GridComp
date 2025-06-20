      subroutine rotate_f ( q_geo,q_cmp,im,jm,lm,
     .                      dlam,dphi,rotation,tilt,precession,
     .                      sgn,norder,check,lattice )
C***********************************************************************
C
C  PURPOSE:
C  ========
C    Compute a forward transformation from the geophysical grid to the
C    computational grid.
C
C  INPUT:
C  ======
C    q_geo ...... Field q_geo(im,jm,lm) on the geophysical grid
C    im ......... Longitudinal dimension of q_geo
C    jm ......... Latitudinal  dimension of q_geo
C    lm ......... Vertical     dimension of q_geo
C    dlam ....... Computational Grid Delta Lambda
C    dphi ....... Computational Grid Delta Phi
C    rotation ... Rotation parameter lam_np (Degrees)
C    tilt ....... Rotation parameter phi_np (Degrees)
C    precession . Rotation parameter lam_0  (Degrees)
C    sgn ........ Flag for scalar field  ( sgn =  1 )
C                    or vector component ( sgn = -1 )
C    norder ..... Order of Interpolation:  Bi-Linear => abs(norder) = 1
C                                          Bi-Cubic  => abs(norder) = 3
C                 Note: If norder < 0, then check for positive definite
C    check ...... Logical Flag to check for Undefined values
C
C  OUTPUT:
C  =======
C    q_cmp ...... Field q_cmp(im,jm,lm) on the computational grid
C
C
C***********************************************************************
C*                  GODDARD LABORATORY FOR ATMOSPHERES                 *
C***********************************************************************

      use g3_dynamics_lattice_module
      implicit none
      type ( dynamics_lattice_type ) lattice

      include 'mpif.h'
      integer  status(mpi_status_size)
      integer   stats(mpi_status_size,0:lattice%nx*lattice%ny-1)
      integer   statr(mpi_status_size,0:lattice%nx*lattice%ny-1)

      integer  sendquest(0:lattice%nx*lattice%ny-1)
      integer  recvquest(0:lattice%nx*lattice%ny-1)
      integer  ierror

      integer, allocatable, save :: sendlen(:)    ! Number of    unique           grid-points needed to SEND to   other PEs
      integer, allocatable, save :: recvlen(:)    ! Number of    unique           grid-points needed to RECV from other PEs
      integer, allocatable, save :: sendloc(:,:)  ! Global index locations of the grid-points needed to SEND to   other PEs
      integer, allocatable, save :: recvloc(:,:)  ! Global index locations of the grid-points needed to RECV from other PEs

      integer, allocatable       :: bucket (:)    ! Dummy bucket for array reduction
      real(kind=8),    allocatable       :: sendbuf(:,:)  ! Dummy buffer to SEND data
      real(kind=8),    allocatable       :: recvbuf(:,:)  ! Dummy buffer to RECV data

c Input Variables
c ---------------
      integer im,jm,lm,norder,sgn
      real(kind=8)    q_geo(im,jm,lm)
      real(kind=8)    q_cmp(im,jm,lm)
      real(kind=8)     dlam(lattice%imglobal)
      real(kind=8)     dphi(lattice%jmglobal)
      real(kind=8)     rotation, tilt, precession
      logical  check

c Local Variables
c ---------------
      integer  i,j,l,m,n,sgnjp1
      integer  peid, num,sgnjm2
      integer  pei_im2
      integer  pei_im1
      integer  pei_ip0
      integer  pei_ip1
      integer  pej_jm2
      integer  pej_jm1
      integer  pej_jp0
      integer  pej_jp1

c Index Locations
c ---------------
      integer, allocatable ::         len(:)  ! Total number of grid-points needed from other PEs
      integer, allocatable ::        nn(:,:)  ! Mapping from total to unique grid-points
      integer, allocatable ::       pos(:,:)  ! Position index of data to RECV from other PEs
      integer, allocatable ::       loc(:,:)

      integer, allocatable ::       ip1(:,:),       ip0(:,:),       im1(:,:),       im2(:,:)
      integer, allocatable ::       jp1(:,:),       jp0(:,:),       jm1(:,:),       jm2(:,:)
      integer, allocatable ::    im2jm2(:,:),    im2jm1(:,:),    im2jp0(:,:),    im2jp1(:,:)
      integer, allocatable ::    im1jm2(:,:),    im1jm1(:,:),    im1jp0(:,:),    im1jp1(:,:)
      integer, allocatable ::    ip0jm2(:,:),    ip0jm1(:,:),    ip0jp0(:,:),    ip0jp1(:,:)
      integer, allocatable ::    ip1jm2(:,:),    ip1jm1(:,:),    ip1jp0(:,:),    ip1jp1(:,:)
      integer, allocatable :: pe_im2jm2(:,:), pe_im2jm1(:,:), pe_im2jp0(:,:), pe_im2jp1(:,:)
      integer, allocatable :: pe_im1jm2(:,:), pe_im1jm1(:,:), pe_im1jp0(:,:), pe_im1jp1(:,:)
      integer, allocatable :: pe_ip0jm2(:,:), pe_ip0jm1(:,:), pe_ip0jp0(:,:), pe_ip0jp1(:,:)
      integer, allocatable :: pe_ip1jm2(:,:), pe_ip1jm1(:,:), pe_ip1jp0(:,:), pe_ip1jp1(:,:)

      integer, allocatable, save ::     msgn(:,:)
      integer, allocatable, save ::     nsgn(:,:)
      integer, allocatable, save :: n_im2jm2(:,:), n_im2jm1(:,:), n_im2jp0(:,:), n_im2jp1(:,:)
      integer, allocatable, save :: n_im1jm2(:,:), n_im1jm1(:,:), n_im1jp0(:,:), n_im1jp1(:,:)
      integer, allocatable, save :: n_ip0jm2(:,:), n_ip0jm1(:,:), n_ip0jp0(:,:), n_ip0jp1(:,:)
      integer, allocatable, save :: n_ip1jm2(:,:), n_ip1jm1(:,:), n_ip1jp0(:,:), n_ip1jp1(:,:)
      integer, allocatable, save :: p_im2jm2(:,:), p_im2jm1(:,:), p_im2jp0(:,:), p_im2jp1(:,:)
      integer, allocatable, save :: p_im1jm2(:,:), p_im1jm1(:,:), p_im1jp0(:,:), p_im1jp1(:,:)
      integer, allocatable, save :: p_ip0jm2(:,:), p_ip0jm1(:,:), p_ip0jp0(:,:), p_ip0jp1(:,:)
      integer, allocatable, save :: p_ip1jm2(:,:), p_ip1jm1(:,:), p_ip1jp0(:,:), p_ip1jp1(:,:)

c Bi-Linear Weights
c -----------------
      real(kind=8), allocatable, save :: wl_ip0jp0 (:,:)
      real(kind=8), allocatable, save :: wl_im1jp0 (:,:)
      real(kind=8), allocatable, save :: wl_ip0jm1 (:,:)
      real(kind=8), allocatable, save :: wl_im1jm1 (:,:)

c Bi-Cubic Weights
c ----------------
      real(kind=8), allocatable, save :: wc_ip1jp1 (:,:)
      real(kind=8), allocatable, save :: wc_ip0jp1 (:,:)
      real(kind=8), allocatable, save :: wc_im1jp1 (:,:)
      real(kind=8), allocatable, save :: wc_im2jp1 (:,:)
      real(kind=8), allocatable, save :: wc_ip1jp0 (:,:)
      real(kind=8), allocatable, save :: wc_ip0jp0 (:,:)
      real(kind=8), allocatable, save :: wc_im1jp0 (:,:)
      real(kind=8), allocatable, save :: wc_im2jp0 (:,:)
      real(kind=8), allocatable, save :: wc_ip1jm1 (:,:)
      real(kind=8), allocatable, save :: wc_ip0jm1 (:,:)
      real(kind=8), allocatable, save :: wc_im1jm1 (:,:)
      real(kind=8), allocatable, save :: wc_im2jm1 (:,:)
      real(kind=8), allocatable, save :: wc_ip1jm2 (:,:)
      real(kind=8), allocatable, save :: wc_ip0jm2 (:,:)
      real(kind=8), allocatable, save :: wc_im1jm2 (:,:)
      real(kind=8), allocatable, save :: wc_im2jm2 (:,:)

      real(kind=8), allocatable, save ::  old_dlam (:)
      real(kind=8), allocatable, save ::  old_dphi (:)

      real(kind=8)    ux, ap1, ap0, am1, am2
      real(kind=8)    uy, bp1, bp0, bm1, bm2

      real(kind=8)    lon_cmp(lattice%imglobal), lon_geo(lattice%imglobal)
      real(kind=8)    lat_cmp(lattice%jmglobal), lat_geo(lattice%jmglobal)
      real(kind=8)    pi,cosnp,sinnp,p1,p2,p3,eps,lam,phi,d
      real(kind=8)    lam1,lam2,phi1,phi2
      real(kind=8)    dl,dp,lam_np,phi_np,lam_0,eps_np
      real(kind=8)    lam_geo, lam_cmp, dlam_max, dlam_min
      real(kind=8)    phi_geo, phi_cmp, dphi_max, dphi_min
      real(kind=8)    undef, getcon
      integer im1_geo,igeo
      integer jm1_geo,jgeo

      logical compute_weights
      real(kind=8)    old_rotation
      real(kind=8)    old_tilt
      real(kind=8)    old_precession
      data    old_rotation   /-999.9/
      data    old_tilt       /-999.9/
      data    old_precession /-999.9/

      parameter ( eps = 1.e-10 )

c Initialization
c --------------
      call my_barrier (lattice%comm)  ! Uncomment for more accurate timing

      pi = 4.*atan(1.)
      dl = 2*pi/ lattice%imglobal     ! Uniform Grid Delta Lambda
      dp =   pi/(lattice%jmglobal-1)  ! Uniform Grid Delta Phi

c Allocate Memory for Weights and Index Locations
c -----------------------------------------------
      if(.not.allocated(old_dlam)) then

      allocate (  old_dlam(lattice%imglobal)    ,  old_dphi(lattice%jmglobal)    )
      allocate ( wl_ip0jp0(im,jm) , wl_im1jp0(im,jm) )
      allocate ( wl_ip0jm1(im,jm) , wl_im1jm1(im,jm) )
      allocate ( wc_ip1jp1(im,jm) , wc_ip0jp1(im,jm) , wc_im1jp1(im,jm) , wc_im2jp1(im,jm) )
      allocate ( wc_ip1jp0(im,jm) , wc_ip0jp0(im,jm) , wc_im1jp0(im,jm) , wc_im2jp0(im,jm) )
      allocate ( wc_ip1jm1(im,jm) , wc_ip0jm1(im,jm) , wc_im1jm1(im,jm) , wc_im2jm1(im,jm) )
      allocate ( wc_ip1jm2(im,jm) , wc_ip0jm2(im,jm) , wc_im1jm2(im,jm) , wc_im2jm2(im,jm) )
      allocate (      msgn(im,jm) ,      nsgn(im,jm) )
      allocate (  n_im2jm2(im,jm) ,  n_im1jm2(im,jm) ,  n_ip0jm2(im,jm) ,  n_ip1jm2(im,jm) )
      allocate (  n_im2jm1(im,jm) ,  n_im1jm1(im,jm) ,  n_ip0jm1(im,jm) ,  n_ip1jm1(im,jm) )
      allocate (  n_im2jp0(im,jm) ,  n_im1jp0(im,jm) ,  n_ip0jp0(im,jm) ,  n_ip1jp0(im,jm) )
      allocate (  n_im2jp1(im,jm) ,  n_im1jp1(im,jm) ,  n_ip0jp1(im,jm) ,  n_ip1jp1(im,jm) )
      allocate (  p_im2jm2(im,jm) ,  p_im1jm2(im,jm) ,  p_ip0jm2(im,jm) ,  p_ip1jm2(im,jm) )
      allocate (  p_im2jm1(im,jm) ,  p_im1jm1(im,jm) ,  p_ip0jm1(im,jm) ,  p_ip1jm1(im,jm) )
      allocate (  p_im2jp0(im,jm) ,  p_im1jp0(im,jm) ,  p_ip0jp0(im,jm) ,  p_ip1jp0(im,jm) )
      allocate (  p_im2jp1(im,jm) ,  p_im1jp1(im,jm) ,  p_ip0jp1(im,jm) ,  p_ip1jp1(im,jm) )
      allocate ( recvloc(16*im*jm,0:lattice%nx*lattice%ny-1) )
      allocate ( recvlen(         0:lattice%nx*lattice%ny-1) )
      allocate ( sendlen(         0:lattice%nx*lattice%ny-1) )
      allocate ( sendloc(   im*jm,0:lattice%nx*lattice%ny-1) )
c     allocate ( sendloc(lattice%imglobal*lattice%jmglobal,0:lattice%nx*lattice%ny-1) )

      do i=1,lattice%imglobal
      old_dlam(i) = 0.0
      enddo
      do j=1,lattice%jmglobal
      old_dphi(j) = 0.0
      enddo

      else
             i =  size (old_dlam)
             j =  size (old_dphi)
          if(i.ne.lattice%imglobal .or. j.ne.lattice%jmglobal) then
          deallocate (  old_dlam ,  old_dphi )
          deallocate ( wl_ip0jp0 , wl_im1jp0 )
          deallocate ( wl_ip0jm1 , wl_im1jm1 )
          deallocate ( wc_ip1jp1 , wc_ip0jp1 , wc_im1jp1 , wc_im2jp1 )
          deallocate ( wc_ip1jp0 , wc_ip0jp0 , wc_im1jp0 , wc_im2jp0 )
          deallocate ( wc_ip1jm1 , wc_ip0jm1 , wc_im1jm1 , wc_im2jm1 )
          deallocate ( wc_ip1jm2 , wc_ip0jm2 , wc_im1jm2 , wc_im2jm2 )
          deallocate (      msgn ,      nsgn )
          deallocate (  n_im2jm2 ,  n_im1jm2 ,  n_ip0jm2 ,  n_ip1jm2 )
          deallocate (  n_im2jm1 ,  n_im1jm1 ,  n_ip0jm1 ,  n_ip1jm1 )
          deallocate (  n_im2jp0 ,  n_im1jp0 ,  n_ip0jp0 ,  n_ip1jp0 )
          deallocate (  n_im2jp1 ,  n_im1jp1 ,  n_ip0jp1 ,  n_ip1jp1 )
          deallocate (  p_im2jm2 ,  p_im1jm2 ,  p_ip0jm2 ,  p_ip1jm2 )
          deallocate (  p_im2jm1 ,  p_im1jm1 ,  p_ip0jm1 ,  p_ip1jm1 )
          deallocate (  p_im2jp0 ,  p_im1jp0 ,  p_ip0jp0 ,  p_ip1jp0 )
          deallocate (  p_im2jp1 ,  p_im1jp1 ,  p_ip0jp1 ,  p_ip1jp1 )
          deallocate ( recvloc,recvlen,sendlen,sendloc )

          allocate (  old_dlam(lattice%imglobal)    ,  old_dphi(lattice%jmglobal)    )
          allocate ( wl_ip0jp0(im,jm) , wl_im1jp0(im,jm) )
          allocate ( wl_ip0jm1(im,jm) , wl_im1jm1(im,jm) )
          allocate ( wc_ip1jp1(im,jm) , wc_ip0jp1(im,jm) , wc_im1jp1(im,jm) , wc_im2jp1(im,jm) )
          allocate ( wc_ip1jp0(im,jm) , wc_ip0jp0(im,jm) , wc_im1jp0(im,jm) , wc_im2jp0(im,jm) )
          allocate ( wc_ip1jm1(im,jm) , wc_ip0jm1(im,jm) , wc_im1jm1(im,jm) , wc_im2jm1(im,jm) )
          allocate ( wc_ip1jm2(im,jm) , wc_ip0jm2(im,jm) , wc_im1jm2(im,jm) , wc_im2jm2(im,jm) )
          allocate (      msgn(im,jm) ,      nsgn(im,jm) )
          allocate (  n_im2jm2(im,jm) ,  n_im1jm2(im,jm) ,  n_ip0jm2(im,jm) ,  n_ip1jm2(im,jm) )
          allocate (  n_im2jm1(im,jm) ,  n_im1jm1(im,jm) ,  n_ip0jm1(im,jm) ,  n_ip1jm1(im,jm) )
          allocate (  n_im2jp0(im,jm) ,  n_im1jp0(im,jm) ,  n_ip0jp0(im,jm) ,  n_ip1jp0(im,jm) )
          allocate (  n_im2jp1(im,jm) ,  n_im1jp1(im,jm) ,  n_ip0jp1(im,jm) ,  n_ip1jp1(im,jm) )
          allocate (  p_im2jm2(im,jm) ,  p_im1jm2(im,jm) ,  p_ip0jm2(im,jm) ,  p_ip1jm2(im,jm) )
          allocate (  p_im2jm1(im,jm) ,  p_im1jm1(im,jm) ,  p_ip0jm1(im,jm) ,  p_ip1jm1(im,jm) )
          allocate (  p_im2jp0(im,jm) ,  p_im1jp0(im,jm) ,  p_ip0jp0(im,jm) ,  p_ip1jp0(im,jm) )
          allocate (  p_im2jp1(im,jm) ,  p_im1jp1(im,jm) ,  p_ip0jp1(im,jm) ,  p_ip1jp1(im,jm) )
          allocate ( recvloc(16*im*jm,0:lattice%nx*lattice%ny-1) )
          allocate ( recvlen(         0:lattice%nx*lattice%ny-1) )
          allocate ( sendlen(         0:lattice%nx*lattice%ny-1) )
          allocate ( sendloc(   im*jm,0:lattice%nx*lattice%ny-1) )

          do i=1,lattice%imglobal
          old_dlam(i) = 0.0
          enddo
          do j=1,lattice%jmglobal
          old_dphi(j) = 0.0
          enddo
          endif
      endif


c Compute Geophysical and Computational Lambda's and Phi's
c --------------------------------------------------------
      lon_geo(1) = -pi
      lon_cmp(1) = -pi
      do i=2,lattice%imglobal
      lon_geo(i) = lon_geo(i-1) + dl
      lon_cmp(i) = lon_cmp(i-1) + dlam(i-1)
      enddo
      lat_geo(1) = -pi*0.5
      lat_cmp(1) = -pi*0.5
      do j=2,lattice%jmglobal-1
      lat_geo(j) = lat_geo(j-1) + dp
      lat_cmp(j) = lat_cmp(j-1) + dphi(j-1)
      enddo
      lat_geo(lattice%jmglobal) =  pi*0.5
      lat_cmp(lattice%jmglobal) =  pi*0.5


c Check for Co-incident Grid-Point Latitude and Pole Locations
c ------------------------------------------------------------
                                          eps_np = 0.0
      do j=1,lattice%jmglobal
      phi_cmp = lat_cmp(j)*180./pi
      if( abs( phi_cmp-tilt ).lt.1.e-3 )  eps_np =  1.e-3
      if( tilt+eps_np .gt. 90. )          eps_np = -1.e-3
      enddo

      d      = dl*dp
      lam_np = pi/180.*rotation
      phi_np = pi/180.*(tilt+eps_np)
      lam_0  = pi/180.*precession

      undef = getcon('UNDEF')
      if( tilt.eq.90. ) then 
      cosnp = 0.0
      sinnp = 1.0
      else if(tilt.eq.-90.0) then
      cosnp = 0.0
      sinnp =-1.0
      else
      cosnp = cos(phi_np)
      sinnp = sin(phi_np)
      endif

c Determine if Weights Need to be Updated
c ---------------------------------------
      compute_weights =  rotation.ne.old_rotation .or.
     .                       tilt.ne.old_tilt     .or.
     .                 precession.ne.old_precession

      i = 1
      do while ( .not.compute_weights .and. i.le.lattice%imglobal )
      compute_weights = dlam(i).ne.old_dlam(i)
      i = i+1
      enddo
      j = 1
      do while ( .not.compute_weights .and. j.le.lattice%jmglobal-1 )
      compute_weights = dphi(j).ne.old_dphi(j)
      j = j+1
      enddo


c Computing Weights for Geophysical to Computational Grid Interpolation
c ---------------------------------------------------------------------
      if( compute_weights ) then

      old_rotation   = rotation
      old_tilt       = tilt
      old_precession = precession

      dlam_min = dlam(1)
      dlam_max = dlam(1)
      do i=1,lattice%imglobal
      old_dlam(i) =  dlam(i)
      if( dlam(i).lt.dlam_min ) dlam_min = dlam(i)
      if( dlam(i).gt.dlam_max ) dlam_max = dlam(i)
      enddo
      dphi_min = dphi(1)
      dphi_max = dphi(1)
      do j=1,lattice%jmglobal-1
      old_dphi(j) = dphi(j)
      if( dphi(j).lt.dphi_min ) dphi_min = dphi(j)
      if( dphi(j).gt.dphi_max ) dphi_max = dphi(j)
      enddo
      old_dphi(lattice%jmglobal) = dphi(lattice%jmglobal)

      if( lattice%myid.eq.0 ) then
      print *, 'Computing Weights for Geophysical to Computational Grid Interpolation'
      print *, '---------------------------------------------------------------------'
      print *, 'lam_np = ',rotation
      print *, 'phi_np = ',tilt,'   (eps_np = ',eps_np,')'
      print *, 'lam_0  = ',precession
      print *
      print *, 'Longitudinal Stretching Factor: ',dlam_max/dlam_min
      print *, ' Latitudinal Stretching Factor: ',dphi_max/dphi_min
      print *
      print *, ' Rotate_F Communication Matrix'
      print *, ' -----------------------------'
      endif

      allocate (      nn(16*im*jm,0:lattice%nx*lattice%ny-1) )
      allocate (     pos(16*im*jm,0:lattice%nx*lattice%ny-1) )
      allocate (     loc(16*im*jm,0:lattice%nx*lattice%ny-1) )
      allocate (     len(         0:lattice%nx*lattice%ny-1) )

      allocate (       ip1(im,jm) ,       ip0(im,jm) ,       im1(im,jm) ,       im2(im,jm) )
      allocate (       jp1(im,jm) ,       jp0(im,jm) ,       jm1(im,jm) ,       jm2(im,jm) )
      allocate (    im2jm2(im,jm) ,    im1jm2(im,jm) ,    ip0jm2(im,jm) ,    ip1jm2(im,jm) )
      allocate (    im2jm1(im,jm) ,    im1jm1(im,jm) ,    ip0jm1(im,jm) ,    ip1jm1(im,jm) )
      allocate (    im2jp0(im,jm) ,    im1jp0(im,jm) ,    ip0jp0(im,jm) ,    ip1jp0(im,jm) )
      allocate (    im2jp1(im,jm) ,    im1jp1(im,jm) ,    ip0jp1(im,jm) ,    ip1jp1(im,jm) )
      allocate ( pe_im2jm2(im,jm) , pe_im1jm2(im,jm) , pe_ip0jm2(im,jm) , pe_ip1jm2(im,jm) )
      allocate ( pe_im2jm1(im,jm) , pe_im1jm1(im,jm) , pe_ip0jm1(im,jm) , pe_ip1jm1(im,jm) )
      allocate ( pe_im2jp0(im,jm) , pe_im1jp0(im,jm) , pe_ip0jp0(im,jm) , pe_ip1jp0(im,jm) )
      allocate ( pe_im2jp1(im,jm) , pe_im1jp1(im,jm) , pe_ip0jp1(im,jm) , pe_ip1jp1(im,jm) )

      do j=1,jm
      do i=1,im
      lam_cmp = lon_cmp( lattice%iglobal(i) )
      phi_cmp = lat_cmp( lattice%jglobal(j) )

      p1 = cosnp*cos(phi_cmp)*cos(lam_cmp-lam_np)
     .   + sin(phi_cmp)*sinnp
      p1 = min(p1, 1.0)
      p1 = max(p1,-1.0)
      phi_geo = asin( p1 )

      if( tilt.eq.90.0 .or. tilt.eq.-90.0 ) then
      p2 = sinnp*cos(lam_cmp-lam_np)
      else
      p2 = sinnp*cos(phi_cmp)*cos(lam_cmp-lam_np)
     .   - sin(phi_cmp)*cosnp
      p2 = p2 / max( cos(phi_geo),eps )
      p2 = min(p2, 1.0)
      p2 = max(p2,-1.0)
      endif
      p2 = acos( p2 )

      p3 = cos(phi_cmp)*sin(lam_cmp-lam_np)
      if( p3.lt.0.0 ) p2 = -p2 
      p2 = p2 - lam_0
      lam_geo = mod( p2+3.0*pi,2.0*pi ) - pi

c Determine (i,j) Indexing Based on Geophysical Grid
c --------------------------------------------------
      im1_geo = 1
      do igeo = 2,lattice%imglobal
      if( lon_geo(igeo).lt.lam_geo ) im1_geo = igeo
      enddo
      jm1_geo = 1
      do jgeo = 2,lattice%jmglobal
      if( lat_geo(jgeo).lt.phi_geo ) jm1_geo = jgeo
      enddo

      im1(i,j) = im1_geo
      ip0(i,j) = im1(i,j) + 1
      ip1(i,j) = ip0(i,j) + 1
      im2(i,j) = im1(i,j) - 1

      jm1(i,j) = jm1_geo
      jp0(i,j) = jm1(i,j) + 1
      jp1(i,j) = jp0(i,j) + 1
      jm2(i,j) = jm1(i,j) - 1

c Fix Longitude Index Boundaries
c ------------------------------
      if(im1(i,j).eq.lattice%imglobal) then
      ip0(i,j) = 1
      ip1(i,j) = 2
      endif
      if(im1(i,j).eq.1) then
      im2(i,j) = lattice%imglobal
      endif
      if(ip0(i,j).eq.lattice%imglobal) then
      ip1(i,j) = 1
      endif

c Compute ij index
c ----------------
         msgn(i,j) = 2
      if( jp0(i,j).eq.lattice%jmglobal ) then
          jp1(i,j) =  lattice%jmglobal-1
       ip1jp1(i,j) = 1+mod(ip1(i,j)+lattice%imglobal/2-1,lattice%imglobal) + ( jp1(i,j)-1 )*lattice%imglobal
       ip0jp1(i,j) = 1+mod(ip0(i,j)+lattice%imglobal/2-1,lattice%imglobal) + ( jp1(i,j)-1 )*lattice%imglobal
       im1jp1(i,j) = 1+mod(im1(i,j)+lattice%imglobal/2-1,lattice%imglobal) + ( jp1(i,j)-1 )*lattice%imglobal
       im2jp1(i,j) = 1+mod(im2(i,j)+lattice%imglobal/2-1,lattice%imglobal) + ( jp1(i,j)-1 )*lattice%imglobal
         msgn(i,j) = 1
      else
       ip1jp1(i,j) = ip1(i,j) + ( jp1(i,j)-1 )*lattice%imglobal
       ip0jp1(i,j) = ip0(i,j) + ( jp1(i,j)-1 )*lattice%imglobal
       im1jp1(i,j) = im1(i,j) + ( jp1(i,j)-1 )*lattice%imglobal
       im2jp1(i,j) = im2(i,j) + ( jp1(i,j)-1 )*lattice%imglobal
      endif

      ip1jp0(i,j) = ip1(i,j) + ( jp0(i,j)-1 )*lattice%imglobal
      ip0jp0(i,j) = ip0(i,j) + ( jp0(i,j)-1 )*lattice%imglobal
      im1jp0(i,j) = im1(i,j) + ( jp0(i,j)-1 )*lattice%imglobal
      im2jp0(i,j) = im2(i,j) + ( jp0(i,j)-1 )*lattice%imglobal

      ip1jm1(i,j) = ip1(i,j) + ( jm1(i,j)-1 )*lattice%imglobal
      ip0jm1(i,j) = ip0(i,j) + ( jm1(i,j)-1 )*lattice%imglobal
      im1jm1(i,j) = im1(i,j) + ( jm1(i,j)-1 )*lattice%imglobal
      im2jm1(i,j) = im2(i,j) + ( jm1(i,j)-1 )*lattice%imglobal

         nsgn(i,j) = 2
      if( jm1(i,j).eq.1 ) then
          jm2(i,j) =  2
       ip1jm2(i,j) = 1+mod(ip1(i,j)+lattice%imglobal/2-1,lattice%imglobal) + ( jm2(i,j)-1 )*lattice%imglobal
       ip0jm2(i,j) = 1+mod(ip0(i,j)+lattice%imglobal/2-1,lattice%imglobal) + ( jm2(i,j)-1 )*lattice%imglobal
       im1jm2(i,j) = 1+mod(im1(i,j)+lattice%imglobal/2-1,lattice%imglobal) + ( jm2(i,j)-1 )*lattice%imglobal
       im2jm2(i,j) = 1+mod(im2(i,j)+lattice%imglobal/2-1,lattice%imglobal) + ( jm2(i,j)-1 )*lattice%imglobal
         nsgn(i,j) = 1
      else
       ip1jm2(i,j) = ip1(i,j) + ( jm2(i,j)-1 )*lattice%imglobal
       ip0jm2(i,j) = ip0(i,j) + ( jm2(i,j)-1 )*lattice%imglobal
       im1jm2(i,j) = im1(i,j) + ( jm2(i,j)-1 )*lattice%imglobal
       im2jm2(i,j) = im2(i,j) + ( jm2(i,j)-1 )*lattice%imglobal
      endif

c Determine PE location of Source Data
c ------------------------------------
      pej_jp1 = lattice%pejglobal( jp1(i,j) )
      pej_jp0 = lattice%pejglobal( jp0(i,j) )
      pej_jm1 = lattice%pejglobal( jm1(i,j) )
      pej_jm2 = lattice%pejglobal( jm2(i,j) )

      pei_ip1 = lattice%peiglobal( ip1(i,j) )
      pei_ip0 = lattice%peiglobal( ip0(i,j) )
      pei_im1 = lattice%peiglobal( im1(i,j) )
      pei_im2 = lattice%peiglobal( im2(i,j) )

      pe_ip1jp0(i,j) = pei_ip1 + pej_jp0*lattice%nx
      pe_ip0jp0(i,j) = pei_ip0 + pej_jp0*lattice%nx
      pe_im1jp0(i,j) = pei_im1 + pej_jp0*lattice%nx
      pe_im2jp0(i,j) = pei_im2 + pej_jp0*lattice%nx

      pe_ip1jm1(i,j) = pei_ip1 + pej_jm1*lattice%nx
      pe_ip0jm1(i,j) = pei_ip0 + pej_jm1*lattice%nx
      pe_im1jm1(i,j) = pei_im1 + pej_jm1*lattice%nx
      pe_im2jm1(i,j) = pei_im2 + pej_jm1*lattice%nx

            n = 1 + (ip1jp1(i,j)-1)/lattice%imglobal 
            m = ip1jp1(i,j) - (n-1)*lattice%imglobal 
      pei_ip1 = lattice%peiglobal( m )
            n = 1 + (ip0jp1(i,j)-1)/lattice%imglobal 
            m = ip0jp1(i,j) - (n-1)*lattice%imglobal 
      pei_ip0 = lattice%peiglobal( m )
            n = 1 + (im1jp1(i,j)-1)/lattice%imglobal 
            m = im1jp1(i,j) - (n-1)*lattice%imglobal 
      pei_im1 = lattice%peiglobal( m )
            n = 1 + (im2jp1(i,j)-1)/lattice%imglobal 
            m = im2jp1(i,j) - (n-1)*lattice%imglobal 
      pei_im2 = lattice%peiglobal( m )

      pe_ip1jp1(i,j) = pei_ip1 + pej_jp1*lattice%nx
      pe_ip0jp1(i,j) = pei_ip0 + pej_jp1*lattice%nx
      pe_im1jp1(i,j) = pei_im1 + pej_jp1*lattice%nx
      pe_im2jp1(i,j) = pei_im2 + pej_jp1*lattice%nx

            n = 1 + (ip1jm2(i,j)-1)/lattice%imglobal 
            m = ip1jm2(i,j) - (n-1)*lattice%imglobal 
      pei_ip1 = lattice%peiglobal( m )
            n = 1 + (ip0jm2(i,j)-1)/lattice%imglobal 
            m = ip0jm2(i,j) - (n-1)*lattice%imglobal 
      pei_ip0 = lattice%peiglobal( m )
            n = 1 + (im1jm2(i,j)-1)/lattice%imglobal 
            m = im1jm2(i,j) - (n-1)*lattice%imglobal 
      pei_im1 = lattice%peiglobal( m )
            n = 1 + (im2jm2(i,j)-1)/lattice%imglobal 
            m = im2jm2(i,j) - (n-1)*lattice%imglobal 
      pei_im2 = lattice%peiglobal( m )

      pe_ip1jm2(i,j) = pei_ip1 + pej_jm2*lattice%nx
      pe_ip0jm2(i,j) = pei_ip0 + pej_jm2*lattice%nx
      pe_im1jm2(i,j) = pei_im1 + pej_jm2*lattice%nx
      pe_im2jm2(i,j) = pei_im2 + pej_jm2*lattice%nx

c Compute Immediate Surrounding Coordinates on Uniform Grid
c ---------------------------------------------------------
      lam  =  lam_geo
      phi  =  phi_geo
      lam1 = -pi     + (im1(i,j)-1)*dl
      phi1 = -pi*0.5 + (jm1(i,j)-1)*dp
      lam2 =  lam1 + dl
      phi2 =  phi1 + dp

c Bi-Linear Weights
c -----------------
      wl_im1jm1(i,j) = (lam2-lam)*(phi2-phi)/d
      wl_ip0jm1(i,j) = (lam-lam1)*(phi2-phi)/d
      wl_im1jp0(i,j) = (lam2-lam)*(phi-phi1)/d
      wl_ip0jp0(i,j) = (lam-lam1)*(phi-phi1)/d

c Bi-Cubic Weights
c ----------------
      ux = (lam2-lam)/dl
      uy = (phi2-phi)/dp

      ap1 = -ux*(ux-1)*(ux-2)/6.
      ap0 = (ux-1)*(ux-2)*(ux+1)/2.
      am1 = -ux*(ux+1)*(ux-2)/2.
      am2 = ux*(ux-1)*(ux+1)/6.

      bp1 = -uy*(uy-1)*(uy-2)/6.
      bp0 = (uy-1)*(uy-2)*(uy+1)/2.
      bm1 = -uy*(uy+1)*(uy-2)/2.
      bm2 = uy*(uy-1)*(uy+1)/6.

      wc_ip1jp1(i,j) = bp1*ap1
      wc_ip0jp1(i,j) = bp1*ap0
      wc_im1jp1(i,j) = bp1*am1
      wc_im2jp1(i,j) = bp1*am2

      wc_ip1jp0(i,j) = bp0*ap1
      wc_ip0jp0(i,j) = bp0*ap0
      wc_im1jp0(i,j) = bp0*am1
      wc_im2jp0(i,j) = bp0*am2

      wc_ip1jm1(i,j) = bm1*ap1
      wc_ip0jm1(i,j) = bm1*ap0
      wc_im1jm1(i,j) = bm1*am1
      wc_im2jm1(i,j) = bm1*am2

      wc_ip1jm2(i,j) = bm2*ap1
      wc_ip0jm2(i,j) = bm2*ap0
      wc_im1jm2(i,j) = bm2*am1
      wc_im2jm2(i,j) = bm2*am2

      enddo
      enddo

c Create Data Index Array for Total and Non-Redundant Data
c Note:  recvlen(peid):   is the total number of unique grid-points needed from a given peid
c -----  recvloc(m,peid): is a list of the global locations of those points
c -------------------------------------------------------------------------
      allocate ( bucket(16*lattice%imglobal*lattice%jmglobal) )
      do peid = 0,lattice%nx*lattice%ny-1
          num = 0
         do j=1,jm
         do i=1,im
            if ( pe_im2jm2(i,j).eq.peid ) then
                                          num  = num + 1
                                          loc(num,peid) = im2jm2(i,j)
                                          pos(num,peid) = i+(j-1)*im
            endif
            if ( pe_im1jm2(i,j).eq.peid ) then
                                          num  = num + 1
                                          loc(num,peid) = im1jm2(i,j)
                                          pos(num,peid) = i+(j-1)*im
            endif
            if ( pe_ip0jm2(i,j).eq.peid ) then 
                                          num  = num + 1
                                          loc(num,peid) = ip0jm2(i,j)
                                          pos(num,peid) = i+(j-1)*im
            endif
            if ( pe_ip1jm2(i,j).eq.peid ) then 
                                          num  = num + 1
                                          loc(num,peid) = ip1jm2(i,j)
                                          pos(num,peid) = i+(j-1)*im
            endif

            if ( pe_im2jm1(i,j).eq.peid ) then
                                          num  = num + 1
                                          loc(num,peid) = im2jm1(i,j)
                                          pos(num,peid) = i+(j-1)*im
            endif
            if ( pe_im1jm1(i,j).eq.peid ) then
                                          num  = num + 1
                                          loc(num,peid) = im1jm1(i,j)
                                          pos(num,peid) = i+(j-1)*im
            endif
            if ( pe_ip0jm1(i,j).eq.peid ) then 
                                          num  = num + 1
                                          loc(num,peid) = ip0jm1(i,j)
                                          pos(num,peid) = i+(j-1)*im
            endif
            if ( pe_ip1jm1(i,j).eq.peid ) then 
                                          num  = num + 1
                                          loc(num,peid) = ip1jm1(i,j)
                                          pos(num,peid) = i+(j-1)*im
            endif

            if ( pe_im2jp0(i,j).eq.peid ) then 
                                          num  = num + 1
                                          loc(num,peid) = im2jp0(i,j)
                                          pos(num,peid) = i+(j-1)*im
            endif
            if ( pe_im1jp0(i,j).eq.peid ) then 
                                          num  = num + 1
                                          loc(num,peid) = im1jp0(i,j)
                                          pos(num,peid) = i+(j-1)*im
            endif
            if ( pe_ip0jp0(i,j).eq.peid ) then 
                                          num  = num + 1
                                          loc(num,peid) = ip0jp0(i,j)
                                          pos(num,peid) = i+(j-1)*im
            endif
            if ( pe_ip1jp0(i,j).eq.peid ) then 
                                          num  = num + 1
                                          loc(num,peid) = ip1jp0(i,j)
                                          pos(num,peid) = i+(j-1)*im
            endif

            if ( pe_im2jp1(i,j).eq.peid ) then 
                                          num  = num + 1
                                          loc(num,peid) = im2jp1(i,j)
                                          pos(num,peid) = i+(j-1)*im
            endif
            if ( pe_im1jp1(i,j).eq.peid ) then 
                                          num  = num + 1
                                          loc(num,peid) = im1jp1(i,j)
                                          pos(num,peid) = i+(j-1)*im
            endif
            if ( pe_ip0jp1(i,j).eq.peid ) then 
                                          num  = num + 1
                                          loc(num,peid) = ip0jp1(i,j)
                                          pos(num,peid) = i+(j-1)*im
            endif
            if ( pe_ip1jp1(i,j).eq.peid ) then 
                                          num  = num + 1
                                          loc(num,peid) = ip1jp1(i,j)
                                          pos(num,peid) = i+(j-1)*im
            endif
         enddo
         enddo
              len(peid) = num

c Strip redundancies from buffer
c ------------------------------
                        m = 0
         if( num.ne.0 ) then
                bucket(:) = 0
                do i=1,num
                bucket( loc(i,peid) ) = 1
                enddo
                do i=1,16*lattice%imglobal*lattice%jmglobal
                   if( bucket(i).eq.1 ) then
                                     m = m+1
                       recvloc(m,peid) = i
                            bucket (i) = m
                   endif
                enddo
                do i=1,num
                nn(i,peid) = bucket( loc(i,peid) )
                enddo
         endif

         recvlen(peid) = m
      enddo
      deallocate ( bucket )

c Store buffer locations
c ----------------------
      do peid = 0,lattice%nx*lattice%ny-1
         do m = 1,len(peid)
            n = nn (m,peid)
            j = 1 + (pos(m,peid)-1)/im
            i = pos(m,peid) - (j-1)*im
            if( recvloc(n,peid).eq.ip1jp1(i,j) ) then
                  n_ip1jp1(i,j) = n
                  p_ip1jp1(i,j) = peid
            endif
            if( recvloc(n,peid).eq.ip0jp1(i,j) ) then
                  n_ip0jp1(i,j) = n
                  p_ip0jp1(i,j) = peid
            endif
            if( recvloc(n,peid).eq.im1jp1(i,j) ) then
                  n_im1jp1(i,j) = n
                  p_im1jp1(i,j) = peid
            endif
            if( recvloc(n,peid).eq.im2jp1(i,j) ) then
                  n_im2jp1(i,j) = n
                  p_im2jp1(i,j) = peid
            endif

            if( recvloc(n,peid).eq.ip1jp0(i,j) ) then
                  n_ip1jp0(i,j) = n
                  p_ip1jp0(i,j) = peid
            endif
            if( recvloc(n,peid).eq.ip0jp0(i,j) ) then
                  n_ip0jp0(i,j) = n
                  p_ip0jp0(i,j) = peid
            endif
            if( recvloc(n,peid).eq.im1jp0(i,j) ) then
                  n_im1jp0(i,j) = n
                  p_im1jp0(i,j) = peid
            endif
            if( recvloc(n,peid).eq.im2jp0(i,j) ) then
                  n_im2jp0(i,j) = n
                  p_im2jp0(i,j) = peid
            endif

            if( recvloc(n,peid).eq.ip1jm1(i,j) ) then
                  n_ip1jm1(i,j) = n
                  p_ip1jm1(i,j) = peid
            endif
            if( recvloc(n,peid).eq.ip0jm1(i,j) ) then
                  n_ip0jm1(i,j) = n
                  p_ip0jm1(i,j) = peid
            endif
            if( recvloc(n,peid).eq.im1jm1(i,j) ) then
                  n_im1jm1(i,j) = n
                  p_im1jm1(i,j) = peid
            endif
            if( recvloc(n,peid).eq.im2jm1(i,j) ) then
                  n_im2jm1(i,j) = n
                  p_im2jm1(i,j) = peid
            endif

            if( recvloc(n,peid).eq.ip1jm2(i,j) ) then
                  n_ip1jm2(i,j) = n
                  p_ip1jm2(i,j) = peid
            endif
            if( recvloc(n,peid).eq.ip0jm2(i,j) ) then
                  n_ip0jm2(i,j) = n
                  p_ip0jm2(i,j) = peid
            endif
            if( recvloc(n,peid).eq.im1jm2(i,j) ) then
                  n_im1jm2(i,j) = n
                  p_im1jm2(i,j) = peid
            endif
            if( recvloc(n,peid).eq.im2jm2(i,j) ) then
                  n_im2jm2(i,j) = n
                  p_im2jm2(i,j) = peid
            endif
         enddo
      enddo

c Send and Recv the data length needed for interpolation
c ------------------------------------------------------
      do peid = 0,lattice%nx*lattice%ny-1
      if( lattice%myid.ne.peid ) then

          call mpi_isend ( recvlen(peid),1,mpi_integer,peid,peid,lattice%comm,sendquest(peid),ierror )

      else
          sendquest(peid) = mpi_request_null
      endif
      enddo

      do peid = 0,lattice%nx*lattice%ny-1
      if( lattice%myid.ne.peid ) then

          call mpi_irecv ( sendlen(peid),1,mpi_integer,peid,lattice%myid,lattice%comm,recvquest(peid),ierror )

      else
            sendlen(peid) = recvlen(peid)
          recvquest(peid) = mpi_request_null
      endif
      enddo

      call mpi_waitall ( lattice%nx*lattice%ny,sendquest(0:lattice%nx*lattice%ny-1),stats(1,0),ierror )
      call mpi_waitall ( lattice%nx*lattice%ny,recvquest(0:lattice%nx*lattice%ny-1),statr(1,0),ierror )

c Send and Recv the data locations needed for interpolation
c ---------------------------------------------------------
      do peid = 0,lattice%nx*lattice%ny-1
          sendquest(peid) = mpi_request_null
      if( lattice%myid.ne.peid ) then
          if( recvlen(peid).ne.0 ) then

          call mpi_isend ( recvloc(1,peid),recvlen(peid),mpi_integer,peid,peid,lattice%comm,sendquest(peid),ierror )

          endif
      endif
      enddo

      do peid = 0,lattice%nx*lattice%ny-1
         recvquest(peid) = mpi_request_null
         if( sendlen(peid).ne.0 ) then
             if( lattice%myid.ne.peid ) then

             call mpi_irecv ( sendloc(1,peid),sendlen(peid),mpi_integer,peid,lattice%myid,lattice%comm,recvquest(peid),ierror )

             else
             do n=1,sendlen(peid)
             sendloc(n,peid) = recvloc(n,peid)
             enddo
             endif
         endif
      enddo

      call mpi_waitall ( lattice%nx*lattice%ny,sendquest(0:lattice%nx*lattice%ny-1),stats(1,0),ierror )
      call mpi_waitall ( lattice%nx*lattice%ny,recvquest(0:lattice%nx*lattice%ny-1),statr(1,0),ierror )

c Print Communication Matrix
c --------------------------
      do n=0,lattice%nx*lattice%ny-1
      if( n.eq.lattice%myid ) then
          write(6,1001) (sendlen(m),m=0,lattice%nx*lattice%ny-1)
      endif
      call my_barrier (lattice%comm)
      enddo
      if( lattice%myid.eq.lattice%nx*lattice%ny-1 ) print *
 1001 format(128(1x,i4))

      deallocate ( nn,pos,loc,len )
      deallocate (       ip1 ,       ip0 ,       im1 ,       im2 )
      deallocate (       jp1 ,       jp0 ,       jm1 ,       jm2 )
      deallocate (    im2jm2 ,    im1jm2 ,    ip0jm2 ,    ip1jm2 )
      deallocate (    im2jm1 ,    im1jm1 ,    ip0jm1 ,    ip1jm1 )
      deallocate (    im2jp0 ,    im1jp0 ,    ip0jp0 ,    ip1jp0 )
      deallocate (    im2jp1 ,    im1jp1 ,    ip0jp1 ,    ip1jp1 )
      deallocate ( pe_im2jm2 , pe_im1jm2 , pe_ip0jm2 , pe_ip1jm2 )
      deallocate ( pe_im2jm1 , pe_im1jm1 , pe_ip0jm1 , pe_ip1jm1 )
      deallocate ( pe_im2jp0 , pe_im1jp0 , pe_ip0jp0 , pe_ip1jp0 )
      deallocate ( pe_im2jp1 , pe_im1jp1 , pe_ip0jp1 , pe_ip1jp1 )

      endif   ! End compute_weights check


C***********************************************************************
C*                     Begin Rotate Interpolation                      *
C***********************************************************************

c Send and Recv the data needed for interpolation
c -----------------------------------------------
      allocate ( sendbuf(im*jm*lm                    ,0:lattice%nx*lattice%ny-1) )
      allocate ( recvbuf(lattice%imax*lattice%jmax*lm,0:lattice%nx*lattice%ny-1) )

      do peid = 0,lattice%nx*lattice%ny-1
          sendquest(peid) = mpi_request_null
          if( sendlen(peid).ne.0 ) then
              if( lattice%myid.ne.peid ) then

                 m = 0
              do n = 1,sendlen(peid)
                 j = 1 + (sendloc(n,peid)-1)/lattice%imglobal
                 i = sendloc(n,peid) - (j-1)*lattice%imglobal
              do L = 1,lm
                 m = m+1
              sendbuf(m,peid) = q_geo( lattice%ilocal(i),lattice%jlocal(j),L )
              enddo
              enddo
              call mpi_isend ( sendbuf(1,peid),m,mpi_double_precision,peid,peid,lattice%comm,
     .                         sendquest(peid),ierror )

              else
                 m = 0
              do n = 1,sendlen(peid)
                 j = 1 + (sendloc(n,peid)-1)/lattice%imglobal
                 i = sendloc(n,peid) - (j-1)*lattice%imglobal
              do L = 1,lm
                 m = m+1
              recvbuf(m,peid) = q_geo( lattice%ilocal(i),lattice%jlocal(j),L )
              enddo
              enddo
          endif
      endif
      enddo

      do peid = 0,lattice%nx*lattice%ny-1
          recvquest(peid) = mpi_request_null
      if( lattice%myid.ne.peid ) then
          if( recvlen(peid).ne.0 ) then

          call mpi_irecv ( recvbuf(1,peid),recvlen(peid)*lm,mpi_double_precision,peid,lattice%myid,lattice%comm,
     .                     recvquest(peid),ierror )

          endif
      endif
      enddo

      call mpi_waitall ( lattice%nx*lattice%ny,sendquest(0:lattice%nx*lattice%ny-1),stats(1,0),ierror )
      call mpi_waitall ( lattice%nx*lattice%ny,recvquest(0:lattice%nx*lattice%ny-1),statr(1,0),ierror )

c Interpolate Geophysical Quantities to Computational Grid Using Bi-Linear
c ------------------------------------------------------------------------
      if( abs(norder).eq.1 ) then

      if( check  ) then
      do j=1,jm
      do i=1,im
      do L=1,lm
         if( recvbuf( L+(n_im1jm1(i,j)-1)*lm,p_im1jm1(i,j) ).ne.undef  .and.
     .       recvbuf( L+(n_ip0jm1(i,j)-1)*lm,p_ip0jm1(i,j) ).ne.undef  .and.
     .       recvbuf( L+(n_im1jp0(i,j)-1)*lm,p_im1jp0(i,j) ).ne.undef  .and.
     .       recvbuf( L+(n_ip0jp0(i,j)-1)*lm,p_ip0jp0(i,j) ).ne.undef ) then
         q_cmp(i,j,L) = wl_im1jm1(i,j) * recvbuf( L+(n_im1jm1(i,j)-1)*lm,p_im1jm1(i,j) )
     .                + wl_ip0jm1(i,j) * recvbuf( L+(n_ip0jm1(i,j)-1)*lm,p_ip0jm1(i,j) )
     .                + wl_im1jp0(i,j) * recvbuf( L+(n_im1jp0(i,j)-1)*lm,p_im1jp0(i,j) )
     .                + wl_ip0jp0(i,j) * recvbuf( L+(n_ip0jp0(i,j)-1)*lm,p_ip0jp0(i,j) )
         else
         q_cmp(i,j,L) = undef
         endif
      enddo
      enddo
      enddo
      endif      ! End UNDEF Check

      if( .not.check  ) then
      do j=1,jm
      do i=1,im
      do L=1,lm
         q_cmp(i,j,L) = wl_im1jm1(i,j) * recvbuf( L+(n_im1jm1(i,j)-1)*lm,p_im1jm1(i,j) )
     .                + wl_ip0jm1(i,j) * recvbuf( L+(n_ip0jm1(i,j)-1)*lm,p_ip0jm1(i,j) )
     .                + wl_im1jp0(i,j) * recvbuf( L+(n_im1jp0(i,j)-1)*lm,p_im1jp0(i,j) )
     .                + wl_ip0jp0(i,j) * recvbuf( L+(n_ip0jp0(i,j)-1)*lm,p_ip0jp0(i,j) )
      enddo
      enddo
      enddo
      endif      ! End UNDEF Check

      endif      ! End Check for Bi-Linear Interpolation

c Interpolate Geophysical Quantities to Computational Grid Using Bi-Cubic
c -----------------------------------------------------------------------
      if( abs(norder).eq.3 ) then

      if( check  ) then
      do j=1,jm
      do i=1,im
         sgnjp1 = sgn**msgn(i,j)
         sgnjm2 = sgn**nsgn(i,j)
      do L=1,lm

      if( recvbuf( L+(n_ip1jp1(i,j)-1)*lm,p_ip1jp1(i,j) ).ne.undef  .and.
     .    recvbuf( L+(n_ip0jp1(i,j)-1)*lm,p_ip0jp1(i,j) ).ne.undef  .and.
     .    recvbuf( L+(n_im1jp1(i,j)-1)*lm,p_im1jp1(i,j) ).ne.undef  .and.
     .    recvbuf( L+(n_im2jp1(i,j)-1)*lm,p_im2jp1(i,j) ).ne.undef  .and.

     .    recvbuf( L+(n_ip1jp0(i,j)-1)*lm,p_ip1jp0(i,j) ).ne.undef  .and.
     .    recvbuf( L+(n_ip0jp0(i,j)-1)*lm,p_ip0jp0(i,j) ).ne.undef  .and.
     .    recvbuf( L+(n_im1jp0(i,j)-1)*lm,p_im1jp0(i,j) ).ne.undef  .and.
     .    recvbuf( L+(n_im2jp0(i,j)-1)*lm,p_im2jp0(i,j) ).ne.undef  .and.

     .    recvbuf( L+(n_ip1jm1(i,j)-1)*lm,p_ip1jm1(i,j) ).ne.undef  .and.
     .    recvbuf( L+(n_ip0jm1(i,j)-1)*lm,p_ip0jm1(i,j) ).ne.undef  .and.
     .    recvbuf( L+(n_im1jm1(i,j)-1)*lm,p_im1jm1(i,j) ).ne.undef  .and.
     .    recvbuf( L+(n_im2jm1(i,j)-1)*lm,p_im2jm1(i,j) ).ne.undef  .and.

     .    recvbuf( L+(n_ip1jm2(i,j)-1)*lm,p_ip1jm2(i,j) ).ne.undef  .and.
     .    recvbuf( L+(n_ip0jm2(i,j)-1)*lm,p_ip0jm2(i,j) ).ne.undef  .and.
     .    recvbuf( L+(n_im1jm2(i,j)-1)*lm,p_im1jm2(i,j) ).ne.undef  .and.
     .    recvbuf( L+(n_im2jm2(i,j)-1)*lm,p_im2jm2(i,j) ).ne.undef ) then

      q_cmp(i,j,L) = wc_ip1jp1(i,j) * recvbuf( L+(n_ip1jp1(i,j)-1)*lm,p_ip1jp1(i,j) )*sgnjp1
     .             + wc_ip0jp1(i,j) * recvbuf( L+(n_ip0jp1(i,j)-1)*lm,p_ip0jp1(i,j) )*sgnjp1
     .             + wc_im1jp1(i,j) * recvbuf( L+(n_im1jp1(i,j)-1)*lm,p_im1jp1(i,j) )*sgnjp1
     .             + wc_im2jp1(i,j) * recvbuf( L+(n_im2jp1(i,j)-1)*lm,p_im2jp1(i,j) )*sgnjp1

     .             + wc_ip1jp0(i,j) * recvbuf( L+(n_ip1jp0(i,j)-1)*lm,p_ip1jp0(i,j) )
     .             + wc_ip0jp0(i,j) * recvbuf( L+(n_ip0jp0(i,j)-1)*lm,p_ip0jp0(i,j) )
     .             + wc_im1jp0(i,j) * recvbuf( L+(n_im1jp0(i,j)-1)*lm,p_im1jp0(i,j) )
     .             + wc_im2jp0(i,j) * recvbuf( L+(n_im2jp0(i,j)-1)*lm,p_im2jp0(i,j) )

     .             + wc_ip1jm1(i,j) * recvbuf( L+(n_ip1jm1(i,j)-1)*lm,p_ip1jm1(i,j) )
     .             + wc_ip0jm1(i,j) * recvbuf( L+(n_ip0jm1(i,j)-1)*lm,p_ip0jm1(i,j) )
     .             + wc_im1jm1(i,j) * recvbuf( L+(n_im1jm1(i,j)-1)*lm,p_im1jm1(i,j) )
     .             + wc_im2jm1(i,j) * recvbuf( L+(n_im2jm1(i,j)-1)*lm,p_im2jm1(i,j) )

     .             + wc_ip1jm2(i,j) * recvbuf( L+(n_ip1jm2(i,j)-1)*lm,p_ip1jm2(i,j) )*sgnjm2
     .             + wc_ip0jm2(i,j) * recvbuf( L+(n_ip0jm2(i,j)-1)*lm,p_ip0jm2(i,j) )*sgnjm2
     .             + wc_im1jm2(i,j) * recvbuf( L+(n_im1jm2(i,j)-1)*lm,p_im1jm2(i,j) )*sgnjm2
     .             + wc_im2jm2(i,j) * recvbuf( L+(n_im2jm2(i,j)-1)*lm,p_im2jm2(i,j) )*sgnjm2

      elseif( recvbuf( L+(n_im1jm1(i,j)-1)*lm,p_im1jm1(i,j) ).ne.undef  .and.
     .        recvbuf( L+(n_ip0jm1(i,j)-1)*lm,p_ip0jm1(i,j) ).ne.undef  .and.
     .        recvbuf( L+(n_im1jp0(i,j)-1)*lm,p_im1jp0(i,j) ).ne.undef  .and.
     .        recvbuf( L+(n_ip0jp0(i,j)-1)*lm,p_ip0jp0(i,j) ).ne.undef ) then

      q_cmp(i,j,L) = wl_im1jm1(i,j) * recvbuf( L+(n_im1jm1(i,j)-1)*lm,p_im1jm1(i,j) )
     .             + wl_ip0jm1(i,j) * recvbuf( L+(n_ip0jm1(i,j)-1)*lm,p_ip0jm1(i,j) )
     .             + wl_im1jp0(i,j) * recvbuf( L+(n_im1jp0(i,j)-1)*lm,p_im1jp0(i,j) )
     .             + wl_ip0jp0(i,j) * recvbuf( L+(n_ip0jp0(i,j)-1)*lm,p_ip0jp0(i,j) )

      else
      q_cmp(i,j,L) = undef
      endif

c Check for Positive Definite
c ---------------------------
      if( norder.lt.0 ) then
          if( q_cmp(i,j,L).ne.undef .and. q_cmp(i,j,L).lt.0.0 ) q_cmp(i,j,L) = 0.0
      endif

      enddo      ! End L Loop
      enddo      ! End i Loop
      enddo      ! End j Loop

      endif      ! End UNDEF Check

      if( .not.check  ) then
      do j=1,jm
      do i=1,im
         sgnjp1 = sgn**msgn(i,j)
         sgnjm2 = sgn**nsgn(i,j)
      do L=1,lm

      q_cmp(i,j,L) = wc_ip1jp1(i,j) * recvbuf( L+(n_ip1jp1(i,j)-1)*lm,p_ip1jp1(i,j) )*sgnjp1
     .             + wc_ip0jp1(i,j) * recvbuf( L+(n_ip0jp1(i,j)-1)*lm,p_ip0jp1(i,j) )*sgnjp1
     .             + wc_im1jp1(i,j) * recvbuf( L+(n_im1jp1(i,j)-1)*lm,p_im1jp1(i,j) )*sgnjp1
     .             + wc_im2jp1(i,j) * recvbuf( L+(n_im2jp1(i,j)-1)*lm,p_im2jp1(i,j) )*sgnjp1

     .             + wc_ip1jp0(i,j) * recvbuf( L+(n_ip1jp0(i,j)-1)*lm,p_ip1jp0(i,j) )
     .             + wc_ip0jp0(i,j) * recvbuf( L+(n_ip0jp0(i,j)-1)*lm,p_ip0jp0(i,j) )
     .             + wc_im1jp0(i,j) * recvbuf( L+(n_im1jp0(i,j)-1)*lm,p_im1jp0(i,j) )
     .             + wc_im2jp0(i,j) * recvbuf( L+(n_im2jp0(i,j)-1)*lm,p_im2jp0(i,j) )

     .             + wc_ip1jm1(i,j) * recvbuf( L+(n_ip1jm1(i,j)-1)*lm,p_ip1jm1(i,j) )
     .             + wc_ip0jm1(i,j) * recvbuf( L+(n_ip0jm1(i,j)-1)*lm,p_ip0jm1(i,j) )
     .             + wc_im1jm1(i,j) * recvbuf( L+(n_im1jm1(i,j)-1)*lm,p_im1jm1(i,j) )
     .             + wc_im2jm1(i,j) * recvbuf( L+(n_im2jm1(i,j)-1)*lm,p_im2jm1(i,j) )

     .             + wc_ip1jm2(i,j) * recvbuf( L+(n_ip1jm2(i,j)-1)*lm,p_ip1jm2(i,j) )*sgnjm2
     .             + wc_ip0jm2(i,j) * recvbuf( L+(n_ip0jm2(i,j)-1)*lm,p_ip0jm2(i,j) )*sgnjm2
     .             + wc_im1jm2(i,j) * recvbuf( L+(n_im1jm2(i,j)-1)*lm,p_im1jm2(i,j) )*sgnjm2
     .             + wc_im2jm2(i,j) * recvbuf( L+(n_im2jm2(i,j)-1)*lm,p_im2jm2(i,j) )*sgnjm2

c Check for Positive Definite
c ---------------------------
      if( norder.lt.0 ) then
          if( q_cmp(i,j,L).ne.undef .and. q_cmp(i,j,L).lt.0.0 )  q_cmp(i,j,L) = 0.0
      endif

      enddo      ! End L Loop
      enddo      ! End i Loop
      enddo      ! End j Loop

      endif      ! End UNDEF Check
      endif      ! End Check for Bi-Cubic Interpolation

      deallocate ( sendbuf )
      deallocate ( recvbuf )

      return
      end

      subroutine rotate_b ( q_cmp,q_geo,im,jm,lm,
     .                      dlam,dphi,rotation,tilt,precession,
     .                      sgn,norder,check,lattice )
C***********************************************************************
C
C  PURPOSE:
C  ========
C    Compute a backward transformation from the computational grid to the
C    geophysical grid.
C
C  INPUT:
C  ======
C    q_cmp ...... Field q_cmp(im,jm,lm) on the computational grid
C    im ......... Longitudinal dimension of q_cmp
C    jm ......... Latitudinal  dimension of q_cmp
C    lm ......... Vertical     dimension of q_cmp
C    dlam ....... Computational Grid Delta Lambda
C    dphi ....... Computational Grid Delta Phi
C    rotation ... Rotation parameter lam_np (Degrees)
C    tilt ....... Rotation parameter phi_np (Degrees)
C    precession . Rotation parameter lam_0  (Degrees)
C    sgn ........ Flag for scalar field  ( msgn =  1 )
C                    or vector component ( msgn = -1 )
C    norder ..... Order of Interpolation:  Bi-Linear => abs(norder) = 1
C                                          Bi-Cubic  => abs(norder) = 3
C                 Note: If norder < 0, then check for positive definite
C    check ...... Logical Flag to check for Undefined values
C
C  OUTPUT:
C  =======
C    q_geo ...... Field q_geo(im,jm,lm) on the geophysical grid
C
C
C***********************************************************************
C*                  GODDARD LABORATORY FOR ATMOSPHERES                 *
C***********************************************************************

      use g3_dynamics_lattice_module
      implicit none
      type ( dynamics_lattice_type ) lattice

      include 'mpif.h'
      integer  status(mpi_status_size)
      integer   stats(mpi_status_size,0:lattice%nx*lattice%ny-1)
      integer   statr(mpi_status_size,0:lattice%nx*lattice%ny-1)

      integer  sendquest(0:lattice%nx*lattice%ny-1)
      integer  recvquest(0:lattice%nx*lattice%ny-1)

      integer  ierror
      integer, allocatable, save :: sendlen(:)    ! Number of    unique           grid-points needed to SEND to   other PEs
      integer, allocatable, save :: recvlen(:)    ! Number of    unique           grid-points needed to RECV from other PEs
      integer, allocatable, save :: sendloc(:,:)  ! Global index locations of the grid-points needed to SEND to   other PEs
      integer, allocatable, save :: recvloc(:,:)  ! Global index locations of the grid-points needed to RECV from other PEs

      integer, allocatable       :: bucket (:)    ! Dummy bucket for array reduction
      real(kind=8),    allocatable       :: sendbuf(:,:)  ! Dummy buffer to SEND data
      real(kind=8),    allocatable       :: recvbuf(:,:)  ! Dummy buffer to RECV data

c Input Variables
c ---------------
      integer im,jm,lm,norder,sgn
      real(kind=8)    q_geo(im,jm,lm)
      real(kind=8)    q_cmp(im,jm,lm)
      real(kind=8)     dlam(lattice%imglobal)
      real(kind=8)     dphi(lattice%jmglobal)
      real(kind=8)     rotation, tilt, precession
      logical  check

c Local Variables
c ---------------
      integer  i,j,l,m,n,sgnjp1
      integer  peid, num,sgnjm2
      integer  pei_im2
      integer  pei_im1
      integer  pei_ip0
      integer  pei_ip1
      integer  pej_jm2
      integer  pej_jm1
      integer  pej_jp0
      integer  pej_jp1

c Index Locations
c ---------------
      integer, allocatable ::         len(:)  ! Total number of grid-points needed from other PEs
      integer, allocatable ::        nn(:,:)  ! Mapping from total to unique grid-points
      integer, allocatable ::       pos(:,:)  ! Position index of data to RECV from other PEs
      integer, allocatable ::       loc(:,:)

      integer, allocatable ::       ip1(:,:),       ip0(:,:),       im1(:,:),       im2(:,:)
      integer, allocatable ::       jp1(:,:),       jp0(:,:),       jm1(:,:),       jm2(:,:)
      integer, allocatable ::    im2jm2(:,:),    im2jm1(:,:),    im2jp0(:,:),    im2jp1(:,:)
      integer, allocatable ::    im1jm2(:,:),    im1jm1(:,:),    im1jp0(:,:),    im1jp1(:,:)
      integer, allocatable ::    ip0jm2(:,:),    ip0jm1(:,:),    ip0jp0(:,:),    ip0jp1(:,:)
      integer, allocatable ::    ip1jm2(:,:),    ip1jm1(:,:),    ip1jp0(:,:),    ip1jp1(:,:)
      integer, allocatable :: pe_im2jm2(:,:), pe_im2jm1(:,:), pe_im2jp0(:,:), pe_im2jp1(:,:)
      integer, allocatable :: pe_im1jm2(:,:), pe_im1jm1(:,:), pe_im1jp0(:,:), pe_im1jp1(:,:)
      integer, allocatable :: pe_ip0jm2(:,:), pe_ip0jm1(:,:), pe_ip0jp0(:,:), pe_ip0jp1(:,:)
      integer, allocatable :: pe_ip1jm2(:,:), pe_ip1jm1(:,:), pe_ip1jp0(:,:), pe_ip1jp1(:,:)

      integer, allocatable, save ::     msgn(:,:)
      integer, allocatable, save ::     nsgn(:,:)
      integer, allocatable, save :: n_im2jm2(:,:), n_im2jm1(:,:), n_im2jp0(:,:), n_im2jp1(:,:)
      integer, allocatable, save :: n_im1jm2(:,:), n_im1jm1(:,:), n_im1jp0(:,:), n_im1jp1(:,:)
      integer, allocatable, save :: n_ip0jm2(:,:), n_ip0jm1(:,:), n_ip0jp0(:,:), n_ip0jp1(:,:)
      integer, allocatable, save :: n_ip1jm2(:,:), n_ip1jm1(:,:), n_ip1jp0(:,:), n_ip1jp1(:,:)
      integer, allocatable, save :: p_im2jm2(:,:), p_im2jm1(:,:), p_im2jp0(:,:), p_im2jp1(:,:)
      integer, allocatable, save :: p_im1jm2(:,:), p_im1jm1(:,:), p_im1jp0(:,:), p_im1jp1(:,:)
      integer, allocatable, save :: p_ip0jm2(:,:), p_ip0jm1(:,:), p_ip0jp0(:,:), p_ip0jp1(:,:)
      integer, allocatable, save :: p_ip1jm2(:,:), p_ip1jm1(:,:), p_ip1jp0(:,:), p_ip1jp1(:,:)

c Bi-Linear Weights
c -----------------
      real(kind=8), allocatable, save :: wl_ip0jp0 (:,:)
      real(kind=8), allocatable, save :: wl_im1jp0 (:,:)
      real(kind=8), allocatable, save :: wl_ip0jm1 (:,:)
      real(kind=8), allocatable, save :: wl_im1jm1 (:,:)

c Bi-Cubic Weights
c ----------------
      real(kind=8), allocatable, save :: wc_ip1jp1 (:,:)
      real(kind=8), allocatable, save :: wc_ip0jp1 (:,:)
      real(kind=8), allocatable, save :: wc_im1jp1 (:,:)
      real(kind=8), allocatable, save :: wc_im2jp1 (:,:)
      real(kind=8), allocatable, save :: wc_ip1jp0 (:,:)
      real(kind=8), allocatable, save :: wc_ip0jp0 (:,:)
      real(kind=8), allocatable, save :: wc_im1jp0 (:,:)
      real(kind=8), allocatable, save :: wc_im2jp0 (:,:)
      real(kind=8), allocatable, save :: wc_ip1jm1 (:,:)
      real(kind=8), allocatable, save :: wc_ip0jm1 (:,:)
      real(kind=8), allocatable, save :: wc_im1jm1 (:,:)
      real(kind=8), allocatable, save :: wc_im2jm1 (:,:)
      real(kind=8), allocatable, save :: wc_ip1jm2 (:,:)
      real(kind=8), allocatable, save :: wc_ip0jm2 (:,:)
      real(kind=8), allocatable, save :: wc_im1jm2 (:,:)
      real(kind=8), allocatable, save :: wc_im2jm2 (:,:)

      real(kind=8), allocatable, save ::    old_dlam (:)
      real(kind=8), allocatable, save ::    old_dphi (:)

      real(kind=8)    ux, ap1, ap0, am1, am2
      real(kind=8)    uy, bp1, bp0, bm1, bm2

      real(kind=8)    lon_cmp(lattice%imglobal), lon_geo(lattice%imglobal)
      real(kind=8)    lat_cmp(lattice%jmglobal), lat_geo(lattice%jmglobal)
      real(kind=8)    pi,cosnp,sinnp,p1,p2,p3,eps,d
      real(kind=8)    lam,lam_ip1,lam_ip0,lam_im1,lam_im2
      real(kind=8)    phi,phi_jp1,phi_jp0,phi_jm1,phi_jm2
      real(kind=8)    dl,dp,lam_np,phi_np,lam_0,eps_np
      real(kind=8)    lam_geo, lam_cmp, dlam_max, dlam_min
      real(kind=8)    phi_geo, phi_cmp, dphi_max, dphi_min
      real(kind=8)    undef, getcon
      integer im1_cmp,icmp
      integer jm1_cmp,jcmp

      logical compute_weights
      real(kind=8)    old_rotation
      real(kind=8)    old_tilt
      real(kind=8)    old_precession
      data    old_rotation   /-999.9/
      data    old_tilt       /-999.9/
      data    old_precession /-999.9/

      parameter ( eps = 1.e-10 )

c Initialization
c --------------
      call my_barrier (lattice%comm)  ! Uncomment for more accurate timing

      pi = 4.*atan(1.)
      dl = 2*pi/ lattice%imglobal     ! Uniform Grid Delta Lambda
      dp =   pi/(lattice%jmglobal-1)  ! Uniform Grid Delta Phi

c Allocate Memory for Weights and Index Locations
c -----------------------------------------------
      if(.not.allocated(old_dlam)) then

      allocate (  old_dlam(lattice%imglobal)    ,  old_dphi(lattice%jmglobal)    )
      allocate ( wl_ip0jp0(im,jm) , wl_im1jp0(im,jm) )
      allocate ( wl_ip0jm1(im,jm) , wl_im1jm1(im,jm) )
      allocate ( wc_ip1jp1(im,jm) , wc_ip0jp1(im,jm) , wc_im1jp1(im,jm) , wc_im2jp1(im,jm) )
      allocate ( wc_ip1jp0(im,jm) , wc_ip0jp0(im,jm) , wc_im1jp0(im,jm) , wc_im2jp0(im,jm) )
      allocate ( wc_ip1jm1(im,jm) , wc_ip0jm1(im,jm) , wc_im1jm1(im,jm) , wc_im2jm1(im,jm) )
      allocate ( wc_ip1jm2(im,jm) , wc_ip0jm2(im,jm) , wc_im1jm2(im,jm) , wc_im2jm2(im,jm) )
      allocate (      msgn(im,jm) ,      nsgn(im,jm) )
      allocate (  n_im2jm2(im,jm) ,  n_im1jm2(im,jm) ,  n_ip0jm2(im,jm) ,  n_ip1jm2(im,jm) )
      allocate (  n_im2jm1(im,jm) ,  n_im1jm1(im,jm) ,  n_ip0jm1(im,jm) ,  n_ip1jm1(im,jm) )
      allocate (  n_im2jp0(im,jm) ,  n_im1jp0(im,jm) ,  n_ip0jp0(im,jm) ,  n_ip1jp0(im,jm) )
      allocate (  n_im2jp1(im,jm) ,  n_im1jp1(im,jm) ,  n_ip0jp1(im,jm) ,  n_ip1jp1(im,jm) )
      allocate (  p_im2jm2(im,jm) ,  p_im1jm2(im,jm) ,  p_ip0jm2(im,jm) ,  p_ip1jm2(im,jm) )
      allocate (  p_im2jm1(im,jm) ,  p_im1jm1(im,jm) ,  p_ip0jm1(im,jm) ,  p_ip1jm1(im,jm) )
      allocate (  p_im2jp0(im,jm) ,  p_im1jp0(im,jm) ,  p_ip0jp0(im,jm) ,  p_ip1jp0(im,jm) )
      allocate (  p_im2jp1(im,jm) ,  p_im1jp1(im,jm) ,  p_ip0jp1(im,jm) ,  p_ip1jp1(im,jm) )
      allocate ( recvloc(16*im*jm,0:lattice%nx*lattice%ny-1) )
      allocate ( recvlen(         0:lattice%nx*lattice%ny-1) )
      allocate ( sendlen(         0:lattice%nx*lattice%ny-1) )
      allocate ( sendloc(   im*jm,0:lattice%nx*lattice%ny-1) )
c     allocate ( sendloc(lattice%imglobal*lattice%jmglobal,0:lattice%nx*lattice%ny-1) )

      do i=1,lattice%imglobal
      old_dlam(i) = 0.0
      enddo
      do j=1,lattice%jmglobal
      old_dphi(j) = 0.0
      enddo

      else
             i =  size (old_dlam)
             j =  size (old_dphi)
          if(i.ne.lattice%imglobal .or. j.ne.lattice%jmglobal) then
          deallocate (  old_dlam ,  old_dphi )
          deallocate ( wl_ip0jp0 , wl_im1jp0 )
          deallocate ( wl_ip0jm1 , wl_im1jm1 )
          deallocate ( wc_ip1jp1 , wc_ip0jp1 , wc_im1jp1 , wc_im2jp1 )
          deallocate ( wc_ip1jp0 , wc_ip0jp0 , wc_im1jp0 , wc_im2jp0 )
          deallocate ( wc_ip1jm1 , wc_ip0jm1 , wc_im1jm1 , wc_im2jm1 )
          deallocate ( wc_ip1jm2 , wc_ip0jm2 , wc_im1jm2 , wc_im2jm2 )
          deallocate (      msgn ,      nsgn )
          deallocate (  n_im2jm2 ,  n_im1jm2 ,  n_ip0jm2 ,  n_ip1jm2 )
          deallocate (  n_im2jm1 ,  n_im1jm1 ,  n_ip0jm1 ,  n_ip1jm1 )
          deallocate (  n_im2jp0 ,  n_im1jp0 ,  n_ip0jp0 ,  n_ip1jp0 )
          deallocate (  n_im2jp1 ,  n_im1jp1 ,  n_ip0jp1 ,  n_ip1jp1 )
          deallocate (  p_im2jm2 ,  p_im1jm2 ,  p_ip0jm2 ,  p_ip1jm2 )
          deallocate (  p_im2jm1 ,  p_im1jm1 ,  p_ip0jm1 ,  p_ip1jm1 )
          deallocate (  p_im2jp0 ,  p_im1jp0 ,  p_ip0jp0 ,  p_ip1jp0 )
          deallocate (  p_im2jp1 ,  p_im1jp1 ,  p_ip0jp1 ,  p_ip1jp1 )
          deallocate ( recvloc,recvlen,sendlen,sendloc )

          allocate (  old_dlam(lattice%imglobal)    ,  old_dphi(lattice%jmglobal)    )
          allocate ( wl_ip0jp0(im,jm) , wl_im1jp0(im,jm) )
          allocate ( wl_ip0jm1(im,jm) , wl_im1jm1(im,jm) )
          allocate ( wc_ip1jp1(im,jm) , wc_ip0jp1(im,jm) , wc_im1jp1(im,jm) , wc_im2jp1(im,jm) )
          allocate ( wc_ip1jp0(im,jm) , wc_ip0jp0(im,jm) , wc_im1jp0(im,jm) , wc_im2jp0(im,jm) )
          allocate ( wc_ip1jm1(im,jm) , wc_ip0jm1(im,jm) , wc_im1jm1(im,jm) , wc_im2jm1(im,jm) )
          allocate ( wc_ip1jm2(im,jm) , wc_ip0jm2(im,jm) , wc_im1jm2(im,jm) , wc_im2jm2(im,jm) )
          allocate (      msgn(im,jm) ,      nsgn(im,jm) )
          allocate (  n_im2jm2(im,jm) ,  n_im1jm2(im,jm) ,  n_ip0jm2(im,jm) ,  n_ip1jm2(im,jm) )
          allocate (  n_im2jm1(im,jm) ,  n_im1jm1(im,jm) ,  n_ip0jm1(im,jm) ,  n_ip1jm1(im,jm) )
          allocate (  n_im2jp0(im,jm) ,  n_im1jp0(im,jm) ,  n_ip0jp0(im,jm) ,  n_ip1jp0(im,jm) )
          allocate (  n_im2jp1(im,jm) ,  n_im1jp1(im,jm) ,  n_ip0jp1(im,jm) ,  n_ip1jp1(im,jm) )
          allocate (  p_im2jm2(im,jm) ,  p_im1jm2(im,jm) ,  p_ip0jm2(im,jm) ,  p_ip1jm2(im,jm) )
          allocate (  p_im2jm1(im,jm) ,  p_im1jm1(im,jm) ,  p_ip0jm1(im,jm) ,  p_ip1jm1(im,jm) )
          allocate (  p_im2jp0(im,jm) ,  p_im1jp0(im,jm) ,  p_ip0jp0(im,jm) ,  p_ip1jp0(im,jm) )
          allocate (  p_im2jp1(im,jm) ,  p_im1jp1(im,jm) ,  p_ip0jp1(im,jm) ,  p_ip1jp1(im,jm) )
          allocate ( recvloc(16*im*jm,0:lattice%nx*lattice%ny-1) )
          allocate ( recvlen(         0:lattice%nx*lattice%ny-1) )
          allocate ( sendlen(         0:lattice%nx*lattice%ny-1) )
          allocate ( sendloc(   im*jm,0:lattice%nx*lattice%ny-1) )

          do i=1,lattice%imglobal
          old_dlam(i) = 0.0
          enddo
          do j=1,lattice%jmglobal
          old_dphi(j) = 0.0
          enddo
          endif
      endif


c Compute Geophysical and Computational Lambda's and Phi's
c --------------------------------------------------------
      lon_geo(1) = -pi
      lon_cmp(1) = -pi
      do i=2,lattice%imglobal
      lon_geo(i) = lon_geo(i-1) + dl
      lon_cmp(i) = lon_cmp(i-1) + dlam(i-1)
      enddo
      lat_geo(1) = -pi*0.5
      lat_cmp(1) = -pi*0.5
      do j=2,lattice%jmglobal-1
      lat_geo(j) = lat_geo(j-1) + dp
      lat_cmp(j) = lat_cmp(j-1) + dphi(j-1)
      enddo
      lat_geo(lattice%jmglobal) =  pi*0.5
      lat_cmp(lattice%jmglobal) =  pi*0.5

c Check for Co-incident Grid-Point Latitude and Pole Locations
c ------------------------------------------------------------
                                          eps_np = 0.0
      do j=1,lattice%jmglobal
      phi_cmp = lat_cmp(j)*180./pi
      if( abs( phi_cmp-tilt ).lt.1.e-3 )  eps_np =  1.e-3
      if( tilt+eps_np .gt. 90. )          eps_np = -1.e-3
      enddo

      lam_np = pi/180.*rotation
      phi_np = pi/180.*(tilt+eps_np)
      lam_0  = pi/180.*precession

      undef = getcon('UNDEF')
      if( tilt.eq.90. ) then 
      cosnp = 0.0
      sinnp = 1.0
      else if(tilt.eq.-90.0) then
      cosnp = 0.0
      sinnp =-1.0
      else
      cosnp = cos(phi_np)
      sinnp = sin(phi_np)
      endif

c Determine if Weights Need to be Updated
c ---------------------------------------
      compute_weights =  rotation.ne.old_rotation .or.
     .                       tilt.ne.old_tilt     .or.
     .                 precession.ne.old_precession

      i = 1
      do while ( .not.compute_weights .and. i.le.lattice%imglobal )
      compute_weights = dlam(i).ne.old_dlam(i)
      i = i+1
      enddo
      j = 1
      do while ( .not.compute_weights .and. j.le.lattice%jmglobal-1 )
      compute_weights = dphi(j).ne.old_dphi(j)
      j = j+1
      enddo

c Compute Weights for Computational to Geophysical Grid Interpolation
c -------------------------------------------------------------------
      if( compute_weights ) then

      old_rotation   = rotation
      old_tilt       = tilt
      old_precession = precession

      dlam_min = dlam(1)
      dlam_max = dlam(1)
      do i=1,lattice%imglobal
      old_dlam(i) =  dlam(i)
      if( dlam(i).lt.dlam_min ) dlam_min = dlam(i)
      if( dlam(i).gt.dlam_max ) dlam_max = dlam(i)
      enddo
      dphi_min = dphi(1)
      dphi_max = dphi(1)
      do j=1,lattice%jmglobal-1
      old_dphi(j) = dphi(j)
      if( dphi(j).lt.dphi_min ) dphi_min = dphi(j)
      if( dphi(j).gt.dphi_max ) dphi_max = dphi(j)
      enddo
      old_dphi(lattice%jmglobal) = dphi(lattice%jmglobal)

      if( lattice%myid.eq.0 ) then
      print *, 'Computing Weights for Computational to Geophysical Grid Interpolation'
      print *, '---------------------------------------------------------------------'
      print *, 'lam_np = ',rotation
      print *, 'phi_np = ',tilt,'   (eps_np = ',eps_np,')'
      print *, 'lam_0  = ',precession
      print *
      print *, 'Longitudinal Stretching Factor: ',dlam_max/dlam_min
      print *, ' Latitudinal Stretching Factor: ',dphi_max/dphi_min
      print *
      print *, ' Rotate_B Communication Matrix'
      print *, ' -----------------------------'
      endif

      allocate (      nn(16*im*jm,0:lattice%nx*lattice%ny-1) )
      allocate (     pos(16*im*jm,0:lattice%nx*lattice%ny-1) )
      allocate (     loc(16*im*jm,0:lattice%nx*lattice%ny-1) )
      allocate (     len(         0:lattice%nx*lattice%ny-1) )

      allocate (       ip1(im,jm) ,       ip0(im,jm) ,       im1(im,jm) ,       im2(im,jm) )
      allocate (       jp1(im,jm) ,       jp0(im,jm) ,       jm1(im,jm) ,       jm2(im,jm) )
      allocate (    im2jm2(im,jm) ,    im1jm2(im,jm) ,    ip0jm2(im,jm) ,    ip1jm2(im,jm) )
      allocate (    im2jm1(im,jm) ,    im1jm1(im,jm) ,    ip0jm1(im,jm) ,    ip1jm1(im,jm) )
      allocate (    im2jp0(im,jm) ,    im1jp0(im,jm) ,    ip0jp0(im,jm) ,    ip1jp0(im,jm) )
      allocate (    im2jp1(im,jm) ,    im1jp1(im,jm) ,    ip0jp1(im,jm) ,    ip1jp1(im,jm) )
      allocate ( pe_im2jm2(im,jm) , pe_im1jm2(im,jm) , pe_ip0jm2(im,jm) , pe_ip1jm2(im,jm) )
      allocate ( pe_im2jm1(im,jm) , pe_im1jm1(im,jm) , pe_ip0jm1(im,jm) , pe_ip1jm1(im,jm) )
      allocate ( pe_im2jp0(im,jm) , pe_im1jp0(im,jm) , pe_ip0jp0(im,jm) , pe_ip1jp0(im,jm) )
      allocate ( pe_im2jp1(im,jm) , pe_im1jp1(im,jm) , pe_ip0jp1(im,jm) , pe_ip1jp1(im,jm) )

      do j=1,jm
      do i=1,im
      lam_geo = lon_geo( lattice%iglobal(i) )
      phi_geo = lat_geo( lattice%jglobal(j) )

      p1 = cosnp*cos(phi_geo)*cos(lam_geo+lam_0-pi)
     .   + sin(phi_geo)*sinnp
      p1 = min(p1, 1.0)
      p1 = max(p1,-1.0)
      phi_cmp = asin( p1 )

      if( tilt.eq.90.0 .or. tilt.eq.-90.0 ) then
      p2 = sinnp*cos(lam_geo+lam_0-pi)
      else
      p2 = sinnp*cos(phi_geo)*cos(lam_geo+lam_0-pi)
     .   - sin(phi_geo)*cosnp
      p2 = p2 / max( cos(phi_cmp),eps )
      p2 = min(p2, 1.0)
      p2 = max(p2,-1.0)
      endif
      p2 = acos( p2 )

      p3 = cos(phi_geo)*sin(lam_geo+lam_0-pi)
      if( p3.lt.0.0 ) p2 = -p2
      p2 = p2 + lam_np - pi
      lam_cmp = mod( p2+3.0*pi,2.0*pi ) - pi

c Determine (i,j) Indexing Based on Computational Grid
c ----------------------------------------------------
      im1_cmp = 1
      do icmp = 2,lattice%imglobal
      if( lon_cmp(icmp).lt.lam_cmp ) im1_cmp = icmp
      enddo
      jm1_cmp = 1
      do jcmp = 2,lattice%jmglobal
      if( lat_cmp(jcmp).lt.phi_cmp ) jm1_cmp = jcmp
      enddo

      im1(i,j) = im1_cmp
      ip0(i,j) = im1(i,j) + 1
      ip1(i,j) = ip0(i,j) + 1
      im2(i,j) = im1(i,j) - 1

      jm1(i,j) = jm1_cmp
      jp0(i,j) = jm1(i,j) + 1
      jp1(i,j) = jp0(i,j) + 1
      jm2(i,j) = jm1(i,j) - 1

c Fix Longitude Index Boundaries
c ------------------------------
      if(im1(i,j).eq.lattice%imglobal) then
      ip0(i,j) = 1
      ip1(i,j) = 2
      endif
      if(im1(i,j).eq.1) then
      im2(i,j) = lattice%imglobal
      endif
      if(ip0(i,j).eq.lattice%imglobal) then
      ip1(i,j) = 1
      endif

c Compume Immediate Surrounding Coordinates
c -----------------------------------------
      lam     =  lam_cmp
      phi     =  phi_cmp

c Compute and Adjust Longitude Weights
c ------------------------------------
      lam_im2 =  lon_cmp(im2(i,j))
      lam_im1 =  lon_cmp(im1(i,j))
      lam_ip0 =  lon_cmp(ip0(i,j))
      lam_ip1 =  lon_cmp(ip1(i,j))

      if( lam_im2.gt.lam_im1 ) lam_im2 = lam_im2 - 2*pi
      if( lam_im1.gt.lam_ip0 ) lam_ip0 = lam_ip0 + 2*pi
      if( lam_im1.gt.lam_ip1 ) lam_ip1 = lam_ip1 + 2*pi
      if( lam_ip0.gt.lam_ip1 ) lam_ip1 = lam_ip1 + 2*pi

c Compute and Adjust Latitude Weights   
c ------------------------------------------------------------------
      phi_jm2 =  lat_cmp(jm2(i,j))
      phi_jm1 =  lat_cmp(jm1(i,j))
      phi_jp0 =  lat_cmp(jp0(i,j))
      phi_jp1 =  lat_cmp(jp1(i,j))

      if( jm2(i,j).eq.0 ) phi_jm2 = phi_jm1 - dphi(1)
      if( jm1(i,j).eq.lattice%jmglobal ) then
                          phi_jp0 = phi_jm1 + dphi(lattice%jmglobal-1)
                          phi_jp1 = phi_jp0 + dphi(lattice%jmglobal-2)
      endif
      if( jp1(i,j).eq.lattice%jmglobal+1 ) phi_jp1 = phi_jp0 + dphi(lattice%jmglobal-1)

c Compute ij index
c ----------------
         msgn(i,j) = 2
      if( jp0(i,j).eq.lattice%jmglobal ) then
          jp1(i,j) =  lattice%jmglobal-1
       ip1jp1(i,j) = 1+mod(ip1(i,j)+lattice%imglobal/2-1,lattice%imglobal) + ( jp1(i,j)-1 )*lattice%imglobal
       ip0jp1(i,j) = 1+mod(ip0(i,j)+lattice%imglobal/2-1,lattice%imglobal) + ( jp1(i,j)-1 )*lattice%imglobal
       im1jp1(i,j) = 1+mod(im1(i,j)+lattice%imglobal/2-1,lattice%imglobal) + ( jp1(i,j)-1 )*lattice%imglobal
       im2jp1(i,j) = 1+mod(im2(i,j)+lattice%imglobal/2-1,lattice%imglobal) + ( jp1(i,j)-1 )*lattice%imglobal
         msgn(i,j) = 1
      else
       ip1jp1(i,j) = ip1(i,j) + ( jp1(i,j)-1 )*lattice%imglobal
       ip0jp1(i,j) = ip0(i,j) + ( jp1(i,j)-1 )*lattice%imglobal
       im1jp1(i,j) = im1(i,j) + ( jp1(i,j)-1 )*lattice%imglobal
       im2jp1(i,j) = im2(i,j) + ( jp1(i,j)-1 )*lattice%imglobal
      endif

      ip1jp0(i,j) = ip1(i,j) + ( jp0(i,j)-1 )*lattice%imglobal
      ip0jp0(i,j) = ip0(i,j) + ( jp0(i,j)-1 )*lattice%imglobal
      im1jp0(i,j) = im1(i,j) + ( jp0(i,j)-1 )*lattice%imglobal
      im2jp0(i,j) = im2(i,j) + ( jp0(i,j)-1 )*lattice%imglobal

      ip1jm1(i,j) = ip1(i,j) + ( jm1(i,j)-1 )*lattice%imglobal
      ip0jm1(i,j) = ip0(i,j) + ( jm1(i,j)-1 )*lattice%imglobal
      im1jm1(i,j) = im1(i,j) + ( jm1(i,j)-1 )*lattice%imglobal
      im2jm1(i,j) = im2(i,j) + ( jm1(i,j)-1 )*lattice%imglobal

         nsgn(i,j) = 2
      if( jm1(i,j).eq.1 ) then
          jm2(i,j) =  2
       ip1jm2(i,j) = 1+mod(ip1(i,j)+lattice%imglobal/2-1,lattice%imglobal) + ( jm2(i,j)-1 )*lattice%imglobal
       ip0jm2(i,j) = 1+mod(ip0(i,j)+lattice%imglobal/2-1,lattice%imglobal) + ( jm2(i,j)-1 )*lattice%imglobal
       im1jm2(i,j) = 1+mod(im1(i,j)+lattice%imglobal/2-1,lattice%imglobal) + ( jm2(i,j)-1 )*lattice%imglobal
       im2jm2(i,j) = 1+mod(im2(i,j)+lattice%imglobal/2-1,lattice%imglobal) + ( jm2(i,j)-1 )*lattice%imglobal
         nsgn(i,j) = 1
      else
       ip1jm2(i,j) = ip1(i,j) + ( jm2(i,j)-1 )*lattice%imglobal
       ip0jm2(i,j) = ip0(i,j) + ( jm2(i,j)-1 )*lattice%imglobal
       im1jm2(i,j) = im1(i,j) + ( jm2(i,j)-1 )*lattice%imglobal
       im2jm2(i,j) = im2(i,j) + ( jm2(i,j)-1 )*lattice%imglobal
      endif

c Determine PE location of Source Data
c ------------------------------------
      pej_jp1 = lattice%pejglobal( jp1(i,j) )
      pej_jp0 = lattice%pejglobal( jp0(i,j) )
      pej_jm1 = lattice%pejglobal( jm1(i,j) )
      pej_jm2 = lattice%pejglobal( jm2(i,j) )

      pei_ip1 = lattice%peiglobal( ip1(i,j) )
      pei_ip0 = lattice%peiglobal( ip0(i,j) )
      pei_im1 = lattice%peiglobal( im1(i,j) )
      pei_im2 = lattice%peiglobal( im2(i,j) )

      pe_ip1jp0(i,j) = pei_ip1 + pej_jp0*lattice%nx
      pe_ip0jp0(i,j) = pei_ip0 + pej_jp0*lattice%nx
      pe_im1jp0(i,j) = pei_im1 + pej_jp0*lattice%nx
      pe_im2jp0(i,j) = pei_im2 + pej_jp0*lattice%nx

      pe_ip1jm1(i,j) = pei_ip1 + pej_jm1*lattice%nx
      pe_ip0jm1(i,j) = pei_ip0 + pej_jm1*lattice%nx
      pe_im1jm1(i,j) = pei_im1 + pej_jm1*lattice%nx
      pe_im2jm1(i,j) = pei_im2 + pej_jm1*lattice%nx

            n = 1 + (ip1jp1(i,j)-1)/lattice%imglobal
            m = ip1jp1(i,j) - (n-1)*lattice%imglobal
      pei_ip1 = lattice%peiglobal( m )
            n = 1 + (ip0jp1(i,j)-1)/lattice%imglobal
            m = ip0jp1(i,j) - (n-1)*lattice%imglobal
      pei_ip0 = lattice%peiglobal( m )
            n = 1 + (im1jp1(i,j)-1)/lattice%imglobal
            m = im1jp1(i,j) - (n-1)*lattice%imglobal
      pei_im1 = lattice%peiglobal( m )
            n = 1 + (im2jp1(i,j)-1)/lattice%imglobal
            m = im2jp1(i,j) - (n-1)*lattice%imglobal
      pei_im2 = lattice%peiglobal( m )

      pe_ip1jp1(i,j) = pei_ip1 + pej_jp1*lattice%nx
      pe_ip0jp1(i,j) = pei_ip0 + pej_jp1*lattice%nx
      pe_im1jp1(i,j) = pei_im1 + pej_jp1*lattice%nx
      pe_im2jp1(i,j) = pei_im2 + pej_jp1*lattice%nx

            n = 1 + (ip1jm2(i,j)-1)/lattice%imglobal
            m = ip1jm2(i,j) - (n-1)*lattice%imglobal
      pei_ip1 = lattice%peiglobal( m )
            n = 1 + (ip0jm2(i,j)-1)/lattice%imglobal
            m = ip0jm2(i,j) - (n-1)*lattice%imglobal
      pei_ip0 = lattice%peiglobal( m )
            n = 1 + (im1jm2(i,j)-1)/lattice%imglobal
            m = im1jm2(i,j) - (n-1)*lattice%imglobal
      pei_im1 = lattice%peiglobal( m )
            n = 1 + (im2jm2(i,j)-1)/lattice%imglobal
            m = im2jm2(i,j) - (n-1)*lattice%imglobal
      pei_im2 = lattice%peiglobal( m )

      pe_ip1jm2(i,j) = pei_ip1 + pej_jm2*lattice%nx
      pe_ip0jm2(i,j) = pei_ip0 + pej_jm2*lattice%nx
      pe_im1jm2(i,j) = pei_im1 + pej_jm2*lattice%nx
      pe_im2jm2(i,j) = pei_im2 + pej_jm2*lattice%nx

c Bi-Linear Weights
c -----------------
              d      = (lam_ip0-lam_im1)*(phi_jp0-phi_jm1)
      wl_im1jm1(i,j) = (lam_ip0-lam    )*(phi_jp0-phi    )/d
      wl_ip0jm1(i,j) = (lam    -lam_im1)*(phi_jp0-phi    )/d
      wl_im1jp0(i,j) = (lam_ip0-lam    )*(phi    -phi_jm1)/d
      wl_ip0jp0(i,j) = (lam    -lam_im1)*(phi    -phi_jm1)/d

c Bi-Cubic Weights
c ----------------
      ap1 = ( (lam    -lam_ip0)*(lam    -lam_im1)*(lam    -lam_im2) )
     .    / ( (lam_ip1-lam_ip0)*(lam_ip1-lam_im1)*(lam_ip1-lam_im2) )
      ap0 = ( (lam_ip1-lam    )*(lam    -lam_im1)*(lam    -lam_im2) )
     .    / ( (lam_ip1-lam_ip0)*(lam_ip0-lam_im1)*(lam_ip0-lam_im2) )
      am1 = ( (lam_ip1-lam    )*(lam_ip0-lam    )*(lam    -lam_im2) )
     .    / ( (lam_ip1-lam_im1)*(lam_ip0-lam_im1)*(lam_im1-lam_im2) )
      am2 = ( (lam_ip1-lam    )*(lam_ip0-lam    )*(lam_im1-lam    ) )
     .    / ( (lam_ip1-lam_im2)*(lam_ip0-lam_im2)*(lam_im1-lam_im2) )

      bp1 = ( (phi    -phi_jp0)*(phi    -phi_jm1)*(phi    -phi_jm2) )
     .    / ( (phi_jp1-phi_jp0)*(phi_jp1-phi_jm1)*(phi_jp1-phi_jm2) )
      bp0 = ( (phi_jp1-phi    )*(phi    -phi_jm1)*(phi    -phi_jm2) )
     .    / ( (phi_jp1-phi_jp0)*(phi_jp0-phi_jm1)*(phi_jp0-phi_jm2) )
      bm1 = ( (phi_jp1-phi    )*(phi_jp0-phi    )*(phi    -phi_jm2) )
     .    / ( (phi_jp1-phi_jm1)*(phi_jp0-phi_jm1)*(phi_jm1-phi_jm2) )
      bm2 = ( (phi_jp1-phi    )*(phi_jp0-phi    )*(phi_jm1-phi    ) )
     .    / ( (phi_jp1-phi_jm2)*(phi_jp0-phi_jm2)*(phi_jm1-phi_jm2) )

      wc_ip1jp1(i,j) = bp1*ap1
      wc_ip0jp1(i,j) = bp1*ap0
      wc_im1jp1(i,j) = bp1*am1
      wc_im2jp1(i,j) = bp1*am2

      wc_ip1jp0(i,j) = bp0*ap1
      wc_ip0jp0(i,j) = bp0*ap0
      wc_im1jp0(i,j) = bp0*am1
      wc_im2jp0(i,j) = bp0*am2

      wc_ip1jm1(i,j) = bm1*ap1
      wc_ip0jm1(i,j) = bm1*ap0
      wc_im1jm1(i,j) = bm1*am1
      wc_im2jm1(i,j) = bm1*am2

      wc_ip1jm2(i,j) = bm2*ap1
      wc_ip0jm2(i,j) = bm2*ap0
      wc_im1jm2(i,j) = bm2*am1
      wc_im2jm2(i,j) = bm2*am2

      enddo
      enddo

c Create Data Index Array for Total and Non-Redundant Data
c Note:  recvlen(peid):   is the total number of unique grid-points needed from a given peid
c -----  recvloc(m,peid): is a list of the global locations of those points
c -------------------------------------------------------------------------
      allocate ( bucket(16*lattice%imglobal*lattice%jmglobal) )
      do peid = 0,lattice%nx*lattice%ny-1
          num = 0
         do j=1,jm
         do i=1,im
            if ( pe_im2jm2(i,j).eq.peid ) then
                                          num  = num + 1
                                          loc(num,peid) = im2jm2(i,j)
                                          pos(num,peid) = i+(j-1)*im
            endif
            if ( pe_im1jm2(i,j).eq.peid ) then
                                          num  = num + 1
                                          loc(num,peid) = im1jm2(i,j)
                                          pos(num,peid) = i+(j-1)*im
            endif
            if ( pe_ip0jm2(i,j).eq.peid ) then
                                          num  = num + 1
                                          loc(num,peid) = ip0jm2(i,j)
                                          pos(num,peid) = i+(j-1)*im
            endif
            if ( pe_ip1jm2(i,j).eq.peid ) then
                                          num  = num + 1
                                          loc(num,peid) = ip1jm2(i,j)
                                          pos(num,peid) = i+(j-1)*im
            endif

            if ( pe_im2jm1(i,j).eq.peid ) then
                                          num  = num + 1
                                          loc(num,peid) = im2jm1(i,j)
                                          pos(num,peid) = i+(j-1)*im
            endif
            if ( pe_im1jm1(i,j).eq.peid ) then
                                          num  = num + 1
                                          loc(num,peid) = im1jm1(i,j)
                                          pos(num,peid) = i+(j-1)*im
            endif
            if ( pe_ip0jm1(i,j).eq.peid ) then
                                          num  = num + 1
                                          loc(num,peid) = ip0jm1(i,j)
                                          pos(num,peid) = i+(j-1)*im
            endif
            if ( pe_ip1jm1(i,j).eq.peid ) then
                                          num  = num + 1
                                          loc(num,peid) = ip1jm1(i,j)
                                          pos(num,peid) = i+(j-1)*im
            endif

            if ( pe_im2jp0(i,j).eq.peid ) then
                                          num  = num + 1
                                          loc(num,peid) = im2jp0(i,j)
                                          pos(num,peid) = i+(j-1)*im
            endif
            if ( pe_im1jp0(i,j).eq.peid ) then
                                          num  = num + 1
                                          loc(num,peid) = im1jp0(i,j)
                                          pos(num,peid) = i+(j-1)*im
            endif
            if ( pe_ip0jp0(i,j).eq.peid ) then
                                          num  = num + 1
                                          loc(num,peid) = ip0jp0(i,j)
                                          pos(num,peid) = i+(j-1)*im
            endif
            if ( pe_ip1jp0(i,j).eq.peid ) then
                                          num  = num + 1
                                          loc(num,peid) = ip1jp0(i,j)
                                          pos(num,peid) = i+(j-1)*im
            endif

            if ( pe_im2jp1(i,j).eq.peid ) then
                                          num  = num + 1
                                          loc(num,peid) = im2jp1(i,j)
                                          pos(num,peid) = i+(j-1)*im
            endif
            if ( pe_im1jp1(i,j).eq.peid ) then
                                          num  = num + 1
                                          loc(num,peid) = im1jp1(i,j)
                                          pos(num,peid) = i+(j-1)*im
            endif
            if ( pe_ip0jp1(i,j).eq.peid ) then
                                          num  = num + 1
                                          loc(num,peid) = ip0jp1(i,j)
                                          pos(num,peid) = i+(j-1)*im
            endif
            if ( pe_ip1jp1(i,j).eq.peid ) then
                                          num  = num + 1
                                          loc(num,peid) = ip1jp1(i,j)
                                          pos(num,peid) = i+(j-1)*im
            endif
         enddo
         enddo
              len(peid) = num

c Strip redundancies from buffer
c ------------------------------
                        m = 0
         if( num.ne.0 ) then
                bucket(:) = 0
                do i=1,num
                bucket( loc(i,peid) ) = 1
                enddo
                do i=1,16*lattice%imglobal*lattice%jmglobal
                   if( bucket(i).eq.1 ) then
                                     m = m+1
                       recvloc(m,peid) = i
                            bucket (i) = m
                   endif
                enddo
                do i=1,num
                nn(i,peid) = bucket( loc(i,peid) )
                enddo
         endif

         recvlen(peid) = m
      enddo
      deallocate ( bucket )

c Store buffer locations
c ----------------------
      do peid = 0,lattice%nx*lattice%ny-1
         do m = 1,len(peid)
            n = nn (m,peid)
            j = 1 + (pos(m,peid)-1)/im
            i = pos(m,peid) - (j-1)*im
            if( recvloc(n,peid).eq.ip1jp1(i,j) ) then
                  n_ip1jp1(i,j) = n
                  p_ip1jp1(i,j) = peid
            endif
            if( recvloc(n,peid).eq.ip0jp1(i,j) ) then
                  n_ip0jp1(i,j) = n
                  p_ip0jp1(i,j) = peid
            endif
            if( recvloc(n,peid).eq.im1jp1(i,j) ) then
                  n_im1jp1(i,j) = n
                  p_im1jp1(i,j) = peid
            endif
            if( recvloc(n,peid).eq.im2jp1(i,j) ) then
                  n_im2jp1(i,j) = n
                  p_im2jp1(i,j) = peid
            endif

            if( recvloc(n,peid).eq.ip1jp0(i,j) ) then
                  n_ip1jp0(i,j) = n
                  p_ip1jp0(i,j) = peid
            endif
            if( recvloc(n,peid).eq.ip0jp0(i,j) ) then
                  n_ip0jp0(i,j) = n
                  p_ip0jp0(i,j) = peid
            endif
            if( recvloc(n,peid).eq.im1jp0(i,j) ) then
                  n_im1jp0(i,j) = n
                  p_im1jp0(i,j) = peid
            endif
            if( recvloc(n,peid).eq.im2jp0(i,j) ) then
                  n_im2jp0(i,j) = n
                  p_im2jp0(i,j) = peid
            endif

            if( recvloc(n,peid).eq.ip1jm1(i,j) ) then
                  n_ip1jm1(i,j) = n
                  p_ip1jm1(i,j) = peid
            endif
            if( recvloc(n,peid).eq.ip0jm1(i,j) ) then
                  n_ip0jm1(i,j) = n
                  p_ip0jm1(i,j) = peid
            endif
            if( recvloc(n,peid).eq.im1jm1(i,j) ) then
                  n_im1jm1(i,j) = n
                  p_im1jm1(i,j) = peid
            endif
            if( recvloc(n,peid).eq.im2jm1(i,j) ) then
                  n_im2jm1(i,j) = n
                  p_im2jm1(i,j) = peid
            endif

            if( recvloc(n,peid).eq.ip1jm2(i,j) ) then
                  n_ip1jm2(i,j) = n
                  p_ip1jm2(i,j) = peid
            endif
            if( recvloc(n,peid).eq.ip0jm2(i,j) ) then
                  n_ip0jm2(i,j) = n
                  p_ip0jm2(i,j) = peid
            endif
            if( recvloc(n,peid).eq.im1jm2(i,j) ) then
                  n_im1jm2(i,j) = n
                  p_im1jm2(i,j) = peid
            endif
            if( recvloc(n,peid).eq.im2jm2(i,j) ) then
                  n_im2jm2(i,j) = n
                  p_im2jm2(i,j) = peid
            endif
         enddo
      enddo


c Send and Recv the data length needed for interpolation
c ------------------------------------------------------
      do peid = 0,lattice%nx*lattice%ny-1
      if( lattice%myid.ne.peid ) then

          call mpi_isend ( recvlen(peid),1,mpi_integer,peid,peid,lattice%comm,sendquest(peid),ierror )

      else
          sendquest(peid) = mpi_request_null
      endif
      enddo

      do peid = 0,lattice%nx*lattice%ny-1
      if( lattice%myid.ne.peid ) then

          call mpi_irecv ( sendlen(peid),1,mpi_integer,peid,lattice%myid,lattice%comm,recvquest(peid),ierror )

      else
            sendlen(peid) = recvlen(peid)
          recvquest(peid) = mpi_request_null
      endif
      enddo

      call mpi_waitall ( lattice%nx*lattice%ny,sendquest(0:lattice%nx*lattice%ny-1),stats(1,0),ierror )
      call mpi_waitall ( lattice%nx*lattice%ny,recvquest(0:lattice%nx*lattice%ny-1),statr(1,0),ierror )

c Send and Recv the data locations needed for interpolation
c ---------------------------------------------------------
      do peid = 0,lattice%nx*lattice%ny-1
          sendquest(peid) = mpi_request_null
      if( lattice%myid.ne.peid ) then
          if( recvlen(peid).ne.0 ) then

          call mpi_isend ( recvloc(1,peid),recvlen(peid),mpi_integer,peid,peid,lattice%comm,sendquest(peid),ierror )

          endif
      endif
      enddo

      do peid = 0,lattice%nx*lattice%ny-1
         recvquest(peid) = mpi_request_null
         if( sendlen(peid).ne.0 ) then
             if( lattice%myid.ne.peid ) then

             call mpi_irecv ( sendloc(1,peid),sendlen(peid),mpi_integer,peid,lattice%myid,lattice%comm,recvquest(peid),ierror )

             else
             do n=1,sendlen(peid)
             sendloc(n,peid) = recvloc(n,peid)
             enddo
             endif
         endif
      enddo

      call mpi_waitall ( lattice%nx*lattice%ny,sendquest(0:lattice%nx*lattice%ny-1),stats(1,0),ierror )
      call mpi_waitall ( lattice%nx*lattice%ny,recvquest(0:lattice%nx*lattice%ny-1),statr(1,0),ierror )

c Print Communication Matrix
c --------------------------
      do n=0,lattice%nx*lattice%ny-1
      if( n.eq.lattice%myid ) then
          write(6,1001) (sendlen(m),m=0,lattice%nx*lattice%ny-1)
      endif
      call my_barrier (lattice%comm)
      enddo
      if( lattice%myid.eq.lattice%nx*lattice%ny-1 ) print *
 1001 format(128(1x,i4))

      deallocate ( nn,pos,loc,len )
      deallocate (       ip1 ,       ip0 ,       im1 ,       im2 )
      deallocate (       jp1 ,       jp0 ,       jm1 ,       jm2 )
      deallocate (    im2jm2 ,    im1jm2 ,    ip0jm2 ,    ip1jm2 )
      deallocate (    im2jm1 ,    im1jm1 ,    ip0jm1 ,    ip1jm1 )
      deallocate (    im2jp0 ,    im1jp0 ,    ip0jp0 ,    ip1jp0 )
      deallocate (    im2jp1 ,    im1jp1 ,    ip0jp1 ,    ip1jp1 )
      deallocate ( pe_im2jm2 , pe_im1jm2 , pe_ip0jm2 , pe_ip1jm2 )
      deallocate ( pe_im2jm1 , pe_im1jm1 , pe_ip0jm1 , pe_ip1jm1 )
      deallocate ( pe_im2jp0 , pe_im1jp0 , pe_ip0jp0 , pe_ip1jp0 )
      deallocate ( pe_im2jp1 , pe_im1jp1 , pe_ip0jp1 , pe_ip1jp1 )

      endif   ! End compute_weights check


C***********************************************************************
C*                     Begin Rotate Interpolation                      *
C***********************************************************************

c Send and Recv the data needed for interpolation
c -----------------------------------------------
      allocate ( sendbuf(im*jm*lm                    ,0:lattice%nx*lattice%ny-1) )
      allocate ( recvbuf(lattice%imax*lattice%jmax*lm,0:lattice%nx*lattice%ny-1) )

      do peid = 0,lattice%nx*lattice%ny-1
          sendquest(peid) = mpi_request_null
          if( sendlen(peid).ne.0 ) then
              if( lattice%myid.ne.peid ) then

                 m = 0
              do n = 1,sendlen(peid)
                 j = 1 + (sendloc(n,peid)-1)/lattice%imglobal
                 i = sendloc(n,peid) - (j-1)*lattice%imglobal
              do L = 1,lm
                 m = m+1
              sendbuf(m,peid) = q_cmp( lattice%ilocal(i),lattice%jlocal(j),L )
              enddo
              enddo
              call mpi_isend ( sendbuf(1,peid),m,mpi_double_precision,peid,peid,lattice%comm,
     .                         sendquest(peid),ierror )

              else
                 m = 0
              do n = 1,sendlen(peid)
                 j = 1 + (sendloc(n,peid)-1)/lattice%imglobal
                 i = sendloc(n,peid) - (j-1)*lattice%imglobal
              do L = 1,lm
                 m = m+1
              recvbuf(m,peid) = q_cmp( lattice%ilocal(i),lattice%jlocal(j),L )
              enddo
              enddo
          endif
      endif
      enddo

      do peid = 0,lattice%nx*lattice%ny-1
          recvquest(peid) = mpi_request_null
      if( lattice%myid.ne.peid ) then
          if( recvlen(peid).ne.0 ) then

          call mpi_irecv ( recvbuf(1,peid),recvlen(peid)*lm,mpi_double_precision,peid,lattice%myid,lattice%comm,
     .                     recvquest(peid),ierror )

          endif
      endif
      enddo

      call mpi_waitall ( lattice%nx*lattice%ny,sendquest(0:lattice%nx*lattice%ny-1),stats(1,0),ierror )
      call mpi_waitall ( lattice%nx*lattice%ny,recvquest(0:lattice%nx*lattice%ny-1),statr(1,0),ierror )

c Interpolate Geophysical Quantities to Computational Grid Using Bi-Linear
c ------------------------------------------------------------------------
      if( abs(norder).eq.1 ) then

      if( check  ) then
      do j=1,jm
      do i=1,im
      do L=1,lm
         if( recvbuf( L+(n_im1jm1(i,j)-1)*lm,p_im1jm1(i,j) ).ne.undef  .and.
     .       recvbuf( L+(n_ip0jm1(i,j)-1)*lm,p_ip0jm1(i,j) ).ne.undef  .and.
     .       recvbuf( L+(n_im1jp0(i,j)-1)*lm,p_im1jp0(i,j) ).ne.undef  .and.
     .       recvbuf( L+(n_ip0jp0(i,j)-1)*lm,p_ip0jp0(i,j) ).ne.undef ) then
         q_geo(i,j,L) = wl_im1jm1(i,j) * recvbuf( L+(n_im1jm1(i,j)-1)*lm,p_im1jm1(i,j) )
     .                + wl_ip0jm1(i,j) * recvbuf( L+(n_ip0jm1(i,j)-1)*lm,p_ip0jm1(i,j) )
     .                + wl_im1jp0(i,j) * recvbuf( L+(n_im1jp0(i,j)-1)*lm,p_im1jp0(i,j) )
     .                + wl_ip0jp0(i,j) * recvbuf( L+(n_ip0jp0(i,j)-1)*lm,p_ip0jp0(i,j) )
         else
         q_geo(i,j,L) = undef
         endif
      enddo
      enddo
      enddo
      endif      ! End UNDEF Check

      if( .not.check  ) then
      do j=1,jm
      do i=1,im
      do L=1,lm
         q_geo(i,j,L) = wl_im1jm1(i,j) * recvbuf( L+(n_im1jm1(i,j)-1)*lm,p_im1jm1(i,j) )
     .                + wl_ip0jm1(i,j) * recvbuf( L+(n_ip0jm1(i,j)-1)*lm,p_ip0jm1(i,j) )
     .                + wl_im1jp0(i,j) * recvbuf( L+(n_im1jp0(i,j)-1)*lm,p_im1jp0(i,j) )
     .                + wl_ip0jp0(i,j) * recvbuf( L+(n_ip0jp0(i,j)-1)*lm,p_ip0jp0(i,j) )
      enddo
      enddo
      enddo
      endif      ! End UNDEF Check

      endif      ! End Check for Bi-Linear Interpolation

c Interpolate Geophysical Quantities to Computational Grid Using Bi-Cubic
c -----------------------------------------------------------------------
      if( abs(norder).eq.3 ) then

      if( check  ) then
      do j=1,jm
      do i=1,im
         sgnjp1 = sgn**msgn(i,j)
         sgnjm2 = sgn**nsgn(i,j)
      do L=1,lm

      if( recvbuf( L+(n_ip1jp1(i,j)-1)*lm,p_ip1jp1(i,j) ).ne.undef  .and.
     .    recvbuf( L+(n_ip0jp1(i,j)-1)*lm,p_ip0jp1(i,j) ).ne.undef  .and.
     .    recvbuf( L+(n_im1jp1(i,j)-1)*lm,p_im1jp1(i,j) ).ne.undef  .and.
     .    recvbuf( L+(n_im2jp1(i,j)-1)*lm,p_im2jp1(i,j) ).ne.undef  .and.

     .    recvbuf( L+(n_ip1jp0(i,j)-1)*lm,p_ip1jp0(i,j) ).ne.undef  .and.
     .    recvbuf( L+(n_ip0jp0(i,j)-1)*lm,p_ip0jp0(i,j) ).ne.undef  .and.
     .    recvbuf( L+(n_im1jp0(i,j)-1)*lm,p_im1jp0(i,j) ).ne.undef  .and.
     .    recvbuf( L+(n_im2jp0(i,j)-1)*lm,p_im2jp0(i,j) ).ne.undef  .and.

     .    recvbuf( L+(n_ip1jm1(i,j)-1)*lm,p_ip1jm1(i,j) ).ne.undef  .and.
     .    recvbuf( L+(n_ip0jm1(i,j)-1)*lm,p_ip0jm1(i,j) ).ne.undef  .and.
     .    recvbuf( L+(n_im1jm1(i,j)-1)*lm,p_im1jm1(i,j) ).ne.undef  .and.
     .    recvbuf( L+(n_im2jm1(i,j)-1)*lm,p_im2jm1(i,j) ).ne.undef  .and.

     .    recvbuf( L+(n_ip1jm2(i,j)-1)*lm,p_ip1jm2(i,j) ).ne.undef  .and.
     .    recvbuf( L+(n_ip0jm2(i,j)-1)*lm,p_ip0jm2(i,j) ).ne.undef  .and.
     .    recvbuf( L+(n_im1jm2(i,j)-1)*lm,p_im1jm2(i,j) ).ne.undef  .and.
     .    recvbuf( L+(n_im2jm2(i,j)-1)*lm,p_im2jm2(i,j) ).ne.undef ) then

      q_geo(i,j,L) = wc_ip1jp1(i,j) * recvbuf( L+(n_ip1jp1(i,j)-1)*lm,p_ip1jp1(i,j) )*sgnjp1
     .             + wc_ip0jp1(i,j) * recvbuf( L+(n_ip0jp1(i,j)-1)*lm,p_ip0jp1(i,j) )*sgnjp1
     .             + wc_im1jp1(i,j) * recvbuf( L+(n_im1jp1(i,j)-1)*lm,p_im1jp1(i,j) )*sgnjp1
     .             + wc_im2jp1(i,j) * recvbuf( L+(n_im2jp1(i,j)-1)*lm,p_im2jp1(i,j) )*sgnjp1

     .             + wc_ip1jp0(i,j) * recvbuf( L+(n_ip1jp0(i,j)-1)*lm,p_ip1jp0(i,j) )
     .             + wc_ip0jp0(i,j) * recvbuf( L+(n_ip0jp0(i,j)-1)*lm,p_ip0jp0(i,j) )
     .             + wc_im1jp0(i,j) * recvbuf( L+(n_im1jp0(i,j)-1)*lm,p_im1jp0(i,j) )
     .             + wc_im2jp0(i,j) * recvbuf( L+(n_im2jp0(i,j)-1)*lm,p_im2jp0(i,j) )

     .             + wc_ip1jm1(i,j) * recvbuf( L+(n_ip1jm1(i,j)-1)*lm,p_ip1jm1(i,j) )
     .             + wc_ip0jm1(i,j) * recvbuf( L+(n_ip0jm1(i,j)-1)*lm,p_ip0jm1(i,j) )
     .             + wc_im1jm1(i,j) * recvbuf( L+(n_im1jm1(i,j)-1)*lm,p_im1jm1(i,j) )
     .             + wc_im2jm1(i,j) * recvbuf( L+(n_im2jm1(i,j)-1)*lm,p_im2jm1(i,j) )

     .             + wc_ip1jm2(i,j) * recvbuf( L+(n_ip1jm2(i,j)-1)*lm,p_ip1jm2(i,j) )*sgnjm2
     .             + wc_ip0jm2(i,j) * recvbuf( L+(n_ip0jm2(i,j)-1)*lm,p_ip0jm2(i,j) )*sgnjm2
     .             + wc_im1jm2(i,j) * recvbuf( L+(n_im1jm2(i,j)-1)*lm,p_im1jm2(i,j) )*sgnjm2
     .             + wc_im2jm2(i,j) * recvbuf( L+(n_im2jm2(i,j)-1)*lm,p_im2jm2(i,j) )*sgnjm2

      elseif( recvbuf( L+(n_im1jm1(i,j)-1)*lm,p_im1jm1(i,j) ).ne.undef  .and.
     .        recvbuf( L+(n_ip0jm1(i,j)-1)*lm,p_ip0jm1(i,j) ).ne.undef  .and.
     .        recvbuf( L+(n_im1jp0(i,j)-1)*lm,p_im1jp0(i,j) ).ne.undef  .and.
     .        recvbuf( L+(n_ip0jp0(i,j)-1)*lm,p_ip0jp0(i,j) ).ne.undef ) then

      q_geo(i,j,L) = wl_im1jm1(i,j) * recvbuf( L+(n_im1jm1(i,j)-1)*lm,p_im1jm1(i,j) )
     .             + wl_ip0jm1(i,j) * recvbuf( L+(n_ip0jm1(i,j)-1)*lm,p_ip0jm1(i,j) )
     .             + wl_im1jp0(i,j) * recvbuf( L+(n_im1jp0(i,j)-1)*lm,p_im1jp0(i,j) )
     .             + wl_ip0jp0(i,j) * recvbuf( L+(n_ip0jp0(i,j)-1)*lm,p_ip0jp0(i,j) )

      else
      q_geo(i,j,L) = undef
      endif

c Check for Positive Definite
c ---------------------------
      if( norder.lt.0 ) then
          if( q_geo(i,j,L).ne.undef .and. q_geo(i,j,L).lt.0.0 ) q_geo(i,j,L) = 0.0
      endif

      enddo      ! End L Loop
      enddo      ! End i Loop
      enddo      ! End j Loop

      endif      ! End UNDEF Check

      if( .not.check  ) then
      do j=1,jm
      do i=1,im
         sgnjp1 = sgn**msgn(i,j)
         sgnjm2 = sgn**nsgn(i,j)
      do L=1,lm

      q_geo(i,j,L) = wc_ip1jp1(i,j) * recvbuf( L+(n_ip1jp1(i,j)-1)*lm,p_ip1jp1(i,j) )*sgnjp1
     .             + wc_ip0jp1(i,j) * recvbuf( L+(n_ip0jp1(i,j)-1)*lm,p_ip0jp1(i,j) )*sgnjp1
     .             + wc_im1jp1(i,j) * recvbuf( L+(n_im1jp1(i,j)-1)*lm,p_im1jp1(i,j) )*sgnjp1
     .             + wc_im2jp1(i,j) * recvbuf( L+(n_im2jp1(i,j)-1)*lm,p_im2jp1(i,j) )*sgnjp1

     .             + wc_ip1jp0(i,j) * recvbuf( L+(n_ip1jp0(i,j)-1)*lm,p_ip1jp0(i,j) )
     .             + wc_ip0jp0(i,j) * recvbuf( L+(n_ip0jp0(i,j)-1)*lm,p_ip0jp0(i,j) )
     .             + wc_im1jp0(i,j) * recvbuf( L+(n_im1jp0(i,j)-1)*lm,p_im1jp0(i,j) )
     .             + wc_im2jp0(i,j) * recvbuf( L+(n_im2jp0(i,j)-1)*lm,p_im2jp0(i,j) )

     .             + wc_ip1jm1(i,j) * recvbuf( L+(n_ip1jm1(i,j)-1)*lm,p_ip1jm1(i,j) )
     .             + wc_ip0jm1(i,j) * recvbuf( L+(n_ip0jm1(i,j)-1)*lm,p_ip0jm1(i,j) )
     .             + wc_im1jm1(i,j) * recvbuf( L+(n_im1jm1(i,j)-1)*lm,p_im1jm1(i,j) )
     .             + wc_im2jm1(i,j) * recvbuf( L+(n_im2jm1(i,j)-1)*lm,p_im2jm1(i,j) )

     .             + wc_ip1jm2(i,j) * recvbuf( L+(n_ip1jm2(i,j)-1)*lm,p_ip1jm2(i,j) )*sgnjm2
     .             + wc_ip0jm2(i,j) * recvbuf( L+(n_ip0jm2(i,j)-1)*lm,p_ip0jm2(i,j) )*sgnjm2
     .             + wc_im1jm2(i,j) * recvbuf( L+(n_im1jm2(i,j)-1)*lm,p_im1jm2(i,j) )*sgnjm2
     .             + wc_im2jm2(i,j) * recvbuf( L+(n_im2jm2(i,j)-1)*lm,p_im2jm2(i,j) )*sgnjm2

c Check for Positive Definite
c ---------------------------
      if( norder.lt.0 ) then
          if( q_geo(i,j,L).ne.undef .and. q_geo(i,j,L).lt.0.0 )  q_geo(i,j,L) = 0.0
      endif

      enddo      ! End L Loop
      enddo      ! End i Loop
      enddo      ! End j Loop

      endif      ! End UNDEF Check
      endif      ! End Check for Bi-Cubic Interpolation

      deallocate ( sendbuf )
      deallocate ( recvbuf )

      return
      end

      subroutine vector ( u,v,im,jm,lm,msgn,
     .                    dlam,dphi,rotation,tilt,precession,lattice )
C***********************************************************************
C
C  PURPOSE:
C  ========
C    Transform wind vectors after rotation to new coordinate system.
C
C  INPUT:
C  ======
C    u .......... U-Wind u(im,jm,lm) after rotation from old system
C    v .......... V-Wind v(im,jm,lm) after rotation from old system
C    im ......... Longitudinal dimension of wind
C    jm ......... Latitudinal  dimension of wind
C    lm ......... Vertical     dimension of wind
C    msgn ....... Flag for a forward  transformation  ( msgn =  1 )
C                       or a backward transformation  ( msgn = -1 )
C    dlam ....... Computational Grid Delta Lambda
C    dphi ....... Computational Grid Delta Phi
C    rotation ... Rotation parameter lam_np (Degrees)
C    tilt ....... Rotation parameter phi_np (Degrees)
C    precession . Rotation parameter lam_0  (Degrees)
C
C  OUTPUT:
C  =======
C    u .......... U-Wind u(im,jm,lm) after rotation in new system
C    v .......... V-Wind v(im,jm,lm) after rotation in new system
C
C
C***********************************************************************
C*                  GODDARD LABORATORY FOR ATMOSPHERES                 *
C***********************************************************************

      use g3_dynamics_lattice_module
      implicit none
      type ( dynamics_lattice_type ) lattice

c Input Variables
c ---------------
      integer  im,jm,lm,msgn
      real(kind=8)     u(im,jm,lm)
      real(kind=8)     v(im,jm,lm)
      real(kind=8)     dlam(lattice%imglobal)
      real(kind=8)     dphi(lattice%jmglobal)
      real(kind=8)     rotation, tilt, precession

c Local Variables
c ---------------
      integer i,j,l,m,n

c Tranformation Functions
c -----------------------
      real(kind=8), allocatable, save ::    cos_chi_f (:,:)
      real(kind=8), allocatable, save ::    sin_chi_f (:,:)
      real(kind=8), allocatable, save ::    cos_chi_b (:,:)
      real(kind=8), allocatable, save ::    sin_chi_b (:,:)

      real(kind=8), allocatable, save ::    old_dlam (:)
      real(kind=8), allocatable, save ::    old_dphi (:)

      real(kind=8)    lon_cmp(lattice%imglobal), lon_geo(lattice%imglobal)
      real(kind=8)    lat_cmp(lattice%jmglobal), lat_geo(lattice%jmglobal)
      real(kind=8)    pi,cosnp,sinnp,p1,p2,p3,eps
      real(kind=8)    dl,dp,lam_np,phi_np,lam_0,eps_np
      real(kind=8)    utmp,vtmp
      real(kind=8)    lam_geo, lam_cmp
      real(kind=8)    phi_geo, phi_cmp
      real(kind=8)    undef, getcon

      logical compute_weights
      real(kind=8)    old_rotation
      real(kind=8)    old_tilt
      real(kind=8)    old_precession
      data    old_rotation   /-999.9/
      data    old_tilt       /-999.9/
      data    old_precession /-999.9/

      parameter ( eps = 1.e-10 )

c Initialization
c --------------
      pi = 4.*atan(1.)
      dl = 2*pi/ lattice%imglobal     ! Uniform Grid Delta Lambda
      dp =   pi/(lattice%jmglobal-1)  ! Uniform Grid Delta Phi

c Allocate Memory for Weights and Index Locations
c -----------------------------------------------
      if(.not.allocated(old_dlam)) then

      allocate (  old_dlam(lattice%imglobal)    ,  old_dphi(lattice%jmglobal)    )
      allocate ( cos_chi_f(im,jm) , sin_chi_f(im,jm) )
      allocate ( cos_chi_b(im,jm) , sin_chi_b(im,jm) )
      do i=1,lattice%imglobal
      old_dlam(i) = 0.0
      enddo
      do j=1,lattice%jmglobal
      old_dphi(j) = 0.0
      enddo

      else
             i =  size (old_dlam)
             j =  size (old_dphi)
          if(i.ne.lattice%imglobal .or. j.ne.lattice%jmglobal) then
          deallocate (  old_dlam ,  old_dphi )
          deallocate ( cos_chi_f , sin_chi_f )
          deallocate ( cos_chi_b , sin_chi_b )
          allocate (  old_dlam(lattice%imglobal)    ,  old_dphi(lattice%jmglobal)    )
          allocate ( cos_chi_f(im,jm) , sin_chi_f(im,jm) )
          allocate ( cos_chi_b(im,jm) , sin_chi_b(im,jm) )
          do i=1,lattice%imglobal
          old_dlam(i) = 0.0
          enddo
          do j=1,lattice%jmglobal
          old_dphi(j) = 0.0
          enddo
          endif
      endif


c Compute Geophysical and Computational Lambda's and Phi's
c --------------------------------------------------------
      lon_geo(1) = -pi
      lon_cmp(1) = -pi
      do i=2,lattice%imglobal
      lon_geo(i) = lon_geo(i-1) + dl
      lon_cmp(i) = lon_cmp(i-1) + dlam(i-1)
      enddo
      lat_geo(1) = -pi*0.5
      lat_cmp(1) = -pi*0.5
      do j=2,lattice%jmglobal-1
      lat_geo(j) = lat_geo(j-1) + dp
      lat_cmp(j) = lat_cmp(j-1) + dphi(j-1)
      enddo
      lat_geo(lattice%jmglobal) =  pi*0.5
      lat_cmp(lattice%jmglobal) =  pi*0.5


c Check for Co-incident Grid-Point Latitude and Pole Locations
c ------------------------------------------------------------
                                          eps_np = 0.0
      do j=1,lattice%jmglobal
      phi_cmp = lat_cmp(j)*180./pi
      if( abs( phi_cmp-tilt ).lt.1.e-3 )  eps_np =  1.e-3
      if( tilt+eps_np .gt. 90. )          eps_np = -1.e-3
      enddo

      lam_np = pi/180.*rotation
      phi_np = pi/180.*(tilt+eps_np)
      lam_0  = pi/180.*precession

      undef = getcon('UNDEF')
      if( tilt.eq.90. ) then 
      cosnp = 0.0
      sinnp = 1.0
      else if(tilt.eq.-90.0) then
      cosnp = 0.0
      sinnp =-1.0
      else
      cosnp = cos(phi_np)
      sinnp = sin(phi_np)
      endif

c Determine if Weights Need to be Updated
c ---------------------------------------
      compute_weights =  rotation.ne.old_rotation .or.
     .                       tilt.ne.old_tilt     .or.
     .                 precession.ne.old_precession

      i = 1
      do while ( .not.compute_weights .and. i.le.lattice%imglobal )
      compute_weights = dlam(i).ne.old_dlam(i)
      i = i+1
      enddo
      j = 1
      do while ( .not.compute_weights .and. j.le.lattice%jmglobal-1 )
      compute_weights = dphi(j).ne.old_dphi(j)
      j = j+1
      enddo

c Computing Rotational Transformations for Vector Wind Field
c ----------------------------------------------------------
      if( compute_weights ) then

      old_rotation   = rotation
      old_tilt       = tilt
      old_precession = precession

      do i=1,lattice%imglobal
      old_dlam(i) =  dlam(i)
      enddo
      do j=1,lattice%jmglobal-1
      old_dphi(j) = dphi(j)
      enddo
      old_dphi(lattice%jmglobal) = dphi(lattice%jmglobal)

      do j=1,jm
      do i=1,im

      if( tilt.eq.90 .or. tilt.eq.-90 ) then

      cos_chi_f(i,j) = sinnp
      sin_chi_f(i,j) = 0.0
      cos_chi_b(i,j) = sinnp
      sin_chi_b(i,j) = 0.0

      else 

c Forward Rotation
c ----------------
      lam_cmp = lon_cmp( lattice%iglobal(i) )
      phi_cmp = lat_cmp( lattice%jglobal(j) )

      p1 = cosnp*cos(phi_cmp)*cos(lam_cmp-lam_np)
     .   + sin(phi_cmp)*sinnp
      p1 = min(p1, 1.0)
      p1 = max(p1,-1.0)
      phi_geo = asin( p1 )

      p2 = sinnp*cos(phi_cmp)*cos(lam_cmp-lam_np)
     .   - sin(phi_cmp)*cosnp
      p2 = p2 / max( cos(phi_geo),eps )
      p2 = min(p2, 1.0)
      p2 = max(p2,-1.0)
      p2 = acos( p2 )

      p3 = cos(phi_cmp)*sin(lam_cmp-lam_np)
      if( p3.lt.0.0 ) p2 = -p2 
      p2 = p2 - lam_0
      lam_geo = mod( p2+3.0*pi,2.0*pi ) - pi

      cos_chi_f(i,j) = (sinnp*cos(phi_cmp)
     .                 -cosnp*sin(phi_cmp)*cos(lam_cmp-lam_np))
     .               /  max( cos(phi_geo),eps )
      sin_chi_f(i,j) = (cosnp*sin(lam_cmp-lam_np))
     .               /  max( cos(phi_geo),eps )


c Backward Rotation
c -----------------
      lam_geo = lon_geo( lattice%iglobal(i) )
      phi_geo = lat_geo( lattice%jglobal(j) )

      p1 = cosnp*cos(phi_geo)*cos(lam_geo+lam_0-pi)
     .   + sin(phi_geo)*sinnp
      p1 = min(p1, 1.0)
      p1 = max(p1,-1.0)
      phi_cmp = asin( p1 )

      p2 = sinnp*cos(phi_geo)*cos(lam_geo+lam_0-pi)
     .   - sin(phi_geo)*cosnp
      p2 = p2 / max( cos(phi_cmp),eps )
      p2 = min(p2, 1.0)
      p2 = max(p2,-1.0)
      p2 = acos( p2 )

      p3 = cos(phi_geo)*sin(lam_geo+lam_0-pi)
      if( p3.lt.0.0 ) p2 = -p2
      p2 = p2 + lam_np - pi
      lam_cmp = mod( p2+3.0*pi,2.0*pi ) - pi

      cos_chi_b(i,j) = (sinnp*cos(phi_geo)
     .                 -cosnp*sin(phi_geo)*cos(lam_geo-pi+lam_0))
     .               /  max( cos(phi_cmp),eps )
      sin_chi_b(i,j) = (cosnp*sin(lam_geo-pi+lam_0))
     .               /  max( cos(phi_cmp),eps )

      endif
      enddo
      enddo

      endif

c Compute New Vector Wind Field
c -----------------------------
      if( msgn.eq.1 ) then
      do L=1,lm
         j=1
      do i=1,im*jm
      if( u(i,j,L).ne.undef .and. v(i,j,L).ne.undef ) then
      utmp = cos_chi_f(i,j)*u(i,j,L) - sin_chi_f(i,j)*v(i,j,L)
      vtmp = cos_chi_f(i,j)*v(i,j,L) + sin_chi_f(i,j)*u(i,j,L)
      else
      utmp = undef
      vtmp = undef
      endif
      u(i,j,L) = utmp
      v(i,j,L) = vtmp
      enddo
      enddo
      endif

      if( msgn.eq.-1 ) then
      do L=1,lm
         j=1
      do i=1,im*jm
      if( u(i,j,L).ne.undef .and. v(i,j,L).ne.undef ) then
      utmp = cos_chi_b(i,j)*u(i,j,L) - sin_chi_b(i,j)*v(i,j,L)
      vtmp = cos_chi_b(i,j)*v(i,j,L) + sin_chi_b(i,j)*u(i,j,L)
      else
      utmp = undef
      vtmp = undef
      endif
      u(i,j,L) = utmp
      v(i,j,L) = vtmp
      enddo
      enddo
      endif

      return
      end

      subroutine atoc ( qa,qc,dlam,dphi,im,jm,km,itype,lattice )

C ******************************************************************
C ****                                                          ****
C ****  This program converts 'A' gridded data                  ****
C ****                     to 'C' gridded data.                 ****
C ****                                                          ****
C ****      A direct transfer   is made      for itype = 0      ****
C ****      Cubic Interpolation is made in x for itype = 1      ****
C ****      Cubic Interpolation is made in y for itype = 2      ****
C ****                                                          ****
C ******************************************************************

      use g3_dynamics_lattice_module
      implicit none
      type ( dynamics_lattice_type ) lattice
      integer  im,jm,km,itype
      integer  im1,im2,ip1,ip2
      integer  jm1,jm2,jp1,jp2,jmx

      real(kind=8) dlam(lattice%imglobal)
      real(kind=8) dphi(lattice%jmglobal)
      real(kind=8) qa  (im,jm,km)
      real(kind=8) qc  (im,jm,km)

      integer    i,j,k
      real(kind=8) qxg (-1:im+2,1:jm), dxg(-1:im+2)
      real(kind=8) qyg (1:im,-1:jm+2), dyg(-1:jm+2)

      real(kind=8)   ap2(im), ap1(im), ap0(im), am1(im)
      real(kind=8)   bp2(jm), bp1(jm), bp0(jm), bm1(jm)

C *********************************************************
C ****                mass points                      ****
C *********************************************************

      if (itype.eq.0) then

      do k = 1,km
      do j = 1,jm
      do i = 1,im
      qc(i,j,k) = qa(i,j,k)
      enddo
      enddo
      enddo

      endif

C *********************************************************
C ****           stagger in x-direction                ****
C *********************************************************

      if (itype.eq.1) then

c Ghost the variable lon-grid increments
c --------------------------------------
      do i=1,im
      dxg(i) = dlam( lattice%iglobal(i) )
      enddo

      im1 = lattice%iglobal(1) -1
      im2 = lattice%iglobal(1) -2
      ip1 = lattice%iglobal(im)+1
      ip2 = lattice%iglobal(im)+2

      if( im1.eq.0 ) then
          im1 = lattice%imglobal
          im2 = lattice%imglobal-1
      endif
      if( im2.eq.0 ) then
          im2 = lattice%imglobal
      endif
      if( ip1.eq.lattice%imglobal+1 ) then
          ip1 = 1
          ip2 = 2
      endif
      if( ip2.eq.lattice%imglobal+1 ) then
          ip2 = 1
      endif

      dxg( 0)   = dlam( im1 )
      dxg(-1)   = dlam( im2 )
      dxg(im+1) = dlam( ip1 )
      dxg(im+2) = dlam( ip2 )


      do i=1,im
      ap2(i) = ( -0.5*dxg(i)*0.5*dxg(i)*(0.5*dxg(i)+dxg(i-1)) )
     .       / ( dxg(i+1)*(dxg(i+1)+dxg(i))*(dxg(i+1)+dxg(i)+dxg(i-1)) )
      ap1(i) = ( (dxg(i+1)+0.5*dxg(i))*0.5*dxg(i)*(0.5*dxg(i)+dxg(i-1)) )
     .       / ( dxg(i+1)*dxg(i)*(dxg(i)+dxg(i-1)) )
      ap0(i) = ( (dxg(i+1)+0.5*dxg(i))*0.5*dxg(i)*(0.5*dxg(i)+dxg(i-1)) )
     .       / ( dxg(i-1)*dxg(i)*(dxg(i)+dxg(i+1)) )
      am1(i) = ( -0.5*dxg(i)*0.5*dxg(i)*(0.5*dxg(i)+dxg(i+1)) )
     .       / ( dxg(i-1)*(dxg(i-1)+dxg(i))*(dxg(i+1)+dxg(i)+dxg(i-1)) )
      enddo

      do k=1,km
      call ghostx ( qa(1,1,k),qxg,im,jm,1,2,lattice,'both' )
      do j=1,jm
      do i=1,im
      qc(i,j,k) = ap2(i)*qxg(i+2,j) + ap1(i)*qxg(i+1,j)
     .          + ap0(i)*qxg(i  ,j) + am1(i)*qxg(i-1,j)
      enddo
      enddo
      enddo

      endif

C *********************************************************
C ****           stagger in y-direction                ****
C *********************************************************

      if (itype.eq.2) then

c Ghost the variable lat-grid increments
c --------------------------------------
                                        jmx = jm
      if( lattice%pej.eq.lattice%ny-1 ) jmx = jm-1

      do j=1,jmx
      dyg(j) = dphi( lattice%jglobal(j) )
      enddo

      jm1 = lattice%jglobal(1)  -1
      jm2 = lattice%jglobal(1)  -2
      jp1 = lattice%jglobal(jmx)+1
      jp2 = lattice%jglobal(jmx)+2

      if( jm1.eq.0 ) then
          jm1 = 1
          jm2 = 2
      endif
      if( jm2.eq.0 ) then
          jm2 = 1
      endif
      if( jp1.eq.lattice%jmglobal ) then
          jp1 =  lattice%jmglobal-1
          jp2 =  lattice%jmglobal-2
      endif
      if( jp2.eq.lattice%jmglobal ) then
          jp2 =  lattice%jmglobal-1
      endif

      dyg( 0)    = dphi( jm1 )
      dyg(-1)    = dphi( jm2 )
      dyg(jmx+1) = dphi( jp1 )
      dyg(jmx+2) = dphi( jp2 )

      do j=1,jmx
      bp2(j) = ( -0.5*dyg(j)*0.5*dyg(j)*(0.5*dyg(j)+dyg(j-1)) )
     .       / ( dyg(j+1)*(dyg(j+1)+dyg(j))*(dyg(j+1)+dyg(j)+dyg(j-1)) )
      bp1(j) = ( (dyg(j+1)+0.5*dyg(j))*0.5*dyg(j)*(0.5*dyg(j)+dyg(j-1)) )
     .       / ( dyg(j+1)*dyg(j)*(dyg(j)+dyg(j-1)) )
      bp0(j) = ( (dyg(j+1)+0.5*dyg(j))*0.5*dyg(j)*(0.5*dyg(j)+dyg(j-1)) )
     .       / ( dyg(j-1)*dyg(j)*(dyg(j)+dyg(j+1)) )
      bm1(j) = ( -0.5*dyg(j)*0.5*dyg(j)*(0.5*dyg(j)+dyg(j+1)) )
     .       / ( dyg(j-1)*(dyg(j-1)+dyg(j))*(dyg(j+1)+dyg(j)+dyg(j-1)) )
      enddo

      do k=1,km
      call ghosty ( qa(1,1,k),qyg,im,jm,1,0,-1,2,lattice,'both' )
      do j=1,jmx
      do i=1,im
      qc(i,j,k) = bp2(j)*qyg(i,j+2) + bp1(j)*qyg(i,j+1)
     .          + bp0(j)*qyg(i,j  ) + bm1(j)*qyg(i,j-1)
      enddo
      enddo
      enddo

      endif

      return
      end

      subroutine ctoa ( qc,qa,dlam,dphi,im,jm,km,itype,lattice )
                                                                        
C ******************************************************************    
C ****                                                          ****    
C ****  This program converts 'C' gridded data                  ****    
C ****                     to 'A' gridded data.                 ****    
C ****                                                          ****    
C ****  A pos.definite transfer     is made      for itype = 0  ****
C ****  A Cubic Interpolation       is made in x for itype = 1  ****
C ****  A Cubic Interpolation       is made in y for itype = 2  ****
C ****  A direct transfer           is made      for itype = 3  ****
C ****  A direct transfer with undef check       for itype = 4  ****
C ****                                                          ****    
C ******************************************************************    

      use g3_dynamics_lattice_module
      implicit none
      type ( dynamics_lattice_type ) lattice
      integer  im,jm,km,itype
      integer  im1,im2,ip1,ip2
      integer  jm1,jm2,jp1,jp2,jmx

      real(kind=8) dlam(lattice%imglobal)
      real(kind=8) dphi(lattice%jmglobal)
      real(kind=8) qa  (im,jm,km)                                         
      real(kind=8) qc  (im,jm,km)                                         

      real(kind=8)        cnp1,cnp2
      real(kind=8)        csp1,csp2

      integer  i,j,k
      real(kind=8)     qnp,qsp,getcon,undef
      real(kind=8)     qxg (-1:im+2,1:jm), dxg(-1:im+2)
      real(kind=8)     qyg (1:im,-1:jm+2), dyg(-1:jm+2)
                                                                        
      real(kind=8)   ap1(im), ap0(im), am1(im), am2(im)
      real(kind=8)   bp1(jm), bp0(jm), bm1(jm), bm2(jm)

      undef = getcon('UNDEF')

C *********************************************************             
C ****                mass points                      ****             
C *********************************************************             
                                                                        
      if (itype.eq.0 .or. itype.eq.3 .or. itype.eq.4) then

      cnp1 = ( dphi(lattice%jmglobal-1)+ dphi(lattice%jmglobal-2) )**2
     .     / ( dphi(lattice%jmglobal-2)*(dphi(lattice%jmglobal-2)+2*dphi(lattice%jmglobal-1)) )
      cnp2 = 1.0-cnp1
      cnp1 = 0.5*cnp1
      cnp2 = 0.5*cnp2

      csp1 = ( dphi(1)+ dphi(2) )**2
     .     / ( dphi(2)*(dphi(2)+2*dphi(1)) )
      csp2 = 1.0-csp1
      csp1 = 0.5*csp1
      csp2 = 0.5*csp2


      do k = 1,km

c Direct Transfer of Interior Points
c ----------------------------------
      do j = 1,jm
      do i = 1,im
      qa(i,j,k) = qc(i,j,k)                                             
      enddo
      enddo

c North Pole Point
c ----------------
      if( lattice%pej.eq.lattice%ny-1 ) then

      call ghosty ( qc(1,1,k),qyg,im,jm,1,0,1,2,lattice,'pole' )

c Interpolate to Pole Value (4th-order) w/o UNDEF check
c -----------------------------------------------------
      if (itype.eq.0 .or. itype.eq.3) then
      do i=1,im
      qa(i,jm,k) = cnp1*( qc(i,jm-1,k) + qyg(i,jm+1) )
     .           + cnp2*( qc(i,jm-2,k) + qyg(i,jm+2) )
      enddo
      call zmean ( qa(1,jm,k),qnp,dlam,im,1,0,.false.,lattice )
      endif

c Interpolate to Pole Value (4th-order) w/UNDEF check
c ---------------------------------------------------
      if (itype.eq.4) then
      do i=1,im
         if( qc(i,jm-1,k).eq.undef .or. qyg(i,jm+1).eq.undef  .or.
     .       qc(i,jm-2,k).eq.undef .or. qyg(i,jm+2).eq.undef ) then
         qa(i,jm,k) = undef
         else
         qa(i,jm,k) = cnp1*( qc(i,jm-1,k) + qyg(i,jm+1) )
     .              + cnp2*( qc(i,jm-2,k) + qyg(i,jm+2) )
         endif
      enddo
      call zmean ( qa(1,jm,k),qnp,dlam,im,1,0,.true.,lattice )
      endif

c Unify Pole Values
c -----------------
      if( itype.eq.0 ) qnp = max( qnp,0.0 ) ! Do not allow Constituents < 0
      do i=1,im
      qa(i,jm,k) = qnp
      enddo

      endif   ! End North Pole Check

c South Pole Point
c ----------------
      if( lattice%pej.eq.0 ) then

      call ghosty ( qc(1,1,k),qyg,im,jm,1,0,1,2,lattice,'pole' )

c Interpolate to Pole Value (4th-order) w/o UNDEF check
c -----------------------------------------------------
      if (itype.eq.0 .or. itype.eq.3) then
      do i=1,im
      qa(i,1,k) = cnp1*( qc(i,2,k) + qyg(i, 0) )
     .          + cnp2*( qc(i,3,k) + qyg(i,-1) )
      enddo
      call zmean ( qa(1,1,k),qsp,dlam,im,1,0,.false.,lattice )
      endif

c Interpolate to Pole Value (4th-order) w/UNDEF check
c ---------------------------------------------------
      if (itype.eq.4) then
      do i=1,im
         if( qc(i,2,k).eq.undef .or. qyg(i, 0).eq.undef  .or.
     .       qc(i,3,k).eq.undef .or. qyg(i,-1).eq.undef ) then
         qa(i,1,k) = undef
         else
         qa(i,1,k) = cnp1*( qc(i,2,k) + qyg(i, 0) )
     .             + cnp2*( qc(i,3,k) + qyg(i,-1) )
         endif
      enddo
      call zmean ( qa(1,1,k),qsp,dlam,im,1,0,.true.,lattice )
      endif

c Unify Pole Values
c -----------------
      if( itype.eq.0 ) qsp = max( qsp,0.0 ) ! Do not allow Constituents < 0
      do i=1,im
      qa(i,1 ,k) = qsp
      enddo

      endif   ! End South Pole Check

      enddo   ! End K-Loop
      endif

C *********************************************************             
C ****           stagger in x-direction                ****             
C *********************************************************             
                                                                        
      if (itype.eq.1) then                                              

c Ghost the variable lon-grid increments
c --------------------------------------
      do i=1,im
      dxg(i) = dlam( lattice%iglobal(i) )
      enddo

      im1 = lattice%iglobal(1) -1
      im2 = lattice%iglobal(1) -2
      ip1 = lattice%iglobal(im)+1
      ip2 = lattice%iglobal(im)+2

      if( im1.eq.0 ) then
          im1 = lattice%imglobal
          im2 = lattice%imglobal-1
      endif
      if( im2.eq.0 ) then
          im2 = lattice%imglobal
      endif
      if( ip1.eq.lattice%imglobal+1 ) then
          ip1 = 1
          ip2 = 2
      endif
      if( ip2.eq.lattice%imglobal+1 ) then
          ip2 = 1
      endif

      dxg( 0)   = dlam( im1 )
      dxg(-1)   = dlam( im2 )
      dxg(im+1) = dlam( ip1 )
      dxg(im+2) = dlam( ip2 )


      do i=1,im
      ap1(i) = ( -0.5*dxg(i-1)*0.5*dxg(i-1)*(0.5*dxg(i-1)+dxg(i-2)) )
     .       / ( dxg(i)*(dxg(i)+dxg(i-1))*(dxg(i)+dxg(i-1)+dxg(i-2)) )
      ap0(i) = ( (dxg(i)+0.5*dxg(i-1))*0.5*dxg(i-1)*(0.5*dxg(i-1)+dxg(i-2)) )
     .       / ( dxg(i)*dxg(i-1)*(dxg(i-1)+dxg(i-2)) )
      am1(i) = ( (dxg(i)+0.5*dxg(i-1))*0.5*dxg(i-1)*(0.5*dxg(i-1)+dxg(i-2)) )
     .       / ( dxg(i-2)*dxg(i-1)*(dxg(i-1)+dxg(i)) )
      am2(i) = ( -0.5*dxg(i-1)*0.5*dxg(i-1)*(0.5*dxg(i-1)+dxg(i)) )
     .       / ( dxg(i-2)*(dxg(i-2)+dxg(i-1))*(dxg(i)+dxg(i-1)+dxg(i-2)) )
      enddo


      do k=1,km                                                         
      call ghostx ( qc(1,1,k),qxg,im,jm,1,2,lattice,'both' )
      do j=1,jm
      do i=1,im
      qa(i,j,k) = ap1(i)*qxg(i+1,j) + ap0(i)*qxg(i  ,j)
     .          + am1(i)*qxg(i-1,j) + am2(i)*qxg(i-2,j)
      enddo
      enddo
      enddo

      endif                                                             
                                                                        
C *********************************************************             
C ****           stagger in y-direction                ****             
C *********************************************************             
                                                                        
      if (itype.eq.2) then                                              
                                                                        
c Ghost the variable lat-grid increments
c --------------------------------------
                                        jmx = jm
      if( lattice%pej.eq.lattice%ny-1 ) jmx = jm-1

      do j=1,jmx
      dyg(j) = dphi( lattice%jglobal(j) )
      enddo

      jm1 = lattice%jglobal(1)  -1
      jm2 = lattice%jglobal(1)  -2
      jp1 = lattice%jglobal(jmx)+1
      jp2 = lattice%jglobal(jmx)+2

      if( jm1.eq.0 ) then
          jm1 = 1
          jm2 = 2
      endif
      if( jm2.eq.0 ) then
          jm2 = 1
      endif
      if( jp1.eq.lattice%jmglobal ) then
          jp1 =  lattice%jmglobal-1
          jp2 =  lattice%jmglobal-2
      endif
      if( jp2.eq.lattice%jmglobal ) then
          jp2 =  lattice%jmglobal-1
      endif

      dyg( 0)    = dphi( jm1 )
      dyg(-1)    = dphi( jm2 )
      dyg(jmx+1) = dphi( jp1 )
      dyg(jmx+2) = dphi( jp2 )


      do j=1,jm
      bp1(j) = ( -0.5*dyg(j-1)*0.5*dyg(j-1)*(0.5*dyg(j-1)+dyg(j-2)) )
     .       / ( dyg(j)*(dyg(j)+dyg(j-1))*(dyg(j)+dyg(j-1)+dyg(j-2)) )
      bp0(j) = ( (dyg(j)+0.5*dyg(j-1))*0.5*dyg(j-1)*(0.5*dyg(j-1)+dyg(j-2)) )
     .       / ( dyg(j)*dyg(j-1)*(dyg(j-1)+dyg(j-2)) )
      bm1(j) = ( (dyg(j)+0.5*dyg(j-1))*0.5*dyg(j-1)*(0.5*dyg(j-1)+dyg(j-2)) )
     .       / ( dyg(j-2)*dyg(j-1)*(dyg(j-1)+dyg(j)) )
      bm2(j) = ( -0.5*dyg(j-1)*0.5*dyg(j-1)*(0.5*dyg(j-1)+dyg(j)) )
     .       / ( dyg(j-2)*(dyg(j-2)+dyg(j-1))*(dyg(j)+dyg(j-1)+dyg(j-2)) )
      enddo


      do k=1,km                                                         
      call ghosty ( qc(1,1,k),qyg,im,jm,1,1,-1,2,lattice,'both' )
      do j=1,jm 
      do i=1,im
      qa(i,j,k) = bp1(j)*qyg(i,j+1) + bp0(j)*qyg(i  ,j)
     .          + bm1(j)*qyg(i,j-1) + bm2(j)*qyg(i,j-2)
      enddo
      enddo
      enddo

      endif

      return
      end

      subroutine polewnd (uz,vz,ustr,vstr,dlam,dphi,im,jm,lattice)
C***********************************************************************
C  PURPOSE     To compute fictitious winds at and near the pole
C  ARGUMENTS   DESCRIPTION
C     Uz ..... U-WIND on C-Grid
C     Vz ..... V-WIND on C-Grid
C     Ustr ... U-WIND on C-Grid at   pole
C     Vstr ... V-WIND on C-Grid near pole
C
C***********************************************************************
C*                  GODDARD LABORATORY FOR ATMOSPHERES                 *
C***********************************************************************

      use g3_dynamics_lattice_module
      implicit none
      type ( dynamics_lattice_type ) lattice

      integer    im,jm

      real(kind=8) Uz   (im,jm)
      real(kind=8) Vz   (im,jm)
      real(kind=8) ustr (im,jm)
      real(kind=8) vstr (im,jm)

      real(kind=8) dlam (lattice%imglobal)
      real(kind=8) dphi (lattice%jmglobal)

c Local Dynamic Space
c -------------------
      real(kind=8) cosu (lattice%imglobal)
      real(kind=8) cosv (lattice%imglobal)
      real(kind=8) sinu (lattice%imglobal)
      real(kind=8) sinv (lattice%imglobal)
      real(kind=8) lamu (lattice%imglobal)
      real(kind=8) lamv (lattice%imglobal)
      real(kind=8)  up  (im)
      real(kind=8)  upc (im)
      real(kind=8)  ups (im)
      real(kind=8)  vp  (im)

      real(kind=8)    pi,a1,a2,a3,uc,us,vc,vs
      integer i,j,jmg

      real(kind=8), allocatable :: ug(:,:)
      real(kind=8), allocatable :: vg(:,:)

      pi  = 4.0*atan(1.0)
      jmg = lattice%jmglobal

C **********************************************************************
C ****                    Initialize Winds                          ****
C **********************************************************************
       
      do j=1,jm
      do i=1,im
      ustr(i,j) = uz(i,j)
      vstr(i,j) = vz(i,j)
      enddo
      enddo
       
C **********************************************************************
C ****            Compute Winds at and near the Pole                ****
C **********************************************************************
       
      if( lattice%pej.eq.0 .or. lattice%pej.eq.lattice%ny-1 ) then

c Compute Computational Lambda's at U & V points
c ----------------------------------------------
      lamv(1) = -pi
      lamu(1) = lamv(1)   + dlam(1)*0.5
      do i=2,lattice%imglobal
      lamv(i) = lamv(i-1) + dlam(i-1)
      lamu(i) = lamv(i)   + dlam(i)*0.5
      enddo

      do i=1,lattice%imglobal
      cosu(i) = cos( lamu(i) )
      sinu(i) = sin( lamu(i) )
      cosv(i) = cos( lamv(i) )
      sinv(i) = sin( lamv(i) )
      enddo

      endif

c North Pole
c ----------
      if( lattice%pej.eq.lattice%ny-1 ) then

      allocate ( ug(1:im,-1:jm+2) )
      allocate ( vg(1:im,-1:jm+2) )
      call ghosty ( uz,ug,im,jm,1,0,-1,2,lattice,'pole' )
      call ghosty ( vz,vg,im,jm,1,1,-1,2,lattice,'pole' )

      a1 = ( dphi(jmg-1)+ dphi(jmg-2) )**2
     .   / ( dphi(jmg-2)*(dphi(jmg-2)+2*dphi(jmg-1)) )
      a2 = 1.0-a1
      a1 = a1*0.5
      a2 = a2*0.5

      do i=1,im
      up (i) = a1*( ug(i,jm-1)+ug(i,jm+1) )
     .       + a2*( ug(i,jm-2)+ug(i,jm+2) )
      upc(i) = up(i)*cosu( lattice%iglobal(i) )
      ups(i) = up(i)*sinu( lattice%iglobal(i) )
      enddo

      call zmean ( upc,uc,dlam,im,1,1,.false.,lattice )
      call zmean ( ups,us,dlam,im,1,1,.false.,lattice )

      uc = 2.0*uc
      us = 2.0*us
      vc =     us
      vs =    -uc

      a3 = - dphi(jmg-1)* (dphi(jmg-1)+dphi(jmg-2))*0.5
     .   / ( dphi(jmg-2)+2*dphi(jmg-1) )**2
      a2 = ( dphi(jmg-1) + a3*(2*dphi(jmg-1)+dphi(jmg-2)) )
     .   / ( dphi(jmg-2) + 2*dphi(jmg-1) )
      a1 = 1.0-a2-a3

      do i=1,im
      ustr(i,jm) = uc*cosu( lattice%iglobal(i) ) + us*sinu( lattice%iglobal(i) )
        vp(i)    = vc*cosv( lattice%iglobal(i) ) + vs*sinv( lattice%iglobal(i) )
      enddo

      do i=1,im
      vstr(i,jm-1) = a1*vp(i) + a2*vg(i,jm-2) + a3*vg(i,jm+1)
      enddo

      deallocate ( ug,vg )
      endif


c South Pole
c ----------
      if( lattice%pej.eq.0 ) then

      allocate ( ug(1:im,-1:jm+2) )
      allocate ( vg(1:im,-1:jm+2) )
      call ghosty ( uz,ug,im,jm,1,0,-1,2,lattice,'pole' )
      call ghosty ( vz,vg,im,jm,1,1,-1,2,lattice,'pole' )

      a1 = ( dphi(1)+ dphi(2) )**2
     .   / ( dphi(2)*(dphi(2)+2*dphi(1)) )
      a2 = 1.0-a1
      a1 = a1*0.5
      a2 = a2*0.5

      do i=1,im
      up (i) = a1*( ug(i,2)+ug(i, 0) )
     .       + a2*( ug(i,3)+ug(i,-1) )
      upc(i) = up(i)*cosu( lattice%iglobal(i) )
      ups(i) = up(i)*sinu( lattice%iglobal(i) )
      enddo

      call zmean ( upc,uc,dlam,im,1,1,.false.,lattice )
      call zmean ( ups,us,dlam,im,1,1,.false.,lattice )

      uc = 2.0*uc
      us = 2.0*us
      vc =    -us
      vs =     uc

      a3 = - dphi(1)* (dphi(1)+dphi(2))*0.5
     .   / ( dphi(2)+2*dphi(1) )**2
      a2 = ( dphi(1) + a3*(2*dphi(1)+dphi(2)) )
     .   / ( dphi(2) + 2*dphi(1) )
      a1 = 1.0-a2-a3

      do i=1,im
      ustr(i,1) = uc*cosu( lattice%iglobal(i) ) + us*sinu( lattice%iglobal(i) )
        vp(i)   = vc*cosv( lattice%iglobal(i) ) + vs*sinv( lattice%iglobal(i) )
      enddo

      do i=1,im
      vstr(i,1) = a1*vp(i) + a2*vg(i,2) + a3*vg(i,-1)
      enddo

      deallocate ( ug,vg )
      endif

      return
      end

      subroutine ctoa_winds ( uc,vc,ua,va,dlam,dphi,im,jm,lm,lattice )
C***********************************************************************
C  PURPOSE                                                              
C     Driver to convert C-Grid wind data to A-Grid wind data                               
C  ARGUMENTS   DESCRIPTION                                              
C     uc ..... U-WIND on C-Grid
C     vc ..... V-WIND on C-Grid
C     ua ..... U-WIND on A-Grid
C     va ..... V-WIND on A-Grid
C     im ..... x-dimension
C     jm ..... y-dimension
C     lm ..... z-dimension
C
C***********************************************************************
C*                  GODDARD LABORATORY FOR ATMOSPHERES                 *
C***********************************************************************

      use g3_dynamics_lattice_module
      implicit none
      type ( dynamics_lattice_type ) lattice

c Input Variables
c ---------------
      integer  im,jm,lm
      real(kind=8)  uc(im,jm,lm)
      real(kind=8)  vc(im,jm,lm)
      real(kind=8)  ua(im,jm,lm)
      real(kind=8)  va(im,jm,lm)

      real(kind=8) dlam(im)
      real(kind=8) dphi(jm)

c Local Variables
c ---------------
      integer L
      real(kind=8)    ut(im,jm,lm)
      real(kind=8)    vt(im,jm,lm)

      do L=1,lm
      call polewnd ( uc(1,1,L),vc(1,1,L),                           
     .               ut(1,1,L),vt(1,1,L),dlam,dphi,im,jm,lattice )
      enddo                                                             
      call ctoa ( ut(1,1,1),ua(1,1,1),dlam,dphi,im,jm,lm,1,lattice )
      call ctoa ( vt(1,1,1),va(1,1,1),dlam,dphi,im,jm,lm,2,lattice )

      return
      end

      subroutine ctoaset ( vars,grid )
C***********************************************************************
C  PURPOSE                                                              
C     Converts C-Grid Dynamics Variables to A-Grid Dynamics Variables                 
C                                                                       
C***********************************************************************
C*                  GODDARD LABORATORY FOR ATMOSPHERES                 *
C***********************************************************************
                                                                        
c Declare Modules and Data Structures
c -----------------------------------
      use    g3_dynamics_state_module
      type ( dynamics_vars_type ) vars
      type ( dynamics_grid_type ) grid

      call ctoa_winds ( vars%u,vars%v, vars%u,vars%v,
     .                  grid%dlam,grid%dphi,grid%im,grid%jm,grid%lm,grid%lattice )

      call ctoa ( vars%p   ,vars%p   ,grid%dlam,grid%dphi,grid%im,grid%jm,1      ,0,grid%lattice )         
      call ctoa ( vars%t   ,vars%t   ,grid%dlam,grid%dphi,grid%im,grid%jm,grid%lm,0,grid%lattice )         

      do m=1,grid%ntracer
      call ctoa ( vars%q(1,1,1,m),vars%q(1,1,1,m),grid%dlam,grid%dphi,grid%im,grid%jm,grid%lm,0,grid%lattice )
      enddo

      return                                                            
      end                                                               

      subroutine atocset ( vars,grid )
C***********************************************************************
C  PURPOSE                                                              
C     Converts A-Grid Dynamics Variables to C-Grid Dynamics Variables    
C                                                                       
C***********************************************************************
C*                  GODDARD LABORATORY FOR ATMOSPHERES                 *
C***********************************************************************

c Declare Modules and Data Structures
c -----------------------------------
      use g3_dynamics_state_module
      type ( dynamics_vars_type ) vars
      type ( dynamics_grid_type ) grid

      call atoc ( vars%p,   vars%p   ,grid%dlam,grid%dphi,grid%im,grid%jm,1      ,0,grid%lattice )
      call atoc ( vars%u,   vars%u   ,grid%dlam,grid%dphi,grid%im,grid%jm,grid%lm,1,grid%lattice )
      call atoc ( vars%v,   vars%v   ,grid%dlam,grid%dphi,grid%im,grid%jm,grid%lm,2,grid%lattice )
      call atoc ( vars%t,   vars%t   ,grid%dlam,grid%dphi,grid%im,grid%jm,grid%lm,0,grid%lattice )

      do m=1,grid%ntracer
      call atoc ( vars%q(1,1,1,m),vars%q(1,1,1,m),grid%dlam,grid%dphi,grid%im,grid%jm,grid%lm,0,grid%lattice )
      enddo

      return
      end
