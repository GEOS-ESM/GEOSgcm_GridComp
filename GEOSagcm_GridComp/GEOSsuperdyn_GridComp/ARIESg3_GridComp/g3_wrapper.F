      subroutine g3_wrapper (dynamics,phis,scheme,dtphy,nsplit,alpha,omega)

C***********************************************************************        
C  Purpose                                                                      
C  -------                                                                      
C     Driver for the C-Grid Atmospheric GEOS GCM (Aires Dynamical Core)          
C
C  Arguments  Description                                                      
C  ----------------------
C     dynamics.. Dynamics State  Data Structure
C     scheme ... Character Description of Time-Scheme (MATS or LEAP)            
C     dtphy .... Time-Step in seconds for PHYSICS
C     nsplit ... Ratio of Time-Step in seconds for PHYSICS/DYNAMICS
C     alpha .... Time Filter Coefficient                                        
C                                                                               
C***********************************************************************        
C*                  GODDARD LABORATORY FOR ATMOSPHERES                 *        
C***********************************************************************        

c Declare Modules and Data Structures
c -----------------------------------
      use g3_dynamics_state_module
      implicit none

c Declare Data Structures
c -----------------------
      type ( dynamics_state_type ), target :: dynamics

c Input Parameters
c ----------------
      real(kind=8)   phis(dynamics%grid%im,dynamics%grid%jm)
      real(kind=8)  omega(dynamics%grid%im,dynamics%grid%jm,dynamics%grid%lm)
      real(kind=8)  dtphy
      real*4        alpha
      integer       nsplit
      character*4   scheme

c Local Variables
c ---------------
      real(kind=8)  dtdyn
      integer       im,jm,lm,nq

      character*4 leapfrog, matsuno
      logical     alarm
      logical     first
      data        first /.true./                                                
      data        leapfrog /'LEAP'/                                                
      data        matsuno  /'MATS'/                                                
                                                                                
      logical  brown_campana
      real(kind=8)     tcoef

c Local Reference to Dynamics Grid
c --------------------------------
c     integer, pointer :: n
      integer             n, nm1

c Temporary Variables for NP1 State
c ---------------------------------
      type (  dynamics_vars_type ) dynamics_vars_np1

c Temporary Variables for Total Tendency
c --------------------------------------
      real(kind=8), allocatable ::  dpdt(:,:)
      real(kind=8), allocatable ::  dudt(:,:,:)
      real(kind=8), allocatable ::  dvdt(:,:,:)
      real(kind=8), allocatable ::  dtdt(:,:,:)
      real(kind=8), allocatable ::  dqdt(:,:,:,:)
      
c Miscellanious Temporary Variables
c ---------------------------------
      integer  i,k
      real(kind=8) getcon, timestep

C *********************************************************************
C *****            Allocate Temporary Work Space                   ****
C *********************************************************************

      im = dynamics%grid%im
      jm = dynamics%grid%jm
      lm = dynamics%grid%lm
      nq = dynamics%grid%ntracer

      call create_dynamics_vars ( dynamics_vars_np1,im,jm,lm,nq )

      allocate ( dpdt(im,jm)       )
      allocate ( dudt(im,jm,lm)    )
      allocate ( dvdt(im,jm,lm)    )
      allocate ( dtdt(im,jm,lm)    )
      allocate ( dqdt(im,jm,lm,nq) )

C **********************************************************************        
C ****                  Perform Inner Dynamics Loop                 ****        
C **********************************************************************        

      dtdyn = dtphy/nsplit

      do i=1,nsplit
                                                                                
      n   = dynamics%grid%n
      nm1 = mod(n,2)+1

C **********************************************************************
C ****                    Initialize Tendencies                     ****
C **********************************************************************
    
      dpdt(:,:)     = 0.0
      dudt(:,:,:)   = 0.0
      dvdt(:,:,:)   = 0.0
      dtdt(:,:,:)   = 0.0
      dqdt(:,:,:,:) = 0.0
 
C **********************************************************************
C ****                 Matsuno Predictor Time Scheme                ****
C **********************************************************************

      if(i.eq.1 .or. scheme.eq.matsuno) then
      tcoef         =  0.0
      timestep      =  dtdyn
      brown_campana = .false.

c Predictor
c ---------
      call g3_dyndrv( dynamics%vars(n),dynamics%vars(n),dynamics%grid,omega,
     .                dpdt,dudt,dvdt,dtdt,dqdt,
     .                phis,tcoef,timestep,brown_campana )

      call shapij   ( dynamics%vars(n)%p,dynamics%vars(n)%u,dynamics%vars(n)%v,
     .                dynamics%vars(n)%t,dynamics%vars(n)%q,
     .                dudt,dvdt,dtdt,dqdt,
     .                im,jm,lm,nq,dynamics%grid )

      call g3_step  ( dynamics%vars(n),dynamics%vars(n),dynamics_vars_np1,
     .                dynamics%grid,dpdt,dudt,dvdt,dtdt,dqdt,
     .                tcoef,timestep )

c Corrector
c ---------
      dpdt(:,:)     = 0.0
      dudt(:,:,:)   = 0.0
      dvdt(:,:,:)   = 0.0
      dtdt(:,:,:)   = 0.0
      dqdt(:,:,:,:) = 0.0

      call g3_dyndrv( dynamics%vars(n),dynamics_vars_np1,dynamics%grid,omega,
     .                dpdt,dudt,dvdt,dtdt,dqdt,
     .                phis,tcoef,timestep,brown_campana )

      call shapij   ( dynamics%vars(n)%p,dynamics%vars(n)%u,dynamics%vars(n)%v,
     .                dynamics%vars(n)%t,dynamics%vars(n)%q,
     .                dudt,dvdt,dtdt,dqdt,
     .                im,jm,lm,nq,dynamics%grid )

      call g3_step  ( dynamics%vars(n),dynamics_vars_np1,dynamics_vars_np1,
     .                dynamics%grid,dpdt,dudt,dvdt,dtdt,dqdt,
     .                tcoef,timestep )

      endif

C **********************************************************************
C ****                     Leapfrog Time Scheme                     ****
C **********************************************************************

      if(i.ne.1 .and. scheme.eq.leapfrog) then
      tcoef         =  alpha
      timestep      =  2*dtdyn
      brown_campana = .true.

      call g3_dyndrv( dynamics%vars(nm1),dynamics%vars(n),dynamics%grid,omega,
     .                dpdt,dudt,dvdt,dtdt,dqdt,
     .                phis,tcoef,timestep,brown_campana )

      call shapij   ( dynamics%vars(nm1)%p,dynamics%vars(nm1)%u,dynamics%vars(nm1)%v,
     .                dynamics%vars(nm1)%t,dynamics%vars(nm1)%q,
     .                dudt,dvdt,dtdt,dqdt,
     .                im,jm,lm,nq,dynamics%grid )

      call g3_step  ( dynamics%vars(nm1),dynamics%vars(n),dynamics_vars_np1,
     .                dynamics%grid,dpdt,dudt,dvdt,dtdt,dqdt,
     .                tcoef,timestep )
      endif

C **********************************************************************        
C ****             Move Updated Fields into Current Arrays          ****        
C **********************************************************************        

      if(scheme.eq.leapfrog) dynamics%grid%n = mod(n,2)+1
                                           n = dynamics%grid%n

      dynamics%vars(n)%p(:,:)     = dynamics_vars_np1%p(:,:)
      dynamics%vars(n)%u(:,:,:)   = dynamics_vars_np1%u(:,:,:)
      dynamics%vars(n)%v(:,:,:)   = dynamics_vars_np1%v(:,:,:)
      dynamics%vars(n)%t(:,:,:)   = dynamics_vars_np1%t(:,:,:)
      dynamics%vars(n)%q(:,:,:,:) = dynamics_vars_np1%q(:,:,:,:)

C **********************************************************************        
C ****                       Update Counter                         ****        
C **********************************************************************        

      enddo
                                                                                
C **********************************************************************        
C ****                Fill Pole Values for Scalars                  ****        
C **********************************************************************        

      call ctoa ( dynamics%vars(n)%p,dynamics%vars(n)%p,
     .            dynamics%grid%dlam,dynamics%grid%dphi, im,jm,1 ,0,dynamics%grid%lattice )
      call ctoa ( dynamics%vars(n)%t,dynamics%vars(n)%t,
     .            dynamics%grid%dlam,dynamics%grid%dphi, im,jm,lm,0,dynamics%grid%lattice )
      call ctoa ( omega,omega,
     .            dynamics%grid%dlam,dynamics%grid%dphi, im,jm,lm,0,dynamics%grid%lattice )
      do k=1,nq
      call ctoa ( dynamics%vars(n)%q(1,1,1,k),dynamics%vars(n)%q(1,1,1,k),
     .            dynamics%grid%dlam,dynamics%grid%dphi, im,jm,lm,0,dynamics%grid%lattice )
      enddo

C **********************************************************************        
C ****                  De-Allocate Workspace                       ****        
C **********************************************************************        
                                                                                
      call destroy_dynamics_vars ( dynamics_vars_np1 )

      deallocate ( dpdt )
      deallocate ( dudt )
      deallocate ( dvdt )
      deallocate ( dtdt )
      deallocate ( dqdt )

      RETURN                                                                    
      END                                                                       

      subroutine g3_dyndrv( prev,curr,grid,omega,
     .                      dpdt,dudt,dvdt,dtdt,dqdt,
     .                      phis,alpha,delt,brown_campana )

C***********************************************************************        
C  PURPOSE                                                                      
C     Driver for the GEOS GCM Dynamics
C                                                                               
C     prev ..... Dynamics State Variable Data Structure (Previous Timestep)     
C     curr ..... Dynamics State Variable Data Structure (Current  Timestep)     
C     grid ..... Dynamics Grid           Data Structure
C     omega..... Dynamics Omega Coupling
C     pke ...... Dynamics PKE   Coupling (Current Timestep)
C     dpdt ..... Dynamics Tendency
C     dudt ..... Dynamics Tendency
C     dvdt ..... Dynamics Tendency
C     dtdt ..... Dynamics Tendency
C     dqdt ..... Dynamics Tendency
C     ALPHA .... Time Filter Coefficient                                        
C     DELT ..... Time-Step in seconds                                           
C
C     BROWN_CAMPANA ... Logical for Brown-Campana Split Explicit Scheme            
C                                                                               
C***********************************************************************        
C*                  GODDARD LABORATORY FOR ATMOSPHERES                 *        
C***********************************************************************        

c Declare Modules and Data Structures
c -----------------------------------
      use   g3_dynamics_state_module
      implicit none
      type   ( dynamics_vars_type ) prev   ! Previous Timestep
      type   ( dynamics_vars_type ) curr   ! Current  Timestep
      type   ( dynamics_grid_type ) grid

c Input Variables
c ---------------
      real(kind=8)  alpha,delt
      logical       brown_campana

      real(kind=8)  dpdt(grid%im,grid%jm)                      ! Dynamics Tendency
      real(kind=8)  dudt(grid%im,grid%jm,grid%lm)              ! Dynamics Tendency
      real(kind=8)  dvdt(grid%im,grid%jm,grid%lm)              ! Dynamics Tendency
      real(kind=8)  dtdt(grid%im,grid%jm,grid%lm)              ! Dynamics Tendency
      real(kind=8)  dqdt(grid%im,grid%jm,grid%lm,grid%ntracer) ! Dynamics Tendency

      real(kind=8)     phis(grid%im,grid%jm)                   ! Earth    Coupling
      real(kind=8)    omega(grid%im,grid%jm,grid%lm)           ! Dynamics Coupling

c Local Copy of Dynamics Grid
c ---------------------------
      real(kind=8) ptop, lam_np,phi_np,lam_0

c Local Variables
c ---------------
      integer     i,j,l,m
      integer     im,jm,lm,nq
      real(kind=8)        rotation
      parameter ( rotation = 86164.09   )

      real(kind=8)        getcon
      real(kind=8)        timestep, akap, pi, om
      real(kind=8)        cp, ae, rgas
      real(kind=8), allocatable ::      pke(:,:,:)
      real(kind=8), allocatable ::  psigdot(:,:,:)
      real(kind=8), allocatable ::     vort(:,:,:)

C *********************************************************************
C *****            Allocate Temporary Work Space                   ****
C *********************************************************************

      im = grid%im
      jm = grid%jm
      lm = grid%lm
      nq = grid%ntracer

      allocate (     pke(im,jm,lm+1) )
      allocate ( psigdot(im,jm,lm)   )
      allocate (    vort(im,jm,lm)   )

      ptop   = grid%ptop
      lam_np = grid%lam_np
      phi_np = grid%phi_np
      lam_0  = grid%lam_0

C **********************************************************************        
C ****                   Begin Executable Code                      ****        
C **********************************************************************        
                                                                                
                          timestep = 0.0
      if( brown_campana ) timestep = delt

      akap    = getcon('KAPPA')
      pi      = 4.*atan(1.)
      om      = 2.*pi/rotation
      cp      = getcon('CP')                                             
      ae      = getcon('EARTH RADIUS')                                             
      rgas    = getcon('RGAS')                                             

      call getpke ( curr%p,pke,grid,im,jm )

      call dycore ( im,jm,lm,grid%sige,ptop,nq,timestep,
     .              om,cp,rgas,ae,lam_np,phi_np,grid%dlam,grid%dphi,
     .              phis,pke,
     .              curr%p, curr%u, curr%v, curr%t, curr%q,
     .              prev%p, prev%u, prev%v, prev%t, prev%q,
     .              dpdt,   dudt,   dvdt,   dtdt,   dqdt, 
     .              omega,  vort,   psigdot,alpha,  grid%lattice )

C **********************************************************************
C ****                  De-Allocate Workspace                       ****
C **********************************************************************

      deallocate (     pke )
      deallocate ( psigdot )
      deallocate (    vort )

      return
      end

      subroutine g3_step ( prev,curr,updt,grid,
     .                     dpdt,dudt,dvdt,dtdt,dqdt,
     .                     alpha,delt )
C***********************************************************************
C  PURPOSE                                                              
C     Updates prognostic fields one time-step
C                                                                       
C***********************************************************************
C*                  GODDARD LABORATORY FOR ATMOSPHERES                 *
C***********************************************************************

c Declare Modules and Data Structures
c -----------------------------------
      use   g3_dynamics_state_module
      implicit none
      type   ( dynamics_vars_type ) prev   ! Previous  Timestep
      type   ( dynamics_vars_type ) curr   ! Current   Timestep
      type   ( dynamics_vars_type ) updt   ! Updated   Timestep
      type   ( dynamics_grid_type ) grid   ! Dynamics  Grid

c Input Variables
c ---------------
      real(kind=8) alpha,delt

      real(kind=8)    dpdt(grid%im,grid%jm)                      ! Dynamics Tendency
      real(kind=8)    dudt(grid%im,grid%jm,grid%lm)              ! Dynamics Tendency
      real(kind=8)    dvdt(grid%im,grid%jm,grid%lm)              ! Dynamics Tendency
      real(kind=8)    dtdt(grid%im,grid%jm,grid%lm)              ! Dynamics Tendency
      real(kind=8)    dqdt(grid%im,grid%jm,grid%lm,grid%ntracer) ! Dynamics Tendency

c Local Variables
c ---------------
      real(kind=8)    pinv (grid%im,grid%jm)
      real(kind=8)     ple (grid%im,grid%jm,grid%lm+1)

      integer  i,j,L,m,im,jm,lm,nq
      real(kind=8)     ptmp1, ptmp2
      real(kind=8)     akap, getcon

      akap  = getcon('KAPPA')

      im = grid%im
      jm = grid%jm
      lm = grid%lm
      nq = grid%ntracer

      call ctoa ( curr%p,curr%p,grid%dlam,grid%dphi,im,jm,1 ,0,grid%lattice )

C **********************************************************************        
C ****       Compute Total Time Tendencies and Update Fields        ****        
C **********************************************************************        
                                                                                
      do j=1,jm
      do i=1,im
      updt%p(i,j) = prev%p(i,j) + delt*dpdt(i,j)
        pinv(i,j) = 1.0 / updt%p(i,j)                                      
      enddo
      enddo

c Mass, Momentum, and Thermodynamic Fields
c ----------------------------------------
      do L=1,lm
      do j=1,jm
      do i=1,im
      updt%u(i,j,L) =             prev%u(i,j,L) + delt*dudt(i,j,L)            
      updt%v(i,j,L) =             prev%v(i,j,L) + delt*dvdt(i,j,L)            
      updt%t(i,j,L) = prev%p(i,j)*prev%t(i,j,L) + delt*dtdt(i,j,L)            
      enddo
      enddo
      enddo

c Active and Passive Tracers
c --------------------------
      do m=1,nq
      do L=1,lm
      do j=1,jm
      do i=1,im
      updt%q(i,j,L,m) = prev%p(i,j)*prev%q(i,j,L,m) + delt*dqdt(i,j,L,m)
      enddo
      enddo
      enddo
      enddo

c Fill Negative Values of Specific Humidity
c -----------------------------------------
      call pqcheck ( updt%q,updt%p,grid,im,jm,lm,delt )

C **********************************************************************
C ****                 Apply Asselin Time Filter                    ****
C **********************************************************************
                                                                        
      if( alpha.ne.0.0 ) then
      call tmfilt ( prev%p,prev%u,prev%v,prev%t,prev%q,
     .              curr%p,curr%u,curr%v,curr%t,curr%q,
     .              updt%p,updt%u,updt%v,updt%t,updt%q,
     .              im,jm,lm,nq,alpha )
      endif

C **********************************************************************
C ****                    Unscale Variables                         ****
C **********************************************************************
                                                                        
      do L=1,lm                                                       
         do j=1,jm
         do i=1,im
         updt%t(i,j,L)   = updt%t(i,j,L)   * pinv(i,j)                             
         enddo
         enddo
         do m=1,nq
         do j=1,jm
         do i=1,im
         updt%q(i,j,L,m) = updt%q(i,j,L,m) * pinv(i,j)            
         enddo
         enddo
         enddo
      enddo
                                                                        
      return
      end

      subroutine tmfilt ( pnm1,unm1,vnm1,tnm1,qnm1,
     .                    pn  ,un  ,vn  ,tn  ,qn  ,
     .                    pnp1,unp1,vnp1,tnp1,qnp1,
     .                    im,jm,lm,nq,alpha )

C***********************************************************************
C  PURPOSE                                                              
C     Perform asselin time filter                                       
C  ARGUMENTS                                                            
C     ALPHA .. Time filter coefficient                                  
C                                                                       
C  Note:  tnp1 and qnp1 variables are mass-weighted                                         
C                                                                       
C***********************************************************************
C*                  GODDARD LABORATORY FOR ATMOSPHERES                 *
C***********************************************************************
                                                                        
      implicit none

c Input Variables
c ---------------
      real(kind=8)    alpha
      integer im,jm,lm,nq     ! Dynamics Grid

      real(kind=8)  pnm1(im,jm)       ! Dynamics State at time-level n-1
      real(kind=8)  unm1(im,jm,lm)    ! Dynamics State at time-level n-1
      real(kind=8)  vnm1(im,jm,lm)    ! Dynamics State at time-level n-1
      real(kind=8)  tnm1(im,jm,lm)    ! Dynamics State at time-level n-1
      real(kind=8)  qnm1(im,jm,lm,nq) ! Dynamics State at time-level n-1

      real(kind=8)  pn  (im,jm)       ! Dynamics State at time-level n
      real(kind=8)  un  (im,jm,lm)    ! Dynamics State at time-level n
      real(kind=8)  vn  (im,jm,lm)    ! Dynamics State at time-level n
      real(kind=8)  tn  (im,jm,lm)    ! Dynamics State at time-level n
      real(kind=8)  qn  (im,jm,lm,nq) ! Dynamics State at time-level n

      real(kind=8)  pnp1(im,jm)       ! Dynamics State at time-level n+1
      real(kind=8)  unp1(im,jm,lm)    ! Dynamics State at time-level n+1
      real(kind=8)  vnp1(im,jm,lm)    ! Dynamics State at time-level n+1
      real(kind=8)  tnp1(im,jm,lm)    ! Dynamics State at time-level n+1
      real(kind=8)  qnp1(im,jm,lm,nq) ! Dynamics State at time-level n+1

                                                                        
c Local Variables
c ---------------
      integer i,j,L,m
      real(kind=8)    one,ahalf, a1,a2
     
      PARAMETER ( ONE    = 1.0 )                                        
      PARAMETER ( AHALF  = 0.5 )                                        
                                                                        
c Local Dynamic Space
c -------------------
      real(kind=8)   pinv(im,jm)                                        
                                                                        
      A1   = ONE   - ALPHA                                              
      A2   = AHALF * ALPHA                                              
                                                                        
C **********************************************************************
C ****           PI-WEIGHT TIME-LEVEL N FOR MASS FIELDS             ****
C **********************************************************************
                                                                        
      DO 10 L=1,lm                                                   
                                                                        
      do j = 1,jm
      do i = 1,im
      tn(i,j,L) = pn(i,j)*tn(i,j,L)                          
      enddo
      enddo
                                                                        
      do m=1,nq
      do j = 1,jm
      do i = 1,im
      qn(i,j,L,m) = pn(i,j)*qn(i,j,L,m)                          
      enddo
      enddo
      enddo
                                                                        
C **********************************************************************
C ****                 APPLY ASSELIN TIME FILTER                    ****
C **********************************************************************
                                                                        
      do j = 1,jm
      do i = 1,im
      un(i,j,L) = A1*un(i,j,L) + A2*( unm1(i,j,L) + unp1(i,j,L) )
      vn(i,j,L) = A1*vn(i,j,L) + A2*( vnm1(i,j,L) + vnp1(i,j,L) )
      tn(i,j,L) = A1*tn(i,j,L) + A2*( tnm1(i,j,L) * pnm1(i,j) + tnp1(i,j,L) )
      enddo
      enddo

      do m=1,nq
      do j = 1,jm
      do i = 1,im
      qn(i,j,L,m) = A1*qn(i,j,L,m) + A2*( qnm1(i,j,L,m)*pnm1(i,j) + qnp1(i,j,L,m) )
      enddo
      enddo
      enddo
                                                                        
   10 CONTINUE                                                          
                                                                        
C **********************************************************************
C ****             APPLY ASSELIN TIME FILTER ON PRESSURE            ****
C ****                  UNSCALE TIME-LEVEL N2                       ****
C **********************************************************************
                                                                        
      do j = 1,jm
      do i = 1,im
      pn  (i,j) = A1*pn(i,j) + A2*( pnm1(i,j) + pnp1(i,j) )           
      pinv(i,j) = ONE / pn(i,j)
      enddo
      enddo

      DO 20 L=1,lm                                                  
      do j = 1,jm
      do i = 1,im
      tn(i,j,L) = tn(i,j,L)*pinv(i,j)
      enddo
      enddo
      do m=1,nq
      do j = 1,jm
      do i = 1,im
      qn(i,j,L,m) = qn(i,j,L,m)*pinv(i,j)
      enddo
      enddo
      enddo
   20 CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               

      subroutine pqcheck ( PQZ,PZ,grid,IM,JM,LM,delt )
C***********************************************************************
C  Purpose                                                              
C     Check Specific Humidity Field for Negative values                 
C
C  Argument Description
C     PQZ ........ (ps-ptop)*Specific Humidity (mb g/g)                                  
C     PZ ......... Pi = ps-ptop (mb)                                  
C     grid ....... Dynamics Grid Structure                                  
C     IM ......... Zonal      Dimension
C     JM ......... Meridional Dimension
C     LM ......... Vertical   Dimension
C     DELT ....... Timestep  (Seconds)
C
C  Note:  The Diagnostics are divided by 2 because:
C         a)  The Timestep = 2*DT for Leapfrog, or 
C         b)  The Routine is called twice for Matsuno
C                                                                       
C***********************************************************************
C*                  GODDARD LABORATORY FOR ATMOSPHERES                 *
C***********************************************************************
                                                                        
      use g3_dynamics_state_module
      implicit none
      type ( dynamics_grid_type ) grid

      integer im,jm,lm,nd
      real(kind=8) delt

      real(kind=8)     PQZ(IM,JM,LM)
      real(kind=8)      PZ(IM,JM)
      real(kind=8)     ple(im,jm,lm+1)

      integer i,j,L
      real(kind=8)    getcon,grav,fact
      real(kind=8)    tmp(im,jm)
      real(kind=8), allocatable :: q(:,:,:)

      grav = getcon('GRAVITY')
      fact = 86400/delt/2
 
      call getple ( pz,ple,grid,im,jm )

      allocate ( q(im,jm,lm) )
      do L=1,lm
      do j=1,jm
      do i=1,im
      q(i,j,L) = pqz(i,j,L)/pz(i,j)
      enddo
      enddo
      enddo
      
c Fill Negative Specific Humidities
c ---------------------------------
      do L=2,lm
      do j=1,jm
      do i=1,im
      if( q(i,j,L-1).lt.0.0 ) THEN                                
          q(i,j,L  ) = q(i,j,L) + q(i,j,L-1)*( ple(i,j,L  )-ple(i,j,L-1) )
     .                                      /( ple(i,j,L+1)-ple(i,j,L  ) )
          q(i,j,L-1) = 0.0
      endif
      enddo
      enddo
      enddo

      do j=1,jm
      do i=1,im
      if( q(i,j,lm).lt.0.0 ) q(i,j,lm) = 0.0
      enddo
      enddo

      do L=1,lm
      do j=1,jm
      do i=1,im
      pqz(i,j,L) = pz(i,j)*q(i,j,L)
      enddo
      enddo
      enddo
      
      deallocate ( q )
      return                                                            
      end                                                               

      subroutine shapij ( PZ,UZ,VZ,TZ,QZ,DUDT,DVDT,DTDT,dqdt,
     .                    im,jm,lm,nq,grid )
C***********************************************************************
C  PURPOSE
C     Used to globally filter U-wind, V-wind, Temperature and Moisture
C     to damp small scale dispersive waves.
C
C  ARGUMENTS   DESCRIPTION
C     PZ ..... SURFACE PRESSURE FIELD (PS-PTOP)
C     UZ ..... U-WIND
C     VZ ..... V-WIND
C     TZ ..... POTENTIAL TEMPERATURE
C     QZ ..... Tracers             
C     DUDT ... U-Wind Tendency
C     DVDT ... V-Wind Tendency
C     DTDT ... Mass weighted Theta  Tendency
C     dqdt ... Mass weighted Tracer Tendencies
C     im ..... X-DIMENSION OF FIELDS
C     jm ..... Y-DIMENSION OF FIELDS
C     lm ..... Z-DIMENSION OF FIELDS
C     nq ..... Number of Tracers
C     GRID ... Dynamics Grid Structure
C
C  REMARKS:
C     1) Full Shapiro Filter is assumed to have a time-scale of TAU
C     2) SHAPT Diagnostic is approximate in Leapfrog Mode since PK 
C        is defined for Timestep N (or N*), while SHAPT is defined for 
C        Timestep NM1 
C
C***********************************************************************
C*                  GODDARD LABORATORY FOR ATMOSPHERES                 *
C***********************************************************************
 
      use g3_dynamics_state_module
      implicit none
      type ( dynamics_grid_type ) grid

      integer im,jm,lm,nd
      real(kind=8)    tau,df
      integer i,j,L,m,nq

      PARAMETER ( TAU   = 5400.0  )
      PARAMETER (  DF   = 1.0/TAU )
 
      real(kind=8) PZ(im,jm)
      real(kind=8) UZ(im,jm,lm)
      real(kind=8) VZ(im,jm,lm)
      real(kind=8) TZ(im,jm,lm)
      real(kind=8) QZ(im,jm,lm,nq)
 
      real(kind=8) DUDT(im,jm,lm)
      real(kind=8) DVDT(im,jm,lm)
      real(kind=8) DTDT(im,jm,lm)
      real(kind=8) dqdt(im,jm,lm,nq)
 
c Local Dynamic Space
c -------------------
      real(kind=8)  SHAPU (im,jm)
      real(kind=8)  SHAPV (im,jm)
      real(kind=8)  SHAPT (im,jm)
      real(kind=8)  SHAPQ (im,jm)
      real(kind=8)   pinv (im,jm)
      real(kind=8)   pk   (im,jm,lm)
      real(kind=8)  pref  (lm)
 
      logical first
      data    first /.true./

      real(kind=8)     dlam_min,pi
      integer, allocatable, save :: nshap(:)

c Determine Order of Shapiro Filter
c ---------------------------------
      if(first) then

      allocate ( nshap(lm) )

            pi = 4.0*atan(1.0)
      dlam_min = 1.0e15
      do i=1,grid%lattice%imglobal
      if( grid%dlam(i).lt.dlam_min ) dlam_min = grid%dlam(i)
      enddo
      dlam_min = dlam_min*180./pi

#if 0
      if( grid%lattice%myid.eq.0 ) then
      print *
      print *, 'DLAM_min = ',dlam_min
      print *
      endif
#endif

      call getpl ( 100000.0-grid%ptop,pref,grid,1,1 )
      pref = pref*0.01
      do L=1,lm
 
c Zonal Dimension >= 5 Degrees
c ----------------------------
      if( dlam_min.ge.4.9 ) then
         if(pref(L).ge.7.) nshap(L) = 16
         if(pref(L).lt.7. .and.
     .      pref(L).ge.1.) nshap(L) =  8
         if(pref(L).lt.1.) nshap(L) =  4
      endif

c Zonal Dimension >= 2.5 Degrees and < 5 Degrees
c ----------------------------------------------
      if( dlam_min.ge.2.4 .and. dlam_min.lt.4.9 ) then
         if(pref(L).ge.7.) nshap(L) =  8
         if(pref(L).lt.7. .and.
     .      pref(L).ge.1.) nshap(L) =  8
         if(pref(L).lt.1.) nshap(L) =  4
      endif
 
c Zonal Dimension < 2.5 Degrees
c -----------------------------
      if( dlam_min.lt.2.4 ) then
          nshap(L) =  4
      endif
 
#if 0
      if( grid%lattice%myid.eq.0 )  write(6,6001)  L,pref(L),nshap(L)
 6001 format(' SHAPIJ: Level ',i3,' Pref: ',f8.3,'  using Shapiro filter of order',i3)
#endif
      enddo
c     if( grid%lattice%myid.eq.0 )  print *
      first = .false.
      endif
 
C **********************************************************************
C ****        Prepare Pole and Near-Pole Values for Shapiro         ****
C **********************************************************************
 
      call ctoa ( tz,tz,grid%dlam,grid%dphi,im,jm,lm,0,grid%lattice )
      do m = 1,nq
      call ctoa ( qz(1,1,1,m),qz(1,1,1,m),grid%dlam,grid%dphi,im,jm,lm,0,grid%lattice )
      enddo
      do L=1,lm
      call polewnd ( uz(1,1,L),vz(1,1,L),
     .               uz(1,1,L),vz(1,1,L),grid%dlam,grid%dphi,im,jm,grid%lattice )
      enddo

C **********************************************************************
C ****        FRACTIONAL FILTER ON FIELDS U, V, T AND Q             ****
C **********************************************************************
 
      do 1000 L=1,lm
      call shap ( uz(1,1,L),shapu,im,jm,1,-1,nshap(L),grid%lattice )
      call shap ( vz(1,1,L),shapv,im,jm,1,-2,nshap(L),grid%lattice )
      call shap ( tz(1,1,L),shapt,im,jm,1, 1,nshap(L),grid%lattice )
 
C **********************************************************************
C ****             MASS-WEIGHT SHAPIRO ADJUSTMENTS                  ****
C **********************************************************************
 
      do j=1,jm
      do I=1,im
      SHAPU(I,j) =           SHAPU(I,j) * DF
      SHAPV(I,j) =           SHAPV(I,j) * DF
      SHAPT(I,j) = PZ(I,j) * SHAPT(I,j) * DF
      enddo
      enddo
 
C **********************************************************************
C ****              UPDATE MASS-WEIGHTED TENDENCIES                 ****
C **********************************************************************
 
      do j=1,jm
      do I=1,im
      DUDT(I,j,L) = DUDT(I,j,L) - SHAPU(I,j)
      DVDT(I,j,L) = DVDT(I,j,L) - SHAPV(I,j)
      DTDT(I,j,L) = DTDT(I,j,L) - SHAPT(I,j)
      enddo
      enddo
 
c Do Tracers
c ----------
      do m = 1,nq
      call shap ( qz(1,1,L,m),shapq,im,jm,1,1,nshap(L),grid%lattice )
      do j=1,jm
      do I=1,im
      shapq(i,j)     =   pz(i,j)     * shapq(i,j)*df
       dqdt(i,j,l,m) = dqdt(i,j,l,m) - shapq(i,j)
      enddo
      enddo
      enddo
 
 1000 continue
 
      return
      end

      subroutine shap ( q,dq,im,jm,lm,msgn,n,lattice )
C***********************************************************************
C  Purpose:
C     Apply Shapiro Filter to a General Field
C  Description of Parameters
C     q ....... Input  Field
C     dq ...... Output Increment: q-filter(q)
C     im ...... Zonal      Dimension
C     jm ...... Meridional Dimension
C     lm ...... Vertical   Dimension
C     msgn .... Set msgn = +1 For Scalar Field
C                   msgn = -1 For U_Wind Field
C                   msgn = -2 For V-Wind Field
C     n ....... Shapiro Filter Order
C
C***********************************************************************
C*                  GODDARD LABORATORY FOR ATMOSPHERES                 *
C***********************************************************************
      use g3_dynamics_lattice_module
      implicit none
      type ( dynamics_lattice_type ) lattice

      integer im,jm,lm,msgn,n
      real(kind=8)    q(im,jm,lm), dq(im,jm,lm)

      integer i,j,k,L,m
      real(kind=8)    fac
      real(kind=8), allocatable :: d1(:,:,:)
      real(kind=8), allocatable :: d2(:,:,:)

      m   =   n/2
      fac = (-0.25)**m

C **********************************************************************
C ****             Filter Mass or U-Wind Points                     ****
C **********************************************************************

      if ( msgn.eq.1 .or. msgn.eq.-1 ) then

c Zonal Direction
c ---------------
      allocate( d1(1-m:im+m,jm,lm) )
      allocate( d2(1-m:im+m,jm,lm) )

      d2 = 0.0
      call ghostx ( q,d1,im,jm,lm,m,lattice,'both' )

      do L=1,lm
      do j=1,jm
         do k=1,m
            do i=2-m,im+m
            d2(i,j,L) = d1(i,j,L)-d1(i-1,j,L)
            enddo
            do i=1-m,im+m-1
            d1(i,j,L) = d2(i+1,j,L)-d2(i,j,L)
            enddo
         enddo
         do i=1,im
         dq(i,j,L) = q(i,j,L)-d1(i,j,L)*fac   ! Filtered Field in X
         enddo
      enddo
      enddo

      deallocate ( d1,d2 )

c Meridional Direction
c --------------------
      allocate( d1(im,1-m:jm+m,lm) )
      allocate( d2(im,1-m:jm+m,lm) )

      d2 = 0.0
      call ghosty ( dq,d1,im,jm,lm,0,msgn,m,lattice,'both' )

      do L=1,lm
      do k=1,m
         do j=2-m,jm+m
         do i=1,im
         d2(i,j,L) = d1(i,j,L)-d1(i,j-1,L)
         enddo
         enddo
         do j=1-m,jm+m-1
         do i=1,im
         d1(i,j,L) = d2(i,j+1,L)-d2(i,j,L)
         enddo
         enddo
      enddo
      enddo

      endif

C **********************************************************************
C ****                 Filter V-Wind Points                         ****
C **********************************************************************

      if ( msgn.eq.-2 ) then

c Zonal Direction
c ---------------
      allocate( d1(1-m:im+m,jm,lm) )
      allocate( d2(1-m:im+m,jm,lm) )

      d2 = 0.0
      call ghostx ( q,d1,im,jm,lm,m,lattice,'both' )

      do L=1,lm
      do j=1,jm
         do k=1,m
            do i=2-m,im+m
            d2(i,j,L) = d1(i,j,L)-d1(i-1,j,L)
            enddo
            do i=1-m,im+m-1
            d1(i,j,L) = d2(i+1,j,L)-d2(i,j,L)
            enddo
         enddo
         do i=1,im
         dq(i,j,L) = q(i,j,L)-d1(i,j,L)*fac   ! Filtered Field in X
         enddo
      enddo
      enddo

      deallocate ( d1,d2 )

c Meridional Direction
c --------------------
      allocate( d1(im,1-m:jm+m,lm) )
      allocate( d2(im,1-m:jm+m,lm) )

      d2 = 0.0
      call ghosty ( dq,d1,im,jm,lm,1,-1,m,lattice,'both' )

      do L=1,lm
      do k=1,m
         do j=2-m,jm+m
         do i=1,im
         d2(i,j,L) = d1(i,j,L)-d1(i,j-1,L)
         enddo
         enddo
         do j=1-m,jm+m-1
         do i=1,im
         d1(i,j,L) = d2(i,j+1,L)-d2(i,j,L)
         enddo
         enddo
      enddo
      enddo

      endif

C **********************************************************************
C ****       Compute Full Field Filter and Increment                ****
C **********************************************************************

      do L=1,lm
      do j=1,jm
      do i=1,im
         dq(i,j,L) = dq(i,j,L)-d1(i,j,L)*fac   ! Filtered Field in X and Y
         dq(i,j,L) =  q(i,j,L)-dq(i,j,L)       ! Shapiro  Filter Increment
      enddo
      enddo
      enddo

      deallocate ( d1,d2 )

      return
      end
