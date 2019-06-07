      subroutine remap ( ple,u,v,thv,qtr,phis_in,phis_out,ak,bk,im,jm,lm,km )

!***********************************************************************
!
!  Purpose
!     Driver for remapping fields to new topography
!
!  Argument Description
!     ple ...... model edge pressure
!     u  ....... model zonal      wind
!     v  ....... model meridional wind
!     thv  ..... model virtual potential  temperature
!     q  ....... model specific   humidity
!     o3  ...... model ozone
!     phis_in... model surface geopotential (input)
!     phis_out.. model surface geopotential (output)
!     ak ....... model vertical   dimension
!     bk ....... model vertical   dimension
!
!     im ....... zonal      dimension
!     jm ....... meridional dimension
!     lm ....... meridional dimension
!
!***********************************************************************

  use MAPL_Mod
  use dynamics_vars, only : T_TRACERS

      implicit none
      integer  im,jm,lm,km

! Input variables
! ---------------
      type(T_TRACERS) qtr(km)
      real*8      ple(im,jm,lm+1)
      real*8        u(im,jm,lm)
      real*8        v(im,jm,lm)
      real*8      thv(im,jm,lm)
      real*8        q(im,jm,lm)
      real*8       o3(im,jm,lm)
      real*8 phis_in (im,jm)
      real*8 phis_out(im,jm)

      real*8       ak(lm+1)
      real*8       bk(lm+1)

! Local variables
! ---------------
      real*8, allocatable ::  ps     (:,:)
      real*8, allocatable ::  phi    (:,:,:)
      real*8, allocatable ::  pke    (:,:,:)
      real*8, allocatable ::  ple_out(:,:,:)
      real*8, allocatable ::  pke_out(:,:,:)

      real*8, allocatable ::     delp(:,:,:)
      real*8, allocatable ::    u_out(:,:,:)
      real*8, allocatable ::    v_out(:,:,:)
      real*8, allocatable ::  thv_out(:,:,:)
      real*8, allocatable ::    q_in (:,:,:,:)
      real*8, allocatable ::    q_out(:,:,:,:)

      real*8    kappa,cp,rgas,eps,rvap
      integer i,j,L,n,k

      kappa = MAPL_KAPPA
      rgas  = MAPL_RGAS
      rvap  = MAPL_RVAP
      eps   = rvap/rgas-1.0
      cp    = rgas/kappa

      allocate(  ps     (im,jm)      )
      allocate(  phi    (im,jm,lm+1) )
      allocate(  pke    (im,jm,lm+1) )
      allocate(  ple_out(im,jm,lm+1) )
      allocate(  pke_out(im,jm,lm+1) )

      allocate(     delp(im,jm,lm)   )
      allocate(    u_out(im,jm,lm)   )
      allocate(    v_out(im,jm,lm)   )
      allocate(  thv_out(im,jm,lm)   )
      allocate(    q_in (im,jm,lm,km))
      allocate(    q_out(im,jm,lm,km))

! Construct Input Heights
! -----------------------
      pke(:,:,:) = ple(:,:,:)**kappa 

      phi(:,:,lm+1) = phis_in(:,:)
      do L=lm,1,-1
      phi(:,:,L) = phi(:,:,L+1) + cp*thv(:,:,L)*( pke(:,:,L+1)-pke(:,:,L) )
      enddo
      
! Compute new surface pressure consistent with output topography
! --------------------------------------------------------------
      do j=1,jm
      do i=1,im
           L = lm
           do while ( phi(i,j,L).lt.phis_out(i,j) )
           L = L-1
           enddo
           ps(i,j) = ple(i,j,L+1)*( 1+(phi(i,j,L+1)-phis_out(i,j))/(cp*thv(i,j,L)*pke(i,j,L+1)) )**(1.0/kappa)
      enddo
      enddo

! Construct fv pressure variables using new surface pressure
! ----------------------------------------------------------
      do L=1,lm+1
      do j=1,jm
      do i=1,im
       ple_out(i,j,L) = ak(L) + bk(L)*ps(i,j)
      enddo
      enddo
      enddo
      pke_out(:,:,:) = ple_out(:,:,:)**kappa 

! Map original fv state onto new eta grid
! ---------------------------------------

      do k=1,size(qtr)
         if(qtr(k)%is_r4) then
            q_in(:,:,:,k) = qtr(k)%content_r4
         else
            q_in(:,:,:,k) = qtr(k)%content
         end if
      enddo

      call gmap ( im,jm,km, kappa,                               &
                  lm, pke    ,ple    ,u    ,v    ,thv    ,q_in , &
                  lm, pke_out,ple_out,u_out,v_out,thv_out,q_out)

      do k=1,size(qtr)
         if(qtr(k)%is_r4) then
            qtr(k)%content_r4 = q_out(:,:,:,k)
         else
            qtr(k)%content    = q_out(:,:,:,k)
         end if
      enddo

      ple(:,:,:) = ple_out(:,:,:)
        u(:,:,:) =   u_out(:,:,:)
        v(:,:,:) =   v_out(:,:,:)
      thv(:,:,:) = thv_out(:,:,:)

      deallocate(  ps      )
      deallocate(  phi     )
      deallocate(  pke     )
      deallocate(  ple_out )
      deallocate(  pke_out )

      deallocate(     delp )
      deallocate(    u_out )
      deallocate(    v_out )
      deallocate(  thv_out )
      deallocate(    q_in  )
      deallocate(    q_out )

      return
      end
