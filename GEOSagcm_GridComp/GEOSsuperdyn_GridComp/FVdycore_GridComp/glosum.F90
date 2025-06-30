!-----------------------------------------------------------------------
!BOP
! !ROUTINE: glosum --- Calculate global sum of an A-Grid Tracer
!
! !INTERFACE:

      subroutine glosum (STATE,NQ,QG)

! !USES:

      use shr_kind_mod,  only : r8 => shr_kind_r8, r4 => shr_kind_r4
      use dynamics_vars, only : T_FVDYCORE_STATE, T_FVDYCORE_GRID
      implicit none

      type(T_FVDYCORE_STATE), intent(IN)  :: STATE
      integer,                intent(IN)  :: NQ
      real(r8),               intent(OUT) :: QG(NQ)

! !DESCRIPTION:
!    Calculate the globally integrated tracer
!
!EOP
!---------------------------------------------------------------------
!BOC

      integer   :: im, jm, km
      integer   :: iam
      integer   :: i1,in,j1,jn,jt,j,k,n

      real (r4), pointer     :: ptr4(:,:,:)
      real (r8), pointer     :: ptr8(:,:,:)
      real (r8), allocatable ::   dp(:,:,:)
      real (r8), allocatable ::   dA(:,:)
      real (r8), allocatable ::   qj(:)
      real (r8), allocatable :: qsum(:,:)
      real (r8), allocatable :: xsum(:)

      real (r8), parameter   :: D0_0 = 0.0_r8

      im  = state%grid%im
      jm  = state%grid%jm
      km  = state%grid%km
      iam = state%grid%iam

      i1  = state%grid%ifirstxy
      in  = state%grid%ilastxy
      j1  = state%grid%jfirstxy
      jn  = state%grid%jlastxy
      jt  = jn - j1 + 1

!-----------------------------------------------------------------------------------------------

      allocate(   dp(i1:in,j1:jn,km) )
      allocate(   dA(i1:in,j1:jn)    )
      allocate( qsum(i1:in,j1:jn)    )
      allocate( xsum(j1:jn)          )
      allocate(   qj(jm)             )

! Compute Grid-Box Area
! ---------------------
      do j=j1,jn
         if ( j == 1  ) then
              dA(:,j) = state%grid%acap    * state%grid%dl*state%grid%dp     ! => 2*pi [1-cos(dp/2)]
         else if ( j == jm ) then
              dA(:,j) = state%grid%acap    * state%grid%dl*state%grid%dp     ! => 2*pi [1-cos(dp/2)]
         else
              dA(:,j) = state%grid%cosp(j) * state%grid%dl*state%grid%dp
         endif
      enddo

! Compute Pressure Thickness
! --------------------------
      do k=1,km
         dp(:,:,k) = state%vars%pe(:,:,k+1)-state%vars%pe(:,:,k)
      enddo

! Loop over Tracers
! -----------------
  do n=1,NQ
         qsum(:,:) = D0_0
     if( STATE%VARS%TRACER(N)%IS_R4 ) then
         do k=1,km
         qsum(:,:) = qsum(:,:) + state%vars%tracer(n)%content_r4(:,:,k)*dp(:,:,k)
         enddo
     else
         do k=1,km
         qsum(:,:) = qsum(:,:) + state%vars%tracer(n)%content   (:,:,k)*dp(:,:,k)
         enddo
     endif
         qsum(:,:) = qsum(:,:) *dA(:,:)

     call par_xsum (state%grid,qsum,jt,xsum)
     do j=j1,jn
         if ( j == 1  ) then
           qj(j) = qsum(i1,j)
         else if ( j == jm ) then
           qj(j) = qsum(i1,j)
         else
           qj(j) = xsum(j)
         endif
     enddo
     call par_vecsum (jm, j1, jn, qj, qg(n), state%grid%commxy_y, state%grid%nprxy_y)

!    if( iam==0 ) print *, 'The global sum for Tracer ',n,' is: ',qg(n)
  enddo

  deallocate( dp )
  deallocate( dA )
  deallocate( qj )
  deallocate( qsum )
  deallocate( xsum )

  return
  end subroutine glosum
