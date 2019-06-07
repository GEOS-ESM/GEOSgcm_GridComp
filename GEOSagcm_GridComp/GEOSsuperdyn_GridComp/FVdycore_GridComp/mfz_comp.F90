!-----------------------------------------------------------------------
!BOP
! !ROUTINE:  mfz_comp --- Calculate vertical mass flux
!
! !INTERFACE:
      subroutine mfz_comp( grid, ae, grav, dt, mfx, mfy, mfz, delp0xy, delpxy )
! !USES:
      use shr_kind_mod, only : r8 => shr_kind_r8
      use dynamics_vars, only : T_FVDYCORE_GRID
#if defined( SPMD )
      use parutilitiesmodule, only: bcstop, sumop,  parcollective
      use mod_comm, only :  commglobal, mp_exit, gid, mp_send3d, mp_recv3d
#endif
      implicit none

#if defined(BILL_DEBUG)
#include "mpif.h"
#endif

! !INPUT PARAMETERS:
      type (T_FVDYCORE_GRID), intent(in) :: grid   ! grid (for XY decomp)

      real(r8), intent(in) :: ae                ! Radius of the Earth (m)
      real(r8), intent(in) :: grav              ! Gravity
      integer,  intent(in) :: dt                ! large time step in seconds
      real(r8), intent(inout) :: mfx(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km) 
      real(r8), intent(inout) :: mfy(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km) 

      real(r8), intent(in) :: delpxy(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km)
      real(r8), intent(in) :: delp0xy(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km)

! !OUTPUT PARAMETERS:
      real(r8), intent(inout) :: mfz(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km+1)

! !DESCRIPTION:
!     Compute cell centered vertical mass flux
!        Vertical integration of the continuity equation using the 
!        divergence of the horizontal mass fluxes
!
! !REVISION HISTORY:
!   WMP  05.12.01   Created
!   WS   09.04.01 : Upgraded to PILGRIM from cam3_6_33
!
! !BUGS:
!   Not yet tested...
!
!EOP
!---------------------------------------------------------------------
!BOC

      real(r8) mfy_north(grid%ifirstxy:grid%ilastxy,grid%km)
      real(r8), allocatable :: mfx_east(:,:)     ! East halo
      real(r8), allocatable :: mfy_sp(:,:), mfy_np(:,:)
      real(r8) delp1xy(grid%im,grid%jm,grid%km)
      real(r8) conv(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km)
      real(r8) pit(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy) ! pressure tendency
      real(r8) area, sum1, sum1a, sum1k, sum2, sum2a, sum2k
      integer  :: im, jm, km
      integer  :: ifirstxy, ilastxy, jfirstxy, jlastxy
      integer  :: ierr, iam, myidxy_y, nprxy_x, nprxy_y, dest, src    ! SPMD related
      integer  :: i,j,k

      im       = grid%im
      jm       = grid%jm
      km       = grid%km

      ifirstxy = grid%ifirstxy
      ilastxy  = grid%ilastxy
      jfirstxy = grid%jfirstxy
      jlastxy  = grid%jlastxy

      iam      = grid%iam
      myidxy_y = grid%myidxy_y
      nprxy_x  = grid%nprxy_x
      nprxy_y  = grid%nprxy_y

      ! Initialize mfy_north
      mfy_north = 0.0_r8

      allocate(mfx_east(jfirstxy:jlastxy,km))     ! East halo
      allocate(mfy_sp(im,km), mfy_np(im,km) )

#if defined( SPMD )
      call mp_send3d( commglobal, iam-nprxy_x, iam+nprxy_x, im, jm, km,    &
                      ifirstxy, ilastxy, jfirstxy, jlastxy, 1, km,         &
                      ifirstxy, ilastxy, jfirstxy, jfirstxy, 1, km, mfy )
      call mp_recv3d( commglobal, iam+nprxy_x, im, jm, km,                 &
                      ifirstxy, ilastxy, jlastxy+1, jlastxy+1, 1, km,      &
                      ifirstxy, ilastxy, jlastxy+1, jlastxy+1, 1, km, mfy_north )
      if (nprxy_x > 1) then
! Nontrivial x decomposition
        dest = myidxy_y*nprxy_x + MOD(iam+nprxy_x-1,nprxy_x)
        src  = myidxy_y*nprxy_x + MOD(iam+1,nprxy_x)
        call mp_send3d( commglobal, dest, src, im, jm, km,                 &
                        ifirstxy, ilastxy, jfirstxy, jlastxy, 1,km,        &
                        ifirstxy, ifirstxy, jfirstxy, jlastxy, 1, km, mfx )
      endif
#endif
!         
! Prepare sum of mfy for divergence at poles
!
      mfy_sp = 0.0_r8
      mfy_np = 0.0_r8
!$omp  parallel do                  &
!$omp  default(shared)              &
!$omp  private(i,k)
      do k=1,km
        if ( jfirstxy == 1 ) then  ! SP
          do i=ifirstxy,ilastxy
            mfy_sp(i,k) = mfy(i,2,k)
          enddo
        endif
        if ( jlastxy == jm ) then  ! NP
          do i=ifirstxy,ilastxy
            mfy_np(i,k) = mfy(i,jm,k)
          enddo
        endif
! Periodic boundary  (for the case of no decomposition in X)
        do j=jfirstxy,jlastxy
           mfx_east(j,k) = mfx(ifirstxy,j,k)
        enddo
      enddo

#if defined( SPMD )
      if ( nprxy_x > 1 ) then
! Non-trivial X decomposition
        call mp_recv3d( commglobal, src, im, jm, km,                             &
                        ilastxy+1, ilastxy+1, jfirstxy, jlastxy, 1, km,          &
                        ilastxy+1, ilastxy+1, jfirstxy, jlastxy, 1, km, mfx_east )
      endif
!                       
! Collect on all PEs the mfy at both poles
!
      if (nprxy_x > 1) then
        call parcollective(grid%commxy_x, sumop, im, km, mfy_sp)
        call parcollective(grid%commxy_x, sumop, im, km, mfy_np)
      endif
#endif

! Compute Convergence of the horizontal Mass flux
!$omp  parallel do                  &
!$omp  default(shared)              &
!$omp  private(i,j,k,sum1,sum2)
      do k=1,km
          do j=jfirstxy,jlastxy-1
             do i=ifirstxy,ilastxy-1
                conv(i,j,k) =  mfx(i,j,k) - mfx(i+1,j,k) +  &
                               mfy(i,j,k) - mfy(i,j+1,k)
             enddo
             conv(ilastxy,j,k) = mfx(ilastxy,j,k) - mfx_east(j,k) +     &
                                 mfy(ilastxy,j,k) - mfy(ilastxy,j+1,k)
          enddo
          j = jlastxy
          do i=ifirstxy,ilastxy-1
             conv(i,j,k) =  mfx(i,j,k) - mfx(i+1,j,k) +  &
                            mfy(i,j,k) - mfy_north(i,k)
          enddo
          conv(ilastxy,j,k) = mfx(ilastxy,j,k) - mfx_east(j,k) +     &
                              mfy(ilastxy,j,k) - mfy_north(ilastxy,k)

! Poles
          if ( jfirstxy == 1  ) then
            sum1 = -SUM(mfy_sp(1:im,k))
            do i=ifirstxy,ilastxy
              conv(i,1,k) = sum1
            enddo
          endif

          if ( jlastxy  == jm ) then
            sum2 =  SUM(mfy_np(1:im,k))
            do i=ifirstxy,ilastxy
              conv(i,jm,k) = sum2
            enddo
          endif
      enddo

! Surface pressure tendency
      pit(:,:) = 0.0
      do k=1,km
         do j=jfirstxy,jlastxy
            do i=ifirstxy,ilastxy
               pit(i,j) = pit(i,j) + conv(i,j,k)
            enddo
         enddo
      enddo

!  Sum over levels
      do k=2,km
         do j=jfirstxy,jlastxy
            do i=ifirstxy,ilastxy
               conv(i,j,k) = conv(i,j,k) + conv(i,j,k-1)
            enddo
         enddo
      enddo

      mfz(:,:,:) = 0.0
      do k=2,km
         do j=MAX(2,jfirstxy),MIN(jlastxy,jm-1)
            area = grid%dl*grid%cosp(j)*ae * grid%dp*ae
            do i=ifirstxy,ilastxy
               mfz(i,j,k) = ( conv(i,j,k-1)  - grid%bk(k)*pit(i,j) )/(grav*area)  ! Kg/m^2/s
            enddo
         enddo
! Poles
          if ( jfirstxy == 1  ) then
            j=1
            area = grid%acap*(grid%dl*ae * grid%dp*ae) 
            do i=ifirstxy,ilastxy
               mfz(i,j,k) = ( conv(i,j,k-1)  - grid%bk(k)*pit(i,j) ) / (grav*area)  ! Kg/m^2/s
            enddo
          endif
          if ( jlastxy  == jm ) then
            j=jm
            area = grid%acap*(grid%dl*ae * grid%dp*ae)                      
            do i=ifirstxy,ilastxy
               mfz(i,j,k) = ( conv(i,j,k-1)  - grid%bk(k)*pit(i,j) ) / (grav*area)  ! Kg/m^2/s
            enddo
          endif
      enddo


#if defined(BILL_DEBUG)
    if (1==1) then
!BMP
    sum1k=0
    sum2k=0
  ! test mass fluxes and d(delp)/dt within each layer
    do k=1,km
      do j=MAX(2,jfirstxy),MIN(jlastxy-1,jm-1)
        area = grid%dl*grid%cosp(j)*ae * grid%dp*ae
        do i=ifirstxy,ilastxy-1

           delp1xy(i,j,k) = dt*( -(mfx(i+1,j,k)-mfx(i,j,k) + mfy(i,j+1,k)-mfy(i,j,k))/area - grav*(mfz(i,j,k+1)-mfz(i,j,k)) )

           if ( (i==100) .and. (j==99) .and. (k==60) ) then
              print*, delpxy(i,j,k)-delp0xy(i,j,k), delp1xy(i,j,k)
           endif

        enddo
        i=ilastxy
        delp1xy(i,j,k) = dt*( -(mfx_east(j,k)-mfx(i,j,k) + mfy(i,j+1,k)-mfy(i,j,k))/area - grav*(mfz(i,j,k+1)-mfz(i,j,k)) )
      enddo
      if (jlastxy /= jm) then
         j=jlastxy
         area = grid%dl*grid%cosp(j)*ae * grid%dp*ae
         do i=ifirstxy,ilastxy-1
            delp1xy(i,j,k) = dt*( -(mfx(i+1,j,k)-mfx(i,j,k) + mfy_north(i,k)-mfy(i,j,k))/area - grav*(mfz(i,j,k+1)-mfz(i,j,k)) )
         enddo
         i=ilastxy
         delp1xy(i,j,k) = dt*( -(mfx_east(j,k)-mfx(i,j,k) + mfy_north(i,k)-mfy(i,j,k))/area - grav*(mfz(i,j,k+1)-mfz(i,j,k)) )
      endif
! Poles
      if ( jfirstxy == 1  ) then
         sum1 = -(SUM(mfy_sp(1:im,k)))*grid%rcap / (grid%dl*ae * grid%dp*ae) 
         j=1
         do i=ifirstxy,ilastxy
            delp1xy(i,j,k) = dt*( sum1 - grav*(mfz(i,j,k+1)-mfz(i,j,k)) )
         enddo
      endif
      if ( jlastxy  == jm ) then
         sum2 =  (SUM(mfy_np(1:im,k)))*grid%rcap / (grid%dl*ae * grid%dp*ae)
         j=jm
         do i=ifirstxy,ilastxy
            delp1xy(i,j,k) = dt*( sum2 - grav*(mfz(i,j,k+1)-mfz(i,j,k)) )
         enddo
      endif

      sum1 = 0
    ! if (jfirstxy == 1) then
    ! do j=1,1
      do j=MAX(1,jfirstxy),MIN(jlastxy,jm)
         do i=ifirstxy,ilastxy
            sum1 = sum1 + delp1xy(i,j,k)
         enddo
      enddo
    ! endif
      call mpi_reduce(sum1, sum1a, 1, MPI_DOUBLE_PRECISION, sumop, 0, grid%comm_y, ierr)
      sum1k=sum1k+sum1a

      sum2 = 0
    ! if (jfirstxy == 1) then
    ! do j=1,1
      do j=MAX(1,jfirstxy),MIN(jlastxy,jm)
         do i=ifirstxy,ilastxy
            sum2 = sum2 + delpxy(i,j,k)-delp0xy(i,j,k)
         enddo
      enddo
    ! endif
      call mpi_reduce(sum2, sum2a, 1, MPI_DOUBLE_PRECISION, sumop, 0, grid%comm_y, ierr)
      sum2k=sum2k+sum2a

      if (gid==0) print*, k, sum1a, sum2a

    enddo

      if (gid==0) print*, sum1k, sum2k

      call mp_exit( commglobal )
      stop

   endif
#endif

      deallocate(mfx_east)
      deallocate(mfy_sp,mfy_np)

      return
      end

