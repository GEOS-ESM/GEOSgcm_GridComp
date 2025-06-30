!-----------------------------------------------------------------------
!BOP
! !ROUTINE: par_xsum --- Calculate x-sum bit-wise consistently
!
! !INTERFACE:
!****6***0*********0*********0*********0*********0*********0**********72
      subroutine par_xsum(grid, a, ltot, xsum)
!****6***0*********0*********0*********0*********0*********0**********72
!
! !USES:
#if defined ( SPMD )
      use parutilitiesmodule, only : parcollective, SUMOP
#endif
      use dynamics_vars, only : T_FVDYCORE_GRID
      use shr_kind_mod, only: r8 => shr_kind_r8

      implicit none

! !INPUT PARAMETERS:
      type (T_FVDYCORE_GRID), intent(in) :: grid
      integer, intent(in) :: ltot       ! number of quantities to be summed
      real (r8) a(grid%ifirstxy:grid%ilastxy,ltot)    ! input vector to be summed

! !OUTPUT PARAMETERS:
      real (r8) xsum(ltot)              ! sum of all vector entries

! !DESCRIPTION:
!     This subroutine calculates the sum of "a" in a reproducible
!     (sequentialized) fashion which should give bit-wise identical
!     results irrespective of the number of MPI processes.
!
! !REVISION HISTORY:
!
!     AAM 00.11.01 : Created
!     WS  03.10.22 : pmgrid removed (now spmd_dyn)
!     WS  04.10.04 : added grid as an argument; removed spmd_dyn
!     WS  05.05.25 : removed ifirst, ilast, im as arguments (in grid)
!     WS  06.12.24 : rewritten to use collective communication call
!
!EOP
!---------------------------------------------------------------------
!BOC
 
! !Local
      real(r8), parameter ::  D0_0                    =  0.0_r8
      real (r8) quan_all(grid%im,ltot)
      integer i,l,ifirst,ilast,im

      quan_all = D0_0  
      ifirst = grid%ifirstxy
      ilast  = grid%ilastxy
      im     = grid%im
      do i=ifirst,ilast
        do l=1,ltot
          quan_all(i,l)=a(i,l)
        enddo
      enddo

#if defined ( SPMD )
      if ( grid%nprxy_x > 1 ) then
        call parcollective( grid%commxy_x, SUMOP, im, ltot, quan_all )
      endif
#endif

      do l=1,ltot
        xsum(l) = D0_0
        do i=1,im
          xsum(l) = xsum(l) + quan_all(i,l)
        enddo
      enddo

      return
!EOC
      end subroutine par_xsum
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: par_xsum_r4 --- Calculate x-sum bit-wise consistently (real4)
!
! !INTERFACE:
!****6***0*********0*********0*********0*********0*********0**********72
      subroutine par_xsum_r4(grid, a, ltot, sum)
!****6***0*********0*********0*********0*********0*********0**********72
!
! !USES:
#if defined ( SPMD )
      use parutilitiesmodule, only : parexchangevector
#endif
      use dynamics_vars, only : T_FVDYCORE_GRID
      use shr_kind_mod, only: r8 => shr_kind_r8, r4 => shr_kind_r4

      implicit none

! !INPUT PARAMETERS:
      type (T_FVDYCORE_GRID), intent(in) :: grid
      integer, intent(in) :: ltot       ! number of quantities to be summed
      real (r4) a(grid%ifirstxy:grid%ilastxy,ltot)    ! input vector to be summed

! !OUTPUT PARAMETERS:
      real (r8) sum(ltot)               ! sum of all vector entries

! !DESCRIPTION:
!     This subroutine calculates the sum of "a" in a reproducible
!     (sequentialized) fashion which should give bit-wise identical
!     results irrespective of the number of MPI processes.
!
! !REVISION HISTORY:
!
!     WS  05.04.08 : Created from par_xsum
!     WS  05.05.25 : removed ifirst, ilast, im as arguments (in grid)
!     WS  06.06.28 : Fixed bug in sequential version
!
!EOP
!---------------------------------------------------------------------
!BOC
 
! !Local
      real(r4), parameter ::  E0_0                    =  0.0_r4
      integer :: nprxy_x, commxy_x
      real (r4) quan_all(grid%im*ltot)
      integer i,l,icount,ipe, ifirst,ilast,im

#if defined ( SPMD )
      real (r4) quan_send(grid%nprxy_x*ltot*(grid%ilastxy-grid%ifirstxy+1))
      integer  sendcount(grid%nprxy_x)
      integer  recvcount(grid%nprxy_x)
#endif

      ifirst = grid%ifirstxy
      ilast  = grid%ilastxy
      im     = grid%im
#if defined ( SPMD )
      nprxy_x = grid%nprxy_x
      commxy_x = grid%commxy_x

      icount=0
      do ipe=1,nprxy_x
        sendcount(ipe) = ltot*(ilast-ifirst+1)
        do i=ifirst,ilast
          do l=1,ltot
            icount=icount+1    
            quan_send(icount)=a(i,l)
          enddo
        enddo
      enddo
      call parexchangevector( commxy_x, sendcount, quan_send,     &
                                          recvcount, quan_all )
#else
      do l=1,ltot
        do i=1,im
          quan_all((i-1)*ltot+l)=a(i,l)
        enddo
      enddo
#endif

      do l=1,ltot
        sum(l) = E0_0
        do i=1,im
          sum(l) = sum(l) + quan_all((i-1)*ltot+l)
        enddo
      enddo

      return
!EOC
      end subroutine par_xsum_r4
!-----------------------------------------------------------------------
