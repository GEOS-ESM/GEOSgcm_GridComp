module BuoyancyMod

  use MAPL_BaseMod, only: MAPL_UNDEF
  use MAPL_Constants, only: MAPL_GRAV
  use ProcessLibraryParameters, only: alhlbcp, gravbcp

  implicit none

  private
  public :: BUOYANCY

contains

  subroutine BUOYANCY( T, Q, QS, DQS, DZ, ZLO, BUOY, CAPE, INHB)

    ! !DESCRIPTION: Computes the buoyancy $ g \frac{T_c-T_e}{T_e} $ at each level
    !  for a parcel raised from the surface. $T_c$ is the virtual temperature of
    !  the parcel and $T_e$ is the virtual temperature of the environment.

    real, dimension(:,:,:),   intent(in)  :: T, Q, QS, DQS, DZ, ZLO
    real, dimension(:,:,:),   intent(out) :: BUOY
    real, dimension(:,:),     intent(out) :: CAPE, INHB

    integer :: I, J, L, IM, JM, LM

    LM = size(T,3)
    IM = size(T,1)
    JM = size(T,2)

!$acc parallel loop gang vector collapse(3)
    do L=LM-1,1,-1
       do J = 1,JM
          do I = 1,IM
             BUOY(I,J,L) = (T(I,J,LM) + gravbcp*ZLO(I,J,LM) + alhlbcp*Q(I,J,LM)) - (T(I,J,L) + gravbcp*ZLO(I,J,L) + alhlbcp*QS(I,J,L))
             BUOY(I,J,L) = MAPL_GRAV*BUOY(I,J,L) / ( (1.+ alhlbcp*DQS(I,J,L))*T(I,J,L) )
          enddo
       enddo
    enddo
!$acc end parallel loop

!$acc parallel loop gang vector collapse(2)
    do J = 1,JM
       do I = 1,IM
          BUOY(I,J,LM) = 0.0
          CAPE(I,J) = 0.
          INHB(I,J) = 0.
       enddo
    enddo
!$acc end parallel loop

!$acc parallel loop gang vector collapse(3)
    do L=1,LM-1
       do J = 1, JM
          do I = 1,IM
             if(BUOY(I,J,L)>0.) then
!$acc atomic update
                CAPE(I,J) = CAPE(I,J) + BUOY(I,J,L)*DZ(I,J,L)
!$acc end atomic
             endif

             if (BUOY(I,J,L)<0.) then
!$acc atomic update
                INHB(I,J) = INHB(I,J) - BUOY(I,J,L)*DZ(I,J,L)
!$acc end atomic
             endif
          enddo
       enddo
    end do
!$acc end parallel loop

!$acc parallel loop gang vector collapse(2)
    do J = 1,JM
       do I = 1,IM
          if(CAPE(I,J) <= 0.0) then
             CAPE(I,J) = mapl_undef
             INHB(I,J) = mapl_undef
          endif
       enddo
    enddo
!$acc end parallel loop

  end subroutine BUOYANCY

end module BuoyancyMod
