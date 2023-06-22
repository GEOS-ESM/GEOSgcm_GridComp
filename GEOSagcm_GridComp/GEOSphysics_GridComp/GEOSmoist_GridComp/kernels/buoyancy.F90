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

    integer :: L, LM

    LM = size(T,3)

    BUOY(:,:,LM) =  T(:,:,LM) + gravbcp*ZLO(:,:,LM) + alhlbcp*Q(:,:,LM)

    do L=LM-1,1,-1
       BUOY(:,:,L) = BUOY(:,:,LM) - (T(:,:,L) + gravbcp*ZLO(:,:,L) + alhlbcp*QS(:,:,L))
       BUOY(:,:,L) = MAPL_GRAV*BUOY(:,:,L) / ( (1.+ alhlbcp*DQS(:,:,L))*T(:,:,L) )
    enddo

    BUOY(:,:,LM) = 0.0

    CAPE = 0.
    INHB = 0.

    do L=1,LM-1
       where(BUOY(:,:,L)>0.)
          CAPE = CAPE + BUOY(:,:,L)*DZ(:,:,L)
       end where
       where(BUOY(:,:,L)<0.)
          INHB = INHB - BUOY(:,:,L)*DZ(:,:,L)
       end where
    end do

    where(CAPE <= 0.0)
       CAPE=MAPL_UNDEF
       INHB=MAPL_UNDEF
    end where

  end subroutine BUOYANCY

end module BuoyancyMod
