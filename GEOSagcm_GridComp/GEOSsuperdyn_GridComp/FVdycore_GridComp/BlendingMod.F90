!-------------------------------------------------------------------------
!      NASA/GSFC, Global Modeling & Assimilation Office, Code 610.1      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  BlendingMod.F90 --- Blend forecast and analysis wind and 
!           virtual potential temperature.
! 
! !INTERFACE:
!

module  BlendingMod
            

! !USES:

    use MAPL_Mod, only : MAPL_CP, MAPL_KAPPA


    implicit none
    
    private

!
! !PUBLIC MEMBER FUNCTIONS:
!
    public  blend_wind_height

!
! !DESCRIPTION: This module implements blending of forecast and analysis 
!               fields.
!
!
! !REVISION HISTORY: 
!
!  12Aug2010  Darmenov    Initial code.
!
!EOP
!-------------------------------------------------------------------------

    integer, parameter :: r8 = 8
    integer, parameter :: r4 = 4

    real(r8), parameter :: c_p   = MAPL_CP
    real(r8), parameter :: kappa = MAPL_KAPPA


contains

!-------------------------------------------------------------------------
!      NASA/GSFC Global Modeling & Assimilation Office, Code 610.1       !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  blend_wind_height --- Blend wind and scaled virtual potential 
!             temperature by tampering wind and geopotential height.
! 
! !INTERFACE:
!
subroutine  blend_wind_height(u_f, v_f, thv_f, pe_f, &
                              u_a, v_a, thv_a, pe_a, &
                              p_below, p_above,      &
                              ims, ime, jms, jme, km)
!
! !USES:
!
    implicit none
!
! !INPUT PARAMETERS: 
!
    real(r8), dimension(ims:ime, jms:jme, km),   intent(in)    :: u_f        ! forecast wind u-component
    real(r8), dimension(ims:ime, jms:jme, km),   intent(in)    :: v_f        ! forecast wind v-component
    real(r8), dimension(ims:ime, jms:jme, km),   intent(in)    :: thv_f      ! forecast scaled virtual potential temperature
    real(r8), dimension(ims:ime, jms:jme, km+1), intent(in)    :: pe_f       ! forecast pressure at edge levels
    real(r8), dimension(ims:ime, jms:jme, km),   intent(inout) :: u_a        ! analysis wind u-component
    real(r8), dimension(ims:ime, jms:jme, km),   intent(inout) :: v_a        ! analysis wind v-component
    real(r8), dimension(ims:ime, jms:jme, km),   intent(inout) :: thv_a      ! analysis scaled virtual potential temperature
    real(r8), dimension(ims:ime, jms:jme, km+1), intent(in)    :: pe_a       ! analysis pressure at edge levels

    real(r8), intent(in) :: p_below                                          ! blending zone - bottom pressure
    real(r8), intent(in) :: p_above                                          ! blending zone - top pressure
    
    integer,  intent(in) :: ims                                              ! tile indices
    integer,  intent(in) :: ime                                              !
    integer,  intent(in) :: jms                                              !
    integer,  intent(in) :: jme                                              !
    integer,  intent(in) :: km                                               ! global vertical dimension

!
! !OUTPUT PARAMETERS:
!

!
! !DESCRIPTION: 
!
! !REVISION HISTORY: 
!
!  12Aug2010 Darmenov    Initial code
!
!EOP
!-------------------------------------------------------------------------

! !LOCAL VARIABLES:

    character(len=*), parameter :: myname = 'blend_wind_height'

    integer, parameter :: k1 = 1
    integer, parameter :: k2 = 2
    
    integer :: k

    real(r8), dimension (:,:,:), allocatable :: pl_a

    real(r8), dimension (:,:,:), allocatable :: dz_pe_f_kappa
    real(r8), dimension (:,:,:), allocatable :: dz_pe_a_kappa

    real(r8), dimension (:,:,:), allocatable :: phi_f
    real(r8), dimension (:,:,:), allocatable :: phi_a

    real(r8), dimension (:,:,:), allocatable :: alpha
    real(r8), dimension (:,:,:), allocatable :: alpha_delta_phi
    real(r8), dimension (:,:,:), allocatable :: dz_alpha_delta_phi

    real(r8), dimension (:,:,:), allocatable :: delta_thv


    allocate(dz_pe_f_kappa(ims:ime,jms:jme,km))
    allocate(dz_pe_a_kappa(ims:ime,jms:jme,km))

    allocate(phi_f(ims:ime,jms:jme,km+1))
    allocate(phi_a(ims:ime,jms:jme,km+1))

    allocate(alpha_delta_phi(ims:ime,jms:jme,km+1))
    allocate(dz_alpha_delta_phi(ims:ime,jms:jme,km))


    ! d/dz(p^kappa) = p^kappa(:,:,k+1) - p^kappa(:,:,k)
    dz_pe_f_kappa(:,:,k1:km) = pe_f(:,:,k2:km+1)**kappa - pe_f(:,:,k1:km)**kappa
    dz_pe_a_kappa(:,:,k1:km) = pe_a(:,:,k2:km+1)**kappa - pe_a(:,:,k1:km)**kappa


    ! forecast and analysis geopotentials
    phi_f(:,:,km+1) = 0.0
    phi_a(:,:,km+1) = 0.0

    do k = km, k1, -1
        phi_f(:,:,k) = phi_f(:,:,k+1) - c_p * thv_f(:,:,k) * dz_pe_f_kappa(:,:,k)
        phi_a(:,:,k) = phi_a(:,:,k+1) - c_p * thv_a(:,:,k) * dz_pe_a_kappa(:,:,k)
    end do


    ! d/dz(blending factor * geopotential increment)
    alpha_delta_phi = (pe_a - p_above)/(p_below - p_above) * (phi_a - phi_f)
    dz_alpha_delta_phi(:,:,k1:km) = alpha_delta_phi(:,:,k2:km+1) - alpha_delta_phi(:,:,k1:km)

    ! not needed past this point
    deallocate(alpha_delta_phi)
    deallocate(phi_f)
    deallocate(phi_a)


    ! blending factor used to tamper wind and geopotential increments
    allocate(alpha(ims:ime,jms:jme,km))
    allocate(delta_thv(ims:ime,jms:jme,km))
    allocate(pl_a(ims:ime,jms:jme,km))

    pl_a(:,:,k1:km) = 0.5 * (pe_a(:,:,k1:km) + pe_a(:,:,k2:km+1))

    ! calculate the blending factor (alpha) and temperature increment
    where (pl_a <= p_above)
        alpha = 0.0
        delta_thv = 0.0
    elsewhere (pl_a < p_below)
        alpha = (pl_a - p_above)/(p_below - p_above)
        delta_thv = (-1/c_p) * dz_alpha_delta_phi * dz_pe_a_kappa
    elsewhere
        alpha = 1.0
        delta_thv = thv_a - thv_f
    end where


    ! increment wind and temperature
    where (pl_a < p_below)        
        u_a = u_f + alpha*(u_a - u_f)
        v_a = v_f + alpha*(v_a - v_f)
        
        thv_a = thv_f + delta_thv
    end where
 
    return
end subroutine blend_wind_height


end module BlendingMod

