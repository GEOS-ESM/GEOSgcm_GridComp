module ion_drag_module

  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  implicit none
  private

  integer, parameter, public :: N_MAJOR_ION_SPECIES = 3
  integer, parameter, public :: ION_OP  = 1
  integer, parameter, public :: ION_O2P = 2
  integer, parameter, public :: ION_NOP = 3

  public :: compute_drag_fields

  real, parameter :: ELEMENTARY_CHARGE_C = 1.602176634e-19
  real, parameter :: ELECTRON_MASS_KG = 9.1093837015e-31
  real, parameter :: AMU_KG = 1.66053906660e-27
  real, parameter :: M3_TO_CM3 = 1.0e-6
  real, parameter :: MIN_TEMPERATURE_K = 1.0
  real, parameter :: MIN_MAGNETIC_FIELD_T = 1.0e-12

  ! Fixed WACCM-X ion-drag physics settings.
  real, parameter :: BURNSIDE_FACTOR = 1.5
  real, parameter :: ELECTRON_COLLISION_FACTOR = 1.0
  real, parameter :: DIP_MIN_RAD = 0.17

  real, parameter :: MASS_OP_KG  = 16.0 * AMU_KG
  real, parameter :: MASS_O2P_KG = 32.0 * AMU_KG
  real, parameter :: MASS_NOP_KG = 30.0 * AMU_KG

contains

  subroutine compute_drag_fields(ui, vi, u_neutral, v_neutral, neutral_temp, &
                                 ion_temp, electron_temp, ion_density_m3, &
                                 n_o_m3, n_o2_m3, n_n2_m3, rho_neutral, &
                                 b_north_t, b_east_t, b_down_t, timestep_s, &
                                 drag_u, drag_v, &
                                 joule_heating_wkg, sigma_pedersen, sigma_hall, &
                                 lxx, lyy, lxy, lyx)
    ! Compute neutral-wind ion-drag tendencies using the WACCM-X/TIE-GCM
    ! conductivity formulation.
    !
    ! The ion-neutral collision coefficients follow the WACCM-X iondrag.F90
    ! and NCAR GLOW conduct.f90 formulations. O+, O2+, and NO+ are treated
    ! explicitly. Neutral O, O2, and N2 number densities and neutral mass
    ! density are supplied by MSIS.
    !
    ! The Pedersen and Hall conductivities are converted to a geographic
    ! horizontal drag tensor. A local implicit 2x2 solve, matching the
    ! WACCM-X treatment, is used to return stable momentum tendencies.
    !
    ! joule_heating_wkg is a diagnostic specific heating rate [W kg-1].
    ! It is not applied to the GEOS temperature by this routine.

    real, intent(in) :: ui(:,:,:), vi(:,:,:)
    real, intent(in) :: u_neutral(:,:,:), v_neutral(:,:,:)
    real, intent(in) :: neutral_temp(:,:,:)
    real, intent(in) :: ion_temp(:,:,:), electron_temp(:,:,:)
    real, intent(in) :: ion_density_m3(:,:,:,:)
    real, intent(in) :: n_o_m3(:,:,:), n_o2_m3(:,:,:), n_n2_m3(:,:,:)
    real, intent(in) :: rho_neutral(:,:,:)
    real, intent(in) :: b_north_t(:,:,:), b_east_t(:,:,:), b_down_t(:,:,:)
    real, intent(in) :: timestep_s

    real, intent(out) :: drag_u(:,:,:), drag_v(:,:,:)
    real, intent(out) :: joule_heating_wkg(:,:,:)
    real, intent(out) :: sigma_pedersen(:,:,:), sigma_hall(:,:,:)
    real, intent(out) :: lxx(:,:,:), lyy(:,:,:), lxy(:,:,:), lyx(:,:,:)

    integer :: i, j, k
    real :: tn, ti, te, tr, sqrt_tr, log10_tr, sqrt_te
    real :: n_o_cm3, n_o2_cm3, n_n2_cm3
    real :: n_op, n_o2p, n_nop, ne_sigma
    real :: nu_o2p_o2, nu_op_o2, nu_nop_o2
    real :: nu_o2p_o, nu_op_o, nu_nop_o
    real :: nu_o2p_n2, nu_op_n2, nu_nop_n2
    real :: nu_o2p, nu_op, nu_nop, nu_e
    real :: bmag, omega_o2p, omega_op, omega_nop, omega_e
    real :: r_o2p, r_op, r_nop, r_e
    real :: q_over_b
    real :: lambda1, lambda2
    real :: dip_angle, dec_angle, sin_dip
    real :: sin_dec, cos_dec, sin2_dec, cos2_dec
    real :: lxx_norot, lyy_norot, lxy_norot, rotation_term
    real :: us, vs
    real :: dti, l11, l12, l21, l22, determinant, detr

    call validate_shapes(ui, vi, u_neutral, v_neutral, neutral_temp, &
                         ion_temp, electron_temp, ion_density_m3, &
                         n_o_m3, n_o2_m3, n_n2_m3, rho_neutral, &
                         b_north_t, b_east_t, b_down_t, &
                         drag_u, drag_v, joule_heating_wkg, &
                         sigma_pedersen, sigma_hall, lxx, lyy, lxy, lyx)

    if (.not. ieee_is_finite(timestep_s) .or. timestep_s <= 0.0) then
      error stop 'IonDrag: timestep_s must be finite and positive'
    end if

    drag_u = 0.0
    drag_v = 0.0
    joule_heating_wkg = 0.0
    sigma_pedersen = 0.0
    sigma_hall = 0.0
    lxx = 0.0
    lyy = 0.0
    lxy = 0.0
    lyx = 0.0

    dti = 1.0 / timestep_s

    do k = 1, size(ui, 3)
      do j = 1, size(ui, 2)
        do i = 1, size(ui, 1)

          if (.not. valid_point(i, j, k, ui, vi, u_neutral, v_neutral, &
                                neutral_temp, ion_temp, electron_temp, &
                                n_o_m3, n_o2_m3, n_n2_m3, rho_neutral, &
                                b_north_t, b_east_t, b_down_t)) cycle

          n_op  = max(ion_density_m3(ION_OP,  i,j,k), 0.0)
          n_o2p = max(ion_density_m3(ION_O2P, i,j,k), 0.0)
          n_nop = max(ion_density_m3(ION_NOP, i,j,k), 0.0)

          if (.not. ieee_is_finite(n_op) .or. &
              .not. ieee_is_finite(n_o2p) .or. &
              .not. ieee_is_finite(n_nop)) cycle

          ne_sigma = n_op + n_o2p + n_nop
          if (ne_sigma <= 0.0) cycle

          tn = max(neutral_temp(i,j,k), MIN_TEMPERATURE_K)
          ti = max(ion_temp(i,j,k), MIN_TEMPERATURE_K)
          te = max(electron_temp(i,j,k), MIN_TEMPERATURE_K)

          tr = 0.5 * (ti + tn)
          sqrt_tr = sqrt(tr)
          log10_tr = log10(tr)

          ! Convert MSIS neutral number densities to cm-3 because the
          ! WACCM-X/GLOW collision coefficients below are in cm3 s-1.
          n_o_cm3  = max(n_o_m3(i,j,k),  0.0) * M3_TO_CM3
          n_o2_cm3 = max(n_o2_m3(i,j,k), 0.0) * M3_TO_CM3
          n_n2_cm3 = max(n_n2_m3(i,j,k), 0.0) * M3_TO_CM3

          ! Ion-neutral momentum-transfer collision coefficients [cm3 s-1].
          nu_o2p_o2 = 2.59e-11 * sqrt_tr * (1.0 - 0.073 * log10_tr)**2
          nu_op_o2  = 6.64e-10
          nu_nop_o2 = 4.27e-10

          nu_o2p_o = 2.31e-10
          nu_op_o  = 3.67e-11 * sqrt_tr * (1.0 - 0.064 * log10_tr)**2 &
                     * BURNSIDE_FACTOR
          nu_nop_o = 2.44e-10

          nu_o2p_n2 = 4.13e-10
          nu_op_n2  = 6.82e-10
          nu_nop_n2 = 4.34e-10

          ! Total ion-neutral collision frequencies [s-1].
          nu_o2p = nu_o2p_o2*n_o2_cm3 + nu_o2p_o*n_o_cm3 + &
                   nu_o2p_n2*n_n2_cm3
          nu_op  = nu_op_o2*n_o2_cm3 + nu_op_o*n_o_cm3 + &
                   nu_op_n2*n_n2_cm3
          nu_nop = nu_nop_o2*n_o2_cm3 + nu_nop_o*n_o_cm3 + &
                   nu_nop_n2*n_n2_cm3

          ! Electron-neutral collision frequency [s-1].
          sqrt_te = sqrt(te)
          nu_e = 2.33e-11*n_n2_cm3*te*(1.0 - 1.21e-4*te) &
               + 1.82e-10*n_o2_cm3*sqrt_te*(1.0 + 3.60e-2*sqrt_te) &
               + 8.90e-11*n_o_cm3*sqrt_te*(1.0 + 5.70e-4*te)
          nu_e = max(nu_e, 0.0) * ELECTRON_COLLISION_FACTOR

          bmag = sqrt(b_north_t(i,j,k)**2 + b_east_t(i,j,k)**2 + &
                      b_down_t(i,j,k)**2)
          if (.not. ieee_is_finite(bmag) .or. bmag <= MIN_MAGNETIC_FIELD_T) cycle

          omega_op  = ELEMENTARY_CHARGE_C * bmag / MASS_OP_KG
          omega_o2p = ELEMENTARY_CHARGE_C * bmag / MASS_O2P_KG
          omega_nop = ELEMENTARY_CHARGE_C * bmag / MASS_NOP_KG
          omega_e   = ELEMENTARY_CHARGE_C * bmag / ELECTRON_MASS_KG

          r_op  = nu_op  / omega_op
          r_o2p = nu_o2p / omega_o2p
          r_nop = nu_nop / omega_nop
          r_e   = nu_e   / omega_e

          q_over_b = ELEMENTARY_CHARGE_C / bmag

          sigma_pedersen(i,j,k) = q_over_b * ( &
               n_op  * r_op  / (1.0 + r_op**2) + &
               n_o2p * r_o2p / (1.0 + r_o2p**2) + &
               n_nop * r_nop / (1.0 + r_nop**2) + &
               ne_sigma * r_e / (1.0 + r_e**2))

          sigma_hall(i,j,k) = q_over_b * ( &
               ne_sigma / (1.0 + r_e**2) - &
               n_op  / (1.0 + r_op**2) - &
               n_o2p / (1.0 + r_o2p**2) - &
               n_nop / (1.0 + r_nop**2))

          ! Convert conductivity to neutral drag-rate coefficients [s-1].
          lambda1 = sigma_pedersen(i,j,k) * bmag*bmag / rho_neutral(i,j,k)
          lambda2 = sigma_hall(i,j,k) * bmag*bmag / rho_neutral(i,j,k)

          ! IGRF local coordinates are north/east/down. Follow the WACCM-X
          ! magnetic geometry and dip-angle floor exactly.
          dip_angle = atan2(b_down_t(i,j,k), &
                            sqrt(b_north_t(i,j,k)**2 + b_east_t(i,j,k)**2))
          dec_angle = -atan2(b_east_t(i,j,k), b_north_t(i,j,k))

          if (abs(dip_angle) >= DIP_MIN_RAD) then
            sin_dip = sin(dip_angle)
          else if (dip_angle >= 0.0) then
            sin_dip = sin(DIP_MIN_RAD)
          else
            sin_dip = -sin(DIP_MIN_RAD)
          end if

          sin_dec = sin(dec_angle)
          cos_dec = cos(dec_angle)
          sin2_dec = sin_dec*sin_dec
          cos2_dec = cos_dec*cos_dec

          lxx_norot = lambda1
          lyy_norot = lambda1 * sin_dip*sin_dip
          lxy_norot = lambda2 * sin_dip

          lxx(i,j,k) = lxx_norot*cos2_dec + lyy_norot*sin2_dec
          lyy(i,j,k) = lyy_norot*cos2_dec + lxx_norot*sin2_dec
          rotation_term = (lyy_norot - lxx_norot) * sin_dec * cos_dec
          lyx(i,j,k) = lxy_norot - rotation_term
          lxy(i,j,k) = lxy_norot + rotation_term

          us = ui(i,j,k) - u_neutral(i,j,k)
          vs = vi(i,j,k) - v_neutral(i,j,k)

          ! Diagnostic ion-neutral friction/Joule heating per unit neutral mass.
          joule_heating_wkg(i,j,k) = us*us*lxx(i,j,k) &
               + us*vs*(lxy(i,j,k) - lyx(i,j,k)) &
               + vs*vs*lyy(i,j,k)

          ! WACCM-X local implicit 2x2 ion-drag solve.
          l11 = dti + lxx(i,j,k)
          l12 = lxy(i,j,k)
          l21 = -lyx(i,j,k)
          l22 = dti + lyy(i,j,k)
          determinant = l11*l22 - l12*l21

          if (.not. ieee_is_finite(determinant) .or. &
              abs(determinant) <= tiny(1.0)) cycle

          detr = dti / determinant
          drag_u(i,j,k) = dti * (detr*(l12*vs - l22*us) + us)
          drag_v(i,j,k) = dti * (detr*(l21*us - l11*vs) + vs)

          if (.not. ieee_is_finite(drag_u(i,j,k))) drag_u(i,j,k) = 0.0
          if (.not. ieee_is_finite(drag_v(i,j,k))) drag_v(i,j,k) = 0.0
          if (.not. ieee_is_finite(joule_heating_wkg(i,j,k))) then
            joule_heating_wkg(i,j,k) = 0.0
          end if

        end do
      end do
    end do

  end subroutine compute_drag_fields


  logical function valid_point(i, j, k, ui, vi, u_neutral, v_neutral, &
                               neutral_temp, ion_temp, electron_temp, &
                               n_o_m3, n_o2_m3, n_n2_m3, rho_neutral, &
                               b_north_t, b_east_t, b_down_t)
    integer, intent(in) :: i, j, k
    real, intent(in) :: ui(:,:,:), vi(:,:,:)
    real, intent(in) :: u_neutral(:,:,:), v_neutral(:,:,:)
    real, intent(in) :: neutral_temp(:,:,:), ion_temp(:,:,:), electron_temp(:,:,:)
    real, intent(in) :: n_o_m3(:,:,:), n_o2_m3(:,:,:), n_n2_m3(:,:,:)
    real, intent(in) :: rho_neutral(:,:,:)
    real, intent(in) :: b_north_t(:,:,:), b_east_t(:,:,:), b_down_t(:,:,:)

    valid_point = &
         ieee_is_finite(ui(i,j,k)) .and. ieee_is_finite(vi(i,j,k)) .and. &
         ieee_is_finite(u_neutral(i,j,k)) .and. &
         ieee_is_finite(v_neutral(i,j,k)) .and. &
         ieee_is_finite(neutral_temp(i,j,k)) .and. neutral_temp(i,j,k) > 0.0 .and. &
         ieee_is_finite(ion_temp(i,j,k)) .and. ion_temp(i,j,k) > 0.0 .and. &
         ieee_is_finite(electron_temp(i,j,k)) .and. electron_temp(i,j,k) > 0.0 .and. &
         ieee_is_finite(n_o_m3(i,j,k)) .and. n_o_m3(i,j,k) >= 0.0 .and. &
         ieee_is_finite(n_o2_m3(i,j,k)) .and. n_o2_m3(i,j,k) >= 0.0 .and. &
         ieee_is_finite(n_n2_m3(i,j,k)) .and. n_n2_m3(i,j,k) >= 0.0 .and. &
         ieee_is_finite(rho_neutral(i,j,k)) .and. rho_neutral(i,j,k) > 0.0 .and. &
         ieee_is_finite(b_north_t(i,j,k)) .and. &
         ieee_is_finite(b_east_t(i,j,k)) .and. &
         ieee_is_finite(b_down_t(i,j,k))

  end function valid_point


  subroutine validate_shapes(ui, vi, u_neutral, v_neutral, neutral_temp, &
                             ion_temp, electron_temp, ion_density_m3, &
                             n_o_m3, n_o2_m3, n_n2_m3, rho_neutral, &
                             b_north_t, b_east_t, b_down_t, &
                             drag_u, drag_v, joule_heating_wkg, &
                             sigma_pedersen, sigma_hall, lxx, lyy, lxy, lyx)
    real, intent(in) :: ui(:,:,:), vi(:,:,:)
    real, intent(in) :: u_neutral(:,:,:), v_neutral(:,:,:)
    real, intent(in) :: neutral_temp(:,:,:), ion_temp(:,:,:), electron_temp(:,:,:)
    real, intent(in) :: ion_density_m3(:,:,:,:)
    real, intent(in) :: n_o_m3(:,:,:), n_o2_m3(:,:,:), n_n2_m3(:,:,:)
    real, intent(in) :: rho_neutral(:,:,:)
    real, intent(in) :: b_north_t(:,:,:), b_east_t(:,:,:), b_down_t(:,:,:)
    real, intent(out) :: drag_u(:,:,:), drag_v(:,:,:), joule_heating_wkg(:,:,:)
    real, intent(out) :: sigma_pedersen(:,:,:), sigma_hall(:,:,:)
    real, intent(out) :: lxx(:,:,:), lyy(:,:,:), lxy(:,:,:), lyx(:,:,:)

    integer :: nx, ny, nz

    nx = size(ui, 1)
    ny = size(ui, 2)
    nz = size(ui, 3)

    if (.not. same_shape_3d(vi, nx, ny, nz)) &
      error stop 'IonDrag: vi shape does not match ui'
    if (.not. same_shape_3d(u_neutral, nx, ny, nz)) &
      error stop 'IonDrag: u_neutral shape does not match ui'
    if (.not. same_shape_3d(v_neutral, nx, ny, nz)) &
      error stop 'IonDrag: v_neutral shape does not match ui'
    if (.not. same_shape_3d(neutral_temp, nx, ny, nz)) &
      error stop 'IonDrag: neutral_temp shape does not match ui'
    if (.not. same_shape_3d(ion_temp, nx, ny, nz)) &
      error stop 'IonDrag: ion_temp shape does not match ui'
    if (.not. same_shape_3d(electron_temp, nx, ny, nz)) &
      error stop 'IonDrag: electron_temp shape does not match ui'

    if (size(ion_density_m3,1) /= N_MAJOR_ION_SPECIES .or. &
        size(ion_density_m3,2) /= nx .or. size(ion_density_m3,3) /= ny .or. &
        size(ion_density_m3,4) /= nz) then
      error stop 'IonDrag: ion_density_m3 shape does not match expected dimensions'
    end if

    if (.not. same_shape_3d(n_o_m3, nx, ny, nz)) &
      error stop 'IonDrag: n_o_m3 shape does not match ui'
    if (.not. same_shape_3d(n_o2_m3, nx, ny, nz)) &
      error stop 'IonDrag: n_o2_m3 shape does not match ui'
    if (.not. same_shape_3d(n_n2_m3, nx, ny, nz)) &
      error stop 'IonDrag: n_n2_m3 shape does not match ui'
    if (.not. same_shape_3d(rho_neutral, nx, ny, nz)) &
      error stop 'IonDrag: rho_neutral shape does not match ui'
    if (.not. same_shape_3d(b_north_t, nx, ny, nz)) &
      error stop 'IonDrag: b_north_t shape does not match ui'
    if (.not. same_shape_3d(b_east_t, nx, ny, nz)) &
      error stop 'IonDrag: b_east_t shape does not match ui'
    if (.not. same_shape_3d(b_down_t, nx, ny, nz)) &
      error stop 'IonDrag: b_down_t shape does not match ui'

    if (.not. same_shape_3d(drag_u, nx, ny, nz)) &
      error stop 'IonDrag: drag_u shape does not match ui'
    if (.not. same_shape_3d(drag_v, nx, ny, nz)) &
      error stop 'IonDrag: drag_v shape does not match ui'
    if (.not. same_shape_3d(joule_heating_wkg, nx, ny, nz)) &
      error stop 'IonDrag: joule_heating_wkg shape does not match ui'
    if (.not. same_shape_3d(sigma_pedersen, nx, ny, nz)) &
      error stop 'IonDrag: sigma_pedersen shape does not match ui'
    if (.not. same_shape_3d(sigma_hall, nx, ny, nz)) &
      error stop 'IonDrag: sigma_hall shape does not match ui'
    if (.not. same_shape_3d(lxx, nx, ny, nz)) &
      error stop 'IonDrag: lxx shape does not match ui'
    if (.not. same_shape_3d(lyy, nx, ny, nz)) &
      error stop 'IonDrag: lyy shape does not match ui'
    if (.not. same_shape_3d(lxy, nx, ny, nz)) &
      error stop 'IonDrag: lxy shape does not match ui'
    if (.not. same_shape_3d(lyx, nx, ny, nz)) &
      error stop 'IonDrag: lyx shape does not match ui'

  end subroutine validate_shapes


  logical function same_shape_3d(field, nx, ny, nz)
    real, intent(in) :: field(:,:,:)
    integer, intent(in) :: nx, ny, nz

    same_shape_3d = size(field,1) == nx .and. &
                    size(field,2) == ny .and. &
                    size(field,3) == nz
  end function same_shape_3d

end module ion_drag_module
