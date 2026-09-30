module ion_drag_module

  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  implicit none
  private

  public :: compute_drag_fields

  real, parameter :: AMU_KG = 1.66053906660e-27
  real, parameter :: TEMP_REF_K = 300.0

  integer, parameter :: N_SUPPORTED_ION_SPECIES = 7

  ! IRI species order from OUTF(5:11):
  ! O+, H+, He+, O2+, NO+, cluster ions, N+.
  real, parameter :: ION_MASS_AMU(N_SUPPORTED_ION_SPECIES) = (/ &
       16.0, 1.0, 4.0, 32.0, 30.0, 30.0, 14.0 /)

  ! Effective ion-neutral momentum-transfer rate coefficients [m3 s-1].
  ! The temperature dependence is applied separately as sqrt(Tn/TEMP_REF_K).
  real, parameter :: ION_NEUTRAL_RATE_COEFF(N_SUPPORTED_ION_SPECIES) = (/ &
       2.5e-15, 2.2e-15, 0.8e-15, 1.0e-15, &
       4.2e-15, 4.0e-15, 4.0e-15 /)

contains

  subroutine compute_drag_fields(ui, vi, ne_m3, species_fraction, &
                                 u_neutral, v_neutral, &
                                 neutral_number_density, neutral_mass_density, &
                                 neutral_temp, drag_u, drag_v, &
                                 frictional_heating)
    ! Compute neutral-wind acceleration and frictional heating due to ion drag.
    !
    ! The formulation follows the structure of the original implementation:
    !
    !   a_n = (rho_i / rho_n) * nu_in * (V_i - V_n)
    !
    ! where rho_i is diagnosed from IRI electron density and ion composition,
    ! rho_n is the MSIS neutral mass density, and nu_in is an effective
    ! ion-neutral collision frequency.
    !
    ! ne_m3 is treated as the total positive-ion number density under
    ! quasi-neutrality. IRI ion-species fractions are normalized locally before
    ! they are used to form mean ion mass and the effective collision rate.
    !
    ! frictional_heating is returned as volumetric heating [W m-3]:
    !
    !   Q_fric = rho_i * nu_in * |V_i - V_n|^2
    !
    ! Inputs and outputs use the GEOS horizontal dimensions followed by level.

    real, intent(in) :: ui(:,:,:), vi(:,:,:)
    real, intent(in) :: ne_m3(:,:,:)
    real, intent(in) :: species_fraction(:,:,:,:)
    real, intent(in) :: u_neutral(:,:,:), v_neutral(:,:,:)
    real, intent(in) :: neutral_number_density(:,:,:)
    real, intent(in) :: neutral_mass_density(:,:,:)
    real, intent(in) :: neutral_temp(:,:,:)

    real, intent(out) :: drag_u(:,:,:), drag_v(:,:,:)
    real, intent(out) :: frictional_heating(:,:,:)

    integer :: i, j, k, species
    integer :: n_species
    real :: fraction_sum
    real :: ion_mean_mass_amu
    real :: rho_i, rho_n
    real :: rate_coeff, nu_in
    real :: thermal_factor
    real :: du, dv
    real :: coupling
    real :: fraction_value

    call validate_shapes( &
         ui, vi, ne_m3, species_fraction, u_neutral, v_neutral, &
         neutral_number_density, neutral_mass_density, neutral_temp, &
         drag_u, drag_v, frictional_heating)

    n_species = size(species_fraction, 1)

    if (n_species /= N_SUPPORTED_ION_SPECIES) then
      error stop 'IonDrag: unsupported number of IRI ion species'
    end if

    drag_u = 0.0
    drag_v = 0.0
    frictional_heating = 0.0

    do k = 1, size(ui, 3)
      do j = 1, size(ui, 2)
        do i = 1, size(ui, 1)

          if (.not. ieee_is_finite(ne_m3(i,j,k)) .or. ne_m3(i,j,k) <= 0.0) cycle
          if (.not. ieee_is_finite(neutral_number_density(i,j,k)) .or. &
              neutral_number_density(i,j,k) <= 0.0) cycle
          if (.not. ieee_is_finite(neutral_mass_density(i,j,k)) .or. &
              neutral_mass_density(i,j,k) <= 0.0) cycle

          fraction_sum = 0.0
          ion_mean_mass_amu = 0.0
          rate_coeff = 0.0

          do species = 1, n_species
            fraction_value = species_fraction(species,i,j,k)

            if (.not. ieee_is_finite(fraction_value)) cycle
            if (fraction_value <= 0.0) cycle

            fraction_sum = fraction_sum + fraction_value
            ion_mean_mass_amu = ion_mean_mass_amu + &
                 fraction_value * ION_MASS_AMU(species)
            rate_coeff = rate_coeff + &
                 fraction_value * ION_NEUTRAL_RATE_COEFF(species)
          end do

          if (fraction_sum <= 0.0) cycle

          ! Normalize because rounded/filtered IRI percentages do not
          ! necessarily sum to exactly one.
          ion_mean_mass_amu = ion_mean_mass_amu / fraction_sum
          rate_coeff = rate_coeff / fraction_sum

          rho_i = ne_m3(i,j,k) * ion_mean_mass_amu * AMU_KG
          rho_n = neutral_mass_density(i,j,k)

          if (.not. ieee_is_finite(rho_i) .or. rho_i <= 0.0) cycle

          if (ieee_is_finite(neutral_temp(i,j,k))) then
            thermal_factor = sqrt(max(neutral_temp(i,j,k), 1.0) / TEMP_REF_K)
          else
            thermal_factor = 1.0
          end if

          nu_in = neutral_number_density(i,j,k) * rate_coeff * thermal_factor

          if (.not. ieee_is_finite(nu_in) .or. nu_in <= 0.0) cycle

          coupling = (rho_i / rho_n) * nu_in

          if (.not. ieee_is_finite(coupling) .or. coupling <= 0.0) cycle

          du = ui(i,j,k) - u_neutral(i,j,k)
          dv = vi(i,j,k) - v_neutral(i,j,k)

          if (.not. ieee_is_finite(du) .or. .not. ieee_is_finite(dv)) cycle

          drag_u(i,j,k) = coupling * du
          drag_v(i,j,k) = coupling * dv

          frictional_heating(i,j,k) = &
               rho_i * nu_in * (du*du + dv*dv)

        end do
      end do
    end do

  end subroutine compute_drag_fields


  subroutine validate_shapes(ui, vi, ne_m3, species_fraction, &
                             u_neutral, v_neutral, &
                             neutral_number_density, neutral_mass_density, &
                             neutral_temp, drag_u, drag_v, frictional_heating)
    ! Verify that all three-dimensional fields share the same grid shape.
    real, intent(in) :: ui(:,:,:), vi(:,:,:)
    real, intent(in) :: ne_m3(:,:,:)
    real, intent(in) :: species_fraction(:,:,:,:)
    real, intent(in) :: u_neutral(:,:,:), v_neutral(:,:,:)
    real, intent(in) :: neutral_number_density(:,:,:)
    real, intent(in) :: neutral_mass_density(:,:,:)
    real, intent(in) :: neutral_temp(:,:,:)
    real, intent(out) :: drag_u(:,:,:), drag_v(:,:,:)
    real, intent(out) :: frictional_heating(:,:,:)

    integer :: nx, ny, nz

    nx = size(ui, 1)
    ny = size(ui, 2)
    nz = size(ui, 3)

    if (size(vi,1) /= nx .or. size(vi,2) /= ny .or. size(vi,3) /= nz) then
      error stop 'IonDrag: vi shape does not match ui'
    end if

    if (size(ne_m3,1) /= nx .or. size(ne_m3,2) /= ny .or. &
        size(ne_m3,3) /= nz) then
      error stop 'IonDrag: ne_m3 shape does not match ui'
    end if

    if (size(species_fraction,2) /= nx .or. &
        size(species_fraction,3) /= ny .or. &
        size(species_fraction,4) /= nz) then
      error stop 'IonDrag: species_fraction shape does not match ui'
    end if

    if (size(u_neutral,1) /= nx .or. size(u_neutral,2) /= ny .or. &
        size(u_neutral,3) /= nz) then
      error stop 'IonDrag: u_neutral shape does not match ui'
    end if

    if (size(v_neutral,1) /= nx .or. size(v_neutral,2) /= ny .or. &
        size(v_neutral,3) /= nz) then
      error stop 'IonDrag: v_neutral shape does not match ui'
    end if

    if (size(neutral_number_density,1) /= nx .or. &
        size(neutral_number_density,2) /= ny .or. &
        size(neutral_number_density,3) /= nz) then
      error stop 'IonDrag: neutral number-density shape does not match ui'
    end if

    if (size(neutral_mass_density,1) /= nx .or. &
        size(neutral_mass_density,2) /= ny .or. &
        size(neutral_mass_density,3) /= nz) then
      error stop 'IonDrag: neutral mass-density shape does not match ui'
    end if

    if (size(neutral_temp,1) /= nx .or. size(neutral_temp,2) /= ny .or. &
        size(neutral_temp,3) /= nz) then
      error stop 'IonDrag: neutral_temp shape does not match ui'
    end if

    if (size(drag_u,1) /= nx .or. size(drag_u,2) /= ny .or. &
        size(drag_u,3) /= nz) then
      error stop 'IonDrag: drag_u shape does not match ui'
    end if

    if (size(drag_v,1) /= nx .or. size(drag_v,2) /= ny .or. &
        size(drag_v,3) /= nz) then
      error stop 'IonDrag: drag_v shape does not match ui'
    end if

    if (size(frictional_heating,1) /= nx .or. &
        size(frictional_heating,2) /= ny .or. &
        size(frictional_heating,3) /= nz) then
      error stop 'IonDrag: frictional_heating shape does not match ui'
    end if

  end subroutine validate_shapes

end module ion_drag_module
