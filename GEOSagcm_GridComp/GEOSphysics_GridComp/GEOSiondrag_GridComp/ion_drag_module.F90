module ion_drag_module
  implicit none
  real, parameter :: AMU_KG = 1.66053906660e-27
  real, parameter :: CM3_TO_M3 = 1.0e6
  real, parameter :: TEMP_REF = 300.0

contains

  subroutine compute_drag_fields(ui, vi, latitudes, ne_m3, species_fraction, &
                                 u_neutral, v_neutral, neutral_density, &
                                 neutral_temp, mean_mass, drag_u, drag_v, &
                                 frictional_heating)
    ! Inputs -- all now (lat, lon, lev) for the top NLEV_IONDRAG levels
    real, intent(in) :: ui(:,:,:), vi(:,:,:)
    real, intent(in) :: latitudes(:)
    real, intent(in) :: ne_m3(:,:,:)
    real, intent(in) :: species_fraction(:,:,:,:)      ! (species, lat, lon, lev)
    real, intent(in) :: u_neutral(:,:,:), v_neutral(:,:,:)
    real, intent(in) :: neutral_density(:,:,:)
    real, intent(in) :: neutral_temp(:,:,:)
    real, intent(in) :: mean_mass

    ! Outputs
    real, intent(out) :: drag_u(:,:,:), drag_v(:,:,:)
    real, intent(out) :: frictional_heating(:,:,:)

    integer :: i, j, k, lev, n_species, nlev
    real :: rho_i, rho_n, nu_in, coupling, du, dv
    real :: thermal_factor, coefficient

    real, parameter :: ION_NEUTRAL_COEFF(7) = (/ 2.5e-15, 2.2e-15, 0.8e-15, 1.0e-15, 4.2e-15, 4.0e-15, 4.0e-15 /)

    n_species = size(species_fraction, 1)
    nlev = size(ui, 3)

    do lev = 1, nlev
      do k = 1, size(ui, 1)   ! Latitude
        do j = 1, size(ui, 2) ! Longitude

          rho_i = 0.0
          do i = 1, n_species
            rho_i = rho_i + species_fraction(i, k, j, lev) * ne_m3(k, j, lev) * AMU_KG
          end do

          rho_n = neutral_density(k, j, lev) * mean_mass * AMU_KG

          coefficient = 0.0
          do i = 1, n_species
            coefficient = coefficient + species_fraction(i, k, j, lev) * ION_NEUTRAL_COEFF(i)
          end do
          thermal_factor = sqrt(max(neutral_temp(k, j, lev), 1.0) / TEMP_REF)
          nu_in = neutral_density(k, j, lev) * coefficient * thermal_factor

          if (rho_n > 0.0) then
            coupling = (rho_i / rho_n) * nu_in
          else
            coupling = 0.0
          end if

          du = ui(k, j, lev) - u_neutral(k, j, lev)
          dv = vi(k, j, lev) - v_neutral(k, j, lev)
          drag_u(k, j, lev) = coupling * du
          drag_v(k, j, lev) = coupling * dv

          frictional_heating(k, j, lev) = rho_i * nu_in * (du**2 + dv**2)

        end do
      end do
    end do

  end subroutine compute_drag_fields

end module ion_drag_module
