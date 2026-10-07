module iri_input_module

  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  implicit none
  private

  integer, parameter, public :: N_MAJOR_ION_SPECIES = 3
  integer, parameter, public :: ION_OP  = 1
  integer, parameter, public :: ION_O2P = 2
  integer, parameter, public :: ION_NOP = 3

  integer, parameter :: IRI_MAX_OUTPUT_LEVELS = 1000
  integer, parameter :: IRI_OUTPUT_INDEX(N_MAJOR_ION_SPECIES) = (/ 5, 8, 9 /)

  ! IRI is evaluated on this internal altitude grid and then interpolated to
  ! the actual GEOS model-level altitudes. These are implementation settings,
  ! not runtime configuration parameters.
  real, parameter :: IRI_PROFILE_BOTTOM_KM = 80.0
  real, parameter :: IRI_PROFILE_TOP_KM = 250.0
  real, parameter :: IRI_PROFILE_STEP_KM = 2.0

  logical, save :: iri_indices_initialized = .false.

  public :: get_iri_state

  interface
    subroutine IRI_SUB(JF, JMAG, ALATI, ALONG, IYYYY, MMDD, DHOUR, &
                       HEIBEG, HEIEND, HEISTP, OUTF, OARR)
      logical, intent(in) :: JF(50)
      integer, intent(in) :: JMAG
      real, intent(in) :: ALATI, ALONG
      integer, intent(in) :: IYYYY, MMDD
      real, intent(in) :: DHOUR, HEIBEG, HEIEND, HEISTP
      real, intent(out) :: OUTF(20, 1000)
      real, intent(inout) :: OARR(100)
    end subroutine IRI_SUB

    subroutine read_ig_rz()
    end subroutine read_ig_rz

    subroutine readapf107()
    end subroutine readapf107
  end interface

contains

  subroutine get_iri_state(latitudes_deg, longitudes_deg, iyear, doy, &
                           ut_hour, f107_daily, f107_81day, alt_km, &
                           ne_m3, ion_density_m3, ion_temp, electron_temp)
    ! Evaluate IRI independently at each GEOS horizontal grid point and
    ! interpolate the requested plasma state to the actual GEOS altitudes.
    !
    ! IRI OUTF fields used here:
    !   OUTF(1,*) : electron density [m-3]
    !   OUTF(3,*) : ion temperature [K]
    !   OUTF(4,*) : electron temperature [K]
    !   OUTF(5,*) : O+ density [m-3] when JF(22)=.false.
    !   OUTF(8,*) : O2+ density [m-3] when JF(22)=.false.
    !   OUTF(9,*) : NO+ density [m-3] when JF(22)=.false.

    real, intent(in) :: latitudes_deg(:,:), longitudes_deg(:,:)
    integer, intent(in) :: iyear, doy
    real, intent(in) :: ut_hour
    real, intent(in) :: f107_daily, f107_81day
    real, intent(in) :: alt_km(:,:,:)

    real, intent(out) :: ne_m3(:,:,:)
    real, intent(out) :: ion_density_m3(:,:,:,:)
    real, intent(out) :: ion_temp(:,:,:)
    real, intent(out) :: electron_temp(:,:,:)

    logical :: jf(50)
    integer :: jmag, mmdd
    real, allocatable :: outf(:,:), oarr(:)
    real :: dhour_ut
    real :: value

    integer :: i, j, lev, species
    integer :: im, jm, nlev
    integer :: n_iri_levels

    im = size(alt_km, 1)
    jm = size(alt_km, 2)
    nlev = size(alt_km, 3)

    call validate_iri_inputs(latitudes_deg, longitudes_deg, alt_km, &
                             f107_daily, f107_81day, ne_m3, &
                             ion_density_m3, ion_temp, electron_temp)

    mmdd = -doy
    dhour_ut = ut_hour + 25.0

    call initialize_iri_indices()

    ! Start from the IRI defaults used by the current GEOS implementation.
    jf = .true.
    jf(4:6) = .false.
    jf(23) = .false.
    jf(30) = .false.
    jf(33) = .false.
    jf(34) = .false.
    jf(35) = .false.
    jf(39:40) = .false.
    jf(47) = .false.

    ! Return individual ion number densities rather than percent composition.
    jf(22) = .false.

    ! Supply daily and 81-day F10.7 from the same GEOS-MLT forcing used by MSIS.
    jf(25) = .false.
    jf(32) = .false.

    ! Preserve the current no-storm configuration for the F2 and E regions.
    jf(26) = .false.
    jf(35) = .false.

    jmag = 0

    n_iri_levels = int((IRI_PROFILE_TOP_KM - IRI_PROFILE_BOTTOM_KM) / &
                       IRI_PROFILE_STEP_KM) + 1

    if (n_iri_levels > IRI_MAX_OUTPUT_LEVELS) then
      error stop 'IRI: internal altitude profile exceeds IRI output capacity'
    end if

    allocate(outf(20, IRI_MAX_OUTPUT_LEVELS), oarr(100))

    ne_m3 = 0.0
    ion_density_m3 = 0.0
    ion_temp = 0.0
    electron_temp = 0.0

    do j = 1, jm
      do i = 1, im

        if (.not. ieee_is_finite(latitudes_deg(i,j)) .or. &
            .not. ieee_is_finite(longitudes_deg(i,j))) cycle

        outf = 0.0
        oarr = -1.0
        oarr(41) = f107_daily
        oarr(46) = f107_81day

        call IRI_SUB(jf, jmag, latitudes_deg(i,j), longitudes_deg(i,j), &
                     iyear, mmdd, dhour_ut, IRI_PROFILE_BOTTOM_KM, &
                     IRI_PROFILE_TOP_KM, IRI_PROFILE_STEP_KM, outf, oarr)

        do lev = 1, nlev
          if (.not. ieee_is_finite(alt_km(i,j,lev))) cycle

          value = interpolate_iri_profile(outf(1,1:n_iri_levels), &
                                          n_iri_levels, alt_km(i,j,lev))
          if (ieee_is_finite(value) .and. value > 0.0) then
            ne_m3(i,j,lev) = value
          end if

          value = interpolate_iri_profile(outf(3,1:n_iri_levels), &
                                          n_iri_levels, alt_km(i,j,lev))
          if (ieee_is_finite(value) .and. value > 0.0) then
            ion_temp(i,j,lev) = value
          end if

          value = interpolate_iri_profile(outf(4,1:n_iri_levels), &
                                          n_iri_levels, alt_km(i,j,lev))
          if (ieee_is_finite(value) .and. value > 0.0) then
            electron_temp(i,j,lev) = value
          end if

          do species = 1, N_MAJOR_ION_SPECIES
            value = interpolate_iri_profile( &
                 outf(IRI_OUTPUT_INDEX(species),1:n_iri_levels), &
                 n_iri_levels, alt_km(i,j,lev))
            if (ieee_is_finite(value) .and. value > 0.0) then
              ion_density_m3(species,i,j,lev) = value
            end if
          end do
        end do
      end do
    end do

    deallocate(outf, oarr)

  end subroutine get_iri_state


  real function interpolate_iri_profile(profile, n_levels, altitude)
    ! Linearly interpolate the fixed 2-km IRI profile to one GEOS altitude.
    ! Values outside 80-250 km are clamped to the nearest profile endpoint.
    real, intent(in) :: profile(:)
    integer, intent(in) :: n_levels
    real, intent(in) :: altitude

    integer :: lower_idx, upper_idx
    real :: position, weight
    real :: lower_value, upper_value

    interpolate_iri_profile = 0.0

    if (n_levels < 1) return
    if (.not. ieee_is_finite(altitude)) return

    if (altitude <= IRI_PROFILE_BOTTOM_KM) then
      if (ieee_is_finite(profile(1))) interpolate_iri_profile = profile(1)
      return
    end if

    position = (altitude - IRI_PROFILE_BOTTOM_KM) / IRI_PROFILE_STEP_KM

    if (position >= real(n_levels - 1)) then
      if (ieee_is_finite(profile(n_levels))) then
        interpolate_iri_profile = profile(n_levels)
      end if
      return
    end if

    lower_idx = floor(position) + 1
    upper_idx = lower_idx + 1
    weight = position - real(lower_idx - 1)

    lower_value = profile(lower_idx)
    upper_value = profile(upper_idx)

    if (ieee_is_finite(lower_value) .and. ieee_is_finite(upper_value)) then
      interpolate_iri_profile = (1.0 - weight)*lower_value + weight*upper_value
    else if (ieee_is_finite(lower_value)) then
      interpolate_iri_profile = lower_value
    else if (ieee_is_finite(upper_value)) then
      interpolate_iri_profile = upper_value
    end if

  end function interpolate_iri_profile


  subroutine initialize_iri_indices()
    ! Read the IRI monthly and daily index tables once per process.
    if (iri_indices_initialized) return

    call read_ig_rz()
    call readapf107()

    iri_indices_initialized = .true.
  end subroutine initialize_iri_indices


  subroutine validate_iri_inputs(latitudes_deg, longitudes_deg, alt_km, &
                                 f107_daily, f107_81day, ne_m3, &
                                 ion_density_m3, ion_temp, electron_temp)
    real, intent(in) :: latitudes_deg(:,:), longitudes_deg(:,:)
    real, intent(in) :: alt_km(:,:,:)
    real, intent(in) :: f107_daily, f107_81day
    real, intent(out) :: ne_m3(:,:,:)
    real, intent(out) :: ion_density_m3(:,:,:,:)
    real, intent(out) :: ion_temp(:,:,:), electron_temp(:,:,:)

    integer :: im, jm, nlev

    im = size(alt_km, 1)
    jm = size(alt_km, 2)
    nlev = size(alt_km, 3)

    if (.not. ieee_is_finite(f107_daily) .or. f107_daily <= 0.0 .or. &
        .not. ieee_is_finite(f107_81day) .or. f107_81day <= 0.0) then
      error stop 'IRI: F10.7 inputs must be finite and positive'
    end if

    if (size(latitudes_deg,1) /= im .or. size(latitudes_deg,2) /= jm) then
      error stop 'IRI: latitude shape does not match alt_km'
    end if
    if (size(longitudes_deg,1) /= im .or. size(longitudes_deg,2) /= jm) then
      error stop 'IRI: longitude shape does not match alt_km'
    end if
    if (size(ne_m3,1) /= im .or. size(ne_m3,2) /= jm .or. &
        size(ne_m3,3) /= nlev) then
      error stop 'IRI: ne_m3 shape does not match alt_km'
    end if
    if (size(ion_density_m3,1) /= N_MAJOR_ION_SPECIES .or. &
        size(ion_density_m3,2) /= im .or. size(ion_density_m3,3) /= jm .or. &
        size(ion_density_m3,4) /= nlev) then
      error stop 'IRI: ion_density_m3 shape does not match expected dimensions'
    end if
    if (size(ion_temp,1) /= im .or. size(ion_temp,2) /= jm .or. &
        size(ion_temp,3) /= nlev) then
      error stop 'IRI: ion_temp shape does not match alt_km'
    end if
    if (size(electron_temp,1) /= im .or. size(electron_temp,2) /= jm .or. &
        size(electron_temp,3) /= nlev) then
      error stop 'IRI: electron_temp shape does not match alt_km'
    end if

  end subroutine validate_iri_inputs

end module iri_input_module
