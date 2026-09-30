module iri_input_module

  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  implicit none
  private

  integer, parameter, public :: N_ION_SPECIES = 7
  integer, parameter :: IRI_MAX_OUTPUT_LEVELS = 1000
  logical, save :: iri_indices_initialized = .false.

  public :: get_iri_densities

  ! Explicit interface for the external IRI routine in irisub.for.
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

  subroutine get_iri_densities(latitudes_deg, longitudes_deg, iyear, doy, &
                               ut_hour, f107_daily, f107_81day, alt_km, &
                               hbeg, hend, hstep, ne_m3, species_fraction)
    ! Evaluate IRI independently at each GEOS horizontal grid point.
    !
    ! latitudes_deg and longitudes_deg are paired two-dimensional cubed-sphere
    ! coordinates, not independent latitude and longitude axes. alt_km contains
    ! the geometric altitude of each requested GEOS model level. IRI itself is
    ! evaluated on its regular hbeg:hend:hstep profile. Electron density and
    ! ion composition are linearly interpolated to each GEOS model altitude.
    !
    ! IRI OUTF units used here:
    !   OUTF(1,*)    electron density [m-3]
    !   OUTF(5:11,*) ion species abundance [%] when JF(22)=.true.
    real, intent(in) :: latitudes_deg(:,:), longitudes_deg(:,:)
    integer, intent(in) :: iyear, doy
    real, intent(in) :: ut_hour
    real, intent(in) :: f107_daily, f107_81day
    real, intent(in) :: alt_km(:,:,:)
    real, intent(in) :: hbeg, hend, hstep

    real, intent(out) :: ne_m3(:,:,:)
    real, intent(out) :: species_fraction(:,:,:,:)

    logical :: jf(50)
    integer :: jmag, mmdd
    real :: outf(20, IRI_MAX_OUTPUT_LEVELS), oarr(100)
    real :: dhour_ut
    real :: value

    integer :: i, j, lev, species
    integer :: im, jm, nlev
    integer :: n_iri_levels

    im = size(alt_km, 1)
    jm = size(alt_km, 2)
    nlev = size(alt_km, 3)

    call validate_iri_inputs( &
         latitudes_deg, longitudes_deg, alt_km, hbeg, hend, hstep, &
         f107_daily, f107_81day, ne_m3, species_fraction)

    ! IRI uses a negative day-of-year value in MMDD to indicate DOY input.
    mmdd = -doy

    ! IRI interprets values above 24 as UT when 25 is added to decimal hour.
    dhour_ut = ut_hour + 25.0

    ! Load the IRI index tables once per process before the first IRI_SUB call.
    ! IRI_SUB expects the caller to initialize these COMMON-block tables.
    call initialize_iri_indices()

    ! Start from the IRI-2020 defaults used by the reference driver.
    jf = .true.
    jf(4:6) = .false.
    jf(23) = .false.
    jf(30) = .false.
    jf(33) = .false.
    jf(34) = .false.
    jf(35) = .false.
    jf(39:40) = .false.
    jf(47) = .false.

    ! Keep ion composition as percent abundance because ion_drag_module uses
    ! electron density plus species fractions to diagnose ion mass density.
    jf(22) = .true.

    ! Supply daily and 81-day F10.7 from the same GEOS-MLT forcing used by MSIS.
    jf(25) = .false.
    jf(32) = .false.

    ! Preserve the current no-storm configuration for the F2 and E regions.
    jf(26) = .false.
    jf(35) = .false.

    ! Geographic coordinates are supplied to IRI.
    jmag = 0

    ! This matches IRI_SUB's internal profile-length calculation for the
    ! positive, ascending altitude grid required by validate_iri_inputs.
    n_iri_levels = int((hend - hbeg) / hstep) + 1
    n_iri_levels = min(n_iri_levels, IRI_MAX_OUTPUT_LEVELS)

    ne_m3 = 0.0
    species_fraction = 0.0

    do j = 1, jm
      do i = 1, im

        ! Leave the pre-zeroed outputs unchanged if a horizontal coordinate
        ! is invalid rather than passing a non-finite value into legacy IRI.
        if (.not. ieee_is_finite(latitudes_deg(i,j)) .or. &
            .not. ieee_is_finite(longitudes_deg(i,j))) cycle

        ! IRI_SUB treats OARR as both input and output. Reset the user-input
        ! values for every column so a previous IRI call cannot alter the
        ! requested F10.7 forcing for the next column.
        outf = 0.0
        oarr = -1.0
        oarr(41) = f107_daily
        oarr(46) = f107_81day

        call IRI_SUB(jf, jmag, latitudes_deg(i,j), longitudes_deg(i,j), &
                     iyear, mmdd, dhour_ut, hbeg, hend, hstep, outf, oarr)

        do lev = 1, nlev
          ! Keep zero output for any invalid altitude. The IonDrag driver also
          ! sanitizes ZLE before reaching this routine, so this is a final guard.
          if (.not. ieee_is_finite(alt_km(i,j,lev))) cycle

          ! Interpolate IRI vertically to the actual GEOS model altitude.
          ! This avoids artificial 10-km steps from nearest-level sampling.
          value = interpolate_iri_profile( &
               outf(1,1:n_iri_levels), n_iri_levels, &
               hbeg, hstep, alt_km(i,j,lev))

          ! IRI documents OUTF(1,*) directly in m-3. Do not apply a cm-3
          ! conversion here.
          if (ieee_is_finite(value) .and. value > 0.0) then
            ne_m3(i,j,lev) = value
          end if

          ! With JF(22)=.true., OUTF(5:11,*) contains percent abundance.
          ! Interpolate each species and convert percent to a dimensionless
          ! fraction for ion_drag_module.
          do species = 1, N_ION_SPECIES
            value = interpolate_iri_profile( &
                 outf(4 + species,1:n_iri_levels), n_iri_levels, &
                 hbeg, hstep, alt_km(i,j,lev))
            if (ieee_is_finite(value) .and. value > 0.0) then
              species_fraction(species,i,j,lev) = value / 100.0
            end if
          end do
        end do
      end do
    end do

  end subroutine get_iri_densities


  real function interpolate_iri_profile(profile, n_levels, hbeg, hstep, altitude)
    ! Linearly interpolate a one-dimensional IRI profile to a GEOS altitude.
    ! Values outside the requested IRI profile are clamped to the nearest
    ! profile endpoint.
    real, intent(in) :: profile(:)
    integer, intent(in) :: n_levels
    real, intent(in) :: hbeg, hstep, altitude

    integer :: lower_idx, upper_idx
    real :: position, weight
    real :: lower_value, upper_value

    interpolate_iri_profile = 0.0

    if (n_levels < 1) return
    if (.not. ieee_is_finite(altitude)) return

    if (altitude <= hbeg) then
      if (ieee_is_finite(profile(1))) then
        interpolate_iri_profile = profile(1)
      end if
      return
    end if

    position = (altitude - hbeg) / hstep

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
      interpolate_iri_profile = &
           (1.0 - weight) * lower_value + weight * upper_value
    else if (ieee_is_finite(lower_value)) then
      interpolate_iri_profile = lower_value
    else if (ieee_is_finite(upper_value)) then
      interpolate_iri_profile = upper_value
    end if

  end function interpolate_iri_profile


  subroutine initialize_iri_indices()
    ! Read the IRI monthly and daily index files once per process.
    ! read_ig_rz populates IG12/Rz12 COMMON-block state from ig_rz.dat.
    ! readapf107 populates Ap/F10.7 COMMON-block state from apf107.dat.
    if (iri_indices_initialized) return

    call read_ig_rz()
    call readapf107()

    iri_indices_initialized = .true.
  end subroutine initialize_iri_indices


  subroutine validate_iri_inputs(latitudes_deg, longitudes_deg, alt_km, &
                                 hbeg, hend, hstep, f107_daily, f107_81day, &
                                 ne_m3, species_fraction)
    ! Validate array shapes and the IRI altitude-profile configuration before
    ! entering the expensive grid-point loop.
    real, intent(in) :: latitudes_deg(:,:), longitudes_deg(:,:)
    real, intent(in) :: alt_km(:,:,:)
    real, intent(in) :: hbeg, hend, hstep
    real, intent(in) :: f107_daily, f107_81day
    real, intent(out) :: ne_m3(:,:,:)
    real, intent(out) :: species_fraction(:,:,:,:)

    integer :: im, jm, nlev
    integer :: n_iri_levels

    im = size(alt_km, 1)
    jm = size(alt_km, 2)
    nlev = size(alt_km, 3)

    if (.not. ieee_is_finite(f107_daily) .or. f107_daily <= 0.0 .or. &
        .not. ieee_is_finite(f107_81day) .or. f107_81day <= 0.0) then
      error stop 'IRI: F10.7 inputs must be finite and positive'
    end if

    if (hstep <= 0.0) then
      error stop 'IRI: IRI_HSTEP must be positive'
    end if

    if (hend < hbeg) then
      error stop 'IRI: IRI_HEND must be greater than or equal to IRI_HBEG'
    end if

    n_iri_levels = int((hend - hbeg) / hstep) + 1
    if (n_iri_levels < 1) then
      error stop 'IRI: invalid altitude profile'
    end if

    if (n_iri_levels > IRI_MAX_OUTPUT_LEVELS) then
      error stop 'IRI: altitude profile exceeds IRI output capacity'
    end if

    if (size(latitudes_deg,1) /= im .or. size(latitudes_deg,2) /= jm) then
      error stop 'IRI: latitude array shape does not match alt_km'
    end if

    if (size(longitudes_deg,1) /= im .or. size(longitudes_deg,2) /= jm) then
      error stop 'IRI: longitude array shape does not match alt_km'
    end if

    if (size(ne_m3,1) /= im .or. size(ne_m3,2) /= jm .or. &
        size(ne_m3,3) /= nlev) then
      error stop 'IRI: ne_m3 shape does not match alt_km'
    end if

    if (size(species_fraction,1) /= N_ION_SPECIES .or. &
        size(species_fraction,2) /= im .or. &
        size(species_fraction,3) /= jm .or. &
        size(species_fraction,4) /= nlev) then
      error stop 'IRI: species_fraction shape does not match expected dimensions'
    end if

  end subroutine validate_iri_inputs

end module iri_input_module
