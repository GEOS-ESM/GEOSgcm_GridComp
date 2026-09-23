module iri_input_module
  implicit none

  integer, parameter :: N_ION_SPECIES = 7

  ! --- Explicit interface for external IRI routine (F77, irisub.for) ---
  interface
    subroutine IRI_SUB(JF,JMAG,ALATI,ALONG,IYYYY,MMDD,DHOUR, &
                        HEIBEG,HEIEND,HEISTP,OUTF,OARR)
      logical, intent(in) :: JF(50)
      integer, intent(in) :: JMAG
      real, intent(in) :: ALATI, ALONG
      integer, intent(in) :: IYYYY, MMDD
      real, intent(in) :: DHOUR, HEIBEG, HEIEND, HEISTP
      real, intent(out) :: OUTF(20,1000)
      real, intent(inout) :: OARR(100)
    end subroutine IRI_SUB
  end interface

contains

  subroutine get_iri_densities(latitudes, longitudes, iyear, doy, ut_hour, &
                                f107_daily, f107_81day, &
                                alt_km, hbeg, hend, hstep, &
                                ne_m3, species_fraction)
    ! alt_km is now per-column, per-level (nlat, nlon, nlev), reflecting
    ! real geometric altitude (e.g. derived from GZ) rather than a single
    ! shared profile. hbeg/hend/hstep still define the fixed range/step
    ! IRI_SUB is asked to compute internally (its own outf(:, alt_idx)
    ! profile spacing); alt_km(k,j,lev) is used only to pick which index
    ! into that profile corresponds to the real altitude of level lev at
    ! that column.
    real, intent(in) :: latitudes(:), longitudes(:)
    integer, intent(in) :: iyear, doy
    real, intent(in) :: ut_hour
    real, intent(in) :: f107_daily, f107_81day
    real, intent(in) :: alt_km(:,:,:)          ! (nlat, nlon, nlev)
    real, intent(in) :: hbeg, hend, hstep

    real, intent(out) :: ne_m3(:,:,:)
    real, intent(out) :: species_fraction(:,:,:,:)

    logical :: jf(50)
    integer :: jmag, mmdd
    real :: outf(20, 1000), oarr(100)
    real :: dhour_ut

    integer :: k, j, lev, nlat, nlon, nlev, alt_idx
    real, parameter :: CM3_TO_M3 = 1.0e6

    nlat = size(latitudes)
    nlon = size(longitudes)
    nlev = size(alt_km, 3)

    mmdd = -doy
    dhour_ut = ut_hour + 25.0

    jf = .true.
    jf(25) = .false.
    jf(32) = .false.
    jf(26) = .false.
    jf(35) = .false.
    jmag = 0

    oarr = -1.0
    oarr(41) = f107_daily
    oarr(46) = f107_81day

    do j = 1, nlon
      do k = 1, nlat
        call IRI_SUB(jf, jmag, latitudes(k), longitudes(j), iyear, mmdd, &
                     dhour_ut, hbeg, hend, hstep, outf, oarr)

        do lev = 1, nlev
          ! Real per-column altitude now indexes into IRI's fixed-step
          ! profile, rather than a shared alt_km(:) applying to every column.
          alt_idx = nint((alt_km(k, j, lev) - hbeg) / hstep) + 1
          alt_idx = max(1, min(alt_idx, size(outf, 2)))

          ne_m3(k, j, lev) = outf(1, alt_idx) * CM3_TO_M3
          species_fraction(1:N_ION_SPECIES, k, j, lev) = &
               max(outf(5:5+N_ION_SPECIES-1, alt_idx), 0.0) / 100.0
        end do
      end do
    end do

  end subroutine get_iri_densities

end module iri_input_module
