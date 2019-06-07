! This code is part of RRTM for GCM Applications - Parallel (RRTMGP)
!
! Contacts: Robert Pincus and Eli Mlawer
! email:  rrtmgp@aer.com
!
! Copyright 2015-2018,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
!
! The cloud optics class used by RRMTGP needs to be initialized with data stored in a netCDF file.
!    RRTMGP itself doesn't include methods for reading the data so we don't conflict with users'
!    local environment. This module provides a straight-forward implementation of reading the data
!    and calling cloud_optics%load().
! Can load either look-up tables or Pade approximates.
! Modified from mo_load_coefficients by Peter Norris 2019.
!
! -------------------------------------------------------------------------------------------------
module mo_load_cloud_optics
  !
  ! Modules for working with rte and rrtmgp
  !
  use mo_rte_kind,     only: wp
  use mo_cloud_optics, only: ty_cloud_optics
  ! --------------------------------------------------
  use mo_simple_netcdf, only: read_field, read_char_vec, read_logical_vec, var_exists, get_dim_size
  use netcdf
  implicit none
  private
  public :: load_and_init_cloud_optics

contains
  subroutine stop_on_err(msg)
    use iso_fortran_env, only : error_unit
    character(len=*), intent(in) :: msg

    if(msg /= "") then
      write(error_unit, *) msg
      stop
    end if
  end subroutine
  !--------------------------------------------------------------------------------------------------------------------
  ! read cloud optical coefficients of specified type from NetCDF file
  subroutine load_and_init_cloud_optics(cloud_optics, filename, co_type, band_lims_wvn)
    class(ty_cloud_optics),   intent(inout) :: cloud_optics
    character(len=*),         intent(in   ) :: filename
    character(len=*),         intent(in   ) :: co_type       ! 'LUT' or 'Pade'
    real(wp), dimension(:,:), intent(in   ) :: band_lims_wvn ! Spectral discretization

    ! Variables that will be passed to cloud_optics%load()
    !
    ! Lookup table interpolation constants
    ! Lower and upper bounds of the tables and constant for calculating interpn indices
    real(wp) :: radliq_lwr, radliq_upr, radliq_fac
    real(wp) :: radice_lwr, radice_upr, radice_fac
    !
    ! Lookup table coefficients
    ! Extinction, single-scattering albedo, and asymmetry parameter for liquid and ice
    real(wp), dimension(:,:),   allocatable :: lut_extliq, lut_ssaliq, lut_asyliq
    real(wp), dimension(:,:,:), allocatable :: lut_extice, lut_ssaice, lut_asyice
    !
    ! Pade coefficients: extinction, single-scattering albedo, and asymmetry factor
    !
    real(wp), dimension(:,:,:),   allocatable :: pade_extliq, pade_ssaliq, pade_asyliq
    real(wp), dimension(:,:,:,:), allocatable :: pade_extice, pade_ssaice, pade_asyice
    !
    ! Boundaries of size regimes. Liquid and ice are separate.
    !   ext is fit to different numbers of size bins than ssa and asy
    !
    real(wp), dimension(:), allocatable :: pade_sizreg_extliq, pade_sizreg_ssaliq, pade_sizreg_asyliq
    real(wp), dimension(:), allocatable :: pade_sizreg_extice, pade_sizreg_ssaice, pade_sizreg_asyice

    ! Band wavenumber limits on file
    real(wp), dimension(:,:), allocatable :: bnd_limits_wavenumber

    ! dimensions
    integer :: nband, nrghice, nsize_liq, nsize_ice, nsizereg, &
               ncoeff_ext, ncoeff_ssa_g, nbound

    ! locals
    integer :: ncid
    character(len=132) :: error_msg

    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("load_and_init_cloud_optics(): can't open file " // trim(fileName))

    ! get shared dimensions
    nband   = get_dim_size(ncid,'nband')
    nrghice = get_dim_size(ncid,'nrghice')

    ! check that the spectral discretization on file is consistent with band_lims_wvn
    if (nband /= size(band_lims_wvn,2)) &
      call stop_on_err("load_and_init_cloud_optics(): inconsistent nbands on file")
    bnd_limits_wavenumber = read_field(ncid, 'bnd_limits_wavenumber', 2, nband)
    if (.not.all(abs(bnd_limits_wavenumber - band_lims_wvn) < 5._wp * spacing(band_lims_wvn))) &
      call stop_on_err("load_and_init_cloud_optics(): inconsistent bands on file")

    select case (co_type)
      case ("LUT")
        ! read data
        radliq_lwr = read_field(ncid, 'radliq_lwr')
        radliq_upr = read_field(ncid, 'radliq_upr')
        radliq_fac = read_field(ncid, 'radliq_fac')
        radice_lwr = read_field(ncid, 'radice_lwr')
        radice_upr = read_field(ncid, 'radice_upr')
        radice_fac = read_field(ncid, 'radice_fac')
        nsize_liq = get_dim_size(ncid, 'nsize_liq')
        lut_extliq = read_field(ncid, 'lut_extliq', nsize_liq, nband)
        lut_ssaliq = read_field(ncid, 'lut_ssaliq', nsize_liq, nband)
        lut_asyliq = read_field(ncid, 'lut_asyliq', nsize_liq, nband)
        nsize_ice = get_dim_size(ncid, 'nsize_ice')
        lut_extice = read_field(ncid, 'lut_extice', nsize_ice, nband, nrghice)
        lut_ssaice = read_field(ncid, 'lut_ssaice', nsize_ice, nband, nrghice)
        lut_asyice = read_field(ncid, 'lut_asyice', nsize_ice, nband, nrghice)
        ! load cloud_optics
        call stop_on_err(cloud_optics%load( &
          band_lims_wvn,  &
          radliq_lwr, radliq_upr, radliq_fac, &
          radice_lwr, radice_upr, radice_fac, &
          lut_extliq, lut_ssaliq, lut_asyliq, &
          lut_extice, lut_ssaice, lut_asyice))
      case ("Pade")
        ! read data
        nsizereg     = get_dim_size(ncid, 'nsizereg')
        ncoeff_ext   = get_dim_size(ncid, 'ncoeff_ext')
        ncoeff_ssa_g = get_dim_size(ncid, 'ncoeff_ssa_g')
        pade_extliq = read_field(ncid, 'pade_extliq', nband, nsizereg, ncoeff_ext)
        pade_ssaliq = read_field(ncid, 'pade_ssaliq', nband, nsizereg, ncoeff_ssa_g)
        pade_asyliq = read_field(ncid, 'pade_asyliq', nband, nsizereg, ncoeff_ssa_g)
        pade_extice = read_field(ncid, 'pade_extice', nband, nsizereg, ncoeff_ext,   nrghice)
        pade_ssaice = read_field(ncid, 'pade_ssaice', nband, nsizereg, ncoeff_ssa_g, nrghice)
        pade_asyice = read_field(ncid, 'pade_asyice', nband, nsizereg, ncoeff_ssa_g, nrghice)
        nbound = get_dim_size(ncid, 'nbound')
        pade_sizreg_extliq = int(read_field(ncid, 'pade_sizreg_extliq', nbound))
        pade_sizreg_ssaliq = int(read_field(ncid, 'pade_sizreg_ssaliq', nbound))
        pade_sizreg_asyliq = int(read_field(ncid, 'pade_sizreg_asyliq', nbound))
        pade_sizreg_extice = int(read_field(ncid, 'pade_sizreg_extice', nbound))
        pade_sizreg_ssaice = int(read_field(ncid, 'pade_sizreg_ssaice', nbound))
        pade_sizreg_asyice = int(read_field(ncid, 'pade_sizreg_asyice', nbound))
        ! load cloud_optics
        call stop_on_err(cloud_optics%load( &
          band_lims_wvn, &
          pade_extliq, pade_ssaliq, pade_asyliq, &
          pade_extice, pade_ssaice, pade_asyice, &
          pade_sizreg_extliq, pade_sizreg_ssaliq, pade_sizreg_asyliq, &
          pade_sizreg_extice, pade_sizreg_ssaice, pade_sizreg_asyice))
      case default
        call stop_on_err("load_and_init_cloud_optics(): bad co_type")
    end select

    ncid = nf90_close(ncid)

  end subroutine load_and_init_cloud_optics
end module
