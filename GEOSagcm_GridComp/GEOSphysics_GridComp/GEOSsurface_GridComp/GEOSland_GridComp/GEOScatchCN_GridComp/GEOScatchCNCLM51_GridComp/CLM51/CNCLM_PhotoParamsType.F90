module CNCLM_PhotoParamsType

  use MAPL_ConstantsMod, ONLY: r8 => MAPL_R4
  use MAPL_ExceptionHandling
  use clm_varctl    , only : use_hydrstress
  use clm_varpar    , only : mxpft, nvegwcs
  use nanMod           , only : nan

  ! !PUBLIC TYPES:
  implicit none
  save

!
! !PUBLIC MEMBER FUNCTIONS:
  public :: init_photo_params_type

  type :: photo_params_type
     real(r8) :: act25  ! Rubisco activity at 25 C (umol CO2/gRubisco/s)
     real(r8) :: fnr  ! Mass ratio of total Rubisco molecular mass to nitrogen in Rubisco (gRubisco/gN in Rubisco)
     real(r8) :: cp25_yr2000  ! CO2 compensation point at 25°C at present day O2 (mol/mol)
     real(r8) :: kc25_coef  ! Michaelis-Menten const. at 25°C for CO2 (unitless)
     real(r8) :: ko25_coef  ! Michaelis-Menten const. at 25°C for O2 (unitless)
     real(r8) :: fnps       ! Fraction of light absorbed by non-photosynthetic pigment (unitless)
     real(r8) :: theta_psii ! Empirical curvature parameter for electron transport rate (unitless)
     real(r8) :: theta_ip   ! Empirical curvature parameter for ap photosynthesis co-limitation (unitless)
     real(r8) :: vcmaxha    ! Activation energy for vcmax (J/mol)
     real(r8) :: jmaxha     ! Activation energy for jmax (J/mol)
     real(r8) :: tpuha      ! Activation energy for tpu (J/mol)
     real(r8) :: lmrha      ! Activation energy for lmr (J/mol)
     real(r8) :: kcha       ! Activation energy for kc (J/mol)
     real(r8) :: koha       ! Activation energy for ko (J/mol)
     real(r8) :: cpha       ! Activation energy for cp (J/mol)
     real(r8) :: vcmaxhd    ! Deactivation energy for vcmax (J/mol)
     real(r8) :: jmaxhd     ! Deactivation energy for jmax (J/mol)
     real(r8) :: tpuhd      ! Deactivation energy for tpu (J/mol)
     real(r8) :: lmrhd      ! Deactivation energy for lmr (J/mol)
     real(r8) :: lmrse      ! Entropy term for lmr (J/mol/K)
     real(r8) :: tpu25ratio ! Ratio of tpu25top to vcmax25top (unitless)
     real(r8) :: kp25ratio  ! Ratio of kp25top to vcmax25top (unitless)
     real(r8) :: vcmaxse_sf ! Scale factor for vcmaxse (unitless)
     real(r8) :: jmaxse_sf  ! Scale factor for jmaxse (unitless)
     real(r8) :: tpuse_sf   ! Scale factor for tpuse (unitless)
     real(r8) :: jmax25top_sf ! Scale factor for jmax25top (unitless)
     real(r8), allocatable, public  :: krmax              (:)
     real(r8), allocatable, private :: kmax               (:,:)
     real(r8), allocatable, private :: psi50              (:,:)
     real(r8), allocatable, private :: ck                 (:,:)
     real(r8), allocatable, private :: lmr_intercept_atkin(:)
     real(r8), allocatable, private :: theta_cj           (:) ! Empirical curvature parameter for ac, aj photosynthesis co-limitation (unitless)

  end type photo_params_type
  type(photo_params_type), public, target, save :: params_inst

contains

!--------------------------------------
 subroutine init_photo_params_type(this)

  ! !DESCRIPTION:
  ! Initialize CTSM photosynthesis parameters needed for calling CTSM routines                                 
  ! jk Oct 2021: type is allocated and initialized to NaN; values are assigned from Catchment states before calls to CLM subroutines are made
  ! this type is only used to be able to pass Catchment states and fluxes to CLM subroutines in the format they expect         
  !                                                                                                                       
  ! !ARGUMENTS:                                                                                                           
    implicit none

    type(photo_params_type), intent(inout):: this

    character(300)                  :: paramfile
    integer                         :: ierr, clm_varid

    real(r8), allocatable, dimension(:)   :: read_tmp_1
    real(r8), allocatable, dimension(:,:) :: read_tmp_2
    real(r8)                              :: read_tmp_3
   !---------------------------------------------------

    allocate( read_tmp_1         (0:mxpft))
    allocate( read_tmp_2         (0:mxpft,nvegwcs))



    allocate( this%krmax       (0:mxpft) )          ; this%krmax(:)        = nan
    allocate( this%theta_cj    (0:mxpft) )          ; this%theta_cj(:)     = nan
    allocate( this%kmax        (0:mxpft,nvegwcs) )  ; this%kmax(:,:)       = nan
    allocate( this%psi50       (0:mxpft,nvegwcs) )  ; this%psi50(:,:)      = nan
    allocate( this%ck          (0:mxpft,nvegwcs) )  ; this%ck(:,:)         = nan

    if ( use_hydrstress .and. nvegwcs /= 4 )then
       _ASSERT(.FALSE.,'Error:: the Plant Hydraulics Stress methodology is for the spacA function is hardcoded for nvegwcs==4')
    end if

    ! jkolassa, Dec 2021: read in parameters from CLM parameter file
    ! TO DO: pass parameter file through rc files rather than hardcoding name here

    paramfile = '/discover/nobackup/jkolassa/CLM/parameter_files/ctsm51_params.c210923.nc'
    ierr = NF90_OPEN(trim(paramfile),NF90_NOWRITE,ncid)
    if (ierr/=0) then
       _ASSERT(.FALSE.,'Error:: the Plant Hydraulics Stress methodology is for the spacA function is hardcoded for nvegwcs==4')
    end if

    ierr = NF90_INQ_VARID(ncid,'krmax',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%krmax(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'lmr_intercept_atkin',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%lmr_intercept_atkin(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'theta_cj',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%theta_cj(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'kmax',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_2)
    this%theta_cj(:,:) = read_tmp_2(0:mxpft,:)

    ierr = NF90_INQ_VARID(ncid,'psi50',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_2)
    this%psi50(:,:) = read_tmp_2(0:mxpft,:)

    ierr = NF90_INQ_VARID(ncid,'ck',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_2)
    this%theta_ck(:,:) = read_tmp_2(0:mxpft,:)

    ierr = NF90_INQ_VARID(ncid,'ko25_coef',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_3)
    this%ko25_coef = read_tmp_3

    ierr = NF90_INQ_VARID(ncid,'kc25_coef',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_3)
    this%kc25_coef = read_tmp_3

    ierr = NF90_INQ_VARID(ncid,'cp25_yr2000',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_3)
    this%cp25_yr2000 = read_tmp_3

    ierr = NF90_INQ_VARID(ncid,'act25',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_3)
    this%act25 = read_tmp_3

    ierr = NF90_INQ_VARID(ncid,'fnr',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_3)
    this%fnr = read_tmp_3

    ierr = NF90_INQ_VARID(ncid,'fnps',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_3)
    this%fnps = read_tmp_3

    ierr = NF90_INQ_VARID(ncid,'theta_psii',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_3)
    this%theta_psii = read_tmp_3

    ierr = NF90_INQ_VARID(ncid,'theta_ip',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_3)
    this%theta_ip = read_tmp_3

    ierr = NF90_INQ_VARID(ncid,'vcmaxha',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_3)
    this%vcmaxha = read_tmp_3

    ierr = NF90_INQ_VARID(ncid,'jmaxha',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_3)
    this%jmaxha = read_tmp_3

    ierr = NF90_INQ_VARID(ncid,'tpuha',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_3)
    this%tpuha = read_tmp_3

    ierr = NF90_INQ_VARID(ncid,'lmrha',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_3)
    this%lmrha = read_tmp_3

    ierr = NF90_INQ_VARID(ncid,'kcha',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_3)
    this%kcha = read_tmp_3

    ierr = NF90_INQ_VARID(ncid,'koha',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_3)
    this%koha = read_tmp_3

    ierr = NF90_INQ_VARID(ncid,'cpha',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_3)
    this%cpha = read_tmp_3

    ierr = NF90_INQ_VARID(ncid,'vcmaxhd',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_3)
    this%vcmaxhd = read_tmp_3

    ierr = NF90_INQ_VARID(ncid,'jmaxhd',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_3)
    this%jmaxhd = read_tmp_3

    ierr = NF90_INQ_VARID(ncid,'tpuhd',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_3)
    this%tpuhd = read_tmp_3

    ierr = NF90_INQ_VARID(ncid,'lmrhd',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_3)
    this%lmrhd = read_tmp_3

    ierr = NF90_INQ_VARID(ncid,'lmrse',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_3)
    this%lmrse = read_tmp_3

    ierr = NF90_INQ_VARID(ncid,'tpu25ratio',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_3)
    this%tpu25ratio = read_tmp_3

    ierr = NF90_INQ_VARID(ncid,'kp25ratio',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_3)
    this%kp25ratio = read_tmp_3

    ierr = NF90_INQ_VARID(ncid,'vcmaxse_sf',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_3)
    this%vcmaxse_sf = read_tmp_3

    ierr = NF90_INQ_VARID(ncid,'jmaxse_sf',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_3)
    this%jmaxse_sf = read_tmp_3

    ierr = NF90_INQ_VARID(ncid,'tpuse_sf',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_3)
    this%tpuse_sf = read_tmp_3

    ierr = NF90_INQ_VARID(ncid,'jmax25top_sf',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_3)
    this%jmax25top_sf = read_tmp_3

    ierr = NF90_CLOSE(ncid)

 end subroutine init_photo_params_type

end module CNCLM_PhotoParamsType
