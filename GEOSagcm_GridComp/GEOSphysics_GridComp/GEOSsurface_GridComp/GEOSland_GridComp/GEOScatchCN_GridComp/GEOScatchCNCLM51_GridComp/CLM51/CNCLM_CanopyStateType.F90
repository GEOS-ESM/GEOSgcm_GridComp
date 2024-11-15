#include "MAPL_Generic.h"

module CanopyStateType

  use MAPL_ConstantsMod, ONLY: r8 => MAPL_R8
  use clm_varpar       , only : nlevcan, nvegwcs, numpft, num_zon, num_veg, &
                                var_col, var_pft
  use clm_varcon       , only : spval
  use nanMod           , only : nan
  use decompMod        , only : bounds_type
  use MAPL_ExceptionHandling

  ! !PUBLIC TYPES:
  implicit none
  save

!
! !PUBLIC MEMBER FUNCTIONS:

  type, public :: canopystate_type

     integer  , pointer :: frac_veg_nosno_patch     (:)   ! patch fraction of vegetation not covered by snow (0 OR 1) [-] 
     integer  , pointer :: frac_veg_nosno_alb_patch (:)   ! patch fraction of vegetation not covered by snow (0 OR 1) [-] 

     real(r8) , pointer :: tlai_patch               (:)   ! patch canopy one-sided leaf area index, no burying by snow
     real(r8) , pointer :: tsai_patch               (:)   ! patch canopy one-sided stem area index, no burying by snow
     real(r8) , pointer :: elai_patch               (:)   ! patch canopy one-sided leaf area index with burying by snow
     real(r8) , pointer :: esai_patch               (:)   ! patch canopy one-sided stem area index with burying by snow
     real(r8) , pointer :: elai240_patch            (:)   ! patch canopy one-sided leaf area index with burying by snow average over 10days 
     real(r8) , pointer :: laisun_patch             (:)   ! patch patch sunlit projected leaf area index  
     real(r8) , pointer :: laisha_patch             (:)   ! patch patch shaded projected leaf area index  
     real(r8) , pointer :: laisun_z_patch           (:,:) ! patch patch sunlit leaf area for canopy layer 
     real(r8) , pointer :: laisha_z_patch           (:,:) ! patch patch shaded leaf area for canopy layer 
     real(r8) , pointer :: mlaidiff_patch           (:)   ! patch difference between lai month one and month two (for dry deposition of chemical tracers)
     real(r8) , pointer :: annlai_patch             (:,:) ! patch 12 months of monthly lai from input data set (for dry deposition of chemical tracers) 
     real(r8) , pointer :: stem_biomass_patch       (:)   ! Aboveground stem biomass (kg/m**2)
     real(r8) , pointer :: leaf_biomass_patch       (:)   ! Aboveground leaf biomass  (kg/m**2)
     real(r8) , pointer :: htop_patch               (:)   ! patch canopy top (m)
     real(r8) , pointer :: hbot_patch               (:)   ! patch canopy bottom (m)
     real(r8) , pointer :: z0m_patch                (:)   ! patch momentum roughness length (m)
     real(r8) , pointer :: displa_patch             (:)   ! patch displacement height (m)
     real(r8) , pointer :: fsun_patch               (:)   ! patch sunlit fraction of canopy         
     real(r8) , pointer :: fsun24_patch             (:)   ! patch 24hr average of sunlit fraction of canopy 
     real(r8) , pointer :: fsun240_patch            (:)   ! patch 240hr average of sunlit fraction of canopy

     real(r8) , pointer :: dleaf_patch              (:)   ! patch characteristic leaf width (diameter) [m]
                                                          ! for non-ED/FATES this is the same as pftcon%dleaf()
     real(r8) , pointer :: rscanopy_patch           (:)   ! patch canopy stomatal resistance (s/m) (ED specific)

     real(r8) , pointer :: vegwp_patch              (:,:) ! patch vegetation water matric potential (mm)
     real(r8) , pointer :: vegwp_ln_patch           (:,:) ! patch vegetation water matric potential at local noon (mm)
     real(r8) , pointer :: vegwp_pd_patch           (:,:) ! patch predawn vegetation water matric potential (mm)

     real(r8)           :: leaf_mr_vcm = spval            ! Scalar constant of leaf respiration with Vcmax

   contains 

    procedure, public :: Init
    procedure, public  :: ReadNML  

  end type canopystate_type
  type(canopystate_type), public, target, save :: canopystate_inst

contains

!--------------------------------------------------------------
  subroutine Init(this, bounds, nch, ityp, fveg, cncol, cnpft, cn5_cold_start, rc)

  ! !DESCRIPTION:
  ! Initialize CTSM canopy state type  needed for calling CTSM routines                                 
  ! jk Oct 2021: type is allocated and initialized to NaN; values are assigned from Catchment states before calls to CLM subroutines are made
  ! this type is only used to be able to pass Catchment states and fluxes to CLM subroutines in the format they expect         
  !                                                                                                                       
  ! !ARGUMENTS:                                                                                                           
    implicit none
    ! INPUT/OUTPUT
    type(bounds_type),                                intent(in) :: bounds
    integer,                                          intent(in) :: nch         ! number of Catchment tiles
    integer, dimension(nch,num_veg,num_zon),          intent(in) :: ityp ! PFT index
    real, dimension(nch,num_veg,num_zon),             intent(in) :: fveg    ! PFT fraction  
    real, dimension(nch,num_zon,var_col),             intent(in) :: cncol         ! column-level restart variable array 
    real, dimension(nch,num_zon,num_veg,var_pft),     intent(in) :: cnpft ! pft-level (patch-level) restart variable array
    logical, optional,                                intent(in) :: cn5_cold_start
    class(canopystate_type)                                      :: this
    integer, optional,                                intent(out) :: rc

    ! LOCAL
    integer :: begp, endp
    integer :: begc, endc
    integer :: begg, endg
    integer :: np, nc, nz, p, nv, nw
    logical :: cold_start = .false.
    !---------------------------------

    begp = bounds%begp ; endp = bounds%endp
    begc = bounds%begc ; endc = bounds%endc
    begg = bounds%begg ; endg = bounds%endg

    ! check whether a cn5_cold_start option was set and change cold_start accordingly
    if (present(cn5_cold_start) .and. (cn5_cold_start.eqv..true.)) then
       cold_start = .true.
    end if

    ! jkolassa: if cold_start is false, check that both CNCOL and CNPFT have the expected size for CNCLM50, else abort 
    if ((cold_start.eqv..false.) .and. ((size(cncol,3).ne.var_col) .or. &
       (size(cnpft,4).ne.var_pft))) then
       _ASSERT(.FALSE.,'option CNCLM50_cold_start = .FALSE. requires a CNCLM50 restart file')
    end if 

    allocate(this%frac_veg_nosno_patch     (begp:endp))           ; this%frac_veg_nosno_patch     (:)   = huge(1)
    allocate(this%frac_veg_nosno_alb_patch (begp:endp))           ; this%frac_veg_nosno_alb_patch (:)   = 0
    allocate(this%tlai_patch               (begp:endp))           ; this%tlai_patch               (:)   = 0.
    allocate(this%tsai_patch               (begp:endp))           ; this%tsai_patch               (:)   = 0.
    allocate(this%elai_patch               (begp:endp))           ; this%elai_patch               (:)   = 0.
    allocate(this%elai240_patch            (begp:endp))           ; this%elai240_patch            (:)   = 0.
    allocate(this%esai_patch               (begp:endp))           ; this%esai_patch               (:)   = 0.
    allocate(this%laisun_patch             (begp:endp))           ; this%laisun_patch             (:)   = 0.
    allocate(this%laisha_patch             (begp:endp))           ; this%laisha_patch             (:)   = 0.
    allocate(this%laisun_z_patch           (begp:endp,1:nlevcan)) ; this%laisun_z_patch           (:,:) = 0.
    allocate(this%laisha_z_patch           (begp:endp,1:nlevcan)) ; this%laisha_z_patch           (:,:) = 0.
    allocate(this%mlaidiff_patch           (begp:endp))           ; this%mlaidiff_patch           (:)   = 0.
    allocate(this%annlai_patch          (12,begp:endp))           ; this%annlai_patch             (:,:) = 0.
    allocate(this%stem_biomass_patch       (begp:endp))           ; this%stem_biomass_patch       (:)   = 0.
    allocate(this%leaf_biomass_patch       (begp:endp))           ; this%leaf_biomass_patch       (:)   = 0.
    allocate(this%htop_patch               (begp:endp))           ; this%htop_patch               (:)   = 0.
    allocate(this%hbot_patch               (begp:endp))           ; this%hbot_patch               (:)   = 0.
    allocate(this%z0m_patch                (begp:endp))           ; this%z0m_patch                (:)   = nan
    allocate(this%displa_patch             (begp:endp))           ; this%displa_patch             (:)   = nan
    allocate(this%fsun_patch               (begp:endp))           ; this%fsun_patch               (:)   = spval
    allocate(this%fsun24_patch             (begp:endp))           ; this%fsun24_patch             (:)   = nan
    allocate(this%fsun240_patch            (begp:endp))           ; this%fsun240_patch            (:)   = nan

    allocate(this%dleaf_patch              (begp:endp))           ; this%dleaf_patch              (:)   = nan
    allocate(this%rscanopy_patch           (begp:endp))           ; this%rscanopy_patch           (:)   = nan
!    allocate(this%gccanopy_patch           (begp:endp))           ; this%gccanopy_patch           (:)   = 0.0_r8     
    allocate(this%vegwp_patch              (begp:endp,1:nvegwcs)) ; this%vegwp_patch              (:,:) = nan
    allocate(this%vegwp_ln_patch           (begp:endp,1:nvegwcs)) ; this%vegwp_ln_patch           (:,:) = nan
    allocate(this%vegwp_pd_patch           (begp:endp,1:nvegwcs)) ; this%vegwp_pd_patch           (:,:) = nan

    ! set parameters to default values or read from parameter file

!    this%leaf_mr_vcm = 0.032 !0.015      ! jkolassa Mar 2022: 0.015 is default value in CTSM5.1, but accoring to ChangeLog 0.032 should be used for Atkin leaf respiration method, which we are using


    ! initialize variables from restart file or set to cold start value

    np = 0
    do nc = 1,nch        ! catchment tile loop
       do nz = 1,num_zon    ! CN zone loop
          do p = 0,numpft  ! PFT index loop
             np = np + 1
             do nv = 1,num_veg ! defined veg loop
                if(ityp(nc,nv,nz)==p .and. fveg(nc,nv,nz)>1.e-4) then

                  ! "old" variables: CNCLM45 and before
                  this%elai_patch  (np) = cnpft(nc,nz,nv, 69)
                  this%esai_patch  (np) = cnpft(nc,nz,nv, 70)
                  this%hbot_patch  (np) = cnpft(nc,nz,nv, 71)
                  this%htop_patch  (np) = cnpft(nc,nz,nv, 72)
                  this%tlai_patch  (np) = cnpft(nc,nz,nv, 73)
                  this%tsai_patch  (np) = cnpft(nc,nz,nv, 74)

                  ! "new" variables: introduced in CNCLM50
                  if (cold_start.eqv..false.) then
                     do nw = 1,nvegwcs
                        this%vegwp_patch(np,nw)    = cnpft(nc,nz,nv, 78+(nw-1))
                     end do
                  elseif (cold_start) then
                     this%vegwp_patch(np,1:nvegwcs)    = -2.5e4_r8
                  else
                     _ASSERT(.FALSE.,'missing CNCLM50_cold_start setting')
                  end if

                  ! jkolassa Mar 2022: these two quantites are computed in Photosynthesis,
                  ! so maybe the do not need to be initialized here
                  this%vegwp_ln_patch(np,1:nvegwcs) = -2.5e4_r8
                  this%vegwp_pd_patch(np,1:nvegwcs) = -2.5e4_r8

                  ! jkolassa May 2022: we do not model vegetation on snow, so the variable below is 1 always
                  this%frac_veg_nosno_patch(np) = 1

                end if ! ityp = p

          end do !nv
       end do ! p
     end do ! nz
  end do ! nc

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine ReadNML( this, NLFilename )
    !
    ! Read in canopy parameter namelist
    !       
    ! USES:
    use shr_mpi_mod   , only : shr_mpi_bcast
    use abortutils    , only : endrun
    use spmdMod       , only : masterproc, mpicom
    use fileutils     , only : getavu, relavu, opnfil
    use shr_nl_mod    , only : shr_nl_find_group_name
    use shr_mpi_mod   , only : shr_mpi_bcast
    use clm_varctl    , only : iulog
    use shr_log_mod   , only : errMsg => shr_log_errMsg
    !
    ! ARGUMENTS:
    implicit none
    class(canopystate_type)      :: this
    character(len=*), intent(IN) :: NLFilename ! Namelist filename
    ! LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file
    real(r8) :: leaf_mr_vcm         ! Scalar of leaf respiration to vcmax
    character(len=32) :: subname = 'CanopyStateType::ReadNML'  ! subroutine name
    !-----------------------------------------------------------------------
    namelist / clm_canopy_inparm / leaf_mr_vcm

    ! ----------------------------------------------------------------------
    ! Read namelist from input namelist filename
    ! ----------------------------------------------------------------------

    if ( masterproc )then

       unitn = getavu()
       write(iulog,*) 'Read in clm_canopy_inparm  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, 'clm_canopy_inparm', status=ierr)
       if (ierr == 0) then
          read(unitn, clm_canopy_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading clm_canopy_inparm namelist"//errmsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg="ERROR finding clm_canopy_inparm namelist"//errmsg(sourcefile, __LINE__))
       end if
       call relavu( unitn )

    end if

    ! Broadcast namelist variables read in
    call shr_mpi_bcast(leaf_mr_vcm, mpicom)
    this%leaf_mr_vcm = leaf_mr_vcm

  end subroutine ReadNML


end module CanopyStateType
