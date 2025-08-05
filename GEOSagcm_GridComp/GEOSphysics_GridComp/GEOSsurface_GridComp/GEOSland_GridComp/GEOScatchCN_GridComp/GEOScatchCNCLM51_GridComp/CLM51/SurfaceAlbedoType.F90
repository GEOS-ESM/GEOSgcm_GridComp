module SurfaceAlbedoType

  use shr_kind_mod     , only : r8 => shr_kind_r8
  use nanMod           , only : nan
  use clm_varpar       , only : numrad, nlevcan, nlevsno, numpft, num_zon, num_veg, &
                                var_col, var_pft
  use clm_varcon       , only : spval, ispval
  use decompMod        , only : bounds_type

  ! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:

 !
  type, public :: surfalb_type

     real(r8), pointer :: coszen_col           (:)   ! col cosine of solar zenith angle
     real(r8), pointer :: albd_patch           (:,:) ! patch surface albedo (direct)   (numrad)                    
     real(r8), pointer :: albi_patch           (:,:) ! patch surface albedo (diffuse)  (numrad)                    
     real(r8), pointer :: albdSF_patch         (:,:) ! patch snow-free surface albedo (direct)   (numrad)
     real(r8), pointer :: albiSF_patch         (:,:) ! patch snow-free surface albedo (diffuse)  (numrad)
     real(r8), pointer :: albgrd_pur_col       (:,:) ! col pure snow ground direct albedo     (numrad)             
     real(r8), pointer :: albgri_pur_col       (:,:) ! col pure snow ground diffuse albedo    (numrad)             
     real(r8), pointer :: albgrd_bc_col        (:,:) ! col ground direct  albedo without BC   (numrad)             
     real(r8), pointer :: albgri_bc_col        (:,:) ! col ground diffuse albedo without BC   (numrad)             
     real(r8), pointer :: albgrd_oc_col        (:,:) ! col ground direct  albedo without OC   (numrad)             
     real(r8), pointer :: albgri_oc_col        (:,:) ! col ground diffuse albedo without OC   (numrad)             
     real(r8), pointer :: albgrd_dst_col       (:,:) ! col ground direct  albedo without dust (numrad)             
     real(r8), pointer :: albgri_dst_col       (:,:) ! col ground diffuse albedo without dust (numrad)             
     real(r8), pointer :: albgrd_col           (:,:) ! col ground albedo (direct)  (numrad)                        
     real(r8), pointer :: albgri_col           (:,:) ! col ground albedo (diffuse) (numrad)                        
     real(r8), pointer :: albsod_col           (:,:) ! col soil albedo: direct  (col,bnd) [frc]                    
     real(r8), pointer :: albsoi_col           (:,:) ! col soil albedo: diffuse (col,bnd) [frc]                    
     real(r8), pointer :: albsnd_hst_col       (:,:) ! col snow albedo, direct , for history files (col,bnd) [frc] 
     real(r8), pointer :: albsni_hst_col       (:,:) ! col snow albedo, diffuse, for history files (col,bnd) [frc] 

     real(r8), pointer :: ftdd_patch           (:,:) ! patch down direct flux below canopy per unit direct flx    (numrad)
     real(r8), pointer :: ftid_patch           (:,:) ! patch down diffuse flux below canopy per unit direct flx   (numrad)
     real(r8), pointer :: ftii_patch           (:,:) ! patch down diffuse flux below canopy per unit diffuse flx  (numrad)
     real(r8), pointer :: fabd_patch           (:,:) ! patch flux absorbed by canopy per unit direct flux         (numrad)
     real(r8), pointer :: fabd_sun_patch       (:,:) ! patch flux absorbed by sunlit canopy per unit direct flux  (numrad)
     real(r8), pointer :: fabd_sha_patch       (:,:) ! patch flux absorbed by shaded canopy per unit direct flux  (numrad)
     real(r8), pointer :: fabi_patch           (:,:) ! patch flux absorbed by canopy per unit diffuse flux        (numrad)
     real(r8), pointer :: fabi_sun_patch       (:,:) ! patch flux absorbed by sunlit canopy per unit diffuse flux (numrad)
     real(r8), pointer :: fabi_sha_patch       (:,:) ! patch flux absorbed by shaded canopy per unit diffuse flux (numrad)
     real(r8), pointer :: fabd_sun_z_patch     (:,:) ! patch absorbed sunlit leaf direct  PAR (per unit lai+sai) for each canopy layer
     real(r8), pointer :: fabd_sha_z_patch     (:,:) ! patch absorbed shaded leaf direct  PAR (per unit lai+sai) for each canopy layer
     real(r8), pointer :: fabi_sun_z_patch     (:,:) ! patch absorbed sunlit leaf diffuse PAR (per unit lai+sai) for each canopy layer
     real(r8), pointer :: fabi_sha_z_patch     (:,:) ! patch absorbed shaded leaf diffuse PAR (per unit lai+sai) for each canopy layer
     real(r8), pointer :: flx_absdv_col        (:,:) ! col absorbed flux per unit incident direct flux:  VIS (col,lyr) [frc]
     real(r8), pointer :: flx_absdn_col        (:,:) ! col absorbed flux per unit incident direct flux:  NIR (col,lyr) [frc]
     real(r8), pointer :: flx_absiv_col        (:,:) ! col absorbed flux per unit incident diffuse flux: VIS (col,lyr) [frc]
     real(r8), pointer :: flx_absin_col        (:,:) ! col absorbed flux per unit incident diffuse flux: NIR (col,lyr) [frc]

     real(r8) , pointer :: fsun_z_patch        (:,:) ! patch patch sunlit fraction of canopy layer
     real(r8) , pointer :: tlai_z_patch        (:,:) ! patch tlai increment for canopy layer                         
     real(r8) , pointer :: tsai_z_patch        (:,:) ! patch tsai increment for canopy layer                         
     integer  , pointer :: ncan_patch          (:)   ! patch number of canopy layers
     integer  , pointer :: nrad_patch          (:)   ! patch number of canopy layers, above snow for radiative transfer
     real(r8) , pointer :: vcmaxcintsun_patch  (:)   ! patch leaf to canopy scaling coefficient, sunlit leaf vcmax   
     real(r8) , pointer :: vcmaxcintsha_patch  (:)   ! patch leaf to canopy scaling coefficient, shaded leaf vcmax   

   contains 

     procedure, public :: Init

end type surfalb_type
type(surfalb_type), public, target, save :: surfalb_inst

contains

!---------------------------------------------------
  subroutine Init(this, bounds, nch, cncol, cnpft)

  ! !DESCRIPTION:
! Initialize CTSM surface albedo needed for calling CTSM routines                                 
! jk Apr 2021: type is allocated and initialized to NaN; values are assigned from Catchment states before calls to CLM subroutines are made
! this type is only used to be able to pass Catchment states and fluxes to CLM subroutines in the format they expect         
!                                                                                                                       
! !ARGUMENTS:                                                                                                           
    implicit none
    !INPUT/OUTPUT
    type(bounds_type),                            intent(in) :: bounds
    integer,                                      intent(in) :: nch    ! number of Catchment tiles
    real, dimension(nch,num_zon,var_col),         intent(in) :: cncol  ! column-level restart variable array 
    real, dimension(nch,num_zon,num_veg,var_pft), intent(in) :: cnpft  ! pft-level (patch-level) restart variable array
    class(surfalb_type)                                      :: this

    ! LOCAL
    integer :: begp, endp
    integer :: begc, endc
    integer :: np, nc, nz, p, nv, n
    !-------------------------------

    begp = bounds%begp ; endp = bounds%endp
    begc = bounds%begc ; endc = bounds%endc    

    allocate(this%coszen_col         (begc:endc))              ; this%coszen_col         (:)   = nan
    allocate(this%albgrd_col         (begc:endc,numrad))       ; this%albgrd_col         (:,:) = nan
    allocate(this%albgri_col         (begc:endc,numrad))       ; this%albgri_col         (:,:) = nan
    allocate(this%albsnd_hst_col     (begc:endc,numrad))       ; this%albsnd_hst_col     (:,:) = spval
    allocate(this%albsni_hst_col     (begc:endc,numrad))       ; this%albsni_hst_col     (:,:) = spval
    allocate(this%albsod_col         (begc:endc,numrad))       ; this%albsod_col         (:,:) = spval
    allocate(this%albsoi_col         (begc:endc,numrad))       ; this%albsoi_col         (:,:) = spval
    allocate(this%albgrd_pur_col     (begc:endc,numrad))       ; this%albgrd_pur_col     (:,:) = nan
    allocate(this%albgri_pur_col     (begc:endc,numrad))       ; this%albgri_pur_col     (:,:) = nan
    allocate(this%albgrd_bc_col      (begc:endc,numrad))       ; this%albgrd_bc_col      (:,:) = nan
    allocate(this%albgri_bc_col      (begc:endc,numrad))       ; this%albgri_bc_col      (:,:) = nan
    allocate(this%albgrd_oc_col      (begc:endc,numrad))       ; this%albgrd_oc_col      (:,:) = nan
    allocate(this%albgri_oc_col      (begc:endc,numrad))       ; this%albgri_oc_col      (:,:) = nan
    allocate(this%albgrd_dst_col     (begc:endc,numrad))       ; this%albgrd_dst_col     (:,:) = nan
    allocate(this%albgri_dst_col     (begc:endc,numrad))       ; this%albgri_dst_col     (:,:) = nan
    allocate(this%albd_patch         (begp:endp,numrad))       ; this%albd_patch         (:,:) = nan
    allocate(this%albi_patch         (begp:endp,numrad))       ; this%albi_patch         (:,:) = nan
    allocate(this%albdSF_patch       (begp:endp,numrad))       ; this%albdSF_patch       (:,:) = nan
    allocate(this%albiSF_patch       (begp:endp,numrad))       ; this%albiSF_patch       (:,:) = nan
    allocate(this%ftdd_patch         (begp:endp,numrad))       ; this%ftdd_patch         (:,:) = nan
    allocate(this%ftid_patch         (begp:endp,numrad))       ; this%ftid_patch         (:,:) = nan
    allocate(this%ftii_patch         (begp:endp,numrad))       ; this%ftii_patch         (:,:) = nan
    allocate(this%fabd_patch         (begp:endp,numrad))       ; this%fabd_patch         (:,:) = nan
    allocate(this%fabd_sun_patch     (begp:endp,numrad))       ; this%fabd_sun_patch     (:,:) = nan
    allocate(this%fabd_sha_patch     (begp:endp,numrad))       ; this%fabd_sha_patch     (:,:) = nan
    allocate(this%fabi_patch         (begp:endp,numrad))       ; this%fabi_patch         (:,:) = nan
    allocate(this%fabi_sun_patch     (begp:endp,numrad))       ; this%fabi_sun_patch     (:,:) = nan
    allocate(this%fabi_sha_patch     (begp:endp,numrad))       ; this%fabi_sha_patch     (:,:) = nan
    allocate(this%fabd_sun_z_patch   (begp:endp,nlevcan))      ; this%fabd_sun_z_patch   (:,:) = 0._r8
    allocate(this%fabd_sha_z_patch   (begp:endp,nlevcan))      ; this%fabd_sha_z_patch   (:,:) = 0._r8
    allocate(this%fabi_sun_z_patch   (begp:endp,nlevcan))      ; this%fabi_sun_z_patch   (:,:) = 0._r8
    allocate(this%fabi_sha_z_patch   (begp:endp,nlevcan))      ; this%fabi_sha_z_patch   (:,:) = 0._r8
    allocate(this%flx_absdv_col      (begc:endc,-nlevsno+1:1)) ; this%flx_absdv_col      (:,:) = spval
    allocate(this%flx_absdn_col      (begc:endc,-nlevsno+1:1)) ; this%flx_absdn_col      (:,:) = spval
    allocate(this%flx_absiv_col      (begc:endc,-nlevsno+1:1)) ; this%flx_absiv_col      (:,:) = spval
    allocate(this%flx_absin_col      (begc:endc,-nlevsno+1:1)) ; this%flx_absin_col      (:,:) = spval

    allocate(this%fsun_z_patch       (begp:endp,nlevcan))      ; this%fsun_z_patch       (:,:) = 0._r8
    allocate(this%tlai_z_patch       (begp:endp,nlevcan))      ; this%tlai_z_patch       (:,:) = 0._r8
    allocate(this%tsai_z_patch       (begp:endp,nlevcan))      ; this%tsai_z_patch       (:,:) = 0._r8
    allocate(this%ncan_patch         (begp:endp))              ; this%ncan_patch         (:)   = 0
    allocate(this%nrad_patch         (begp:endp))              ; this%nrad_patch         (:)   = 0
    allocate(this%vcmaxcintsun_patch (begp:endp))              ; this%vcmaxcintsun_patch (:)   = nan
    allocate(this%vcmaxcintsha_patch (begp:endp))              ; this%vcmaxcintsha_patch (:)   = nan

    ! initialize variables from restart files

    np = 0
    do nc = 1,nch        ! catchment tile loop
       do nz = 1,num_zon    ! CN zone loop
          do p = 0,numpft  ! PFT index loop
             np = np + 1

             this%nrad_patch(np) = 1

             do nv = 1,num_veg ! defined veg loop
                do n = 1,nlevcan
                   this%tlai_z_patch(np,n) = cnpft(nc,nz,nv, 73)
                   if (isnan(this%tlai_z_patch(np,n))) then
                      this%tlai_z_patch(np,n) = 0.
                   end if
                   this%tsai_z_patch(np,n) = cnpft(nc,nz,nv, 74)
                   if (isnan(this%tsai_z_patch(np,n))) then
                      this%tsai_z_patch(np,n) = 0.
                   end if
                   this%vcmaxcintsha_patch(np) = 1._r8
                   this%vcmaxcintsun_patch(np) = 1._r8
                end do
             end do !nv
          end do ! p
       end do ! nz
    end do ! nc

  end subroutine Init

end module SurfaceAlbedoType
