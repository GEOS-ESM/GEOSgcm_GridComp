module PatchType

  use MAPL_ConstantsMod, ONLY: r8 => MAPL_R8
  use nanMod           , only : nan
  use decompMod        , only : bounds_type
  use clm_varcon       , only : ispval
  use clm_varctl       , only : use_fates
  use clm_varpar       , only : numpft, NUM_ZON, NUM_VEG, CN_zone_weight

 !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Patch data type allocation 
  ! -------------------------------------------------------- 
  ! patch types can have values of
  ! -------------------------------------------------------- 
  !   0  => not_vegetated
  !   1  => needleleaf_evergreen_temperate_tree
  !   2  => needleleaf_evergreen_boreal_tree
  !   3  => needleleaf_deciduous_boreal_tree
  !   4  => broadleaf_evergreen_tropical_tree
  !   5  => broadleaf_evergreen_temperate_tree
  !   6  => broadleaf_deciduous_tropical_tree
  !   7  => broadleaf_deciduous_temperate_tree
  !   8  => broadleaf_deciduous_boreal_tree
  !   9  => broadleaf_evergreen_shrub
  !   10 => broadleaf_deciduous_temperate_shrub
  !   11 => broadleaf_deciduous_boreal_shrub
  !   12 => c3_arctic_grass
  !   13 => c3_non-arctic_grass
  !   14 => c4_grass
  !   15 => c3_crop

  ! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:

  type, public :: patch_type

     ! g/l/c/p hierarchy, local g/l/c/p cells only
     integer , pointer :: column   (:) ! index into column level quantities
     real(r8), pointer :: wtcol    (:) ! weight (relative to column) 
     integer , pointer :: landunit (:) ! index into landunit level quantities
     real(r8), pointer :: wtlunit  (:) ! weight (relative to landunit) 
     integer , pointer :: gridcell (:) ! index into gridcell level quantities
     real(r8), pointer :: wtgcell  (:) ! weight (relative to gridcell) 

     ! Non-ED only 
     integer , pointer :: itype    (:) ! patch vegetation 
     integer , pointer :: mxy      (:) ! m index for laixy(i,j,m),etc. (undefined for special landunits)
     logical , pointer :: active   (:) ! true=>do computations on this patch

     ! fates only
     logical , pointer :: is_veg   (:) ! This is an ACTIVE fates patch
     logical , pointer :: is_bareground  (:)
     real(r8), pointer :: wt_ed       (:) !TODO mv ? can this be removed


     logical, pointer  :: is_fates (:) ! true for patch vector space reserved
                                       ! for FATES.
                                       ! this is static and is true for all 
                                       ! patches within fates jurisdiction
                                       ! including patches which are not currently
                                       ! associated with a FATES linked-list patch

    contains

     procedure, public :: Init

  end type patch_type
  type(patch_type), public, target :: patch

 contains

!----------------------------------------------------
  subroutine Init(this, bounds, nch, ityp, fveg)

  ! !ARGUMENTS:                                                                                                           
    implicit none

  ! INPUT:
    type(bounds_type),                       intent(in) :: bounds
    integer,                                 intent(in) :: nch    ! number of Catchment tiles
    integer, dimension(nch,num_veg,num_zon), intent(in) :: ityp   ! PFT index
    real, dimension(nch,num_veg,num_zon),    intent(in) :: fveg   ! PFT fraction
    class(patch_type)                                   :: this

  ! LOCAL: 
    integer :: begp,endp
    integer :: np, nc, nz, p, nv, n
  !-------------------------------

    begp = bounds%begp
    endp = bounds%endp

    allocate(this%gridcell      (begp:endp)); this%gridcell   (:) = ispval
    allocate(this%wtgcell       (begp:endp)); this%wtgcell    (:) = nan

    allocate(this%landunit      (begp:endp)); this%landunit   (:) = ispval
    allocate(this%wtlunit       (begp:endp)); this%wtlunit    (:) = nan

    allocate(this%column        (begp:endp)); this%column     (:) = ispval
    allocate(this%wtcol         (begp:endp)); this%wtcol      (:) = nan

    allocate(this%mxy           (begp:endp)); this%mxy        (:) = ispval
    allocate(this%active        (begp:endp)); this%active     (:) = .false.

    ! TODO (MV, 10-17-14): The following must be commented out because
    ! currently the logic checking if patch%itype(p) is not equal to noveg
    ! is used in RootBiogeophysMod in zeng2001_rootfr- a filter is not used
    ! in that routine - which would elimate this problem

    allocate(this%itype      (begp:endp)); this%itype      (:) = ispval

    allocate(this%is_fates   (begp:endp)); this%is_fates   (:) = .false.

    if (use_fates) then
       allocate(this%is_veg  (begp:endp)); this%is_veg  (:) = .false.
       allocate(this%is_bareground (begp:endp)); this%is_bareground (:) = .false.
       allocate(this%wt_ed      (begp:endp)); this%wt_ed      (:) = nan
    end if

   ! initialize values from restart files 

    np = 0
    n = 0
    do nc = 1,nch        ! catchment tile loop
       do nz = 1,num_zon    ! CN zone loop
          n = n + 1
          do p = 0,numpft  ! PFT index loop
             np = np + 1
             this%itype(np) = p
             this%wtcol(np) = 0.
             this%column(np) = n
             this%gridcell(np) = nc
             this%wtgcell(np)  = 0.
             this%landunit(np) = nc
             this%wtlunit(np)  = 0.
             do nv = 1,num_veg ! defined veg loop
                if(ityp(nc,nv,nz)==p .and. fveg(nc,nv,nz)>1.e-4) then
                   this%active(np) = .true.
                   this%wtcol(np) = this%wtcol(np) + fveg(nc,nv,nz)
                   this%wtgcell(np)  = this%wtgcell(np) + (fveg(nc,nv,nz)*CN_zone_weight(nz))
                   this%wtlunit(np)  = this%wtlunit(np) + (fveg(nc,nv,nz)*CN_zone_weight(nz))
                end if
             end do ! nv
          end do ! p
        end do ! nz
     end do ! nc
  end subroutine Init
end module PatchType
