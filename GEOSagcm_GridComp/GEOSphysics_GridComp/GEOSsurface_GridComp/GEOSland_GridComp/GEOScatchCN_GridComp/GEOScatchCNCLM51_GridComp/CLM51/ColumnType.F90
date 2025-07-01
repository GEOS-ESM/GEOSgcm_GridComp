module ColumnType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Column data type allocation and initialization
  ! -------------------------------------------------------- 
  ! column types can have values of
  ! -------------------------------------------------------- 
  !   1  => (istsoil)          soil (vegetated or bare soil)
  !   2  => (istcrop)          crop (only for crop configuration)
  !   3  => (UNUSED)           (formerly non-multiple elevation class land ice; currently unused)
  !   4  => (istice_mec)       land ice (multiple elevation classes)   
  !   5  => (istdlak)          deep lake
  !   6  => (istwet)           wetland
  !   71 => (icol_roof)        urban roof
  !   72 => (icol_sunwall)     urban sunwall
  !   73 => (icol_shadewall)   urban shadewall
  !   74 => (icol_road_imperv) urban impervious road
  !   75 => (icol_road_perv)   urban pervious road


  use MAPL_ConstantsMod, ONLY: r8 => MAPL_R8
  use nanMod           , only : nan
  use decompMod       , only : bounds_type
  use clm_varcon     , only : zsoi, dzsoi, zisoi, dzsoi_decomp, spval, ispval
  use clm_varctl     , only : use_fates
  use clm_varpar     , only : nlevsno, nlevgrnd, nlevlak, nlevmaxurbgrnd,nlevurb, &
                              CN_zone_weight, numpft, num_zon


  ! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:

  type, public :: column_type
     ! g/l/c/p hierarchy, local g/l/c/p cells only
     integer , pointer :: landunit             (:)   ! index into landunit level quantities
     real(r8), pointer :: wtlunit              (:)   ! weight (relative to landunit)
     integer , pointer :: gridcell             (:)   ! index into gridcell level quantities
     real(r8), pointer :: wtgcell              (:)   ! weight (relative to gridcell)
     integer , pointer :: patchi               (:)   ! beginning patch index for each column
     integer , pointer :: patchf               (:)   ! ending patch index for each column
     integer , pointer :: npatches             (:)   ! number of patches for each column

     ! topological mapping functionality
     integer , pointer :: itype                (:)   ! column type (after init, should only be modified via update_itype routine)
     integer , pointer :: lun_itype            (:)   ! landunit type (col%lun_itype(ci) is the same as lun%itype(col%landunit(ci)), but is often a more convenient way to access this type
     logical , pointer :: active               (:)   ! true=>do computations on this column
     logical , pointer :: type_is_dynamic      (:)   ! true=>itype can change throughout the run

     ! topography
     ! TODO(wjs, 2016-04-05) Probably move these things into topoMod
     real(r8), pointer :: micro_sigma          (:)   ! microtopography pdf sigma (m)
     real(r8), pointer :: topo_slope           (:)   ! gridcell topographic slope
     real(r8), pointer :: topo_std             (:)   ! gridcell elevation standard deviation

     ! vertical levels
     integer , pointer :: snl                  (:)   ! number of snow layers
     real(r8), pointer :: dz                   (:,:) ! layer thickness (m)  (-nlevsno+1:nlevgrnd) 
     real(r8), pointer :: z                    (:,:) ! layer depth (m) (-nlevsno+1:nlevgrnd) 
     real(r8), pointer :: zi                   (:,:) ! interface level below a "z" level (m) (-nlevsno+0:nlevgrnd) 
     real(r8), pointer :: zii                  (:)   ! convective boundary height [m]
     real(r8), pointer :: dz_lake              (:,:) ! lake layer thickness (m)  (1:nlevlak)
     real(r8), pointer :: z_lake               (:,:) ! layer depth for lake (m)
     real(r8), pointer :: lakedepth            (:)   ! variable lake depth (m)                             
     integer , pointer :: nbedrock             (:)   ! variable depth to bedrock index

     ! other column characteristics
     logical , pointer :: hydrologically_active(:)   ! true if this column is a hydrologically active type
     logical , pointer :: urbpoi               (:)   ! true=>urban point


     ! levgrnd_class gives the class in which each layer falls. This is relevant for
     ! columns where there are 2 or more fundamentally different layer types. For
     ! example, this distinguishes between soil and bedrock layers. The particular value
     ! assigned to each class is irrelevant; the important thing is that different
     ! classes (e.g., soil vs. bedrock) have different values of levgrnd_class.
     !
     ! levgrnd_class = ispval indicates that the given layer is completely unused for
     ! this column (i.e., this column doesn't use the full nlevgrnd layers).
     integer , pointer :: levgrnd_class        (:,:) ! class in which each layer falls (1:nlevgrnd)

   contains

     procedure, public :: Init

  end type column_type
  type(column_type), public, target :: col

 contains

!-----------------------------------------------------
  subroutine Init(this, bounds,nch)

  ! !ARGUMENTS:                                                                                                           
    implicit none

  ! INPUT:
    type(bounds_type), intent(in) :: bounds
    integer,           intent(in) :: nch         ! number of Catchment tiles
    class(column_type)            :: this

  ! LOCAL:

    integer :: begc, endc
    integer :: nc, nz, n, c
  !----------------------------

  begc = bounds%begc ; endc = bounds%endc

      ! The following is set in initGridCellsMod
    allocate(this%gridcell    (begc:endc))                     ; this%gridcell    (:)   = ispval
    allocate(this%wtgcell     (begc:endc))                     ; this%wtgcell     (:)   = nan
    allocate(this%landunit    (begc:endc))                     ; this%landunit    (:)   = ispval
    allocate(this%wtlunit     (begc:endc))                     ; this%wtlunit     (:)   = nan
    allocate(this%patchi      (begc:endc))                     ; this%patchi      (:)   = ispval
    allocate(this%patchf      (begc:endc))                     ; this%patchf      (:)   = ispval
    allocate(this%npatches     (begc:endc))                    ; this%npatches     (:)   = ispval
    allocate(this%itype       (begc:endc))                     ; this%itype       (:)   = ispval
    allocate(this%lun_itype   (begc:endc))                     ; this%lun_itype   (:)   = ispval
    allocate(this%active      (begc:endc))                     ; this%active      (:)   = .false.
    allocate(this%type_is_dynamic(begc:endc))                  ; this%type_is_dynamic(:) = .false.

    ! The following is set in initVerticalMod
    allocate(this%snl         (begc:endc))                     ; this%snl         (:)   = ispval  !* cannot be averaged up
    allocate(this%dz          (begc:endc,-nlevsno+1:nlevmaxurbgrnd)) ; this%dz          (:,:) = nan
    allocate(this%z           (begc:endc,-nlevsno+1:nlevmaxurbgrnd)) ; this%z           (:,:) = nan
    allocate(this%zi          (begc:endc,-nlevsno+0:nlevmaxurbgrnd)) ; this%zi          (:,:) = nan
    allocate(this%zii         (begc:endc))                     ; this%zii         (:)   = nan
    allocate(this%lakedepth   (begc:endc))                     ; this%lakedepth   (:)   = spval
    allocate(this%dz_lake     (begc:endc,nlevlak))             ; this%dz_lake     (:,:) = nan
    allocate(this%z_lake      (begc:endc,nlevlak))             ; this%z_lake      (:,:) = nan

    allocate(this%nbedrock   (begc:endc))                      ; this%nbedrock   (:)   = ispval
    allocate(this%levgrnd_class(begc:endc,nlevmaxurbgrnd))     ; this%levgrnd_class(:,:) = ispval
    allocate(this%micro_sigma (begc:endc))                     ; this%micro_sigma (:)   = nan
    allocate(this%topo_slope  (begc:endc))                     ; this%topo_slope  (:)   = nan
    allocate(this%topo_std    (begc:endc))                     ; this%topo_std    (:)   = nan

    allocate(this%hydrologically_active(begc:endc))            ; this%hydrologically_active(:) = .false.
    allocate(this%urbpoi      (begc:endc))                     ; this%urbpoi      (:)   = .false.


    this%nbedrock(:) = 1  !jkolassa: set this to 1, since we only have one soil layer

     do c = bounds%begc,bounds%endc
        this%z(c,1:nlevgrnd)  = zsoi(1:nlevgrnd)
        this%zi(c,0:nlevgrnd) = zisoi(0:nlevgrnd)
        this%dz(c,1:nlevgrnd) = dzsoi(1:nlevgrnd)
        if (nlevgrnd < nlevurb) then
           this%z(c,nlevgrnd+1:nlevurb)  = spval
           this%zi(c,nlevgrnd+1:nlevurb) = spval
           this%dz(c,nlevgrnd+1:nlevurb) = spval
        end if
     end do


  
    n = 0
    do nc = 1,nch        ! catchment tile loop
       do nz = 1,num_zon    ! CN zone loop
          n = n + 1
          this%active(n)   = .true.
          this%gridcell(n) = nc
          this%wtgcell(n)  = CN_zone_weight(nz)
          this%landunit(n) = nc
          this%wtlunit(n)  = CN_zone_weight(nz)
          this%patchi(n)   = (numpft+1)*(n-1) + 1
          this%patchf(n)   = (numpft+1)*n
        end do ! nz
     end do ! nc

 end subroutine Init
end module ColumnType
