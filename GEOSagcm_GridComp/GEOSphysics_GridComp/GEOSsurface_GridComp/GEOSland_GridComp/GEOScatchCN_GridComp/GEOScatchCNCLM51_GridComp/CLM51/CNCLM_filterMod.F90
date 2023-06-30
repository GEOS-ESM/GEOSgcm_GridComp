module filterMod

  use MAPL_ConstantsMod, ONLY: r8 => MAPL_R8
  use nanMod           , only : nan
  use decompMod  , only : bounds_type
  use clm_varpar , only : NUM_ZON, NUM_VEG, numpft
  use pftconMod  , only : npcropmin

  ! !PUBLIC TYPES:
  implicit none
  save

!
! !PUBLIC MEMBER FUNCTIONS:
  public allocFilters         ! allocate memory for filters
  ! PRIVATE
  private init_filter_type



  type clumpfilter
     integer, pointer :: allc(:)         ! all columns
     integer :: num_allc                 ! number of points in allc filter

     integer, pointer :: natvegp(:)      ! CNDV nat-vegetated (present) filter (pfts)
     integer :: num_natvegp              ! number of pfts in nat-vegetated filter

     integer, pointer :: pcropp(:)       ! prognostic crop filter (pfts)
     integer :: num_pcropp               ! number of pfts in prognostic crop filter
     integer, pointer :: soilnopcropp(:) ! soil w/o prog. crops (pfts)
     integer :: num_soilnopcropp         ! number of pfts in soil w/o prog crops

     integer, pointer :: lakep(:)        ! lake filter (pfts)
     integer :: num_lakep                ! number of pfts in lake filter
     integer, pointer :: nolakep(:)      ! non-lake filter (pfts)
     integer :: num_nolakep              ! number of pfts in non-lake filter
     integer, pointer :: lakec(:)        ! lake filter (columns)
     integer :: num_lakec                ! number of columns in lake filter
     integer, pointer :: nolakec(:)      ! non-lake filter (columns)
     integer :: num_nolakec              ! number of columns in non-lake filter

     integer, pointer :: soilc(:)        ! soil filter (columns)
     integer :: num_soilc                ! number of columns in soil filter 
     integer, pointer :: soilp(:)        ! soil filter (pfts)
     integer :: num_soilp                ! number of pfts in soil filter 

     integer, pointer :: snowc(:)        ! snow filter (columns) 
     integer :: num_snowc                ! number of columns in snow filter 
     integer, pointer :: nosnowc(:)      ! non-snow filter (columns) 
     integer :: num_nosnowc              ! number of columns in non-snow filter 

     integer, pointer :: lakesnowc(:)    ! snow filter (columns) 
     integer :: num_lakesnowc            ! number of columns in snow filter 
     integer, pointer :: lakenosnowc(:)  ! non-snow filter (columns) 
     integer :: num_lakenosnowc          ! number of columns in non-snow filter 

     integer, pointer :: exposedvegp(:)  ! patches where frac_veg_nosno is non-zero
     integer :: num_exposedvegp          ! number of patches in exposedvegp filter
     integer, pointer :: noexposedvegp(:)! patches where frac_veg_nosno is 0 (does NOT include lake or urban)
     integer :: num_noexposedvegp        ! number of patches in noexposedvegp filter

     integer, pointer :: hydrologyc(:)   ! hydrology filter (columns)
     integer :: num_hydrologyc           ! number of columns in hydrology filter 

     integer, pointer :: urbanl(:)       ! urban filter (landunits)
     integer :: num_urbanl               ! number of landunits in urban filter 

     integer, pointer :: nourbanl(:)     ! non-urban filter (landunits)
     integer :: num_nourbanl             ! number of landunits in non-urban filter 

     integer, pointer :: urbanc(:)       ! urban filter (columns)
     integer :: num_urbanc               ! number of columns in urban filter
     integer, pointer :: nourbanc(:)     ! non-urban filter (columns)
     integer :: num_nourbanc             ! number of columns in non-urban filter

     integer, pointer :: urbanp(:)       ! urban filter (pfts)
     integer :: num_urbanp               ! number of pfts in urban filter
     integer, pointer :: nourbanp(:)     ! non-urban filter (pfts)
     integer :: num_nourbanp             ! number of pfts in non-urban filter

     integer, pointer :: nolakeurbanp(:) ! non-lake, non-urban filter (pfts)
     integer :: num_nolakeurbanp         ! number of pfts in non-lake, non-urban filter

     integer, pointer :: icemecc(:)      ! glacier mec filter (cols)
     integer :: num_icemecc              ! number of columns in glacier mec filter

     integer, pointer :: do_smb_c(:)     ! glacier+bareland SMB calculations-on filter (cols)
     integer :: num_do_smb_c             ! number of columns in glacier+bareland SMB mec filter         

     integer, pointer :: actfirec(:)     ! glacier+bareland SMB calculations-on filter (cols)
     integer :: num_actfirec             ! number of columns in glacier+bareland SMB mec filter         

     integer, pointer :: actfirep(:)     ! glacier+bareland SMB calculations-on filter (cols)
     integer :: num_actfirep             ! number of columns in glacier+bareland SMB mec filter         

  end type clumpfilter
  public clumpfilter

  ! This is the standard set of filters, which should be used in most places in the code.
  ! These filters only include 'active' points.
    type(clumpfilter), allocatable, public :: filter(:)

contains

  !------------------------------------------------------------------------
  subroutine allocFilters(bounds, nch, ityp, fveg)
    !
    ! !DESCRIPTION:
    ! Allocate CLM filters.
    !
    ! !REVISION HISTORY:
    ! Created by Bill Sacks

  ! !ARGUMENTS:                                                                                                           
    implicit none
    ! INPUT/OUTPUT
    type(bounds_type),                                intent(in) :: bounds
    integer,                                          intent(in) :: nch         ! number of Catchment tiles
    integer, dimension(nch,num_veg,num_zon),          intent(in) :: ityp ! PFT index
    real, dimension(nch,num_veg,num_zon),             intent(in) :: fveg    ! PFT fraction 

    !------------------------------------------------------------------------

    call init_filter_type(bounds, nch, ityp, fveg,  filter)
    

  end subroutine allocFilters

!--------------------------------------------------------------
  subroutine init_filter_type(bounds, nch, ityp, fveg,  this_filter)

  ! !DESCRIPTION:
  ! Initialize CTSM filters                                 
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
    type(clumpfilter), intent(inout), allocatable :: this_filter(:)  ! the filter to allocate
 
    ! LOCAL:
    integer :: n, nc ,nz, p, np, nv

    !--------------------------------------

       if( .not. allocated(this_filter)) then
          allocate(this_filter(1))
       end if

       allocate(this_filter(1)%allc(bounds%endc-bounds%begc+1))

       allocate(this_filter(1)%lakep(bounds%endp-bounds%begp+1))
       allocate(this_filter(1)%nolakep(bounds%endp-bounds%begp+1))
       allocate(this_filter(1)%nolakeurbanp(bounds%endp-bounds%begp+1))

       allocate(this_filter(1)%lakec(bounds%endc-bounds%begc+1))
       allocate(this_filter(1)%nolakec(bounds%endc-bounds%begc+1))

       allocate(this_filter(1)%soilc(bounds%endc-bounds%begc+1))
       allocate(this_filter(1)%soilp(bounds%endp-bounds%begp+1))

       allocate(this_filter(1)%snowc(bounds%endc-bounds%begc+1))
       allocate(this_filter(1)%nosnowc(bounds%endc-bounds%begc+1))

       allocate(this_filter(1)%lakesnowc(bounds%endc-bounds%begc+1))
       allocate(this_filter(1)%lakenosnowc(bounds%endc-bounds%begc+1))

       allocate(this_filter(1)%exposedvegp(bounds%endp-bounds%begp+1))
       allocate(this_filter(1)%noexposedvegp(bounds%endp-bounds%begp+1))

       allocate(this_filter(1)%natvegp(bounds%endp-bounds%begp+1))

       allocate(this_filter(1)%hydrologyc(bounds%endc-bounds%begc+1))

       allocate(this_filter(1)%urbanp(bounds%endp-bounds%begp+1))
       allocate(this_filter(1)%nourbanp(bounds%endp-bounds%begp+1))

       allocate(this_filter(1)%urbanc(bounds%endc-bounds%begc+1))
       allocate(this_filter(1)%nourbanc(bounds%endc-bounds%begc+1))

       allocate(this_filter(1)%urbanl(bounds%endl-bounds%begl+1))
       allocate(this_filter(1)%nourbanl(bounds%endl-bounds%begl+1))

       allocate(this_filter(1)%pcropp(bounds%endp-bounds%begp+1))
       allocate(this_filter(1)%soilnopcropp(bounds%endp-bounds%begp+1))

       allocate(this_filter(1)%icemecc(bounds%endc-bounds%begc+1))
       allocate(this_filter(1)%do_smb_c(bounds%endc-bounds%begc+1))

       allocate(this_filter(1)%actfirec(bounds%endc-bounds%begc+1))
       allocate(this_filter(1)%actfirep(bounds%endp-bounds%begp+1))

       this_filter(1)%num_actfirep = 1
       this_filter(1)%num_actfirec = 1
      
      ! initialize
 
      this_filter(1)%num_soilc = 0
      this_filter(1)%num_soilp = 0
      this_filter(1)%num_pcropp = 0 
      this_filter(1)%num_exposedvegp = 0
      this_filter(1)%num_noexposedvegp = 0
      this_filter(1)%num_nourbanp = 0
      this_filter(1)%num_allc = 0      

      n = 0
      np = 0
      do nc = 1,nch
         do nz = 1,num_zon
            n = n + 1

            this_filter(1)%num_soilc = this_filter(1)%num_soilc + 1
            this_filter(1)%soilc(this_filter(1)%num_soilc) = n 
            this_filter(1)%num_allc = this_filter(1)%num_allc + 1
            this_filter(1)%allc(this_filter(1)%num_allc) = n

            do p = 0,numpft  ! PFT index loop
               np = np + 1
               do nv = 1,num_veg ! defined veg loop
                  if(ityp(nc,nv,nz)==p) then

                    if (fveg(nc,nv,nz)>1.e-4) then
                    
                    this_filter(1)%num_nourbanp = this_filter(1)%num_nourbanp + 1
                    this_filter(1)%nourbanp(this_filter(1)%num_nourbanp) = np

                    this_filter(1)%num_soilp = this_filter(1)%num_soilp + 1
                    this_filter(1)%soilp(this_filter(1)%num_soilp) = np

                    ! jkolassa: not sure this is needed, since we do not use prognostic crop information
                    if(ityp(nc,nv,nz) >= npcropmin) then
                      this_filter(1)%num_pcropp = this_filter(1)%num_pcropp + 1
                      this_filter(1)%pcropp(this_filter(1)%num_pcropp) = np
                    endif


!                    if (fveg(nc,nv,nz)>1.e-4) then

                       this_filter(1)%num_exposedvegp = this_filter(1)%num_exposedvegp + 1
                       this_filter(1)%exposedvegp(this_filter(1)%num_exposedvegp) = np

                    elseif (fveg(nc,nv,nz)<=1.e-4) then
            
                       this_filter(1)%num_noexposedvegp = this_filter(1)%num_noexposedvegp + 1
                       this_filter(1)%noexposedvegp(this_filter(1)%num_noexposedvegp) = np
                       
                    end if
                  end if
               end do ! nv
            end do !p         
         end do !nz
      end do !nc

  end subroutine init_filter_type
end module filterMod
