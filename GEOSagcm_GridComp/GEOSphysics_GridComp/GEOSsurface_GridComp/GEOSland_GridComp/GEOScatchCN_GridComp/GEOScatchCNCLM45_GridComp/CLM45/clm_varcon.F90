module clm_varcon

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clm_varcon
!
! !DESCRIPTION:
! Module containing various model constants
!
! !USES:
  use MAPL_ConstantsMod, ONLY: r8 => MAPL_R4
  use shr_const_mod, only: SHR_CONST_G,     &
                           SHR_CONST_RHOFW, &
                           SHR_CONST_TKFRZ, &
                           SHR_CONST_CDAY,  &
                           SHR_CONST_RGAS,  &
                           SHR_CONST_PI,    &
                           SHR_CONST_PDB
  use clm_varpar   , only: nlevgrnd, nlevdecomp_full
  use clm_varpar   , only: ngases
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! 27 February 2008: Keith Oleson; Add forcing height and aerodynamic parameters
!
!EOP
!-----------------------------------------------------------------------

  !------------------------------------------------------------------
  ! Initialize mathmatical constants
  !------------------------------------------------------------------

  real(r8) :: rpi    = SHR_CONST_PI

  !------------------------------------------------------------------
  ! Initialize physical constants
  !------------------------------------------------------------------

  real(r8) :: grav   = SHR_CONST_G      !gravity constant [m/s2]
  real(r8) :: denh2o = SHR_CONST_RHOFW  !density of liquid water [kg/m3]
  real(r8) :: rgas   = SHR_CONST_RGAS   !universal gas constant [J/K/kmole]
  real(r8) :: tfrz   = SHR_CONST_TKFRZ  !freezing temperature [K]
  real(r8), public, parameter ::  secspday= SHR_CONST_CDAY  ! Seconds per day
  real(r8), public, parameter ::  spval = 1.e36_r8  ! special value for real data
  integer , public, parameter :: ispval = -9999     ! special value for int data

  !!! C13
  real(r8), parameter :: preind_atm_del13c = -6.0   ! preindustrial value for atmospheric del13C
  real(r8), parameter :: preind_atm_ratio = SHR_CONST_PDB + (preind_atm_del13c * SHR_CONST_PDB)/1000.0  ! 13C/12C
  real(r8) :: c13ratio = preind_atm_ratio/(1.0+preind_atm_ratio) ! 13C/(12+13)C preind atmosphere

  !!! C14
  real(r8) :: c14ratio = 1.e-12_r8
  ! real(r8) :: c14ratio = 1._r8  ! debug lets set to 1 to try to avoid numerical errors
 
  integer, private :: i  ! loop index


#ifdef NITRIF_DENITRIF
 !  real(r8), parameter :: nitrif_n2o_loss_frac = 0.02_r8   !fraction of N lost as N2O in nitrification (Parton et al., 2001)
  real(r8), parameter :: nitrif_n2o_loss_frac = 6.e-4_r8   !fraction of N lost as N2O in nitrification (Li et al., 2000)
  real(r8), parameter :: frac_minrlztn_to_no3 = 0.2_r8   !fraction of N mineralized that is dieverted to the nitrification stream (Parton et al., 2001)
#endif


  !------------------------------------------------------------------
  ! Initialize water type constants
  !------------------------------------------------------------------

  ! "land unit " types
  !   1     soil (includes vegetated landunits)
  !   2     land ice (glacier)
  !   3     deep lake
  !  (DEPRECATED: New lake model has variable depth) 4     shallow lake
  !   5     wetland (swamp, marsh, etc.)
  !   6     urban
  !   7     land ice (glacier) with multiple elevation classes
  !   8     crop

  integer, parameter :: istsoil    = 1  !soil         landunit type
  integer, parameter :: istice     = 2  !land ice     landunit type
  integer, parameter :: istdlak    = 3  !deep lake    landunit type
  ! Not used; now 3 is used for all lakes, which have variable depth.
  integer, parameter :: istslak    = 4  !shallow lake landunit type
  integer, parameter :: istwet     = 5  !wetland      landunit type
  integer, parameter :: isturb     = 6  !urban        landunit type
  integer, parameter :: istice_mec = 7  !land ice (multiple elevation classes) landunit type
  integer, parameter :: istcrop    = 8  !crop         landunit type
  integer, parameter :: max_lunit  = 8  !maximum value that lun%itype can have
                             !(i.e., largest value in the above list)

  ! urban column types

  integer, parameter :: icol_roof        = 61
  integer, parameter :: icol_sunwall     = 62
  integer, parameter :: icol_shadewall   = 63
  integer, parameter :: icol_road_imperv = 64
  integer, parameter :: icol_road_perv   = 65

  !------------------------------------------------------------------
  ! Initialize miscellaneous radiation constants
  !------------------------------------------------------------------

  ! Lake Model Constants will be defined in SLakeCon.

  !------------------------------------------------------------------
  ! Soil depths are constants for now; lake depths can vary by gridcell
  ! zlak and dzlak correspond to the default 50 m lake depth.
  ! The values for the following arrays are set in routine iniTimeConst
  !------------------------------------------------------------------

  real(r8), pointer :: zsoi(:)         !soil z  (layers)
  real(r8), pointer :: dzsoi(:)        !soil dz (thickness)
  real(r8), pointer :: zisoi(:)        !soil zi (interfaces)
  real(r8), pointer :: dzsoi_decomp(:) !soil dz (thickness)

  !------------------------------------------------------------------
  ! (Non-tunable) Constants for the CH4 submodel (Tuneable constants in ch4varcon)
  !------------------------------------------------------------------
  ! Note some of these constants are also used in CNNitrifDenitrifMod

  real(r8), parameter :: catomw = 12.011_r8 ! molar mass of C atoms (g/mol)

  real(r8) :: s_con(ngases,4)    ! Schmidt # calculation constants (spp, #)
  data (s_con(1,i),i=1,4) /1898_r8, -110.1_r8, 2.834_r8, -0.02791_r8/ ! CH4
  data (s_con(2,i),i=1,4) /1801_r8, -120.1_r8, 3.7818_r8, -0.047608_r8/ ! O2
  data (s_con(3,i),i=1,4) /1911_r8, -113.7_r8, 2.967_r8, -0.02943_r8/ ! CO2

  real(r8) :: d_con_w(ngases,3)    ! water diffusivity constants (spp, #)  (mult. by 10^-4)
  data (d_con_w(1,i),i=1,3) /0.9798_r8, 0.02986_r8, 0.0004381_r8/ ! CH4
  data (d_con_w(2,i),i=1,3) /1.172_r8, 0.03443_r8, 0.0005048_r8/ ! O2
  data (d_con_w(3,i),i=1,3) /0.939_r8, 0.02671_r8, 0.0004095_r8/ ! CO2

  real(r8) :: d_con_g(ngases,2)    ! gas diffusivity constants (spp, #) (cm^2/s) (mult. by 10^-9)
  data (d_con_g(1,i),i=1,2) /0.1875_r8, 0.0013_r8/ ! CH4
  data (d_con_g(2,i),i=1,2) /0.1759_r8, 0.00117_r8/ ! O2
  data (d_con_g(3,i),i=1,2) /0.1325_r8, 0.0009_r8/ ! CO2

  real(r8) :: c_h_inv(ngases)    ! constant (K) for Henry's law (4.12, Wania)
  data c_h_inv(1:3) /1600._r8, 1500._r8, 2400._r8/ ! CH4, O2, CO2
  real(r8) :: kh_theta(ngases)    ! Henry's constant (L.atm/mol) at standard temperature (298K)
  data kh_theta(1:3) /714.29_r8, 769.23_r8, 29.4_r8/ ! CH4, O2, CO2
  real(r8) :: kh_tbase = 298._r8 ! base temperature for calculation of Henry's constant (K)

! !PUBLIC MEMBER FUNCTIONS:
  public clm_varcon_init          ! Initialze constants that need to be initialized

! !REVISION HISTORY:
! Created by Mariana Vertenstein

!EOP
!-----------------------------------------------------------------------
contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_varcon_init
!
! !INTERFACE:
  subroutine clm_varcon_init()
!
! !DESCRIPTION:
! This subroutine initializes constants in clm_varcon. MUST be called 
! after the clm_varpar_init.
!
! !USES:
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
!   Created by E. Kluzek
!
!
! !LOCAL VARIABLES:
!
!EOP
!------------------------------------------------------------------------------
  allocate( zsoi(1:nlevgrnd) )
  allocate( dzsoi(1:nlevgrnd) )
  allocate( zisoi(0:nlevgrnd) )
  allocate( dzsoi_decomp(1:nlevdecomp_full) )
  
  !! Use this setting temporarily. Need to improve if VERTSOILC is turned on!! fzeng, 13 Mar 2017  
  zsoi(1)  = 0.5
  dzsoi(1) = 1.
  zisoi(0) = 0.
  zisoi(1) = 1.
  dzsoi_decomp(1) = dzsoi(1)

  end subroutine clm_varcon_init

end module clm_varcon
