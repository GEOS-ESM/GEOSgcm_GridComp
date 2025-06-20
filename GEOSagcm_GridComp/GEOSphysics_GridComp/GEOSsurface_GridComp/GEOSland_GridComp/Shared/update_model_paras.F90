MODULE update_model_paras
  
  implicit none

  private

  ! ============================================================================================================
  ! 
  ! The parameter values that need to be tested are defined inside this module. The module uses the fortran
  !    SAVE attribute for those parameters thus lsm_calib_driver.F90 can read desired values from
  !    the LSM_Calib.rc resource file and update them with user specified values.
  ! The module implementation requires a "USE  update_model_paras" line (see below example) in physics modules
  !    where the parameter values are being used. The user can add any parameter here with associated changes 
  !    in corresponding physics subroutines.
  !    USE  update_model_paras, ONLY:              &
  !     RZD_RESET1, RZD_RESET2, PT_CHNG1, SRUN_CHNG1, BEE_ADD, RZD_SFC  
  ! No module uses the update_model_paras module in default LDASsa configuration. 
  !
  ! ============================================================================================================

  REAL,    PUBLIC :: RSWILT     =   1500.  
  REAL,    PUBLIC :: RSSAT      =   25.    
  LOGICAL, PUBLIC :: RZD_RESET1 = .false. ! if .true. attempts to reduce surface soil moisture increases associated with precipitation inputs.
  LOGICAL, PUBLIC :: RZD_RESET2 = .false. ! if .true. attempts to reduce the recovery of surface deficits associated with bare soil evaporation.
  LOGICAL, PUBLIC :: PT_CHNG1   = .false. ! if .true. This changes the near surface equilibrium profile of soil 
                                       !           moisture, to reflect the different equilibrium characteristics of 
                                       !           natural soils vs those in the lab.
  LOGICAL, PUBLIC :: SRUN_CHNG1 = .false. ! if .true. increases runoff fraction
  REAL,    PUBLIC :: BEE_ADD    = 0.      ! modify BEE at the surface
  REAL,    PUBLIC :: RZD_SFC    = 1.      !
  REAL,    PUBLIC :: ASTRFR = 0.3333      ! STRESS TRANSITION POINT
  REAL,    PUBLIC :: STEXP  = 1.          ! STRESS RAMPING

  SAVE RSWILT, RZD_RESET1, RZD_RESET2, PT_CHNG1, SRUN_CHNG1, BEE_ADD, RZD_SFC,  RSSAT, ASTRFR, STEXP  

  public :: upd_params

  contains

    subroutine upd_params (PNAME,value)

      implicit none

      character(*), intent (in) :: PNAME
      real,         intent (in) :: value

      select case (PNAME)

      case ('RSWILT')
         
         RSWILT = value

      case ('RZD_RESET1')
         
         RZD_RESET1 = (value == 1.)

      case ('RZD_RESET2')
         
         RZD_RESET2 = (value == 1.)

      case ('PT_CHNG1')
         
         PT_CHNG1   = (value == 1.)

      case ('SRUN_CHNG1')
         
         SRUN_CHNG1 = (value == 1.)

      case ('BEE_ADD')
         
         BEE_ADD = value

      case ('RZD_SFC')
         
         RZD_SFC = value / 100.

      case ('ENG2_ASTRFR')
         
         ASTRFR = value / 100.

      case ('ENG2_STEXP')
         
         STEXP  = value 

      case ('RSSAT')
         
         RSSAT = value

      end select

    end subroutine upd_params

  end MODULE update_model_paras
