      module CLDPARAMS

      use ESMF

      implicit none

      character(len=ESMF_MAXSTR) :: CLDMICRO
      logical                    :: LCLDMICRO

! Define CLDPARAM
! ---------------
      type CLDPARAM_TYPE
           real               :: CNV_BETA              ! 1
           real               :: ANV_BETA              ! 2
           real               :: LS_BETA               ! 3
           real               :: RH_CRIT               ! 4
           real               :: AUTOC_LS              ! 5
           real               :: QC_CRIT_LS            ! 6
           real               :: ACCRETION             ! 7
           real               :: RAIN_REVAP_FAC        ! 8
           real               :: VOL_TO_FRAC           ! 9
           real               :: SUPERSAT              ! 10
           real               :: SHEAR_EVAP_FAC        ! 11
           real               :: MIN_ALLOW_CCW         ! 12
           real               :: CCW_EVAP_EFF          ! 13
           real               :: NSUB_AUTOCONV         ! 14
           real               :: LS_SUND_INTER         ! 15
           real               :: LS_SUND_COLD          ! 16
           real               :: LS_SUND_TEMP1         ! 17
           real               :: ANV_SUND_INTER        ! 18
           real               :: ANV_SUND_COLD         ! 19
           real               :: ANV_SUND_TEMP1        ! 20
           real               :: ANV_TO_LS_TIME        ! 21
           real               :: CCN_OCEAN             ! 22
           real               :: CCN_LAND              ! 23
           real               :: NCCN_ANVIL_NULL       ! 24
           real               :: NCCN_PBL_NULL         ! 25
           real               :: DISABLE_RAD           ! 26
           real               :: ICE_SETTLE            ! 27
           real               :: ANV_ICEFALL           ! 28
           real               :: LS_ICEFALL            ! 29
           real               :: REVAP_OFF_P           ! 30
           real               :: CNV_ENVF              ! 31
           real               :: ANV_ENVF              ! 31
           real               :: SC_ENVF               ! 31
           real               :: LS_ENVF               ! 31
           real               :: WRHODEP               ! 32
           real               :: ICE_RAMP              ! 33
           real               :: CNV_ICEPARAM          ! 34
           real               :: CNV_ICEFRPWR          ! 35
           real               :: CNV_DDRF              ! 36
           real               :: ANV_DDRF              ! 37
           real               :: LS_DDRF               ! 38
           real               :: AUTOC_ANV             ! 39
           real               :: QC_CRIT_ANV           ! 40
           real               :: TANHRHCRIT            ! 41
           real               :: MINRHCRIT             ! 42
           real               :: MAXRHCRIT             ! 43
           real               :: PRECIPRAD             ! 44
           real               :: TURNRHCRIT            ! 45
           real               :: MAXRHCRITLAND         ! 46
           real               :: FR_LS_WAT             ! 47
           real               :: FR_LS_ICE             ! 48
           real               :: FR_AN_WAT             ! 49
           real               :: FR_AN_ICE             ! 50
           real               :: MIN_RL                ! 51
           real               :: MIN_RI                ! 52
           real               :: MAX_RL                ! 53
           real               :: MAX_RI                ! 54
           real               :: FAC_RL                ! 55
           real               :: FAC_RI                ! 56
           real               :: SNOW_REVAP_FAC        ! 57
           real               :: PDFSHAPE              ! 58
           real               :: TURNRHCRIT_UP         ! 59
           real               :: SLOPERHCRIT           ! 60
           real               :: MIN_LTS               ! 61
           integer            :: CFPBL_EXP             ! 62
           real               :: DISP_FACTOR_LIQ       ! 63
           real               :: DISP_FACTOR_ICE       ! 63
           real               :: SCLM_SHALLOW          ! 63
           real               :: SCLM_DEEP             ! 63
      endtype CLDPARAM_TYPE

  end module CLDPARAMS
