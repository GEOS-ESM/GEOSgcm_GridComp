MODULE ConvParGFHarness
  USE PlumeStateModule
  USE PlumeDiagnosticsModule
  USE ConvParGFCore
  IMPLICIT NONE

CONTAINS

! Public interface routines (modernized)
SUBROUTINE GF2020_INTERFACE(...)
  ! Call convection driver
  CALL GF2020_DRV(...)
  ! Extract diagnostics and profiles
  CALL process_all_plume_profiles(...)
  CALL extract_all_plume_diagnostics(...)
END SUBROUTINE GF2020_INTERFACE

SUBROUTINE GF2020_BEFORE(...)
  ! Setup and initial allocation; handle initialization using PlumeStateModule
END SUBROUTINE GF2020_BEFORE

SUBROUTINE GF2020_AFTER(...)
  ! Post-processing; final extraction of profiles/diagnostics
END SUBROUTINE GF2020_AFTER

END MODULE ConvParGFHarness
