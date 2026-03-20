MODULE PlumeDiagnosticsModule
  USE PlumeStateModule
  IMPLICIT NONE

CONTAINS

! Routine to process all plume profiles (originally: process_all_plume_profiles)
SUBROUTINE process_all_plume_profiles(mxp, myp, mzp, do_this_column, flip, REVSU, PRFIL, REVSU_GF, PRFIL_GF)
  INTEGER, INTENT(IN) :: mxp, myp, mzp
  LOGICAL, DIMENSION(mxp,myp), INTENT(IN) :: do_this_column
  INTEGER, DIMENSION(mzp), INTENT(IN) :: flip
  REAL, DIMENSION(mxp,myp,mzp), INTENT(OUT) :: REVSU, PRFIL
  REAL, DIMENSION(mzp,mxp,myp), INTENT(IN) :: REVSU_GF, PRFIL_GF
  INTEGER :: i, j, k
  DO j = 1, myp
    DO i = 1, mxp
      IF(.NOT. do_this_column(i,j)) CYCLE
      DO k = 1, mzp
        REVSU(i,j,k) = REVSU_GF(flip(k),i,j)
        PRFIL(i,j,k) = PRFIL_GF(flip(k),i,j)
      END DO
    END DO
  END DO
END SUBROUTINE process_all_plume_profiles

! Routine to extract all plume diagnostics (originally: extract_all_plume_diagnostics)
SUBROUTINE extract_all_plume_diagnostics(mxp, myp, maxiens, mzp, do_this_column, icumulus_gf, plume_states, press, 
    CNV_TOPP_DP, CNV_TOPP_MD, CNV_TOPP_SH, MFDP, MFSH, MFMD, SIGMA_DEEP, SIGMA_MID, MUPDP, MDNDP, MUPSH, MUPMD)
  INTEGER, INTENT(IN) :: mxp, myp, maxiens, mzp
  LOGICAL, DIMENSION(mxp,myp), INTENT(IN) :: do_this_column
  LOGICAL, DIMENSION(maxiens), INTENT(IN) :: icumulus_gf
  TYPE(PlumeState), DIMENSION(mxp,myp,maxiens), INTENT(IN) :: plume_states
  REAL, DIMENSION(mzp, mxp, myp), INTENT(IN) :: press
  REAL, DIMENSION(mxp,myp), INTENT(OUT) :: CNV_TOPP_DP, CNV_TOPP_MD, CNV_TOPP_SH, MFDP, MFSH, MFMD, SIGMA_DEEP, SIGMA_MID
  REAL, DIMENSION(mxp,myp,mzp), INTENT(OUT) :: MUPDP, MDNDP, MUPSH, MUPMD
  INTEGER :: i, j, plume

  DO j = 1, myp
    DO i = 1, mxp
      IF(.NOT. do_this_column(i,j)) CYCLE
      DO plume = 1, maxiens
        IF(icumulus_gf(plume)) THEN
          CALL extract_plume_diagnostics(plume, i, j, mzp, plume_states, press, CNV_TOPP_DP, CNV_TOPP_MD, CNV_TOPP_SH, MFDP, MFSH, MFMD, SIGMA_DEEP, SIGMA_MID, MUPDP, MDNDP, MUPSH, MUPMD)
        END IF
      END DO
    END DO
  END DO
END SUBROUTINE extract_all_plume_diagnostics

! Routine to extract one plume type's diagnostics (originally: extract_plume_diagnostics)
SUBROUTINE extract_plume_diagnostics(plume_type, i, j, mzp, plume_states, press, 
    CNV_TOPP_DP, CNV_TOPP_MD, CNV_TOPP_SH, MFDP, MFSH, MFMD, SIGMA_DEEP, SIGMA_MID, MUPDP, MDNDP, MUPSH, MUPMD)
  INTEGER, INTENT(IN) :: plume_type, i, j, mzp
  TYPE(PlumeState), DIMENSION(:,:,:), INTENT(IN) :: plume_states
  REAL, DIMENSION(:,:,:), INTENT(IN) :: press
  REAL, DIMENSION(:,:), INTENT(OUT) :: CNV_TOPP_DP, CNV_TOPP_MD, CNV_TOPP_SH, MFDP, MFSH, MFMD, SIGMA_DEEP, SIGMA_MID
  REAL, DIMENSION(:,:,:), INTENT(OUT) :: MUPDP, MDNDP, MUPSH, MUPMD

  SELECT CASE(plume_type)
  CASE(DEEP)
    CNV_TOPP_DP(i,j)  = press(plume_states(i,j,DEEP)%ktop,i,j)
    MFDP(i,j)         = plume_states(i,j,DEEP)%xmb
    SIGMA_DEEP(i,j)   = plume_states(i,j,DEEP)%sigma
    MUPDP(i,j,1:mzp)  = plume_states(i,j,DEEP)%zup(1:mzp) * plume_states(i,j,DEEP)%xmb
    MDNDP(i,j,1:mzp)  = plume_states(i,j,DEEP)%zdn(1:mzp) * plume_states(i,j,DEEP)%edt * plume_states(i,j,DEEP)%xmb
  CASE(SHAL)
    CNV_TOPP_SH(i,j) = press(plume_states(i,j,SHAL)%ktop,i,j)
    MFSH(i,j)        = plume_states(i,j,SHAL)%xmb
    MUPSH(i,j,1:mzp) = plume_states(i,j,SHAL)%zup(1:mzp) * plume_states(i,j,SHAL)%xmb
  CASE(MID)
    CNV_TOPP_MD(i,j) = press(plume_states(i,j,MID)%ktop,i,j)
    MFMD(i,j)        = plume_states(i,j,MID)%xmb
    SIGMA_MID(i,j)   = plume_states(i,j,MID)%sigma
    MUPMD(i,j,1:mzp) = plume_states(i,j,MID)%zup(1:mzp) * plume_states(i,j,MID)%xmb
  END SELECT
END SUBROUTINE extract_plume_diagnostics

END MODULE PlumeDiagnosticsModule
