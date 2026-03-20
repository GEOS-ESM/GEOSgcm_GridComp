MODULE PlumeStateModule
  IMPLICIT NONE
  TYPE PlumeState
    INTEGER :: ierr
    INTEGER :: jmin, klcl, k224, kbcon, ktop, kstabi, kstabm
    REAL :: cprr, xmb, edt, pwav, sigma
    REAL, DIMENSION(:), ALLOCATABLE :: entr, pcup, up_massentr, up_massdetr, dd_massentr, dd_massdetr, &
         zup, zdn, prup, prdn, clwup, tup, conv_cld_fr, sgs_vvel
  END TYPE PlumeState

CONTAINS

  SUBROUTINE initialize_plume_states(mxp, myp, maxiens, mzp, plume_states)
    INTEGER, INTENT(IN) :: mxp, myp, maxiens, mzp
    TYPE(PlumeState), DIMENSION(mxp,myp,maxiens), INTENT(INOUT) :: plume_states
    INTEGER :: i, j, plume
    DO i = 1, mxp
      DO j = 1, myp
        DO plume = 1, maxiens
          ALLOCATE(plume_states(i,j,plume)%entr(mzp))
          ALLOCATE(plume_states(i,j,plume)%pcup(mzp))
          ALLOCATE(plume_states(i,j,plume)%up_massentr(mzp))
          ALLOCATE(plume_states(i,j,plume)%up_massdetr(mzp))
          ALLOCATE(plume_states(i,j,plume)%dd_massentr(mzp))
          ALLOCATE(plume_states(i,j,plume)%dd_massdetr(mzp))
          ALLOCATE(plume_states(i,j,plume)%zup(mzp))
          ALLOCATE(plume_states(i,j,plume)%zdn(mzp))
          ALLOCATE(plume_states(i,j,plume)%prup(mzp))
          ALLOCATE(plume_states(i,j,plume)%prdn(mzp))
          ALLOCATE(plume_states(i,j,plume)%clwup(mzp))
          ALLOCATE(plume_states(i,j,plume)%tup(mzp))
          ALLOCATE(plume_states(i,j,plume)%conv_cld_fr(mzp))
          ALLOCATE(plume_states(i,j,plume)%sgs_vvel(mzp))
        END DO
      END DO
    END DO
  END SUBROUTINE initialize_plume_states

END MODULE PlumeStateModule
