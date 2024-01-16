module radcoup_loop

    use Process_Library_standalone
    use GEOS_UtilsMod

    implicit none

    public radcoup_loop_standalone


    contains

    subroutine radcoup_loop_standalone(IM, JM, LM, FAC_RL, MIN_RL, MAX_RL, FAC_RI, MIN_RI, MAX_RI, do_qa, &
                                       PLmb, QRAIN, QSNOW, QGRAUPEL, NACTL, NACTI, &
                                       Q, T, QLLS, QILS, CLLS, QLCN, QICN, CLCN, &
                                       RAD_QV, RAD_QL, RAD_QI, RAD_QR, RAD_QS, &
                                       RAD_QG, RAD_CF, CLDREFFL, CLDREFFI, RHX)
        integer, intent(IN) :: IM, JM, LM

        real, intent(IN) :: FAC_RL, MIN_RL, MAX_RL, FAC_RI, MIN_RI, MAX_RI

        logical, intent(IN) :: do_qa

        real, dimension(:,:,:), intent(IN)    :: PLmb, QRAIN, QSNOW, QGRAUPEL, NACTL, NACTI
        real, dimension(:,:,:), intent(INOUT) :: Q, T, QLLS, QILS, CLLS, QLCN, QICN, CLCN
        real, dimension(:,:,:), intent(OUT)   :: RAD_QV, RAD_QL, RAD_QI, RAD_QR, RAD_QS, &
                                                 RAD_QG, RAD_CF, CLDREFFL, CLDREFFI, RHX

        integer :: I, J, L

        ! This calls ESINIT so that subsequent calls to GEOS_QSAT do not result in
        ! a race condition.
        call call_ESINIT

!!$acc parallel loop gang vector collapse(3)
!$omp target teams distribute parallel do collapse(3)
        do L = 1, LM
            do J = 1, JM
                do I = 1, IM
                    ! cleanup clouds
                    call FIX_UP_CLOUDS( Q(I,J,L), T(I,J,L), QLLS(I,J,L), QILS(I,J,L), &
                        CLLS(I,J,L), QLCN(I,J,L), QICN(I,J,L), CLCN(I,J,L) )
                    ! get radiative properties
                    call RADCOUPLE ( T(I,J,L), PLmb(I,J,L), CLLS(I,J,L), CLCN(I,J,L), &
                            Q(I,J,L), QLLS(I,J,L), QILS(I,J,L), QLCN(I,J,L), QICN(I,J,L), &
                            QRAIN(I,J,L), QSNOW(I,J,L), QGRAUPEL(I,J,L), NACTL(I,J,L), NACTI(I,J,L), &
                            RAD_QV(I,J,L), RAD_QL(I,J,L), RAD_QI(I,J,L), RAD_QR(I,J,L), RAD_QS(I,J,L), &
                            RAD_QG(I,J,L), RAD_CF(I,J,L), CLDREFFL(I,J,L), CLDREFFI(I,J,L), &
                            FAC_RL, MIN_RL, MAX_RL, FAC_RI, MIN_RI, MAX_RI)
                    ! Note : GEOS_QSAT has not been tested in standalone
                    if (do_qa) RHX(I,J,L) = Q(I,J,L)/GEOS_QSAT( T(I,J,L), PLmb(I,J,L) )
                enddo
            enddo
        enddo
!$omp end target teams distribute parallel do
!!$acc end parallel loop
    end subroutine

end module

! NASA Docket No. GSC-15,354-1, and identified as "GEOS-5 GCM Modeling Software”
  
! “Copyright © 2008 United States Government as represented by the Administrator
! of the National Aeronautics and Space Administration. All Rights Reserved.”
  
! Licensed under the Apache License, Version 2.0 (the "License"); you may not use
! this file except in compliance with the License. You may obtain a copy of the
! License at
  
! http://www.apache.org/licenses/LICENSE-2.0
  
! Unless required by applicable law or agreed to in writing, software distributed
! under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
! CONDITIONS OF ANY KIND, either express or implied. See the License for the
! specific language governing permissions and limitations under the License.