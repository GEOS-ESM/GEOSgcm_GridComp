module EDMFPARAMS

 implicit none

 type EDMFPARAMS_TYPE
    integer :: DISCRETE
    integer :: IMPLICIT
    integer :: ENTRAIN
    integer :: DOCLASP
    integer :: THERMAL_PLUME
    integer :: TEST
    integer :: DEBUG
    integer :: ET
    real    :: L0
    real    :: L0fac
    real    :: STOCHFRAC
    real    :: ENTWFAC
    real    :: EDFAC
    real    :: ENT0
    real    :: ALPHATH
    real    :: ALPHAQT
    real    :: ALPHAW
    real    :: PWMAX
    real    :: PWMIN
    real    :: WA
    real    :: WB
    real    :: AU0
    real    :: CTH1
    real    :: CTH2
    real    :: RH0_QB
    real    :: C_KH_MF
    real    :: ICE_RAMP
 endtype EDMFPARAMS_TYPE

end module EDMFPARAMS

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