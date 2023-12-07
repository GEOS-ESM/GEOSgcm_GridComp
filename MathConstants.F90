module MAPL_MathConstantsMod

    use, intrinsic :: iso_fortran_env, only: REAL64, REAL32
 
    implicit none
 
 !=============================================================================
 !BOP
 
 ! !MODULE: -- A container module for MAPL mathematical constants
 
 ! !PUBLIC VARIABLES:
    real(kind=REAL64), parameter :: MAPL_PI_R8              = 3.14159265358979323846d0
    real(kind=REAL32), parameter :: MAPL_PI                 = MAPL_PI_R8
    real(kind=REAL64), parameter :: MAPL_DEGREES_TO_RADIANS_R8 = MAPL_PI_R8 / 180.d0
    real(kind=REAL32), parameter :: MAPL_DEGREES_TO_RADIANS = MAPL_PI / 180.
    real(kind=REAL64), parameter :: MAPL_RADIANS_TO_DEGREES = 180.d0 / MAPL_PI_R8
 
 !EOP
 
 end module MAPL_MathConstantsMod

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