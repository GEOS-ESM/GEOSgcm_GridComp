module module_gate
  IMPLICIT NONE

  !-for GEOS-5 runs, set use_gate=.false.
  LOGICAL, PARAMETER :: use_gate = .false. 
  
  !- Here are the place for data related with the GATE soundings
  INTEGER, PARAMETER :: gate = 1   ! flag to turn on/off : 1/0
  INTEGER, PARAMETER :: klon = 161 ! number of soundings for gate
  INTEGER, PARAMETER :: klev = 41  ! number of vertical levels
  INTEGER, PARAMETER :: ktrac= 2   ! number of chemical tracers
  INTEGER, PARAMETER :: levs=klev  
  INTEGER, PARAMETER :: nvar_grads=200
  
  TYPE cupout_vars
    REAL             ,POINTER     :: varp(:,:)
    CHARACTER(LEN=80),ALLOCATABLE :: varn(:)
  END  TYPE cupout_vars
  
  TYPE (cupout_vars), ALLOCATABLE :: cupout(:)
  
  REAL, DIMENSION(klon,klev):: pgeo,ppres,ptemp,pq,pu,pv,pvervel, &
                               zrvten,ztten,zq1,zq2,zqr,zadvt,zadvq

  INTEGER :: JL,KLEV_SOUND
  character (len=128) :: runname, runlabel,rundata="NONE" 

end module module_gate

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