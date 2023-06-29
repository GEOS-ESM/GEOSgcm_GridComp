
select case (TYPE)

case (Starr)
   TT = TI - ZEROC
   if    (TT < TSTARR1                   ) then
      EX = (TT*(TT*(TT*(TT*(TT*(TT*S16+S15)+S14)+S13)+S12)+S11)+S10)
   elseif(TT >= TSTARR1 .and. TT < TSTARR2) then
      W = (TSTARR2 - TT)/(TSTARR2-TSTARR1)
      EX =       W *(TT*(TT*(TT*(TT*(TT*(TT*S16+S15)+S14)+S13)+S12)+S11)+S10) &
           + (1.-W)*(TT*(TT*(TT*(TT*(TT*(TT*S26+S25)+S24)+S23)+S22)+S21)+S20)
   elseif(TT >= TSTARR2 .and. TT < TSTARR3) then
      EX = (TT*(TT*(TT*(TT*(TT*(TT*S26+S25)+S24)+S23)+S22)+S21)+S20)
   elseif(TT >= TSTARR3 .and. TT < TSTARR4) then
      W = (TSTARR4 - TT)/(TSTARR4-TSTARR3)
      EX =       W *(TT*(TT*(TT*(TT*(TT*(TT*S26+S25)+S24)+S23)+S22)+S21)+S20) &
           + (1.-W)*(TT*(TT*(TT*(TT*(TT*(TT*BI6+BI5)+BI4)+BI3)+BI2)+BI1)+BI0)
   else
      EX = (TT*(TT*(TT*(TT*(TT*(TT*BI6+BI5)+BI4)+BI3)+BI2)+BI1)+BI0)
   endif

case (GoffGratch)
   TT = ZEROC/TI
   EX = DI(0) * exp( -(DI(1)/TT + DI(2)*log(TT) + DI(3)*TT) )

case (MurphyKoop)
   EX = exp( CI(0)+ CI(1)/TI + CI(2)*log(TI) + CI(3)*TI )

end select

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