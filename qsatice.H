#define TABLE   ESTBLE
#define FILE_    "esatice.H"
#define PHASE   ICE

if(TX>=TMAX(ICE)) then

#if KIND_==4
#define FUNC Qsat0
#else
#define FUNC Qsatd0
#endif

   if(present(PL)) then
      if(present(DQ)) then
         EX = FUNC(TX,PL=PX,DQ=DX)
      else
         EX = FUNC(TX,PL=PX)
      endif
   else
      if(present(DQ)) then
         EX = FUNC(TX,DQ=DX)
      else
         EX = FUNC(TX)
      endif
   endif

#undef FUNC

else

#include "eqsat.H"

endif

#undef TABLE
#undef FILE_
#undef PHASE

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