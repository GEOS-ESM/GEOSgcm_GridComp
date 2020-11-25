!     path:      $Source$
!     author:    $Author$
!     revision:  $Revision$
!     created:   $Date$

!----------------------------------------------------------------------------
! Copyright (c) 2002-2016, Atmospheric & Environmental Research, Inc. (AER)
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!  * Redistributions of source code must retain the above copyright
!    notice, this list of conditions and the following disclaimer.
!  * Redistributions in binary form must reproduce the above copyright
!    notice, this list of conditions and the following disclaimer in the
!    documentation and/or other materials provided with the distribution.
!  * Neither the name of Atmospheric & Environmental Research, Inc., nor
!    the names of its contributors may be used to endorse or promote products
!    derived from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
! ARE DISCLAIMED. IN NO EVENT SHALL ATMOSPHERIC & ENVIRONMENTAL RESEARCH, INC., 
! BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF 
! THE POSSIBILITY OF SUCH DAMAGE.
!                        (http://www.rtweb.aer.com/)                        
!----------------------------------------------------------------------------

! **************************************************************************
!      subroutine sw_kgbnn
! **************************************************************************
!  RRTM Shortwave Radiative Transfer Model
!  Atmospheric and Environmental Research, Inc., Cambridge, MA
!
!  Original by J.Delamere, Atmospheric & Environmental Research.
!  Reformatted for F90: JJMorcrette, ECMWF
!  Further F90 and GCM revisions:  MJIacono, AER, July 2002
!  Solar variability revisions:  MJIacono, AER, November 2015
!
!  This file contains 14 subroutines that include the 
!  absorption coefficients and other data for each of the 14 shortwave
!  spectral bands used in RRTM_SW.  Here, the data are defined for 16
!  g-points, or sub-intervals, per band.  These data are combined and
!  weighted using a mapping procedure in routine RRTMG_SW_INIT to reduce
!  the total number of g-points from 224 to 112 for use in the GCM.
! **************************************************************************
      subroutine sw_kgb26
! **************************************************************************

      !use parkind, only : im => kind , rb => kind  
      use rrsw_kg26, only : sfluxrefo, raylo, &
                            irradnceo, facbrghto, snsptdrko

      implicit none
      save

! Kurucz solar source function
      sfluxrefo(:) = (/ &
!         &     129.462, 15*0. /)
        &   29.0079,  28.4088,     20.3099,  13.0283 &
        &,  11.8619,  9.95840,     6.68696,  5.38987 &
        &,  3.49829, 0.407693,    0.299027, 0.236827 &
        &, 0.188502, 0.163489, 4.64335e-02, 2.72662e-03 /)

! Solar variability components: time-invariant baseline quiet sun irradiance
      irradnceo(:) = (/ &
        & 3.04845E+01, 2.82835E+01, 2.06520E+01, 1.33606E+01,&
        & 1.19253E+01, 9.46266E+00, 6.49295E+00, 4.95470E+00,&
        & 3.04469E+00, 3.68830E-01, 2.69394E-01, 2.19130E-01,&
        & 1.67527E-01, 1.51647E-01, 3.99119E-02, 0.00000E+00/)
! Solar variability components: facular brightening
      facbrghto(:) = (/ &
        & 3.61383E-02, 3.05195E-02, 2.92002E-02, 3.95202E-02,&
        & 2.30227E-02, 1.34627E-02, 1.21898E-02, 8.33975E-03,&
        & 5.78516E-03, 6.39681E-04, 5.38531E-04, 4.27531E-04,&
        & 3.12727E-04, 1.98617E-04, 7.47491E-05, 1.03928E-05/)
! Solar variability components: sunspot darkening
      snsptdrko(:) = (/ &
        &-1.84373E-02,-1.86059E-02,-1.38321E-02,-8.72988E-03,&
        &-8.16733E-03,-6.16612E-03,-4.27388E-03,-3.57388E-03,&
        &-2.10064E-03,-2.32508E-04,-1.92973E-04,-1.51402E-04,&
        &-1.09697E-04,-6.89804E-05,-2.57998E-05,-3.58079E-06/)

! Rayleigh extinction coefficient at all v 
      raylo(:) = (/ &
        &  1.21263e-06,1.43428e-06,1.67677e-06,1.93255e-06 &
        &, 2.19177e-06,2.44195e-06,2.66926e-06,2.85990e-06 &
        &, 3.00380e-06,3.06996e-06,3.08184e-06,3.09172e-06 &
        &, 3.09938e-06,3.10456e-06,3.10727e-06,3.10818e-06 /)

      end subroutine sw_kgb26
