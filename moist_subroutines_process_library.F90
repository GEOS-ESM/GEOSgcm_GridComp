module GEOSmoist_Process_Library

    use MAPL_ConstantsMod
    ! use GEOS_UtilsMod

    implicit none

    private

    interface ICE_FRACTION
        module procedure ICE_FRACTION_3D
        module procedure ICE_FRACTION_2D
        module procedure ICE_FRACTION_1D
        module procedure ICE_FRACTION_SC
    end interface ICE_FRACTION

    ! In anvil/convective clouds
    real, parameter :: aT_ICE_ALL = 252.16
    real, parameter :: aT_ICE_MAX = 268.16
    real, parameter :: aICEFRPWR  = 2.0
    ! Over snow/ice SRF_TYPE = 2
    real, parameter :: iT_ICE_ALL = 236.16
    real, parameter :: iT_ICE_MAX = 261.16
    real, parameter :: iICEFRPWR  = 6.0
    ! Over Land     SRF_TYPE = 1
    real, parameter :: lT_ICE_ALL = 239.16
    real, parameter :: lT_ICE_MAX = 261.16
    real, parameter :: lICEFRPWR  = 2.0
    ! Over Oceans   SRF_TYPE = 0
    real, parameter :: oT_ICE_ALL = 238.16
    real, parameter :: oT_ICE_MAX = 263.16
    real, parameter :: oICEFRPWR  = 4.0

    ! parameters
    real, parameter :: EPSILON =  MAPL_H2OMW/MAPL_AIRMW
    real, parameter :: K_COND  =  2.4e-2    ! J m**-1 s**-1 K**-1
    real, parameter :: DIFFU   =  2.2e-5    ! m**2 s**-1
    real, parameter :: taufrz  =  450.0
    real, parameter :: dQCmax  =  1.e-4
    ! LDRADIUS4
    ! Liquid  based on DOI 10.1088/1748-9326/3/4/045021
    real, parameter :: RHO_W   = 1000.0  ! Density of liquid water in kg/m^3
    real, parameter :: Ldiss   = 0.07    ! tunable dispersion effect
    real, parameter :: Lk      = 0.75    ! tunable shape effect (0.5:1)
    real, parameter :: Lbe     = 1./3. - 0.14
    real, parameter :: Lbx     = Ldiss*1.e3*(3./(4.*MAPL_PI*Lk*RHO_W*1.e-3))**(1./3.)
                                ! LDRADIUS eqs are in cgs units
    ! Ice
    real, parameter :: RHO_I   =  916.8  ! Density of ice crystal in kg/m^3

    ! combined constantc
    real, parameter :: cpbgrav = MAPL_CP/MAPL_GRAV
    real, parameter :: gravbcp = MAPL_GRAV/MAPL_CP
    real, parameter :: alhlbcp = MAPL_ALHL/MAPL_CP
    real, parameter :: alhfbcp = MAPL_ALHF/MAPL_CP
    real, parameter :: alhsbcp = MAPL_ALHS/MAPL_CP

    real, parameter :: mapl_undef = 1.0e15  ! NOTE : This is the value pulled from MAPL_Mod

    ! Tracer Bundle things for convection
    type CNV_Tracer_Type
    real, pointer              :: Q(:,:,:) => null()
    real                       :: fscav = 0.0
    real                       :: Vect_Hcts(4)
    character(len=100) :: QNAME ! Tracer Name
    character(len=100) :: CNAME ! Component Name
    end type CNV_Tracer_Type
    type(CNV_Tracer_Type), allocatable :: CNV_Tracers(:)

    public :: ICE_FRACTION, CNV_Tracers
!$acc declare create(aT_ICE_ALL, aT_ICE_MAX, aICEFRPWR, iT_ICE_ALL, iT_ICE_MAX, iICEFRPWR, &
!$acc                lT_ICE_ALL, lT_ICE_MAX, lICEFRPWR, oT_ICE_ALL, oT_ICE_MAX, oICEFRPWR, &
!$acc                EPSILON, K_COND, DIFFU, taufrz, dQCmax, RHO_W, Ldiss, Lk, Lbe, Lbx, RHO_I, &
!$acc                cpbgrav, gravbcp, alhlbcp, alhfbcp, alhsbcp, mapl_undef)
    contains    

    function ICE_FRACTION_3D (TEMP,CNV_FRACTION,SRF_TYPE) RESULT(ICEFRCT)
        real, intent(in) :: TEMP(:,:,:),CNV_FRACTION(:,:),SRF_TYPE(:,:)
        real :: ICEFRCT(size(TEMP,1),size(TEMP,2),size(TEMP,3))
        integer :: i,j,l
        do l=1,size(TEMP,3)
            do j=1,size(TEMP,2)
                do i=1,size(TEMP,1)
                    ICEFRCT(i,j,l) = ICE_FRACTION_SC(TEMP(i,j,l),CNV_FRACTION(i,j),SRF_TYPE(i,j))
                enddo
            enddo
        enddo
    end function ICE_FRACTION_3D

    function ICE_FRACTION_2D (TEMP,CNV_FRACTION,SRF_TYPE) RESULT(ICEFRCT)
        real, intent(in) :: TEMP(:,:),CNV_FRACTION(:,:),SRF_TYPE(:,:)
        real :: ICEFRCT(size(TEMP,1),size(TEMP,2))
        integer :: i,j
        do j=1,size(TEMP,2)
            do i=1,size(TEMP,1)
                ICEFRCT(i,j) = ICE_FRACTION_SC(TEMP(i,j),CNV_FRACTION(i,j),SRF_TYPE(i,j))
            enddo
        enddo
    end function ICE_FRACTION_2D

    function ICE_FRACTION_1D (TEMP,CNV_FRACTION,SRF_TYPE) RESULT(ICEFRCT)
        real, intent(in) :: TEMP(:),CNV_FRACTION(:),SRF_TYPE(:)
        real :: ICEFRCT(size(TEMP))
        integer :: i
        do i=1,size(TEMP)
            ICEFRCT(i) = ICE_FRACTION_SC(TEMP(i),CNV_FRACTION(i),SRF_TYPE(i))
        enddo
    end function ICE_FRACTION_1D

    function ICE_FRACTION_SC (TEMP,CNV_FRACTION,SRF_TYPE) RESULT(ICEFRCT)
    !$acc routine seq
        real, intent(in) :: TEMP,CNV_FRACTION,SRF_TYPE
        real             :: ICEFRCT
        real             :: tc, ptc
        real             :: ICEFRCT_C, ICEFRCT_M
  
        ! Anvil clouds
        ! Anvil-Convective sigmoidal function like figure 6(right)
        ! Sigmoidal functions Hu et al 2010, doi:10.1029/2009JD012384
        ICEFRCT_C  = 0.00
        if ( TEMP <= aT_ICE_ALL ) then
            ICEFRCT_C = 1.000
        else if ( (TEMP > aT_ICE_ALL) .AND. (TEMP <= aT_ICE_MAX) ) then
            ICEFRCT_C = SIN( 0.5*MAPL_PI*( 1.00 - ( TEMP - aT_ICE_ALL ) / ( aT_ICE_MAX - aT_ICE_ALL ) ) )
        end if
        ICEFRCT_C = MIN(ICEFRCT_C,1.00)
        ICEFRCT_C = MAX(ICEFRCT_C,0.00)
        ICEFRCT_C = ICEFRCT_C**aICEFRPWR
#ifdef MODIS_ICE_POLY
        ! Use MODIS polynomial from Hu et al, DOI: (10.1029/2009JD012384) 
        tc = MAX(-46.0,MIN(TEMP-MAPL_TICE,46.0)) ! convert to celcius and limit range from -46:46 C
        ptc = 7.6725 + 1.0118*tc + 0.1422*tc**2 + 0.0106*tc**3 + 0.000339*tc**4 + 0.00000395*tc**5
        ICEFRCT_M = 1.0 - (1.0/(1.0 + exp(-1*ptc)))
#else
        ! Sigmoidal functions like figure 6b/6c of Hu et al 2010, doi:10.1029/2009JD012384
        if (SRF_TYPE == 2.0) then
            ! Over snow/ice
            ICEFRCT_M  = 0.00
            if ( TEMP <= iT_ICE_ALL ) then
                ICEFRCT_M = 1.000
            else if ( (TEMP > iT_ICE_ALL) .AND. (TEMP <= iT_ICE_MAX) ) then
                ICEFRCT_M = SIN( 0.5*MAPL_PI*( 1.00 - ( TEMP - iT_ICE_ALL ) / ( iT_ICE_MAX - iT_ICE_ALL ) ) )
            end if
            ICEFRCT_M = MIN(ICEFRCT_M,1.00)
            ICEFRCT_M = MAX(ICEFRCT_M,0.00)
            ICEFRCT_M = ICEFRCT_M**iICEFRPWR
        else if (SRF_TYPE > 1.0) then
            ! Over Land
            ICEFRCT_M  = 0.00
            if ( TEMP <= lT_ICE_ALL ) then
                ICEFRCT_M = 1.000
            else if ( (TEMP > lT_ICE_ALL) .AND. (TEMP <= lT_ICE_MAX) ) then
                ICEFRCT_M = SIN( 0.5*MAPL_PI*( 1.00 - ( TEMP - lT_ICE_ALL ) / ( lT_ICE_MAX - lT_ICE_ALL ) ) )
            end if
            ICEFRCT_M = MIN(ICEFRCT_M,1.00)
            ICEFRCT_M = MAX(ICEFRCT_M,0.00)
            ICEFRCT_M = ICEFRCT_M**lICEFRPWR
        else
            ! Over Oceans
            ICEFRCT_M  = 0.00
            if ( TEMP <= oT_ICE_ALL ) then
                ICEFRCT_M = 1.000
            else if ( (TEMP > oT_ICE_ALL) .AND. (TEMP <= oT_ICE_MAX) ) then
                ICEFRCT_M = SIN( 0.5*MAPL_PI*( 1.00 - ( TEMP - oT_ICE_ALL ) / ( oT_ICE_MAX - oT_ICE_ALL ) ) )
            end if
            ICEFRCT_M = MIN(ICEFRCT_M,1.00)
            ICEFRCT_M = MAX(ICEFRCT_M,0.00)
            ICEFRCT_M = ICEFRCT_M**oICEFRPWR
        endif
#endif
        ! Combine the Convective and MODIS functions
        ICEFRCT  = ICEFRCT_M*(1.0-CNV_FRACTION) + ICEFRCT_C*(CNV_FRACTION)
    end function ICE_FRACTION_SC

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