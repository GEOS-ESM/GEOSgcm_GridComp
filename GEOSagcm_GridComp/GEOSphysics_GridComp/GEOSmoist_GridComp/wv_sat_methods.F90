module wv_sat_methods

! This portable module contains all CAM methods for estimating
! the saturation vapor pressure of water.
!
! wv_saturation provides CAM-specific interfaces and utilities
! based on these formulae.
!
! Typical usage of this module:
!
! Init:
! call wv_sat_methods_init(r8, <constants>, errstring)
!
! Get scheme index from a name string:
! scheme_idx = wv_sat_get_scheme_idx(scheme_name)
! if (.not. wv_sat_valid_idx(scheme_idx)) <throw some error>
!
! Get pressures:
! es = wv_sat_svp_water(t, scheme_idx)
! es = wv_sat_svp_ice(t, scheme_idx)
!
! Use ice/water transition range:
! es = wv_sat_svp_trice(t, ttrice, scheme_idx)
!
! Note that elemental functions cannot be pointed to, nor passed
! as arguments. If you need to do either, it is recommended to
! wrap the function so that it can be given an explicit (non-
! elemental) interface.

implicit none
private
save

integer, parameter :: r8 = selected_real_kind(12) ! 8 byte real

real(r8) :: tmelt   ! Melting point of water at 1 atm (K)
real(r8) :: h2otrip ! Triple point temperature of water (K)
real(r8) :: tboil   ! Boiling point of water at 1 atm (K)

real(r8) :: ttrice  ! Ice-water transition range

real(r8) :: epsilo  ! Ice-water transition range
real(r8) :: omeps   ! 1._r8 - epsilo

! Indices representing individual schemes
integer, parameter :: Invalid_idx = -1
integer, parameter :: OldGoffGratch_idx = 0
integer, parameter :: GoffGratch_idx = 1
integer, parameter :: MurphyKoop_idx = 2
integer, parameter :: Bolton_idx = 3

! Index representing the current default scheme.
integer, parameter :: initial_default_idx = GoffGratch_idx
integer :: default_idx = initial_default_idx

public wv_sat_methods_init
public wv_sat_get_scheme_idx
public wv_sat_valid_idx

public wv_sat_set_default
public wv_sat_reset_default

public wv_sat_svp_water
public wv_sat_svp_ice
public wv_sat_svp_trans

! pressure -> humidity conversion
public wv_sat_svp_to_qsat

! Combined qsat operations
public wv_sat_qsat_water
public wv_sat_qsat_ice
public wv_sat_qsat_trans


public shr_spfn_gamma !DONIF

contains

!---------------------------------------------------------------------
! ADMINISTRATIVE FUNCTIONS
!---------------------------------------------------------------------

! Get physical constants
subroutine wv_sat_methods_init(kind, tmelt_in, h2otrip_in, tboil_in, &
     ttrice_in, epsilo_in, errstring)
  integer, intent(in) :: kind
  real(r8), intent(in) :: tmelt_in
  real(r8), intent(in) :: h2otrip_in
  real(r8), intent(in) :: tboil_in
  real(r8), intent(in) :: ttrice_in
  real(r8), intent(in) :: epsilo_in
  character(len=*), intent(out)  :: errstring

  errstring = ' '

  if (kind /= r8) then
     write(errstring,*) 'wv_sat_methods_init: ERROR: ', &
          kind,' was input kind but ',r8,' is internal kind.'
     return
  end if

  if (ttrice_in < 0._r8) then
     write(errstring,*) 'wv_sat_methods_init: ERROR: ', &
          ttrice_in,' was input for ttrice, but negative range is invalid.'
     return
  end if

  tmelt = tmelt_in
  h2otrip = h2otrip_in
  tboil = tboil_in
  ttrice = ttrice_in
  epsilo = epsilo_in

  omeps = 1._r8 - epsilo

end subroutine wv_sat_methods_init

! Look up index by name.
pure function wv_sat_get_scheme_idx(name) result(idx)
  character(len=*), intent(in) :: name
  integer :: idx
  
  select case (name)
  case("GoffGratch")
     idx = GoffGratch_idx
  case("MurphyKoop")
     idx = MurphyKoop_idx
  case("OldGoffGratch")
     idx = OldGoffGratch_idx
  case("Bolton")
     idx = Bolton_idx
  case default
     idx = Invalid_idx
  end select

end function wv_sat_get_scheme_idx

! Check validity of an index from the above routine.
pure function wv_sat_valid_idx(idx) result(status)
  integer, intent(in) :: idx
  logical :: status

  status = (idx /= Invalid_idx)

end function wv_sat_valid_idx

! Set default scheme (otherwise, Goff & Gratch is default)
! Returns a logical representing success (.true.) or
! failure (.false.).
function wv_sat_set_default(name) result(status)
  character(len=*), intent(in) :: name
  logical :: status

  ! Don't want to overwrite valid default with invalid,
  ! so assign to temporary and check it first.
  integer :: tmp_idx

  tmp_idx = wv_sat_get_scheme_idx(name)

  status = wv_sat_valid_idx(tmp_idx)

  if (status) default_idx = tmp_idx

end function wv_sat_set_default

! Reset default scheme to initial value.
! The same thing can be accomplished with wv_sat_set_default;
! the real reason to provide this routine is to reset the
! module for testing purposes.
subroutine wv_sat_reset_default()

  default_idx = initial_default_idx

end subroutine wv_sat_reset_default

!---------------------------------------------------------------------
! UTILITIES
!---------------------------------------------------------------------

! Get saturation specific humidity given pressure and SVP.
! Specific humidity is limited to range 0-1.
elemental function wv_sat_svp_to_qsat(es, p) result(qs)

  real(r8), intent(in) :: es  ! SVP
  real(r8), intent(in) :: p   ! Current pressure.
  real(r8) :: qs

  ! If pressure is less than SVP, set qs to maximum of 1.
  if ( (p - es) <= 0._r8 ) then
     qs = 1.0_r8
  else
     qs = epsilo*es / (p - omeps*es)
  end if

end function wv_sat_svp_to_qsat

elemental subroutine wv_sat_qsat_water(t, p, es, qs, idx)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate SVP over water at a given temperature, and then      !
  !   calculate and return saturation specific humidity.             !
  !------------------------------------------------------------------!

  ! Inputs
  real(r8), intent(in) :: t    ! Temperature
  real(r8), intent(in) :: p    ! Pressure
  ! Outputs
  real(r8), intent(out) :: es  ! Saturation vapor pressure
  real(r8), intent(out) :: qs  ! Saturation specific humidity

  integer,  intent(in), optional :: idx ! Scheme index

  es = wv_sat_svp_water(t, idx)

  qs = wv_sat_svp_to_qsat(es, p)

  ! Ensures returned es is consistent with limiters on qs.
  es = min(es, p)

end subroutine wv_sat_qsat_water

elemental subroutine wv_sat_qsat_ice(t, p, es, qs, idx)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate SVP over ice at a given temperature, and then        !
  !   calculate and return saturation specific humidity.             !
  !------------------------------------------------------------------!

  ! Inputs
  real(r8), intent(in) :: t    ! Temperature
  real(r8), intent(in) :: p    ! Pressure
  ! Outputs
  real(r8), intent(out) :: es  ! Saturation vapor pressure
  real(r8), intent(out) :: qs  ! Saturation specific humidity

  integer,  intent(in), optional :: idx ! Scheme index

  es = wv_sat_svp_ice(t, idx)

  qs = wv_sat_svp_to_qsat(es, p)

  ! Ensures returned es is consistent with limiters on qs.
  es = min(es, p)

end subroutine wv_sat_qsat_ice

elemental subroutine wv_sat_qsat_trans(t, p, es, qs, idx)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate SVP over ice at a given temperature, and then        !
  !   calculate and return saturation specific humidity.             !
  !------------------------------------------------------------------!

  ! Inputs
  real(r8), intent(in) :: t    ! Temperature
  real(r8), intent(in) :: p    ! Pressure
  ! Outputs
  real(r8), intent(out) :: es  ! Saturation vapor pressure
  real(r8), intent(out) :: qs  ! Saturation specific humidity

  integer,  intent(in), optional :: idx ! Scheme index

  es = wv_sat_svp_trans(t, idx)

  qs = wv_sat_svp_to_qsat(es, p)

  ! Ensures returned es is consistent with limiters on qs.
  es = min(es, p)

end subroutine wv_sat_qsat_trans

!---------------------------------------------------------------------
! SVP INTERFACE FUNCTIONS
!---------------------------------------------------------------------

elemental function wv_sat_svp_water(t, idx) result(es)
  real(r8), intent(in) :: t
  integer,  intent(in), optional :: idx
  real(r8) :: es

  integer :: use_idx

  if (present(idx)) then
     use_idx = idx
  else
     use_idx = default_idx
  end if

  select case (use_idx)
  case(GoffGratch_idx)
     es = GoffGratch_svp_water(t)
  case(MurphyKoop_idx)
     es = MurphyKoop_svp_water(t)
  case(OldGoffGratch_idx)
     es = OldGoffGratch_svp_water(t)
  case(Bolton_idx)
     es = Bolton_svp_water(t)
  end select

end function wv_sat_svp_water

elemental function wv_sat_svp_ice(t, idx) result(es)
  real(r8), intent(in) :: t
  integer,  intent(in), optional :: idx
  real(r8) :: es

  integer :: use_idx

  if (present(idx)) then
     use_idx = idx
  else
     use_idx = default_idx
  end if

  select case (use_idx)
  case(GoffGratch_idx)
     es = GoffGratch_svp_ice(t)
  case(MurphyKoop_idx)
     es = MurphyKoop_svp_ice(t)
  case(OldGoffGratch_idx)
     es = OldGoffGratch_svp_ice(t)
  case(Bolton_idx)
     es = Bolton_svp_water(t)
  end select

end function wv_sat_svp_ice

elemental function wv_sat_svp_trans(t, idx) result (es)

  real(r8), intent(in) :: t
  integer,  intent(in), optional :: idx

  real(r8) :: es

  real(r8) :: esice      ! Saturation vapor pressure over ice
  real(r8) :: weight     ! Intermediate scratch variable for es transition

!
! Water
!
  if (t >= (tmelt - ttrice)) then
     es = wv_sat_svp_water(t,idx)
  else
     es = 0.0_r8
  end if

!
! Ice
!
  if (t < tmelt) then

     esice = wv_sat_svp_ice(t,idx)

     if ( (tmelt - t) > ttrice ) then
        weight = 1.0_r8
     else
        weight = (tmelt - t)/ttrice
     end if

     es = weight*esice + (1.0_r8 - weight)*es
  end if

end function wv_sat_svp_trans

!---------------------------------------------------------------------
! SVP METHODS
!---------------------------------------------------------------------

! Goff & Gratch (1946)

elemental function GoffGratch_svp_water(t) result(es)
  real(r8), intent(in) :: t  ! Temperature in Kelvin
  real(r8) :: es             ! SVP in Pa

  ! uncertain below -70 C
  es = 10._r8**(-7.90298_r8*(tboil/t-1._r8)+ &
       5.02808_r8*log10(tboil/t)- &
       1.3816e-7_r8*(10._r8**(11.344_r8*(1._r8-t/tboil))-1._r8)+ &
       8.1328e-3_r8*(10._r8**(-3.49149_r8*(tboil/t-1._r8))-1._r8)+ &
       log10(1013.246_r8))*100._r8

end function GoffGratch_svp_water

elemental function GoffGratch_svp_ice(t) result(es)
  real(r8), intent(in) :: t  ! Temperature in Kelvin
  real(r8) :: es             ! SVP in Pa

  ! good down to -100 C
  es = 10._r8**(-9.09718_r8*(h2otrip/t-1._r8)-3.56654_r8* &
       log10(h2otrip/t)+0.876793_r8*(1._r8-t/h2otrip)+ &
       log10(6.1071_r8))*100._r8

end function GoffGratch_svp_ice

! Murphy & Koop (2005)

elemental function MurphyKoop_svp_water(t) result(es)
  real(r8), intent(in) :: t  ! Temperature in Kelvin
  real(r8) :: es             ! SVP in Pa

  ! (good for 123 < T < 332 K)
  es = exp(54.842763_r8 - (6763.22_r8 / t) - (4.210_r8 * log(t)) + &
       (0.000367_r8 * t) + (tanh(0.0415_r8 * (t - 218.8_r8)) * &
       (53.878_r8 - (1331.22_r8 / t) - (9.44523_r8 * log(t)) + &
       0.014025_r8 * t)))

end function MurphyKoop_svp_water

elemental function MurphyKoop_svp_ice(t) result(es)
  real(r8), intent(in) :: t  ! Temperature in Kelvin
  real(r8) :: es             ! SVP in Pa

  ! (good down to 110 K)
  es = exp(9.550426_r8 - (5723.265_r8 / t) + (3.53068_r8 * log(t)) &
       - (0.00728332_r8 * t))

end function MurphyKoop_svp_ice

! Old CAM implementation, also labelled Goff & Gratch (1946)

! The water formula differs only due to compiler-dependent order of
! operations, so differences are roundoff level, usually 0.

! The ice formula gives fairly close answers to the current
! implementation, but has been rearranged, and uses the
! 1 atm melting point of water as the triple point.
! Differences are thus small but above roundoff.

! A curious fact: although using the melting point of water was
! probably a mistake, it mildly improves accuracy for ice svp,
! since it compensates for a systematic error in Goff & Gratch.

elemental function OldGoffGratch_svp_water(t) result(es)
  real(r8), intent(in) :: t
  real(r8) :: es
  real(r8) :: ps, e1, e2, f1, f2, f3, f4, f5, f

  ps = 1013.246_r8
  e1 = 11.344_r8*(1.0_r8 - t/tboil)
  e2 = -3.49149_r8*(tboil/t - 1.0_r8)
  f1 = -7.90298_r8*(tboil/t - 1.0_r8)
  f2 = 5.02808_r8*log10(tboil/t)
  f3 = -1.3816_r8*(10.0_r8**e1 - 1.0_r8)/10000000.0_r8
  f4 = 8.1328_r8*(10.0_r8**e2 - 1.0_r8)/1000.0_r8
  f5 = log10(ps)
  f  = f1 + f2 + f3 + f4 + f5

  es = (10.0_r8**f)*100.0_r8
  
end function OldGoffGratch_svp_water

elemental function OldGoffGratch_svp_ice(t) result(es)
  real(r8), intent(in) :: t
  real(r8) :: es
  real(r8) :: term1, term2, term3

  term1 = 2.01889049_r8/(tmelt/t)
  term2 = 3.56654_r8*log(tmelt/t)
  term3 = 20.947031_r8*(tmelt/t)

  es = 575.185606e10_r8*exp(-(term1 + term2 + term3))
  
end function OldGoffGratch_svp_ice

! Bolton (1980)
! zm_conv deep convection scheme contained this SVP calculation.
! It appears to be from D. Bolton, 1980, Monthly Weather Review.
! Unlike the other schemes, no distinct ice formula is associated
! with it. (However, a Bolton ice formula exists in CLUBB.)

! The original formula used degrees C, but this function
! takes Kelvin and internally converts.

elemental function Bolton_svp_water(t) result(es)
  real(r8),parameter :: c1 = 611.2_r8
  real(r8),parameter :: c2 = 17.67_r8
  real(r8),parameter :: c3 = 243.5_r8

  real(r8), intent(in) :: t  ! Temperature in Kelvin
  real(r8) :: es             ! SVP in Pa

  es = c1*exp( (c2*(t - tmelt))/((t - tmelt)+c3) )

end function Bolton_svp_water


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   pure  FUNCTION shr_spfn_gamma(X)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!D    DOUBLE PRECISION FUNCTION DGAMMA(X)
!----------------------------------------------------------------------
!
! THIS ROUTINE CALCULATES THE GAMMA FUNCTION FOR A REAL ARGUMENT X.
!   COMPUTATION IS BASED ON AN ALGORITHM OUTLINED IN REFERENCE 1.
!   THE PROGRAM USES RATIONAL FUNCTIONS THAT APPROXIMATE THE GAMMA
!   FUNCTION TO AT LEAST 20 SIGNIFICANT DECIMAL DIGITS.  COEFFICIENTS
!   FOR THE APPROXIMATION OVER THE INTERVAL (1,2) ARE UNPUBLISHED.
!   THOSE FOR THE APPROXIMATION FOR X .GE. 12 ARE FROM REFERENCE 2.
!   THE ACCURACY ACHIEVED DEPENDS ON THE ARITHMETIC SYSTEM, THE
!   COMPILER, THE INTRINSIC FUNCTIONS, AND PROPER SELECTION OF THE
!   MACHINE-DEPENDENT CONSTANTS.
!
!
!*******************************************************************
!*******************************************************************
!
! EXPLANATION OF MACHINE-DEPENDENT CONSTANTS
!
! BETA   - RADIX FOR THE FLOATING-POINT REPRESENTATION
! MAXEXP - THE SMALLEST POSITIVE POWER OF BETA THAT OVERFLOWS
! XBIG   - THE LARGEST ARGUMENT FOR WHICH GAMMA(X) IS REPRESENTABLE
!          IN THE MACHINE, I.E., THE SOLUTION TO THE EQUATION
!                  GAMMA(XBIG) = BETA**MAXEXP
! XINF   - THE LARGEST MACHINE REPRESENTABLE FLOATING-POINT NUMBER;
!          APPROXIMATELY BETA**MAXEXP
! EPS    - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!          1.0+EPS .GT. 1.0
! XMININ - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!          1/XMININ IS MACHINE REPRESENTABLE
!
!     APPROXIMATE VALUES FOR SOME IMPORTANT MACHINES ARE:
!
!                            BETA       MAXEXP        XBIG
!
! CRAY-1         (S.P.)        2         8191        966.961
! CYBER 180/855
!   UNDER NOS    (S.P.)        2         1070        177.803
! IEEE (IBM/XT,
!   SUN, ETC.)   (S.P.)        2          128        35.040
! IEEE (IBM/XT,
!   SUN, ETC.)   (D.P.)        2         1024        171.624
! IBM 3033       (D.P.)       16           63        57.574
! VAX D-FORMAT   (D.P.)        2          127        34.844
! VAX G-FORMAT   (D.P.)        2         1023        171.489
!
!                            XINF         EPS        XMININ
!
! CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
! CYBER 180/855
!   UNDER NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
! IEEE (IBM/XT,
!   SUN, ETC.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
! IEEE (IBM/XT,
!   SUN, ETC.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
! IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
! VAX D-FORMAT   (D.P.)   1.70D+38     1.39D-17    5.88D-39
! VAX G-FORMAT   (D.P.)   8.98D+307    1.11D-16    1.12D-308
!
!*******************************************************************
!*******************************************************************
!
! ERROR RETURNS
!
!  THE PROGRAM RETURNS THE VALUE XINF FOR SINGULARITIES OR
!     WHEN OVERFLOW WOULD OCCUR.  THE COMPUTATION IS BELIEVED
!     TO BE FREE OF UNDERFLOW AND OVERFLOW.
!
!
!  INTRINSIC FUNCTIONS REQUIRED ARE:
!
!     INT, DBLE, EXP, LOG, REAL, SIN
!
!
! REFERENCES:  AN OVERVIEW OF SOFTWARE DEVELOPMENT FOR SPECIAL
!              FUNCTIONS   W. J. CODY, LECTURE NOTES IN MATHEMATICS,
!              506, NUMERICAL ANALYSIS DUNDEE, 1975, G. A. WATSON
!              (ED.), SPRINGER VERLAG, BERLIN, 1976.
!
!              COMPUTER APPROXIMATIONS, HART, ET. AL., WILEY AND
!              SONS, NEW YORK, 1968.
!
!  LATEST MODIFICATION: OCTOBER 12, 1989
!
!  AUTHORS: W. J. CODY AND L. STOLTZ
!           APPLIED MATHEMATICS DIVISION
!           ARGONNE NATIONAL LABORATORY
!           ARGONNE, IL 60439
!
!----------------------------------------------------------------------
      INTEGER I,N
      LOGICAL PARITY

      real(r8) shr_spfn_gamma
      real(r8), intent(in) :: X
      REAL(r8) &
!D    DOUBLE PRECISION
         C,CONV,EPS,FACT,HALF,ONE,P,PI,Q,RES,SQRTPI,SUM,TWELVE, &
         TWO, XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
      DIMENSION C(7),P(8),Q(8)
            
      !----------------------------------------------------------------------
!  STATEMENT FUNCTIONS FOR CONVERSION BETWEEN INTEGER AND FLOAT
!----------------------------------------------------------------------
      CONV(I) = REAL(I,r8)
!D    CONV(I) = DBLE(I)
      
      
!----------------------------------------------------------------------
!  MATHEMATICAL CONSTANTS
!----------------------------------------------------------------------
    !  DATA ONE,HALF,TWELVE,TWO,ZERO/1.0E0_r8,0.5E0_r8,12.0E0_r8,2.0E0_r8,0.0E0_r8/, &
    !      SQRTPI/0.9189385332046727417803297E0_r8/, &
    !      PI/3.1415926535897932384626434E0_r8/
!D    DATA ONE,HALF,TWELVE,TWO,ZERO/1.0D0,0.5D0,12.0D0,2.0D0,0.0D0/,
!D   1     SQRTPI/0.9189385332046727417803297D0/,
!D   2     PI/3.1415926535897932384626434D0/
!----------------------------------------------------------------------
!  MACHINE DEPENDENT PARAMETERS
!----------------------------------------------------------------------
 !     DATA XBIG,XMININ,EPS/35.040E0_r8,1.18E-38_r8,1.19E-7_r8/, &
 !         XINF/3.4E38_r8/
!D    DATA XBIG,XMININ,EPS/171.624D0,2.23D-308,2.22D-16/,
!D   1     XINF/1.79D308/
!----------------------------------------------------------------------
!  NUMERATOR AND DENOMINATOR COEFFICIENTS FOR RATIONAL MINIMAX
!     APPROXIMATION OVER (1,2).
!----------------------------------------------------------------------
      !DATA P/-1.71618513886549492533811E+0_r8,2.47656508055759199108314E+1_r8,&
      !      -3.79804256470945635097577E+2_r8,6.29331155312818442661052E+2_r8,&
      !      8.66966202790413211295064E+2_r8,-3.14512729688483675254357E+4_r8,&
      !      -3.61444134186911729807069E+4_r8,6.64561438202405440627855E+4_r8/
   !   DATA Q/-3.08402300119738975254353E+1_r8,3.15350626979604161529144E+2_r8,&
   !        -1.01515636749021914166146E+3_r8,-3.10777167157231109440444E+3_r8,&
   !          2.25381184209801510330112E+4_r8,4.75584627752788110767815E+3_r8,&
   !        -1.34659959864969306392456E+5_r8,-1.15132259675553483497211E+5_r8/
!D    DATA P/-1.71618513886549492533811D+0,2.47656508055759199108314D+1,
!D   1       -3.79804256470945635097577D+2,6.29331155312818442661052D+2,
!D   2       8.66966202790413211295064D+2,-3.14512729688483675254357D+4,
!D   3       -3.61444134186911729807069D+4,6.64561438202405440627855D+4/
!D    DATA Q/-3.08402300119738975254353D+1,3.15350626979604161529144D+2,
!D   1      -1.01515636749021914166146D+3,-3.10777167157231109440444D+3,
!D   2        2.25381184209801510330112D+4,4.75584627752788110767815D+3,
!D   3      -1.34659959864969306392456D+5,-1.15132259675553483497211D+5/
!----------------------------------------------------------------------
!  COEFFICIENTS FOR MINIMAX APPROXIMATION OVER (12, INF).
!----------------------------------------------------------------------
   !   DATA C/-1.910444077728E-03_r8,8.4171387781295E-04_r8, &
   !       -5.952379913043012E-04_r8,7.93650793500350248E-04_r8,&
   !       -2.777777777777681622553E-03_r8,8.333333333333333331554247E-02_r8,&
   !        5.7083835261E-03_r8/
!D    DATA C/-1.910444077728D-03,8.4171387781295D-04,
!D   1     -5.952379913043012D-04,7.93650793500350248D-04,
!D   2     -2.777777777777681622553D-03,8.333333333333333331554247D-02,
!D   3      5.7083835261D-03/



! DATA ONE,HALF,TWELVE,TWO,ZERO/1.0E0_r8,0.5E0_r8,12.0E0_r8,2.0E0_r8,0.0E0_r8/, &
 !         SQRTPI/0.9189385332046727417803297E0_r8/, &
 !         PI/3.1415926535897932384626434E0_r8/
	ONE  = 1.0E0_r8
	HALF = 0.5E0_r8
	TWELVE= 12.0E0_r8
	TWO = 2.0E0_r8
	ZERO= 0.0E0_r8
	SQRTPI= 0.9189385332046727417803297E0_r8
	PI= 3.1415926535897932384626434E0_r8
	
	  
	  XBIG =35.040E0_r8
	  XMININ = 1.18E-38_r8
	  EPS = 1.19E-7_r8
	  XINF = 3.4E38_r8 
	
    P = (/-1.71618513886549492533811E+0_r8,2.47656508055759199108314E+1_r8, &
            -3.79804256470945635097577E+2_r8,6.29331155312818442661052E+2_r8, &
            8.66966202790413211295064E+2_r8,-3.14512729688483675254357E+4_r8,  &
            -3.61444134186911729807069E+4_r8,6.64561438202405440627855E+4_r8/)
      
   Q=(/-3.08402300119738975254353E+1_r8,3.15350626979604161529144E+2_r8,&
           -1.01515636749021914166146E+3_r8,-3.10777167157231109440444E+3_r8,&
             2.25381184209801510330112E+4_r8,4.75584627752788110767815E+3_r8,&
           -1.34659959864969306392456E+5_r8,-1.15132259675553483497211E+5_r8/)
	   
   C= (/-1.910444077728E-03_r8,8.4171387781295E-04_r8, &
          -5.952379913043012E-04_r8,7.93650793500350248E-04_r8,&
          -2.777777777777681622553E-03_r8,8.333333333333333331554247E-02_r8,&
           5.7083835261E-03_r8/)




      PARITY=.FALSE.
      FACT=ONE
      N=0
      Y=X
      IF(Y.LE.ZERO)THEN
!----------------------------------------------------------------------
!  ARGUMENT IS NEGATIVE
!----------------------------------------------------------------------
        Y=-X
        Y1=AINT(Y)
        RES=Y-Y1
        IF(RES.NE.ZERO)THEN
          IF(Y1.NE.AINT(Y1*HALF)*TWO)PARITY=.TRUE.
          FACT=-PI/SIN(PI*RES)
          Y=Y+ONE
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ENDIF
!----------------------------------------------------------------------
!  ARGUMENT IS POSITIVE
!----------------------------------------------------------------------
      IF(Y.LT.EPS)THEN
!----------------------------------------------------------------------
!  ARGUMENT .LT. EPS
!----------------------------------------------------------------------
        IF(Y.GE.XMININ)THEN
          RES=ONE/Y
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ELSEIF(Y.LT.TWELVE)THEN
        Y1=Y
        IF(Y.LT.ONE)THEN
!----------------------------------------------------------------------
!  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
          Z=Y
          Y=Y+ONE
        ELSE
!----------------------------------------------------------------------
!  1.0 .LT. ARGUMENT .LT. 12.0, REDUCE ARGUMENT IF NECESSARY
!----------------------------------------------------------------------
          N=INT(Y)-1
          Y=Y-CONV(N)
          Z=Y-ONE
        ENDIF
!----------------------------------------------------------------------
!  EVALUATE APPROXIMATION FOR 1.0 .LT. ARGUMENT .LT. 2.0
!----------------------------------------------------------------------
        XNUM=ZERO
        XDEN=ONE
        DO 260 I=1,8
          XNUM=(XNUM+P(I))*Z
          XDEN=XDEN*Z+Q(I)
  260   CONTINUE
        RES=XNUM/XDEN+ONE
        IF(Y1.LT.Y)THEN
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
          RES=RES/Y1
        ELSEIF(Y1.GT.Y)THEN
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  2.0 .LT. ARGUMENT .LT. 12.0
!----------------------------------------------------------------------
          DO 290 I=1,N
            RES=RES*Y
            Y=Y+ONE
  290     CONTINUE
        ENDIF
      ELSE
!----------------------------------------------------------------------
!  EVALUATE FOR ARGUMENT .GE. 12.0,
!----------------------------------------------------------------------
        IF(Y.LE.XBIG)THEN
          YSQ=Y*Y
          SUM=C(7)
          DO 350 I=1,6
            SUM=SUM/YSQ+C(I)
  350     CONTINUE
          SUM=SUM/Y-Y+SQRTPI
          SUM=SUM+(Y-HALF)*LOG(Y)
          RES=EXP(SUM)
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ENDIF
!----------------------------------------------------------------------
!  FINAL ADJUSTMENTS AND RETURN
!----------------------------------------------------------------------
      IF(PARITY)RES=-RES
      IF(FACT.NE.ONE)RES=FACT/RES
  900 shr_spfn_gamma=RES
!D900 DGAMMA = RES
      RETURN
! ---------- LAST LINE OF GAMMA ----------
      END function shr_spfn_gamma

end module wv_sat_methods
