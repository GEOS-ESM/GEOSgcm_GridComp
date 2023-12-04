!
! Implementation of the Canadian Forest Fire 
! Canadian Forest Fire Weather Index System (CFFWI)
!
! References: 
! 1) Van Wagner, C.E. and  T.L. Picket, 1985.
! 2) Van Wagner, C.E., A method of computing fine fuel 
!    moisture content throughout the diurnal cycle, 1977.
! 3) Lawson, B.D. and Armitage, O.B, Weather guide 
!    for the Canadian Forest Fire Danger Rating System, 2008.
! 4) Kerry Anderson, A comparison of hourly fine fuel moisture 
!    code calculations within Canada, 2009.
! 5) Wotton, B.M., A grass moisture model for the 
!    Canadian Forest Fire Danger Rating System, 2009.
! 6) Wotton, B.M., Alexander, M.E. and Taylor, S.W., Updates and 
!    Revisions to the 1992 Canadian Forest Fire Behavior 
!    Prediction System, GLC-X-10, 2009.
!
! Milestones:
! Anton Darmenov, NASA, 2023 -- add grass FMC
! Anton Darmenov, NASA, 2023 -- add hourly FFMC
! Anton Darmenov, NASA, 2023 -- adjustments for applying the system globally
! Anton Darmenov, NASA, 2022 -- refactor
! Anton Darmenov, NASA, 2010 -- initial implementation
!

module cffwi

implicit none

private

public :: fine_fuel_moisture_code
public :: duff_moisture_code
public :: drought_code
public :: initial_spread_index
public :: buildup_index
public :: fire_weather_index
public :: daily_severity_rating

public :: grass_fuel_moisture_code

public :: drought_code_lf
public :: cffwi_indexes

real, public, parameter :: FFMC_INIT = 85.0
real, public, parameter :: DMC_INIT  =  6.0
real, public, parameter :: DC_INIT   = 15.0

real, public, parameter :: CFFWI_REFERENCE_LATITUDE = 46.0 ! Canada, 46N

real, public, parameter :: NOMINAL_FINE_FUEL_LOAD = 0.35 ! kg m-2

contains


elemental real function ff_scale_mc(ffmc)

    !
    ! Converts Fine Fuel Moisture Code (FFMC) to moisture content (%).
    !
    ! Default values and units:
    !     ffmc = FFMC (default = 85)
    
    implicit none

    real, intent(in) :: ffmc

    ff_scale_mc = 147.27723*(101 - ffmc) / (59.5 + ffmc)

end function ff_scale_mc



elemental real function ff_scale_ffmc(mc)

    !
    ! Converts moisture content (%) to Fine Fuel Moisture Code (FFMC).
    !
    ! Default values and units:
    !     mc = 16.31 (corresponds to FFMC=85 on the FF scale)
    
    implicit none

    real, intent(in) :: mc

    ff_scale_ffmc = 59.5*(250 - mc)/(147.27723 + mc)

end function ff_scale_ffmc



elemental real function fine_fuel_moisture_code(ffmc, T, RH, wind, Pr, dt)
  
    !
    ! Calculates either hourly or daily Fine Fuel Moisture Code (FFMC).
    !
    ! Default values and units:
    !     ffmc = initial FFMC (default = 85)
    !     T    = temperature, C
    !     RH   = relative humidity, %
    !     wind = wind speed, m/s
    !     Pr   = precip, mm;     | dt <= 1hr: precip over the time step
    !                            ! dt  > 1hr: 24-hour precip
    !     dt   = time step, hr;  | dt <= 1hr: trigers the hourly FFMC model
    !                            | dt  > 1hr: trigers the daily  FFMC model
    !
    !
    ! Note:
    !     Weather data are observed either hourly or at local noon.
    !
    
    implicit none

    real, intent(in) :: ffmc, T, RH, wind, Pr
    real, intent(in) :: dt

    ! local
    real :: f_0, w, m_0, h
    real :: r_f, m_r, dm, e_d, e_w, m
    real :: k_0, k_l, k_d, k_w, f_k, f_t
    real :: E_H_term, E_T_term
    real :: result

    real, parameter :: k_factor_hourly = 0.0579
    real, parameter :: k_factor_daily  = 0.581


    ! use the same variable names as in the FFMC equation and 
    ! convert units if necessary
    f_0 = ffmc
    w = 3.6 * wind  ! convert from m/s to km/h


    ! initial fuel moisture content (FF scale)
    m_0 = ff_scale_mc(f_0)
    
    ! current fuel moisture content
    m_r = m_0

    if (dt <= 1.0) then
        f_k = k_factor_hourly
        f_t = dt

        ! canopy effect is not considered in the hourly calculations
        r_f = Pr
    else
        f_k = k_factor_daily
        f_t = 1.0
        
        ! rainfall correction due to canopy effect
        r_f = max(0.0, Pr - 0.5)
    end if

    ! rainfall effect
    if (r_f > tiny(r_f)) then
        m_r = m_0 + 42.5 * r_f * exp(-100.0/(251.0 - m_0)) * (1 - exp(-6.93/r_f))

        if (m_0 > 150) then
            dm = m_0 - 150
            m_r = m_r + 0.0015 * (dm*dm) * sqrt(r_f)
        end if

        ! fuel moisture content has upper limit of 250
        if (m_r > 250) then
            m_r = 250.0
        end if
    end if

    ! equilibrium moisture contents for drying (E_d) and wetting (E_w) conditions
    E_H_term = exp(0.1 * (RH - 100.0))
    E_T_term = 0.18*(21.1 - T) * (1 - exp(-0.115*RH))
    E_d = 0.942*(RH**0.679) + 11*E_H_term + E_T_term
    E_w = 0.618*(RH**0.753) + 10*E_H_term + E_T_term

    h = RH / 100.0

    if (m_r > E_d) then
        ! drying is in effect
        k_0 = 0.424 * (1 - (h**1.7)) + 0.0694 * sqrt(w) * (1 - (h**8))
        k_d = k_0 * f_k * exp(0.0365 * T)

        m = E_d + (m_r - E_d)*(10.0**(-k_d*f_t))
    else
        if (m_r < E_w) then
            ! wetting is in effect
            k_l = 0.424 * (1 - (1 - h)**1.7) + 0.0694 * sqrt(w)* (1 - (1 - h)**8)
            k_w = k_l * f_k * exp(0.0365 * T)

            m = E_w - (E_w - m_r)*(10.0**(-k_w*f_t))
        else
            ! maintain moisture
            m = m_r
        end if
    end if


    ! current FFMC (FF scale)
    result = ff_scale_ffmc(m)
    
    ! clamp FFMC within [0, 101]
    fine_fuel_moisture_code = max(0.0, min(101.0, result))

end function fine_fuel_moisture_code



elemental real function grass_fuel_moisture_code(gfmc, T, RH, wind, Pr, sw_down, ff_load, dt)
  
    !
    ! Calculates hourly Grass Fuel Moisture Code (GFMC).
    !
    ! Default values and units:
    !     gfmc    = initial GFMC (default = 85)
    !     T       = temperature, C
    !     RH      = relative humidity, %
    !     wind    = wind speed, m/s
    !     Pr      = precip, mm;     | dt <= 1hr: precip over the time step
    !                               ! dt  > 1hr: 24-hour precip
    !     SW_down = incident shortwave flux, W m-2
    !     ff_load = fuel load of the fine fuel layer, kg m-2 (default = 0.3 kg m-2)
    !     dt      = time step, hr
    !
    !
    ! Note:
    !     Weather data are observed hourly.
    !
    
    implicit none

    real, intent(in) :: gfmc, T, RH, wind, Pr, sw_down, ff_load
    real, intent(in) :: dt

    ! local
    real :: f_0, w, m_0, h, I_sol
    real :: r_f, m_r, e_d, e_w, m
    real :: k_0, k_l, k_d, k_w, f_k, f_t
    real :: E_T_term
    real :: T_fuel, RH_fuel, svp, svp_fuel
    real :: result

    real, parameter :: k_factor_hourly = 0.389633


    
    ! convert units if necessary
    f_0 = gfmc
    w = 3.6 * wind         ! convert from m/s to km/h
    I_sol = 1e-3 * sw_down ! convert from W m-2 to kW m-2

    ! initial fuel moisture content (FF scale)
    m_0 = ff_scale_mc(f_0)
    
    ! current fuel moisture content
    m_r = m_0

    f_k = k_factor_hourly
    f_t = dt

    ! canopy effect is not considered in the hourly calculations
    r_f = Pr

    ! rainfall effect
    m_r = m_r + 100 * (r_f / ff_load)
    
    ! fuel moisture content has upper limit of 250
    m_r = min(m_r, 250.0)

    ! fuel temperature
    T_fuel = T + 35.07 * I_sol * exp(-0.06215 * w)

    ! saturation vapor pressure
    svp = 6.108 * 10**(7.5 * T / (237.3 + T)) ! Tetens equation

    ! saturation vapor pressure for fuel temperature
    svp_fuel = 6.108 * 10**(7.5 * T_fuel / (237.3 + T_fuel))

    ! fuel level relative humidity
    RH_fuel = RH * (svp / svp_fuel)

    ! equilibrium moisture contents for drying (E_d) and wetting (E_w) conditions
    E_T_term = 0.27 * (26.7 - T_fuel) * (1 - exp(-0.115 * RH_fuel))
    E_d = 1.62*RH_fuel**0.532 + 13.7*exp((RH_fuel - 100) / 13.0) + E_T_term
    E_w = 1.42*RH_fuel**0.512 + 12.0*exp((RH_fuel - 100) / 18.0) + E_T_term

    h = RH_fuel / 100.0

    if (m_r > E_d) then
        ! drying is in effect
        k_0 = 0.424 * (1 - (h**1.7)) + 0.0694 * sqrt(w) * (1 - (h**8))
        k_d = k_0 * f_k * exp(0.0365 * T_fuel)

        m = E_d + (m_r - E_d)*(10.0**(-k_d*f_t))
    else
        if (m_r < E_w) then
            ! wetting is in effect
            k_l = 0.424 * (1 - (1 - h)**1.7) + 0.0694 * sqrt(w)* (1 - (1 - h)**8)
            k_w = k_l * f_k * exp(0.0365 * T_fuel)

            m = E_w - (E_w - m_r)*(10.0**(-k_w*f_t))
        else
            ! maintain moisture
            m = m_r
        end if
    end if


    ! current GFMC (FF scale)
    result = ff_scale_ffmc(m)
    
    ! clamp GFMC within [0, 101]
    grass_fuel_moisture_code = max(0.0, min(101.0, result))

end function grass_fuel_moisture_code



elemental real function duff_moisture_code(dmc, T, RH, Pr, month)

    !
    ! Calculates the Duff Moisture Code (DMC).
    !
    ! Default values and units:
    !     dmc = initial DMC (default = 6)
    !     T   = temperature, C
    !     RH  = relative humidity, %
    !     Pr  = 24-hour precipitation, mm
    !     month = [1, 12]
    !
    !  Note: 
    !     Weather data are measured at local noon.
    !

    implicit none

    real,    intent (in) :: dmc, T, RH, Pr
    integer, intent (in) :: month


    !local
    real, dimension (12), parameter ::                               &
        EFFECTIVE_DAY_LENGTH = (/ 6.5,  7.5,  9.0, 12.8, 13.9, 13.9, &
                                 12.4, 10.9,  9.4,  8.0,  7.0,  6.0 /)

    real :: p_0, L_e
    real :: p, M_0, b, M_r, p_r, k, r_e


    ! use the same variable names as in the DMC equation and 
    ! convert the units if necessary
    p_0 = dmc
    L_e = EFFECTIVE_DAY_LENGTH(month)

    
    p = p_0
    
    ! rainfall correction
    if (Pr > 1.5) then
        r_e = 0.92*Pr - 1.27

        M_0 = 20.0 + exp(5.6348 - p_0/43.43)

        if (p_0 <= 33) then
            b = 100/(0.5 + 0.3*p_0)
        else if (p_0 <= 65) then
            b = 14 - 1.3 * log(p_0)
        else ! p_0 > 65
            b = 6.2 * log(p_0) - 17.2
        end if

        M_r = M_0 + 1000*r_e/(48.77 + b*r_e)
        p_r = 244.72 - 43.43 * log(M_r - 20.0)

        ! p_r should be positive
        p = max(0.0, p_r)
    end if


    ! if T < -1.1C then k is set to 0.0
    k = 1.894e-6 * max(0.0, (T + 1.1)) * (100 - RH) * L_e

    ! current DMC
    duff_moisture_code = p + 100*k

end function duff_moisture_code



elemental real function drought_code(dc, T, Pr, latitude, month)

    !
    ! Calculates the Drought Code (DC).
    ! 
    ! Default values and units:
    !     dc    = initial DC (default = 15.0)
    !     T     = temperature, C
    !     Pr    = 24-hour precipitation, mm
    !     month = [1..12]
    !
    ! Note: 
    !     Weather data are measured at local noon.
    !

    implicit none

    real,    intent(in) :: dc, T, Pr
    real,    intent(in) :: latitude
    integer, intent(in) :: month


    ! local
    real :: d_0, L_f
    real :: d, Q_0, Q_r, r_d, d_r, V, result
    
    ! use the same variable names as in the DC equation and 
    ! convert the units if necessary
    d_0 = dc
    L_f = drought_code_lf(latitude, month)


    d = d_0

    ! rainfall correction
    if (Pr > 2.8) then
        r_d = 0.83 * Pr - 1.27

        Q_0 = 800 * exp(-d_0 / 400)
        Q_r = Q_0 + 3.937 * r_d

        d_r = 400 * log(800 / Q_r)

        ! d_r should be positive
        d = max(0.0, d_r)
    end if

    ! temperature correction
    if (T > -2.8) then
        V = 0.36*(T + 2.8) + L_f
    else
        V = 0.0
    end if


    ! current DC
    result = d + 0.5*V

    drought_code = max(0.0, result)

end function drought_code



elemental real function drought_code_lf(latitude, month)
    
    !
    ! Calculates monthly day length adjustment factors for 
    ! Drought Code (Lf). 
    ! 
    ! Includes latitude considerations in adapting the system 
    ! for global use -- can be applied for northern and southern 
    ! hemispheres. Based on Appendix 3 in Lawson and Armitage (2008)
    !

    implicit none

    real,    intent(in) :: latitude  ! degrees
    integer, intent(in) :: month

   
    ! local
    real, dimension(12), parameter :: &
        DAY_LENGTH_FACTOR = (/ -1.6, -1.6, -1.6, 0.9,  3.8,  5.8, &
                                6.4,  5.0,  2.4, 0.4, -1.6, -1.6 /)
   
    real, parameter :: DAY_LENGTH_FACTOR_EQUATOR = 1.4

    real :: L_f


    if (latitude > 10.0) then
        ! use the reference values (Canada) north of 10N 
        L_f = DAY_LENGTH_FACTOR(month)
    else if (latitude < -10.0) then
        ! reverse the standard values used in Canada for seasons
        ! in the southern hemisphere: NH Jul -> SH Jan, NH Aug -> SH Feb, etc.  
        L_f = DAY_LENGTH_FACTOR(mod(month+5, 12) + 1) 
    else
        ! for locations near the equator, from 10S to 10N, 
        ! use the mean DC day length adjustment value (Lf = 1.4) year-round
        L_f = DAY_LENGTH_FACTOR_EQUATOR
    end if

    drought_code_lf = L_f

end function drought_code_lf



elemental real function initial_spread_index(ffmc, wind)

    !
    ! Calculates the Initial Spread Index (ISI).
    ! 
    ! Default values and units:
    !     ffmc  = initial FFMC (default = 85.0)
    !     wind  = wind speed, m/s
    !
    ! Note:
    !     Weather data are measured at local noon.
    !

    implicit none

    real, intent (in) :: ffmc, wind

    ! local
    real :: w, m, fun_w, fun_f


    ! use the same variable names as in the DC equation and 
    ! convert the units if necessary
    w = 3.6 * wind  ! convert from m/s to km/h

    ! FF scale
    m = ff_scale_mc(ffmc)
   
    if (w < 40) then
        fun_w = exp(0.05039 * w)
    else
        ! modification at the extreme end of winds
        fun_w = 12 * (1 - exp(-0.0818 * (w - 28)))
    end if

    fun_f = 91.9 * exp(-0.1386 * m) * (1 + (m**5.31)/4.93e7)

    ! current ISI
    initial_spread_index = 0.208 * fun_w * fun_f

end function initial_spread_index



elemental real function buildup_index(dmc, dc)

    !
    ! Calculates the Buildup Index (BUI).
    ! 
    ! Default values and units:
    !     dmc = Duff Moisture Code (default = 6.5)
    !     dc  = Drought Code (default = 15.0)
    !

    implicit none

    real, intent (in) :: dmc, dc

    ! local
    real :: result

    if (dmc > 0 .and. dc > 0) then
        if (dmc > 0.4*dc) then
            result = dmc - (1 - 0.8*dc/(dmc + 0.4*dc)) * (0.92 + (0.0114*dmc)**1.7)
        else
            result = 0.8 * dmc * dc / (dmc + 0.4*dc)
        end if
    else
        result = 0.0
    end if

    ! current BUI
    buildup_index = max(0.0, result)

end function buildup_index



elemental real function fire_weather_index(isi, bui)

    !
    ! Calculates the Fire Weather Index (FWI).
    ! 
    ! Default values and units:
    !     isi = Initial Spread Index (default = 0.0)
    !     bui = Buildup Index (default = 0.0)
    !

    implicit none

    real, intent(in) :: isi, bui

    ! local
    real :: f, B, result


    if (bui > 80) then
        f = 1000/(25 + 108.64 * exp(-0.023 * bui))
    else
        f = 0.626 * (bui**0.809) + 2
    end if

    B = 0.1 * isi * f

    if (B > 1) then
        result = exp(2.72 * (0.434 * log(B))**0.647)
    else
        result = B
    end if

    ! current FWI
    fire_weather_index = result

end function fire_weather_index



elemental real function daily_severity_rating(fwi)

    !
    ! Calculates the Daily Severity Rating (DSR).
    ! 
    ! Default values and units:
    !     fwi = Fire Weather Index
    !

    implicit none

    real, intent(in) :: fwi

    ! current DSR
    daily_severity_rating = 0.0272 * fwi**1.77

end function daily_severity_rating



subroutine cffwi_indexes(ffmc_initial, dmc_initial, dc_initial, &
                         T, RH, wind, Pr,                       &
                         latitude, month, time_step,            &
                         ffmc, dmc, dc, isi, bui, fwi, dsr)

    ! Calculates FFMC, DMC, DC, ISI, BUI, FWI and DSR indexes.
    !
    ! Default values and units:
    !     ffmc_initial  = previous hour|day FFMC (default = 85)
    !     dmc_initial   = previous hour|day DMC  (default =  6)
    !     dc_initial    = previous hour|day DC   (default = 15)
    !     T             = temperature, C
    !     RH            = relative humidity, %
    !     wind          = wind speed, m/s
    !     Pr            = 24-hour precipitation, mm
    !     time_step     = time step, hr;  values <= 1, trigger hourly FFMC
    !                                     values  > 1, trigger daily  FFMC
    !
    ! Note:
    !     Weather data are measured either hourly or at local noon.
    !

    implicit none

    real, intent(in)    :: ffmc_initial, dmc_initial, dc_initial
    real, intent(in)    :: T, RH, wind, Pr
    integer, intent(in) :: month
    real, intent(in)    :: latitude
    real, intent(in)    :: time_step
    real, intent(out)   :: ffmc, dmc, dc, isi, bui, fwi, dsr


    ! update fuel moisture codes 
    ffmc = fine_fuel_moisture_code(ffmc_initial, T, RH, wind, Pr, time_step)
    dmc  = duff_moisture_code(dmc_initial, T, RH, Pr, month)
    dc   = drought_code(dc_initial, T, Pr, latitude, month)

    ! update fire behavior indexes
    isi = initial_spread_index(ffmc, wind)
    bui = buildup_index(dmc, dc)

    fwi = fire_weather_index(isi, bui)
    dsr = daily_severity_rating(fwi)

end subroutine cffwi_indexes

end module cffwi
