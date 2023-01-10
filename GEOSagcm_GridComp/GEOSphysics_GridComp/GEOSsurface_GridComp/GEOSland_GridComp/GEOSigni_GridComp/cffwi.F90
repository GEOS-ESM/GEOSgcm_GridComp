!
! cffwi.f - Implementation of the Canadian Forest Fire Weather 
!           Index (CFFWI), Van Wagner, C.E. and  T.L. Picket, 1985
!
!
! 
! Anton Darmenov, NASA GSFC, 2021  -- refactor
! Anton Darmenov, NASA GSFC, 2010  -- initial implementation
!

module cffwi

implicit none

private

public:: fine_fuel_moisture_code
public:: duff_moisture_code
public:: drought_code
public:: initial_spread_index
public:: buildup_index
public:: fire_weather_index
public:: daily_severity_rating

public:: cffwi_indexes

real, public, parameter :: FFMC_INIT = 85.0
real, public, parameter :: DMC_INIT  =  6.0
real, public, parameter :: DC_INIT   = 15.0

contains

elemental real function fine_fuel_moisture_code(ffmc, T, RH, wind, Pr)
  
    !
    ! Calculates the Fine Fuel Moisture Code (FFMC).
    !
    ! Default values and units:
    !     ffmc = previous day FFMC (default = 85)
    !     T    = temperature, C
    !     RH   = relative humidity, %
    !     w    = wind speed, m/s
    !     Pr   = 24-hour precipitation, mm
    !
    ! Note: 
    !     Weather data are measured at local noon.
    !
    
    implicit none

    real, intent(in) :: ffmc, T, RH, wind, Pr

    ! local
    real :: f_0, w, m_0, h
    real :: r_f, m_r, dm, e_d, e_w, m
    real :: k_0, k_l, k_d, k_w
    real :: E_H_term, E_T_term, result


    ! use the same variable names as in the FFMC equation and 
    ! convert the units if necessary
    f_0 = ffmc
    w   = 3.6 * wind  ! convert from m/s to km/h


    ! fuel moisture content from previous day
    m_0 = 147.2*(101 - f_0) / (59.5 + f_0)
    
    ! current fuel moisture content
    m_r = m_0

    ! rainfall correction
    if (Pr > 0.5) then
        r_f = Pr - 0.5
        
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

    E_H_term = exp(0.1 * (RH - 100.0))
    E_T_term = 0.18*(21.1 - T) * (1 - exp(-0.115*RH))
    E_d = 0.942* (RH**0.679) + 11*E_H_term + E_T_term

    h = RH / 100.0

    if (m_r > E_d) then
        k_0 = 0.424 * (1 - (h**1.7)) + 0.0694 * sqrt(w) * (1 - (h**8))
        k_d = k_0 * 0.581 * exp(0.0365 * T)

        m = E_d + (m_r - E_d)*(10.0**(-k_d))
    else
        E_w = 0.618*(RH**0.753) + 10*E_H_term + E_T_term

        if (m_r < E_w) then
            k_l = 0.424 * (1 - (1 - h)**1.7) + 0.0694 * sqrt(w)* (1 - (1 - h)**8)
            k_w = k_l * 0.581 * exp(0.0365 * T)

            m = E_w - (E_w - m_r)*(10.0**(-k_w))
        else
            m = m_r
        end if
    end if    


    ! current FFMC
    result = 59.5*(250 - m)/(147.2 + m)
    
    ! clamp FFMC within [0, 101]
    fine_fuel_moisture_code = max(0.0, min(101.0, result))

end function fine_fuel_moisture_code



elemental real function duff_moisture_code(dmc, T, RH, Pr, month)

    !
    ! Calculates the Duff Moisture Code (DMC).
    !
    ! Default values and units:
    !     dmc = previous day DMC (default = 6)
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



elemental real function drought_code(dc, T, Pr, month)

    !
    ! Calculates the Drought Code (DC).
    ! 
    ! Default values and units:
    !     dc    = previous day DC (default = 15.0)
    !     T     = temperature, C
    !     Pr    = 24-hour precipitation, mm
    !     month = [1..12]
    !
    ! Note: 
    !     Weather data are measured at local noon.
    !

    implicit none

    real,    intent(in) :: dc, T, Pr
    integer, intent(in) :: month

   
    ! local
    real, dimension(12), parameter ::                             &
        DAY_LENGTH_FACTOR = (/ -1.6, -1.6, -1.6, 0.9,  3.8,  5.8, &
                                6.4,  5.0,  2.4, 0.4, -1.6, -1.6 /)

    real :: d_0, L_f
    real :: d, Q_0, Q_r, r_d, d_r, V, result


    ! use the same variable names as in the DC equation and 
    ! convert the units if necessary
    d_0 = dc
    L_f = DAY_LENGTH_FACTOR(month)


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

    ! don't allow T < -2.8C
    if (T < -2.8) then
        V = 0.0
    else
        V = 0.36*(T + 2.8) + L_f
    end if


    ! current DC
    result = d + 0.5*V

    drought_code = max(0.0, result)

end function drought_code



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
    
    m = 147.2 * (101.053 - ffmc)/(59.5 + ffmc)
    
    fun_w = exp(0.05039 * w)
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


subroutine cffwi_indexes(prev_day_ffmc, prev_day_dmc, prev_day_dc, &
                         T, RH, wind, Pr,                          &
                         month,                                    &
                         ffmc, dmc, dc, isi, bui, fwi, dsr)

    ! Calculates FFMC, DMC, DC, ISI, BUI, FWI and DSR indexes.
    !
    ! Default values and units:
    !     prev_day_ffmc = previous day FFMC (default = 85)
    !     prev_day_dmc  = previous day DMC  (default =  6)
    !     prev_day_dc   = previous day DC   (default = 15)
    !     T             = temperature, C
    !     RH            = relative humidity, %
    !     wind          = wind speed, m/s
    !     Pr            = 24-hour precipitation, mm
    !
    ! Note: 
    !     Weather data are measured at local noon.
    !

    implicit none

    real, intent(in)    :: prev_day_ffmc, prev_day_dmc, prev_day_dc
    real, intent(in)    :: T, RH, wind, Pr
    integer, intent(in) :: month
    real, intent(out)   :: ffmc, dmc, dc, isi, bui, fwi, dsr


    ! update fuel moisture codes 
    ffmc = fine_fuel_moisture_code(prev_day_ffmc, T, RH, wind, Pr)
    dmc  = duff_moisture_code(prev_day_dmc, T, RH, Pr, month)
    dc   = drought_code(prev_day_dc, T, Pr, month)

    ! update fire behavior indexes
    isi = initial_spread_index(ffmc, wind)
    bui = buildup_index(dmc, dc)

    fwi = fire_weather_index(isi, bui)
    dsr = daily_severity_rating(fwi)

end subroutine cffwi_indexes

end module cffwi
