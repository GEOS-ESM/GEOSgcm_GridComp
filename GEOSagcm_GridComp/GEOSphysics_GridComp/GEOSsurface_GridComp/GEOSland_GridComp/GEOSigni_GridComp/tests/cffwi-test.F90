!
! Test the implementation of the Canadian Forest Fire 
! Weather Index (CFFWI).
!
!
! Anton Darmenov, NASA, 2022
!


!
! Buld executable: 
!    $ gfortran -g -Wall -o cffwi-test cffwi.F90 cffwi-test.F90
!

program cffwi_test

    use cffwi, only: fine_fuel_moisture_code, duff_moisture_code,       &
                     drought_code, initial_spread_index, buildup_index, &
                     fire_weather_index, daily_severity_rating,         &
                     cffwi_indexes,                                     &
                     FFMC_INIT, DMC_INIT, DC_INIT,                      &
                     CFFWI_REFERENCE_LATITUDE


    implicit none

    logical, parameter :: INDIVIDUAL_INDEXES = .False.

    real, parameter :: FFMC_DEFAULT = FFMC_INIT
    real, parameter :: DMC_DEFAULT  = DMC_INIT
    real, parameter :: DC_DEFAULT   = DC_INIT


    ! local
    real    :: ffmc_pd, dmc_pd, dc_pd
    real    :: ffmc, dmc, dc, fwi, bui, dsr, isi
    real    :: T, RH, wind, rain
    real    :: time_step
    real    :: latitude
    integer :: month
    

    time_step = 24.0 ! values larger than 1 will triger daily FFMC
    latitude  = CFFWI_REFERENCE_LATITUDE
    
    ffmc_pd = FFMC_DEFAULT
    dmc_pd  = DMC_DEFAULT
    dc_pd   = DC_DEFAULT


    !! conditions - heavy rain
    !ffmc_pd = 89.7
    !dmc_pd  = 57.4
    !dc_pd   =108.8
    !
    !T    = 16.0
    !RH   = 50.0
    !wind = 22.0
    !rain = 12.2
    !month= 5

    !! conditions- moderate rain
    !ffmc_pd = 87.7
    !dmc_pd  = 8.5
    !dc_pd   = 19.0
    !
    !T    = 20.0
    !RH   = 21.0
    !wind = 25.0
    !rain =  2.4
    !month= 4

    !! conditions - no rain
    !ffmc_pd = 86.2
    !dmc_pd  = 10.4
    !dc_pd   = 23.6
    !
    !T    = 8.5
    !RH   = 40.0
    !wind = 17.0
    !rain =  0.0
    !month= 4

    !! conditions - no rain and hot weather
    ffmc_pd = 77.6
    dmc_pd  = 18.7
    dc_pd   =117.7
    
    T    = 30.0
    RH   = 38.0
    wind = 22.0
    rain =  0.0
    month= 5
    
    ! convert from km/h to m/s
    wind = wind / 3.6


    if (INDIVIDUAL_INDEXES) then
        write (*, '(A)') "Testing the individual cffwi indexes:"

        ffmc = fine_fuel_moisture_code(ffmc_pd, T, RH, wind, rain, time_step)

        dmc = duff_moisture_code(dmc_pd, T, RH, rain, month)


        dc = drought_code(dc_pd, T, rain, latitude, month) 

        isi = initial_spread_index(ffmc, wind)
        bui = buildup_index(dmc, dc)

        fwi = fire_weather_index(isi, bui)
        dsr = daily_severity_rating(fwi)
    else
        write (*, '(A)') "Testing the cffwi_indexes():"

        call cffwi_indexes(ffmc_pd, dmc_pd, dc_pd,   &
                           T, RH, wind, rain,        &
                           latitude, month,          &
                           time_step,                &
                           ffmc, dmc, dc, isi, bui, fwi, dsr)
    end if

    write (*, '(A, F6.1)') "  FFMC = ", ffmc
    write (*, '(A, F6.1)') "  DMC  = ", dmc
    write (*, '(A, F6.1)') "  DC   = ", dc
    write (*, '(A)')       "  -------------"
    write (*, '(A, F6.1)') "  ISI  = ", isi
    write (*, '(A, F6.1)') "  BUI  = ", bui
    write (*, '(A)')       "  -------------"
    write (*, '(A, F6.1)') "  FWI  = ", fwi
    write (*, '(A, F6.2)') "  DSR  = ", dsr

end program
