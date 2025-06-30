!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !INTERFACE:

 REAL FUNCTION GetCO2(Year, DayOfYear)

! !DESCRIPTION:
!
!  Given the year and day-of-year, this function returns the RCP45 CO2 concentration in 
!  mole fraction (volume mixing ratio).  If Year is less than 1765, the value for 1765
!  is returned.  If Year is greater than 2150, the value for 2150 is returned.  In the
!  original dataset, the value for 2150 is used for all years through 2500.  We choose
!  to truncate the list at 2151 for this application.
!
!  DayOfYear is expected to have a value of 1.00 at 0:00 UTC Jan 1.
!
!  In-line documentation from the source dataset is reproduced below:
! 
! RCP45 Midyear Atmospheric CO2 Concentrations (ppmv)
!
! CONTENT:           CMIP5 recommendations for annual average, global mean concentrations.
! RUN:               RCP4.5, FINAL RELEASE, 26 Nov. 2009
! RCP4.5 CONTACT:    MiniCAM group, Allison Thomson (Allison.Thomson@pnl.gov)
! DATE:              26/11/2009 09:00:37  (updated description, 30 May 2010)
! MAGICC-VERSION:    6.3.09, 25 November 2009
! FILE PRODUCED BY:  RCP Concentration Calculation & Data Group, M. Meinshausen, S. Smith,  
!                     K. Riahi, D. van Vuuren
! DOCUMENTATION:     M. Meinshausen, S. Smith et al. "The RCP GHG concentrations and
!                     their extension from 1765 to 2500", in prep., Climatic Change.
! CMIP5 INFO:        http://cmip-pcmdi.llnl.gov/cmip5/
! RCP DATABASE:      http://www.iiasa.ac.at/web-apps/tnt/RcpDb
! FURTHER INFO:      For data sources, aknowledgements and further information, see 
!                     http://www.pik-potsdam.de/~mmalte/rcps
! NOTES:             RCP4.5 starts 2005; 20th century data and earlier is provided for
!                     convenience.  
!
! !REVISION HISTORY:
!  29 Oct 2010  Nielsen, adapted for GEOS-5.
! 
!EOP
! ---------------------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Year
  INTEGER, INTENT(IN) :: DayOfYear

  REAL :: f,i,n
  INTEGER :: previous,current,next
  
  INTEGER, PARAMETER :: firstYear = 1764
  INTEGER, PARAMETER :: finalYear = 2151
  INTEGER, PARAMETER :: tableLength = finalYear-firstYear+1

  REAL, SAVE :: CO2ppmv(tableLength) = (/                                278.052, &
   278.052,  278.106,  278.220,  278.343,  278.471,  278.600,  278.733,  278.869, &
   279.009,  279.153,  279.302,  279.457,  279.618,  279.782,  279.943,  280.097, &
   280.243,  280.382,  280.518,  280.657,  280.803,  280.957,  281.118,  281.282, &
   281.443,  281.598,  281.747,  281.891,  282.031,  282.167,  282.299,  282.427, &
   282.551,  282.671,  282.787,  282.899,  283.007,  283.111,  283.211,  283.307, &
   283.400,  283.490,  283.578,  283.661,  283.735,  283.797,  283.847,  283.889, &
   283.926,  283.963,  284.001,  284.043,  284.086,  284.129,  284.167,  284.198, &
   284.223,  284.244,  284.263,  284.281,  284.300,  284.320,  284.340,  284.360, &
   284.380,  284.400,  284.385,  284.280,  284.125,  283.975,  283.825,  283.675, &
   283.525,  283.425,  283.400,  283.400,  283.425,  283.500,  283.600,  283.725, &
   283.900,  284.075,  284.225,  284.400,  284.575,  284.725,  284.875,  285.000, &
   285.125,  285.275,  285.425,  285.575,  285.725,  285.900,  286.075,  286.225, &
   286.375,  286.500,  286.625,  286.775,  286.900,  287.000,  287.100,  287.225, &
   287.375,  287.525,  287.700,  287.900,  288.125,  288.400,  288.700,  289.025, &
   289.400,  289.800,  290.225,  290.700,  291.200,  291.675,  292.125,  292.575, &
   292.975,  293.300,  293.575,  293.800,  294.000,  294.175,  294.325,  294.475, &
   294.600,  294.700,  294.800,  294.900,  295.025,  295.225,  295.500,  295.800, &
   296.125,  296.475,  296.825,  297.200,  297.625,  298.075,  298.500,  298.900, &
   299.300,  299.700,  300.075,  300.425,  300.775,  301.100,  301.400,  301.725, &
   302.075,  302.400,  302.700,  303.025,  303.400,  303.775,  304.125,  304.525, &
   304.975,  305.400,  305.825,  306.300,  306.775,  307.225,  307.700,  308.175, &
   308.600,  309.000,  309.400,  309.750,  310.000,  310.175,  310.300,  310.375, &
   310.375,  310.300,  310.200,  310.125,  310.100,  310.125,  310.200,  310.325, &
   310.500,  310.750,  311.100,  311.500,  311.925,  312.425,  313.000,  313.600, &
   314.225,  314.848,  315.500,  316.272,  317.075,  317.795,  318.397,  318.925, &
   319.647,  320.647,  321.605,  322.635,  323.902,  324.985,  325.855,  327.140, &
   328.677,  329.742,  330.585,  331.747,  333.272,  334.848,  336.525,  338.360, &
   339.728,  340.793,  342.198,  343.783,  345.283,  346.797,  348.645,  350.737, &
   352.487,  353.855,  355.017,  355.885,  356.777,  358.128,  359.837,  361.462, &
   363.155,  365.323,  367.348,  368.865,  370.467,  372.522,  374.760,  376.812, &
   378.812,  380.828,  382.777,  384.800,  386.952,  389.128,  391.274,  393.421, &
   395.583,  397.764,  399.966,  402.184,  404.411,  406.643,  408.882,  411.129, &
   413.378,  415.639,  417.936,  420.274,  422.656,  425.080,  427.538,  430.021, &
   432.523,  435.046,  437.589,  440.131,  442.664,  445.207,  447.770,  450.355, &
   452.963,  455.586,  458.215,  460.845,  463.475,  466.093,  468.678,  471.234, &
   473.780,  476.328,  478.881,  481.438,  483.993,  486.535,  489.060,  491.536, &
   493.932,  496.244,  498.474,  500.645,  502.768,  504.847,  506.884,  508.871, &
   510.799,  512.647,  514.401,  516.065,  517.629,  519.096,  520.488,  521.818, &
   523.089,  524.302,  525.451,  526.509,  527.457,  528.296,  529.027,  529.643, &
   530.144,  530.553,  530.883,  531.138,  531.319,  531.490,  531.702,  531.942, &
   532.205,  532.487,  532.776,  533.070,  533.388,  533.741,  534.131,  534.558, &
   535.011,  535.480,  535.955,  536.435,  536.920,  537.399,  537.871,  538.358, &
   538.872,  539.388,  539.884,  540.352,  540.782,  541.168,  541.510,  541.808, &
   542.053,  542.246,  542.408,  542.559,  542.712,  542.866,  543.013,  543.139, &
   543.239,  543.311,  543.355,  543.360,  543.327,  543.288,  543.266,  543.264, &
   543.282,  543.310,  543.337,  543.354,  543.361,  543.356,  543.330,  543.283, &
   543.238,  543.208,  543.197,  543.205,  543.223,  543.239,  543.246,  543.242, &
   543.227,  543.190,  543.131,  543.072,  543.025,  542.993,  542.979,  542.973, &
   542.963,  542.955,  542.955 /)

! Establish location in table for current, previous and next year.
! ----------------------------------------------------------------
  current  = Year-firstYear+1
  current  = MAX(current,1)
  current  = MIN(current,tableLength)

  previous = current-1
  previous = MAX(previous,1)
  previous = MIN(previous,tableLength)
  IF(Year > finalYear) previous = tableLength

  next     = current+1
  next     = MAX(next,1)
  next     = MIN(next,tableLength)
  IF(Year < firstYear) next = 1

! Divide the year into halves.
! ----------------------------
  IF(dayOfYear <= 182) THEN
   i = CO2ppmv(previous)
   f = CO2ppmv(current)
   n = dayOfYear+183
  ELSE
   i = CO2ppmv(current)
   f = CO2ppmv(next)
   n = dayOfYear-183
  END IF

! Linear interpolation to the given day-of-year.
! ----------------------------------------------
  GetCO2 = i + (f-i)*n/365.00

! Convert to mole fraction (volume mixing ratio).
! -----------------------------------------------
  GetCO2 = GetCO2*1.00E-06

 END FUNCTION GetCO2
