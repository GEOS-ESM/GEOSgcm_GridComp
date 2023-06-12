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
!  CCMI Midyear Atmospheric CO2 Concentrations (ppmv)
!
! CONTENT:                 CCMI RefD1 recommendations for annual average, global mean concentrations.
! Sources (through 2014):  UoM-CMIP-1-2-0: Historical GHG mole fractions from NOAA & AGAGE networks and many others contributors (see references in Meinshausen et al.) until 2014, then 
! Source (post-2014):      UoM-MESSAGE-GLOBIOM-ssp245-1-2-1: GHG mole fractions from (historical) NOAA & AGAGE networks and many others contributors (see references in Meinshausen et al.). Future GHG concentrations modelled with MAGICC7.0.      
! References:              Malte Meinshausen, Elisabeth Vogel, Alexander Nauels, Katja Lorbacher, Nicolai Meinshausen, David Etheridge, Paul Fraser, Stephen A. Montzka, Peter Rayner, Cathy Trudinger, Paul Krummel, Urs Beyerle, Josep G. Canadell, John S. Daniel, Ian Enting, Rachel M. Law, Simon O\'Doherty, Ron G. Prinn, Stefan Reimann, Mauro Rubino, Martin K. Vollmer, Ray Weiss (2016), Historical greenhouse gas concentrations, doi:10.5194/gmd-2016-169; Malte Meinshausen, Zebedee Nicholls, et al. Future surface greenhouse gas concentrations, GMD Special Issue (please update for final reference) http://www.geosci-model-dev.net/special_issue590.html
! Input4MIPS Info:         https://esgf-node.llnl.gov/search/input4mips/
! Original Files:          mole-fraction-of-carbon-dioxide-in-air_input4MIPs_GHGConcentrations_CMIP_UoM-CMIP-1-2-0_gr1-GMNHSH_0000-2014.nc and mole-fraction-of-carbon-dioxide-in-air_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-MESSAGE-GLOBIOM-ssp245-1-2-1_gr1-GMNHSH_2015-2500.nc downloaded on Dec. 14, 2020 from input4MIPS
! 
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

  REAL, SAVE :: CO2ppmv(tableLength) = (/                                277.859, &
   277.913,  277.958,  278.017,  278.078,  278.135,  278.194,  278.270,  278.334,  &
   278.412,  278.491,  278.586,  278.683,  278.774,  278.866,  278.959,  279.054,  &
   279.145,  279.243,  279.340,  279.434,  279.532,  279.629,  279.729,  279.831,  &
   279.933,  280.037,  280.141,  280.241,  280.332,  280.415,  280.492,  280.555,  &
   280.623,  280.682,  280.753,  280.834,  280.919,  280.994,  281.073,  281.157,  &
   281.234,  281.322,  281.409,  281.496,  281.566,  281.648,  281.718,  281.790,  &
   281.858,  281.931,  281.992,  282.059,  282.127,  282.200,  282.272,  282.346,  &
   282.415,  282.477,  282.540,  282.618,  282.676,  282.693,  282.693,  282.708,  &
   282.761,  282.825,  282.891,  282.922,  282.966,  283.028,  283.086,  283.149,  &
   283.198,  283.253,  283.311,  283.382,  283.467,  283.528,  283.566,  283.609,  &
   283.725,  283.847,  283.936,  284.068,  284.207,  284.317,  284.451,  284.598,  &
   284.731,  284.846,  284.941,  285.049,  285.204,  285.369,  285.545,  285.739,  &
   285.933,  286.100,  286.271,  286.442,  286.614,  286.781,  286.955,  287.105,  &
   287.225,  287.355,  287.494,  287.664,  287.860,  288.061,  288.291,  288.520,  &
   288.752,  288.993,  289.221,  289.470,  289.737,  290.019,  290.263,  290.512,  &
   290.797,  291.100,  291.414,  291.763,  292.113,  292.458,  292.816,  293.167,  &
   293.477,  293.791,  294.079,  294.365,  294.646,  294.954,  295.300,  295.675,  &
   296.007,  296.325,  296.654,  296.954,  297.289,  297.662,  298.098,  298.518,  &
   298.936,  299.377,  299.829,  300.353,  300.910,  301.419,  301.937,  302.485,  &
   303.011,  303.449,  303.814,  304.246,  304.600,  304.945,  305.271,  305.630,  &
   305.813,  305.954,  306.177,  306.329,  306.495,  306.620,  306.822,  307.093,  &
   307.402,  307.785,  308.227,  309.012,  309.764,  310.294,  310.851,  311.357,  &
   311.811,  312.172,  312.390,  312.413,  312.385,  312.390,  312.486,  312.521,  &
   312.632,  312.821,  313.014,  313.342,  313.730,  314.095,  314.415,  314.698,  &
   314.992,  315.345,  315.807,  316.625,  317.299,  318.044,  318.650,  319.333,  &
   319.816,  320.880,  321.480,  322.389,  323.251,  324.783,  325.400,  327.349,  &
   329.909,  330.756,  330.827,  331.545,  333.353,  335.010,  336.605,  338.705,  &
   340.059,  340.644,  342.266,  344.008,  345.459,  346.903,  348.775,  351.276,  &
   352.894,  354.073,  355.353,  356.229,  356.925,  358.254,  360.239,  362.005,  &
   363.252,  365.933,  367.845,  369.125,  370.673,  372.835,  375.411,  376.987,  &
   378.907,  381.010,  382.603,  384.739,  386.280,  388.717,  390.944,  393.016,  &
   395.725,  397.547,  399.949,  403.117,  405.762,  408.632,  411.506,  414.390,  &
   417.287,  420.198,  423.125,  426.069,  429.030,  432.011,  435.012,  438.033,  &
   441.077,  444.143,  447.232,  450.333,  453.435,  456.541,  459.652,  462.770,  &
   465.897,  469.033,  472.179,  475.336,  478.504,  481.674,  484.838,  487.996,  &
   491.150,  494.300,  497.447,  500.592,  503.734,  506.875,  510.014,  513.137,  &
   516.230,  519.292,  522.326,  525.332,  528.310,  531.261,  534.186,  537.085,  &
   539.958,  542.795,  545.587,  548.335,  551.039,  553.700,  556.319,  558.895,  &
   561.431,  563.925,  566.380,  568.772,  571.083,  573.314,  575.466,  577.540,  &
   579.537,  581.458,  583.304,  585.074,  586.770,  588.378,  589.886,  591.296,  &
   592.609,  593.825,  594.945,  595.970,  596.900,  597.735,  598.477,  599.150,  &
   599.776,  600.355,  600.886,  601.366,  601.797,  602.177,  602.505,  602.782,  &
   603.005,  603.218,  603.459,  603.724,  604.012,  604.322,  604.650,  604.996,  &
   605.357,  605.734,  606.124,  606.527,  606.942,  607.367,  607.803,  608.248,  &
   608.701,  609.163,  609.633,  610.110,  610.593,  611.082,  611.578,  612.079,  &
   612.585,  613.095,  613.611,  614.130,  614.654,  615.181,  615.712,  616.246,  &
   616.784,  617.324,  617.867,  618.413,  618.961,  619.512,  620.064,  620.619,  &
   621.176,  621.734,  622.295,  622.857,  623.420,  623.985,  624.551,  625.119,  &
   625.687,  626.257,  626.828 /)

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
