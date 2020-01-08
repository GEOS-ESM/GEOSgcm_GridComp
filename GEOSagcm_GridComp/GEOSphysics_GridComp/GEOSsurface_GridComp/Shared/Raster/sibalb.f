      SUBROUTINE SIBALB (
     O			 AVISDR, ANIRDR, AVISDF, ANIRDF,
     I			 VLAI, VGRN, ZTH, SNW, ITYP, IRUN
     &       		)

	IMPLICIT NONE

	REAL ALVDRS, ALIDRS
	REAL ALVDRD, ALIDRD
	REAL ALVDRI, ALIDRI

	PARAMETER (
     `         ALVDRS = 0.100
     `,		   ALIDRS = 0.200

     `,        ALVDRD = 0.300
     `,		   ALIDRD = 0.350

     `,        ALVDRI = 0.700
     `,		   ALIDRI = 0.700
     `		  )

* ALVDRS:  Albedo of soil for visible   direct  solar radiation.
* ALIDRS:  Albedo of soil for infra-red direct  solar radiation.
* ALVDFS:  Albedo of soil for visible   diffuse solar radiation.
* ALIDFS:  Albedo of soil for infra-red diffuse solar radiation.


      INTEGER NLAI

      PARAMETER (NLAI = 14 )


      REAL EPSLN, BLAI, DLAI

      PARAMETER (EPSLN = 1.E-6)
      PARAMETER (BLAI = 0.5)
      PARAMETER (DLAI = 0.5)

      REAL ZERO, ONE

      PARAMETER (ZERO=0., ONE=1.0)

      REAL ALATRM

      PARAMETER (ALATRM = BLAI + (NLAI - 1) * DLAI - EPSLN)

      INTEGER NTYPS

      PARAMETER (NTYPS=9)

      INTEGER IRUN

      REAL AVISDR (IRUN), ANIRDR (IRUN), AVISDF (IRUN), ANIRDF (IRUN),
     `	     VLAI   (IRUN),   VGRN (IRUN),   ZTH  (IRUN),    SNW (IRUN)

* OUTPUTS:

* AVISDR:   visible, direct albedo.
* ANIRDR:   near infra-red, direct albedo.
* AVISDF:   visible, diffuse albedo.
* ANIRDF:   near infra-red, diffuse albedo.

* INPUTS:

* VLAI:     the leaf area index.
* VGRN:     the greenness index.
* ZTH:      The cosine of the solar zenith angle.
* SNW:      Snow cover in meters water equivalent.
*



        INTEGER ITYP (IRUN)

C ITYP: Vegetation type as follows:
C                  1:  BROADLEAF EVERGREEN TREES
C                  2:  BROADLEAF DECIDUOUS TREES
C                  3:  NEEDLELEAF TREES
C                  4:  GROUND COVER
C                  5:  BROADLEAF SHRUBS
C                  6:  DWARF TREES (TUNDRA)
C                  7:  BARE SOIL
C                  8:  DESSERT
C                  9:  ICE
C  IRUN: Chip index
C


* [ Definition of Variables: ]
*
	INTEGER I, LAI


	REAL FAC,               GAMMA,          BETA,          ALPHA,
     `	     DX,                DY,             ALA,           GRN (2),
     `	     SNWALB (4, NTYPS), SNWMID (NTYPS)

*

* [ Definition of Functions: ]
*
	REAL COEFF

C --------------------------------------------------



C   Constants used in albedo calculations:

      REAL ALVDR (NLAI, 2, NTYPS)
      REAL BTVDR (NLAI, 2, NTYPS)
      REAL GMVDR (NLAI, 2, NTYPS)
      REAL ALIDR (NLAI, 2, NTYPS)
      REAL BTIDR (NLAI, 2, NTYPS)
      REAL GMIDR (NLAI, 2, NTYPS)
      
C  (Data statements for ALVDR described in full; data statements for
C   other constants follow same framework.)

C    BROADLEAF EVERGREEN (ITYP=4); GREEN=0.33; LAI: .5-7
	DATA (ALVDR (I, 1, 1), I = 1, 14)
     `	  /0.0808, 0.0796, 0.0792, 0.0790, 10*0.0789/

C    BROADLEAF EVERGREEN (ITYP=4); GREEN=0.67; LAI: .5-7
	DATA (ALVDR (I, 2, 1), I = 1, 14)
     `	  /0.0788, 0.0775, 0.0771, 0.0769, 10*0.0768/

C    BROADLEAF DECIDUOUS (ITYP=1); GREEN=0.33; LAI: .5-7
	DATA (ALVDR (I, 1, 2), I = 1, 14)
     `	  /0.0803, 0.0790, 0.0785, 0.0784, 3*0.0783, 7*0.0782/

C    BROADLEAF DECIDUOUS (ITYP=1); GREEN=0.67; LAI: .5-7
	DATA (ALVDR (I, 2, 2), I = 1, 14)
     `	  /0.0782, 0.0770, 0.0765, 0.0763, 10*0.0762/

C    NEEDLELEAF (ITYP=3); GREEN=0.33; LAI=.5-7
	DATA (ALVDR (I, 1, 3), I = 1, 14)
     `	  /0.0758, 0.0746, 0.0742, 0.0740, 10*0.0739/

C    NEEDLELEAF (ITYP=3); GREEN=0.67; LAI=.5-7
	DATA (ALVDR (I, 2, 3), I = 1, 14)
     `	  /0.0683, 0.0672, 0.0667, 2*0.0665, 9*0.0664/

C    GROUNDCOVER (ITYP=2); GREEN=0.33; LAI=.5-7
	DATA (ALVDR (I, 1, 4), I = 1, 14)
     `	  /0.2436, 0.2470, 0.2486, 0.2494, 0.2498, 0.2500, 2*0.2501,
     `		6*0.2502
     `	  /
C    GROUNDCOVER (ITYP=2); GREEN=0.67; LAI=.5-7
	DATA (ALVDR (I, 2, 4), I = 1, 14) /14*0.1637/

C    BROADLEAF SHRUBS (ITYP=5); GREEN=0.33,LAI=.5-7
        DATA (ALVDR (I, 1, 5), I = 1, 14)
     &    /0.0807, 0.0798, 0.0794, 0.0792, 0.0792, 9*0.0791/

C    BROADLEAF SHRUBS (ITYP=5); GREEN=0.67,LAI=.5-7
        DATA (ALVDR (I, 2, 5), I = 1, 14)
     &    /0.0787, 0.0777, 0.0772, 0.0771, 10*0.0770/

C    DWARF TREES, OR TUNDRA (ITYP=6); GREEN=0.33,LAI=.5-7
        DATA (ALVDR (I, 1, 6), I = 1, 14)
     &    /0.0802, 0.0791, 0.0787, 0.0786, 10*0.0785/

C    DWARF TREES, OR TUNDRA (ITYP=6); GREEN=0.67,LAI=.5-7
        DATA (ALVDR (I, 2, 6), I = 1, 14)
     &    /0.0781, 0.0771, 0.0767, 0.0765, 0.0765, 9*0.0764/


C    BARE SOIL
	DATA (ALVDR (I, 1, 7), I = 1, 14) /14*ALVDRS/
	DATA (ALVDR (I, 2, 7), I = 1, 14) /14*ALVDRS/

C    DESERT
	DATA (ALVDR (I, 1, 8), I = 1, 14) /14*ALVDRD/
	DATA (ALVDR (I, 2, 8), I = 1, 14) /14*ALVDRD/

C    ICE
	DATA (ALVDR (I, 1, 9), I = 1, 14) /14*ALVDRI/
	DATA (ALVDR (I, 2, 9), I = 1, 14) /14*ALVDRI/
C****
C**** -------------------------------------------------
	DATA (BTVDR (I, 1, 1), I = 1, 14)
     `	  /0.0153, 0.0372, 0.0506, 0.0587, 0.0630, 0.0652, 0.0663,
     `		0.0668, 0.0671, 0.0672, 4*0.0673
     `	  /
	DATA (BTVDR (I, 2, 1), I = 1, 14)
     *	  /0.0135, 0.0354, 0.0487, 0.0568, 0.0611, 0.0633, 0.0644,
     `		0.0650, 0.0652, 0.0654, 0.0654, 3*0.0655
     `	  /
	DATA (BTVDR (I, 1, 2), I = 1, 14)
     *	  /0.0148, 0.0357, 0.0462, 0.0524, 0.0554, 0.0569, 0.0576,
     `		0.0579, 0.0580, 0.0581, 0.0581, 3*0.0582
     `	  /
	DATA (BTVDR (I, 2, 2), I = 1, 14)
     *	  /0.0131, 0.0342, 0.0446, 0.0508, 0.0539, 0.0554, 0.0560,
     `		0.0564, 0.0565, 5*0.0566
     `	  /
	DATA (BTVDR (I, 1, 3), I = 1, 14)
     *	  /0.0108, 0.0334, 0.0478, 0.0571, 0.0624, 0.0652, 0.0666,
     `		0.0673, 0.0677, 0.0679, 4*0.0680
     `	  /
	DATA (BTVDR (I, 2, 3), I = 1, 14)
     *	  /0.0034, 0.0272, 0.0408, 0.0501, 0.0554, 0.0582, 0.0597,
     *		0.0604, 0.0608, 0.0610, 4*0.0611
     `	  /
	DATA (BTVDR (I, 1, 4), I = 1, 14)
     *	  /0.2050, 0.2524, 0.2799, 0.2947, 0.3022, 0.3059, 0.3076,
     *		0.3085, 0.3088, 0.3090, 4*0.3091
     `	  /
	DATA (BTVDR (I, 2, 4), I = 1, 14)
     *	  /0.1084, 0.1404, 0.1617, 0.1754, 0.1837, 0.1887, 0.1915,
     *		0.1931, 0.1940, 0.1946, 0.1948, 0.1950, 2*0.1951
     `	  /
        DATA (BTVDR (I, 1, 5), I = 1, 14)
     &    /0.0203, 0.0406, 0.0548, 0.0632, 0.0679, 0.0703, 0.0716,
     &     0.0722, 0.0726, 0.0727, 0.0728, 0.0728, 0.0728, 0.0729
     `	  /

        DATA (BTVDR (I, 2, 5), I = 1, 14)
     &    /0.0184, 0.0385, 0.0526, 0.0611,  0.0658, 0.0683, 0.0696,
     &     0.0702, 0.0705, 0.0707, 4*0.0708
     `	  /

        DATA (BTVDR (I, 1, 6), I = 1, 14)
     &    /0.0199, 0.0388, 0.0494,  0.0554, 0.0584, 0.0599, 0.0606,
     &     0.0609, 0.0611, 5*0.0612
     `	  /

        DATA (BTVDR (I, 2, 6), I = 1, 14)
     &    /0.0181, 0.0371, 0.0476, 0.0537,  0.0568, 0.0583, 0.0590,
     &     0.0593, 0.0595, 0.0595, 4*0.0596
     `	  /

	DATA (BTVDR (I, 1, 7), I = 1, 14) /14*0./
	DATA (BTVDR (I, 2, 7), I = 1, 14) /14*0./

	DATA (BTVDR (I, 1, 8), I = 1, 14) /14*0./
	DATA (BTVDR (I, 2, 8), I = 1, 14) /14*0./

	DATA (BTVDR (I, 1, 9), I = 1, 14) /14*0./
	DATA (BTVDR (I, 2, 9), I = 1, 14) /14*0./

C****
C**** -----------------------------------------------------------
	DATA (GMVDR (I, 1, 1), I = 1, 14)
     `	  /0.0814, 0.1361, 0.2078, 0.2650, 0.2986, 0.3169,  0.3265,
     * 	   0.3313, 0.3337, 0.3348, 0.3354, 0.3357, 2*0.3358
     `	  /
	DATA (GMVDR (I, 2, 1), I = 1, 14)
     *	  /0.0760, 0.1336, 0.2034, 0.2622, 0.2969, 0.3159,  0.3259,
     *	   0.3309, 0.3333, 0.3346, 0.3352, 0.3354, 2*0.3356
     `	  /
	DATA (GMVDR (I, 1, 2), I = 1, 14)
     *	  /0.0834, 0.1252, 0.1558, 0.1927, 0.2131,   0.2237, 0.2290,
     *	   0.2315, 0.2327, 0.2332, 0.2335, 2*0.2336, 0.2337
     `	  /
	DATA (GMVDR (I, 2, 2), I = 1, 14)
     *	  /0.0789, 0.1235, 0.1531, 0.1912, 0.2122, 0.2232,  0.2286,
     *	   0.2312, 0.2324, 0.2330, 0.2333, 0.2334, 2*0.2335
     `	  /
	DATA (GMVDR (I, 1, 3), I = 1, 14)
     *	  /0.0647, 0.1342, 0.2215, 0.2968, 0.3432, 0.3696, 0.3838,
     *	   0.3912, 0.3950, 0.3968, 0.3978, 0.3982, 0.3984, 0.3985
     `	  /
	DATA (GMVDR (I, 2, 3), I = 1, 14)
     *	  /0.0258, 0.1227, 0.1999, 0.2825, 0.3339, 0.3634, 0.3794,
     *	   0.3877, 0.3919, 0.3940, 0.3950, 0.3956, 0.3958, 0.3959
     `	  /
	DATA (GMVDR (I, 1, 4), I = 1, 14)
     *	  /0.3371, 0.5762, 0.7159, 0.7927, 0.8324, 0.8526,  0.8624,
     *	   0.8671, 0.8693, 0.8704, 0.8709, 0.8710, 2*0.8712
     `	  /
	DATA (GMVDR (I, 2, 4), I = 1, 14)
     *	  /0.2634, 0.4375, 0.5532, 0.6291, 0.6763, 0.7048, 0.7213,
     *	   0.7310, 0.7363, 0.7395, 0.7411, 0.7420, 0.7426, 0.7428
     `	  /
        DATA (GMVDR (I, 1, 5), I = 1, 14)
     &    /0.0971, 0.1544, 0.2511, 0.3157, 0.3548, 0.3768, 0.3886,
     &     0.3948, 0.3978, 0.3994, 0.4001, 0.4006, 0.4007, 0.4008
     `	  /

        DATA (GMVDR (I, 2, 5), I = 1, 14)
     &    /0.0924, 0.1470, 0.2458, 0.3123, 0.3527, 0.3756, 0.3877,
     &     0.3942, 0.3974, 0.3990, 0.3998, 0.4002, 0.4004, 0.4005
     `	  /

        DATA (GMVDR (I, 1, 6), I = 1, 14)
     &    /0.0970, 0.1355, 0.1841, 0.2230, 0.2447,  0.2561, 0.2617,
     &     0.2645, 0.2658, 0.2664, 0.2667, 3*0.2669
     `	  /

        DATA (GMVDR (I, 2, 6), I = 1, 14)
     &    /0.0934, 0.1337, 0.1812, 0.2213, 0.2437, 0.2554, 0.2613,
     &     0.2642, 0.2656, 0.2662, 0.2665, 0.2667, 0.2667, 0.2668
     `	  /

	DATA (GMVDR (I, 1, 7), I = 1, 14) /14*1./
	DATA (GMVDR (I, 2, 7), I = 1, 14) /14*1./

	DATA (GMVDR (I, 1, 8), I = 1, 14) /14*1./
	DATA (GMVDR (I, 2, 8), I = 1, 14) /14*1./

	DATA (GMVDR (I, 1, 9), I = 1, 14) /14*1./
	DATA (GMVDR (I, 2, 9), I = 1, 14) /14*1./

C****
C****  -----------------------------------------------------------

	DATA (ALIDR (I, 1, 1), I = 1, 14)
     *	  /0.2867,  0.2840, 0.2828, 0.2822, 0.2819, 0.2818, 2*0.2817,
     *	   6*0.2816
     `	  /
	DATA (ALIDR (I, 2, 1), I = 1, 14)
     *	  /0.3564, 0.3573, 0.3577, 0.3580, 2*0.3581, 8*0.3582/
	DATA (ALIDR (I, 1, 2), I = 1, 14)
     *	  /0.2848, 0.2819, 0.2804, 0.2798, 0.2795, 2*0.2793, 7*0.2792/
	DATA (ALIDR (I, 2, 2), I = 1, 14)
     *	  /0.3544, 0.3550, 0.3553, 2*0.3555, 9*0.3556/
	DATA (ALIDR (I, 1, 3), I = 1, 14)
     *	  /0.2350, 0.2311, 0.2293, 0.2285, 0.2281, 0.2280, 8*0.2279/
	DATA (ALIDR (I, 2, 3), I = 1, 14)
     *	  /0.2474, 0.2436, 0.2418, 0.2410, 0.2406, 0.2405, 3*0.2404,
     *	   5*0.2403
     `    /
	DATA (ALIDR (I, 1, 4), I = 1, 14)
     *	  /0.5816, 0.6157, 0.6391, 0.6556, 0.6673, 0.6758, 0.6820,
     *	   0.6866, 0.6899, 0.6924, 0.6943, 0.6956, 0.6966, 0.6974
     `	  /
	DATA (ALIDR (I, 2, 4), I = 1, 14)
     *	  /0.5489, 0.5770, 0.5955, 0.6079, 0.6163, 0.6221, 0.6261,
     *	   0.6288, 0.6308, 0.6321, 0.6330, 0.6337, 0.6341, 0.6344
     `	  /
        DATA (ALIDR (I, 1, 5), I = 1, 14)
     &    /0.2845, 0.2837, 0.2832, 0.2831, 0.2830, 9*0.2829/
        DATA (ALIDR (I, 2, 5), I = 1, 14)
     &    /0.3532, 0.3562, 0.3578,  0.3586, 0.3590, 0.3592, 0.3594,
     &     0.3594, 0.3594, 5*0.3595
     `	  /
        DATA (ALIDR (I, 1, 6), I = 1, 14)
     &    /0.2825, 0.2812, 0.2806, 0.2803, 0.2802, 9*0.2801/
        DATA (ALIDR (I, 2, 6), I = 1, 14)
     &    /0.3512, 0.3538,  0.3552, 0.3559, 0.3562, 0.3564, 0.3565,
     &     0.3565, 6*0.3566
     `	  /

	DATA (ALIDR (I, 1, 7), I = 1, 14) /14*ALIDRS/
	DATA (ALIDR (I, 2, 7), I = 1, 14) /14*ALIDRS/

	DATA (ALIDR (I, 1, 8), I = 1, 14) /14*ALIDRD/
	DATA (ALIDR (I, 2, 8), I = 1, 14) /14*ALIDRD/

	DATA (ALIDR (I, 1, 9), I = 1, 14) /14*ALIDRI/
	DATA (ALIDR (I, 2, 9), I = 1, 14) /14*ALIDRI/

C****
C**** -----------------------------------------------------------
	DATA (BTIDR (I, 1, 1), I = 1, 14)
     *	  /0.1291, 0.1707, 0.1969, 0.2125, 0.2216,   0.2267, 0.2295,
     *	   0.2311, 0.2319, 0.2323, 0.2326, 2*0.2327, 0.2328
     `	  /
	DATA (BTIDR (I, 2, 1), I = 1, 14)
     *	  /0.1939, 0.2357, 0.2598, 0.2735, 0.2810,  0.2851, 0.2874,
     *	   0.2885, 0.2892, 0.2895, 0.2897, 3*0.2898
     `	  /
	DATA (BTIDR (I, 1, 2), I = 1, 14)
     *	  /0.1217, 0.1522, 0.1713, 0.1820,   0.1879,  0.1910, 0.1926,
     *	   0.1935, 0.1939, 0.1942, 2*0.1943, 2*0.1944
     `	  /
	DATA (BTIDR (I, 2, 2), I = 1, 14)
     *	  /0.1781, 0.2067, 0.2221, 0.2301,   0.2342,  0.2363, 0.2374,
     *	   0.2379, 0.2382, 0.2383, 2*0.2384, 2*0.2385
     `	  /
	DATA (BTIDR (I, 1, 3), I = 1, 14)
     *	  /0.0846, 0.1299, 0.1614, 0.1814, 0.1935,   0.2004, 0.2043,
     *	   0.2064, 0.2076, 0.2082, 0.2085, 2*0.2087, 0.2088
     `	  /
	DATA (BTIDR (I, 2, 3), I = 1, 14)
     *	  /0.0950, 0.1410, 0.1722, 0.1921, 0.2042, 0.2111,  0.2151,
     *	   0.2172, 0.2184, 0.2191, 0.2194, 0.2196, 2*0.2197
     `	  /
	DATA (BTIDR (I, 1, 4), I = 1, 14)
     *	  /0.5256, 0.7444, 0.9908, 1.2700, 1.5680, 1.8505, 2.0767,
     *	   2.2211, 2.2808, 2.2774, 2.2362, 2.1779, 2.1160, 2.0564
     `	  /
	DATA (BTIDR (I, 2, 4), I = 1, 14)
     *	  /0.4843, 0.6714, 0.8577, 1.0335, 1.1812, 1.2858, 1.3458,
     *	   1.3688, 1.3685, 1.3546, 1.3360, 1.3168, 1.2989, 1.2838
     `	  /
	DATA (BTIDR (I, 1, 5), I = 1, 14)
     &    /0.1498, 0.1930, 0.2201, 0.2364, 0.2460, 0.2514, 0.2544,
     &     0.2560, 0.2569, 0.2574, 0.2577, 0.2578, 0.2579, 0.2579
     `	  /

        DATA (BTIDR (I, 2, 5), I = 1, 14)
     &    /0.2184, 0.2656, 0.2927, 0.3078, 0.3159,  0.3202, 0.3224,
     &     0.3235, 0.3241, 0.3244, 0.3245, 3*0.3246
     `	  /

        DATA (BTIDR (I, 1, 6), I = 1, 14)
     &    /0.1369, 0.1681, 0.1860, 0.1958, 0.2010,  0.2038, 0.2053,
     &     0.2060, 0.2064, 0.2066, 0.2067, 3*0.2068
     `	  /

        DATA (BTIDR (I, 2, 6), I = 1, 14)
     &    /0.1969, 0.2268, 0.2416,  0.2488, 0.2521, 0.2537, 0.2544,
     &     0.2547, 0.2548, 5*0.2549
     `	  / 


	DATA (BTIDR (I, 1, 7), I = 1, 14) /14*0./
	DATA (BTIDR (I, 2, 7), I = 1, 14) /14*0./

	DATA (BTIDR (I, 1, 8), I = 1, 14) /14*0./
	DATA (BTIDR (I, 2, 8), I = 1, 14) /14*0./

	DATA (BTIDR (I, 1, 9), I = 1, 14) /14*0./
	DATA (BTIDR (I, 2, 9), I = 1, 14) /14*0./

C****
C**** --------------------------------------------------------------
	DATA (GMIDR (I, 1, 1), I = 1, 14)
     *	  /0.1582, 0.2581, 0.3227, 0.3635, 0.3882, 0.4026, 0.4108,
     *	   0.4154, 0.4179, 0.4193, 0.4200, 0.4204, 0.4206, 0.4207
     `	  /
	DATA (GMIDR (I, 2, 1), I = 1, 14)
     *	  /0.1934, 0.3141, 0.3818, 0.4200, 0.4415, 0.4533, 0.4598,
     *	   0.4633, 0.4651, 0.4662, 0.4667, 0.4671, 2*0.4672
     `	  /
	DATA (GMIDR (I, 1, 2), I = 1, 14)
     *	  /0.1347, 0.1871, 0.2277, 0.2515, 0.2651, 0.2727, 0.2768,
     *	   0.2790, 0.2801, 0.2808, 0.2811, 0.2812, 0.2813, 0.2814
     `	  /
	DATA (GMIDR (I, 2, 2), I = 1, 14)
     *	  /0.1440, 0.2217, 0.2629, 0.2839, 0.2947, 0.3003, 0.3031,
     *	   0.3046, 0.3054, 0.3058, 0.3060, 2*0.3061, 0.3062
     `	  /
	DATA (GMIDR (I, 1, 3), I = 1, 14)
     *	  /0.1372, 0.2368, 0.3235, 0.3839, 0.4229, 0.4465, 0.4602,
     *	   0.4679, 0.4722, 0.4745, 0.4758, 0.4764, 0.4768, 0.4770
     `	  /
	DATA (GMIDR (I, 2, 3), I = 1, 14)
     *	  /0.1435, 0.2524, 0.3370, 0.3955, 0.4332, 0.4563, 0.4697,
     *	   0.4773, 0.4815, 0.4839, 0.4851, 0.4858, 0.4861, 0.4863
     `	  /
	DATA (GMIDR (I, 1, 4), I = 1, 14)
     *	  /0.4298, 0.9651, 1.6189, 2.4084, 3.2992, 4.1928, 4.9611,
     *	   5.5095, 5.8085, 5.9069, 5.8726, 5.7674, 5.6346, 5.4944
     `	  /
	DATA (GMIDR (I, 2, 4), I = 1, 14)
     *	  /0.4167, 0.8974, 1.4160, 1.9414, 2.4147, 2.7803, 3.0202,
     *	   3.1468, 3.1954, 3.1932, 3.1676, 3.1328, 3.0958, 3.0625
     `	  /
        DATA (GMIDR (I, 1, 5), I = 1, 14)
     &    /0.1959, 0.3203, 0.3985, 0.4472, 0.4766, 0.4937, 0.5034,
     &     0.5088, 0.5117, 0.5134, 0.5143, 0.5147, 0.5150, 0.5152
     `	  /

        DATA (GMIDR (I, 2, 5), I = 1, 14)
     &    /0.2328, 0.3859, 0.4734, 0.5227, 0.5498, 0.5644, 0.5720,
     &     0.5761, 0.5781, 0.5792, 0.5797, 0.5800, 0.5802, 0.5802
     `	  /

        DATA (GMIDR (I, 1, 6), I = 1, 14)
     &    /0.1447, 0.2244, 0.2698, 0.2953, 0.3094, 0.3170, 0.3211,
     &     0.3233, 0.3244, 0.3250, 0.3253, 0.3255, 0.3256, 0.3256
     `	  /

        DATA (GMIDR (I, 2, 6), I = 1, 14)
     &    /0.1643, 0.2624, 0.3110, 0.3347, 0.3461, 0.3517, 0.3543,
     &     0.3556, 0.3562, 0.3564, 0.3565, 0.3566, 0.3566, 0.3566
     `	  /

	DATA (GMIDR (I, 1, 7), I = 1, 14) /14*1./
	DATA (GMIDR (I, 2, 7), I = 1, 14) /14*1./

	DATA (GMIDR (I, 1, 8), I = 1, 14) /14*1./
	DATA (GMIDR (I, 2, 8), I = 1, 14) /14*1./

	DATA (GMIDR (I, 1, 9), I = 1, 14) /14*1./
	DATA (GMIDR (I, 2, 9), I = 1, 14) /14*1./


C**** -----------------------------------------------------------

      DATA GRN /0.33, 0.67/

      DATA SNWALB /.85, .50, .85, .50,
     *             .85, .50, .85, .50,
     *             .85, .50, .85, .50,
     *             .85, .50, .85, .50,
     *             .85, .50, .85, .50,
     &             .85, .50, .85, .50,
     &             .85, .50, .85, .50,
     &             .85, .50, .85, .50,
     &             .85, .50, .85, .50
     `		  /

      DATA SNWMID /50.,50.,50.,2.,50.,2.,2.,2.,2./


CFPP$ EXPAND (COEFF)

      DO 100 I=1,IRUN
C         print*,'1 ok'
	  ALA = AMIN1 (AMAX1 (ZERO, VLAI(I)), ALATRM)
	  LAI = 1 + MAX(0, INT((ALA-BLAI)/DLAI) )
	  DX = (ALA - (BLAI+(LAI-1)*DLAI)) * (ONE/DLAI)
	  DY = (VGRN(I)- GRN(1)) * (ONE/(GRN(2) - GRN(1)))

C         print*,'2 ok'
	  ALPHA = COEFF (ALVDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)
	  BETA  = COEFF (BTVDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)
	  GAMMA = COEFF (GMVDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)

          GAMMA = MAX(GAMMA,0.01)

C         print*,'3 ok'
!	  AVISDR(I) = ALPHA - ZTH(I)*BETA / (GAMMA+ZTH(I))
	  AVISDF(I) = ALPHA-BETA
     *          + 2.*BETA*GAMMA*(1.-GAMMA*ALOG((1.+GAMMA)/GAMMA))

C         print*,'4 ok'
	  ALPHA = COEFF (ALIDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)
	  BETA  = COEFF (BTIDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)
	  GAMMA = COEFF (GMIDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)

          GAMMA = MAX(GAMMA,0.01)

C         print*,'5 ok',zth(i),gamma,ityp(i),lai,dx,dy,vgrn(i),vlai(i)
!	  ANIRDR(I) = ALPHA - ZTH(I)*BETA / (GAMMA+ZTH(I))
C         print*,'5b ok'
	  ANIRDF(I) = ALPHA-BETA
     *          + 2.*BETA*GAMMA*(1.-GAMMA*ALOG((1.+GAMMA)/GAMMA))

C         print*,'6 ok'
	  IF (SNW (I) .GT. ZERO) THEN
	   FAC = SNW(I) / (SNW(I) + SNWMID(ITYP(I)))

C         print*,'7 ok'
!	   AVISDR(I) = AVISDR(I) + (SNWALB(1,ITYP(I)) - AVISDR(I)) * FAC
!	   ANIRDR(I) = ANIRDR(I) + (SNWALB(2,ITYP(I)) - ANIRDR(I)) * FAC
	   AVISDF(I) = AVISDF(I) + (SNWALB(3,ITYP(I)) - AVISDF(I)) * FAC
	   ANIRDF(I) = ANIRDF(I) + (SNWALB(4,ITYP(I)) - ANIRDF(I)) * FAC
	  ENDIF
C         print*,'8 ok'

 100  CONTINUE

      RETURN
      END


      FUNCTION COEFF(TABLE, NTABL, LAI ,DX, DY)

      INTEGER NTABL, LAI

      REAL TABLE (NTABL, 2), DX, DY

      COEFF = (TABLE(LAI,  1)
     *      + (TABLE(LAI  ,2) - TABLE(LAI  ,1)) * DY ) * (1.0-DX)
     *      + (TABLE(LAI+1,1)
     *      + (TABLE(LAI+1,2) - TABLE(LAI+1,1)) * DY ) * DX

      RETURN
      END
