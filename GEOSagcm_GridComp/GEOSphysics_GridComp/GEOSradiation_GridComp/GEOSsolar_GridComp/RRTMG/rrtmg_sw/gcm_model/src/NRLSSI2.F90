module NRLSSI2

   ! -----------------------------------------------------------------
   ! Explanation of the NRLSSI2 model of solar irradiance in RRTMG SW:
   ! -----------------------------------------------------------------
   ! The Sun has an approximately 11-year magnetic cycle. Faculae brightening
   ! is measured by the Bremen Mg II index. Sunspot dimming is measured by the
   ! SPOT67 index. Both have different wavelength signatures than the dominant
   ! quiet sun.
   !    The RRTMG_SW NRLSSI2 model averages over 12 of these solar cycles
   ! (cycles 13-24, late 1800’s to near present) to come up with an 11-year
   ! “average solar cycle”. This is the idealized 11-year solar cycle (call
   ! it “AvgCyc11”) used in RRTMG_SW.
   !    The spectral solar irradiance, SSI(𝜆,t), and its spectral integral,
   ! TSI(t), a.k.a. the Solar “Constant”, are defined for the mean Sun-Earth
   ! distance of 1AU. The spectral and time dependence are assumed independent:
   !    SSI(𝜆,t) = SVAR_I    * Irradnce(𝜆) +
   !               SVAR_F(t) * Facbrght(𝜆) +
   !               SVAR_S(t) * Snsptdrk(𝜆)
   ! The spectral terms are valid for the time mean of AvgCyc11 (denoted <>),
   ! and are hardcoded into RRTMG_SW. The time-dependent SVAR_(t) terms are
   ! used to scale these spectral responses to model solar variability. The
   ! spectral terms integrated over 𝜆 are:
   !    Iint = 1360.37, Fint = 0.996047, Sint = -0.511590 [in Wm-2]
   ! These sum to <TSI> = 1360.85 Wm-2, the NRLSSI2 mean solar constant, cf.
   ! the Kurucz value (1368.22 Wm-2).
   !    The SVAR_F/S are modeled as linear functions of the Mg and SB indices,
   ! and must have a <SVAR_F|S> = 1. SVAR_I is 1 for AvgCyc11. This gives the
   ! correct <TSI>.
   ! -----------------------------------------------------------------

   use iso_fortran_env, only : error_unit

   implicit none

   ! state is saved outside of scope
   save

   ! all fields and methods hidden unless declared public below
   private

   ! Mean quiet sun, facular brightening, and sunspot dimming contributions to TSI
   ! for NRLSSI2, 100-50000 cm-1 (spectrally integrated from hi-res values after
   ! mapping to g-point space). Mean is time-average over Solar cycles 13-24.
   ! -----------------------------------------------------------------------------

   real, parameter :: Iint = 1360.37        ! mean quiet sun TSI contribution
   real, parameter :: Fint = 0.996047       ! mean facular brightening TSI contribn
   real, parameter :: Sint = -0.511590      ! mean sunspot dimming TSI contribution

   ! "Average 11-year solar cycle" of Mg and SB indices developed from
   ! Solar cycles 13-24. There are 132 (11*12) mid-month values (elements
   ! 2-132) and identical cycle start/end values (elements 1 and 134).
   ! -----------------------------------------------------------------

   ! number of array elements
   integer, parameter :: nsolfrac = 134

   ! Facular index from NRLSSI2 Mg "Bremen" index
   real, parameter :: mgavgcyc (nsolfrac) = (/ &
     &   0.150737,  0.150746,  0.150733,  0.150718,  0.150725,  0.150762, &
     &   0.150828,  0.150918,  0.151017,  0.151113,  0.151201,  0.151292, &
     &   0.151403,  0.151557,  0.151766,  0.152023,  0.152322,  0.152646, &
     &   0.152969,  0.153277,  0.153579,  0.153899,  0.154252,  0.154651, &
     &   0.155104,  0.155608,  0.156144,  0.156681,  0.157178,  0.157605, &
     &   0.157971,  0.158320,  0.158702,  0.159133,  0.159583,  0.160018, &
     &   0.160408,  0.160725,  0.160960,  0.161131,  0.161280,  0.161454, &
     &   0.161701,  0.162034,  0.162411,  0.162801,  0.163186,  0.163545, &
     &   0.163844,  0.164029,  0.164054,  0.163910,  0.163621,  0.163239, &
     &   0.162842,  0.162525,  0.162344,  0.162275,  0.162288,  0.162369, &
     &   0.162500,  0.162671,  0.162878,  0.163091,  0.163251,  0.163320, &
     &   0.163287,  0.163153,  0.162927,  0.162630,  0.162328,  0.162083, &
     &   0.161906,  0.161766,  0.161622,  0.161458,  0.161266,  0.161014, &
     &   0.160666,  0.160213,  0.159690,  0.159190,  0.158831,  0.158664, &
     &   0.158634,  0.158605,  0.158460,  0.158152,  0.157691,  0.157152, &
     &   0.156631,  0.156180,  0.155827,  0.155575,  0.155406,  0.155280, &
     &   0.155145,  0.154972,  0.154762,  0.154554,  0.154388,  0.154267, &
     &   0.154152,  0.154002,  0.153800,  0.153567,  0.153348,  0.153175, &
     &   0.153044,  0.152923,  0.152793,  0.152652,  0.152510,  0.152384, &
     &   0.152282,  0.152194,  0.152099,  0.151980,  0.151844,  0.151706, &
     &   0.151585,  0.151496,  0.151437,  0.151390,  0.151347,  0.151295, &
     &   0.151220,  0.151115,  0.150993,  0.150883,  0.150802,  0.150752, &
     &   0.150729,  0.150737/)

   ! Sunspot index from NRLSSI2 SB "SPOT67" index
   real, parameter :: sbavgcyc (nsolfrac) = (/ &
     &    50.3550,   44.1322,   52.0179,   59.2231,   66.3702,   71.7545, &
     &    76.8671,   83.4723,   91.1574,   98.4915,  105.3173,  115.1791, &
     &   130.9432,  155.0483,  186.5379,  221.5456,  256.9212,  291.5276, &
     &   325.2953,  356.4789,  387.2470,  422.8557,  466.1698,  521.5139, &
     &   593.2833,  676.6234,  763.6930,  849.1200,  928.4259,  994.9705, &
     &  1044.2605, 1087.5703, 1145.0623, 1224.3491, 1320.6497, 1413.0979, &
     &  1472.1591, 1485.7531, 1464.1610, 1439.1617, 1446.2449, 1496.4323, &
     &  1577.8394, 1669.5933, 1753.0408, 1821.9296, 1873.2789, 1906.5240, &
     &  1920.4482, 1904.6881, 1861.8397, 1802.7661, 1734.0215, 1665.0562, &
     &  1608.8999, 1584.8208, 1594.0162, 1616.1486, 1646.6031, 1687.1962, &
     &  1736.4778, 1787.2419, 1824.9084, 1835.5236, 1810.2161, 1768.6124, &
     &  1745.1085, 1748.7762, 1756.1239, 1738.9929, 1700.0656, 1658.2209, &
     &  1629.2925, 1620.9709, 1622.5157, 1623.4703, 1612.3083, 1577.3031, &
     &  1516.7953, 1430.0403, 1331.5112, 1255.5171, 1226.7653, 1241.4419, &
     &  1264.6549, 1255.5559, 1203.0286, 1120.2747, 1025.5101,  935.4602, &
     &   855.0434,  781.0189,  718.0328,  678.5850,  670.4219,  684.1906, &
     &   697.0376,  694.8083,  674.1456,  638.8199,  602.3454,  577.6292, &
     &   565.6213,  553.7846,  531.7452,  503.9732,  476.9708,  452.4296, &
     &   426.2826,  394.6636,  360.1086,  324.9731,  297.2957,  286.1536, &
     &   287.4195,  288.9029,  282.7594,  267.7211,  246.6594,  224.7318, &
     &   209.2318,  204.5217,  204.1653,  200.0440,  191.0689,  175.7699, &
     &   153.9869,  128.4389,  103.8445,   85.6083,   73.6264,   64.4393, &
     &    56.5779,   50.3550/)

   ! Fractional interval length of mgavgcyc and sbavgcyc and its half
   real, parameter :: intrvl_len = 1.0 / (nsolfrac-2)
   real, parameter :: intrvl_len_hf = 0.5 * intrvl_len

   ! Specification of linear dependence of SVAR_F/S on the Mg and SB indices
   ! -----------------------------------------------------------------------
   ! SVAR_F = (Mg - Mg_0) / (Mg_avg - Mg_0), where Mg_avg = <Mg>, and
   ! SVAR_S = (SB - SB_0) / (SB_avg - SB_0), where SB_avg = <SB>, where
   ! <> is the time average over AvgCyc11 (i.e., over the 132 inner values
   ! of the respective mgavgcyc or sbavgcyc). Clearly, e.g., <SVAR_F> = 1.

   ! means over AvgCyc11
   real, parameter :: Mg_avg = 0.1567652  ! NRLSSI2 Mg "Bremen" index
   real, parameter :: SB_avg = 909.71260  ! NRLSSI2 SB "SPOT67" index

   ! offsets where SVAR is zero
   real, parameter :: Mg_0 = 0.14959542   ! facular Mg offset
   real, parameter :: SB_0 = 0.00066696   ! sunspot SB offset

   ! initialization status used to avoid repeating once-only calculations
   ! --------------------------------------------------------------------
   logical :: NRLSSI2_initialized = .false.

   ! and related variables
   integer :: isolvar_old = -999
   real :: indsolvar_old(2) = -999.

   ! module level saved variables for isolvar==1 case (calculated only once)
   ! ------------------------------------------------------------------------
   real :: isolvar_1_mean_svar_f
   real :: isolvar_1_mean_svar_s

   ! public interface
   ! ----------------
   public :: initialize_NRLSSI2
   public :: adjust_solcyc_amplitudes
   public :: interpolate_indices
   public :: Iint, Fint, Sint
   public :: Mg_avg, Mg_0
   public :: SB_avg, SB_0
   public :: isolvar_1_mean_svar_f
   public :: isolvar_1_mean_svar_s

contains

   ! Initialization (once-only calculations).
   ! Should be called each time, but will only do calculations as necessary,
   ! that is, when either its never been called or the argument list changes.

   subroutine initialize_NRLSSI2 (isolvar, indsolvar_opt)

      integer, intent(in) :: isolvar
      real, intent(in), optional :: indsolvar_opt (2)

      real :: indsolvar (2)
      real :: solcycfr, indsolvar_scl (2)
      real :: iscl1_mean, iscl1_Mg_mean
      real :: iscl2_mean, iscl2_SB_mean
      logical :: rerun, scl1, scl2
      integer :: n

      ! actualize optional arguments
      if (present(indsolvar_opt)) then
         indsolvar = indsolvar_opt
      else
         indsolvar = 1.
      endif

      ! decide if need to rerun initialization
      rerun = (.not. NRLSSI2_initialized)
      if (.not. rerun) rerun = (isolvar .ne. isolvar_old)
      if (.not. rerun) rerun = (.not. all(indsolvar == indsolvar_old))
      if (rerun) then

         ! means needed for isolvar == 1
         if (isolvar == 1) then

            ! no scaling defaults
            isolvar_1_mean_svar_f = 1.
            isolvar_1_mean_svar_s = 1.

            scl1 = (indsolvar(1) .ne. 1.)
            scl2 = (indsolvar(2) .ne. 1.)
            if (scl1 .or. scl2) then

               ! mean <indsolvar_scl> over a cycle
               if (scl1) iscl1_mean = (1. + indsolvar(1)) / 2.
               if (scl2) iscl2_mean = (1. + indsolvar(2)) / 2.

               ! mean <iscl{1,2} * {Mg,SB}> over a cycle
               iscl1_Mg_mean = 0.
               iscl2_SB_mean = 0.
               solcycfr = intrvl_len_hf
               do n = 2,nsolfrac-1
                  call adjust_solcyc_amplitudes(solcycfr, indsolvar, indsolvar_scl)
                  if (scl1) iscl1_Mg_mean = iscl1_Mg_mean + indsolvar_scl(1) * mgavgcyc(n)
                  if (scl2) iscl2_SB_mean = iscl2_SB_mean + indsolvar_scl(2) * sbavgcyc(n)
                  solcycfr = solcycfr + intrvl_len
               enddo
               if (scl1) iscl1_Mg_mean = iscl1_Mg_mean / (nsolfrac-2)
               if (scl2) iscl2_SB_mean = iscl2_SB_mean / (nsolfrac-2)

               ! means with indsolvar scaling
               if (scl1) isolvar_1_mean_svar_f = (iscl1_Mg_mean - iscl1_mean * Mg_0) / (Mg_avg - Mg_0)
               if (scl2) isolvar_1_mean_svar_s = (iscl2_SB_mean - iscl2_mean * SB_0) / (SB_avg - SB_0)

            end if

         end if

         ! "once-only" calculations done
         isolvar_old = isolvar
         indsolvar_old = indsolvar
         NRLSSI2_initialized = .true.

      end if

   end subroutine initialize_NRLSSI2


   ! Adjust amplitude scaling of mean solar cycle to be unity at
   ! solar minimum (solcycfrac_min), to be the requested indsolvar
   ! at solar maximum (solcycfrac_max), and to vary linearly with
   ! solcycfrac between those values.

   subroutine adjust_solcyc_amplitudes( &
      solcycfr, indsolvar, indsolvar_scl)

      real, intent(in) :: solcycfr  ! position inn cycle [0,1]

      ! see subroutine comment
      ! 1=faculae scaling, 2=sunspot scaling
      real, intent(in)  :: indsolvar     (2)
      real, intent(out) :: indsolvar_scl (2)

      ! timing of minima and maxima within solar cycle
      real, parameter :: solcycfrac_min = 0.0189    ! Solar cycle frac at solar minimum
      real, parameter :: solcycfrac_max = 0.3750    ! Solar cycle frac at solar maximum
      real, parameter :: fracdiff_min2max = solcycfrac_max - solcycfrac_min
      real, parameter :: fracdiff_max2min = 1. - fracdiff_min2max

      real :: wgt

      if (solcycfr >= 0. .and. solcycfr < solcycfrac_min) then
         wgt = (solcycfr+1.-solcycfrac_max) / fracdiff_max2min
         indsolvar_scl(1) = indsolvar(1) + wgt * (1.-indsolvar(1))
         indsolvar_scl(2) = indsolvar(2) + wgt * (1.-indsolvar(2))
      elseif (solcycfr >= solcycfrac_min .and. solcycfr <= solcycfrac_max) then
         wgt = (solcycfr-solcycfrac_min) / fracdiff_min2max
         indsolvar_scl(1) = 1. + wgt * (indsolvar(1)-1.)
         indsolvar_scl(2) = 1. + wgt * (indsolvar(2)-1.)
      elseif (solcycfr > solcycfrac_max .and. solcycfr <= 1.) then
         wgt = (solcycfr-solcycfrac_max) / fracdiff_max2min
         indsolvar_scl(1) = indsolvar(1) + wgt * (1.-indsolvar(1))
         indsolvar_scl(2) = indsolvar(2) + wgt * (1.-indsolvar(2))
      else
         write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
         error stop 'RRTMG_SW: solcycfr must be in [0,1]'
      endif

   end subroutine adjust_solcyc_amplitudes


   ! Evaluate Mg and SB indices from mean solar cycle AvgCyc11
   !   at solcycfr in [0,1].

   subroutine interpolate_indices (solcycfr, Mg, SB)

      real, intent(in) :: solcycfr  ! position inn cycle [0,1]
      real, intent(out) :: Mg, SB   ! faculae and sunspot indices

      integer :: sfid
      real :: fraclo, frachi, intfrac

      if (solcycfr > 0. .and. solcycfr < 1.) then

         if (solcycfr <= intrvl_len_hf) then

            ! Initial half interval
            sfid = 1
            fraclo = 0.
            frachi = intrvl_len_hf

         elseif (solcycfr > intrvl_len_hf .and. solcycfr < 1. - intrvl_len_hf) then

            ! Intervening whole intervals
            sfid = floor((solcycfr-intrvl_len_hf) * (nsolfrac-2)) + 2
            fraclo = (sfid-2) * intrvl_len + intrvl_len_hf
            frachi = fraclo + intrvl_len

         elseif (solcycfr >= 1. - intrvl_len_hf) then

            ! Final half interval
            sfid = (nsolfrac-2) + 1
            fraclo = 1. - intrvl_len_hf
            frachi = 1.

         endif

         ! interpolate to Mg and SB at solcycfr
         intfrac = (solcycfr - fraclo) / (frachi - fraclo)
         Mg = mgavgcyc(sfid) + intfrac * (mgavgcyc(sfid+1) - mgavgcyc(sfid))
         SB = sbavgcyc(sfid) + intfrac * (sbavgcyc(sfid+1) - sbavgcyc(sfid))

      elseif (solcycfr == 0.) then

         Mg = mgavgcyc(1)
         SB = sbavgcyc(1)

      elseif (solcycfr == 1.) then

         Mg = mgavgcyc(nsolfrac)
         SB = sbavgcyc(nsolfrac)

      else

         write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
         error stop 'RRTMG_SW: solcycfr must be in [0,1]'

      endif

   end subroutine interpolate_indices

end module NRLSSI2
