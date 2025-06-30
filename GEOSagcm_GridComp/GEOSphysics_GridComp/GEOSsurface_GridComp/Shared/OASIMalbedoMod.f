      module OCalbedoMod

      implicit none
      private

      public ocalbedo

      integer, parameter :: nlt=33
      integer, parameter :: nltgmao=8

!  Spectral seawater absorption and total scattering coefficients in 
!  units of /m.  Derived from Smith and Baker 1981 (200-370 nm),
!  and (730-800 nm), Pope and Fry 1997 (380-720), Circio and Petty 
!  1951 (800nm-2.5um), and Maul 1985 (2.5-4um).  (see Gregg and Casey 
!  2008 for complete references)

      real    aw(nlt),bw(nlt)
      integer lam(nlt)

      data lam /   250,    325,    350,    375,    400,    425,    450,
     *     475,    500,    525,    550,    575,    600,    625,    650, 
     *     675,    700,    725,    775,    850, 
     *     950,   1050,   1150,   1250,    1350,     1450,    1550,
     *    1650,   1750,    1900,  2200,    2900,     3700/

      data aw / 0.6112, 0.0762, 0.0461, 0.0182, 0.0063, 0.0051, 0.0083,
     *  0.0119, 0.0215, 0.0407, 0.0550, 0.0849, 0.1995, 0.2850, 0.3512,  
     *  0.4559, 0.6433, 1.4449, 2.3900, 3.7382, 
     *  27.481, 19.347, 67.180, 94.998, 363.126, 1118.607, 944.876, 
     *  519.60, 646.72, 3768.56, 2628.08, 437623.0, 1338404.0/  

      data bw/ 0.0567, 0.0187, 0.0135, 0.0100, 0.0076, 0.0058, 0.0045,
     * 0.0036, 0.0029, 0.0023, 0.0019, 0.0016, 0.0014, 0.0012, 0.0009,
     * 0.0007, 0.0007, 0.0006, 0.0004, 0.0002,
     * 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

!  Coefficients and factors for spectral reflectance
!  (from Frouin et al., 1996 Spectral reflectance of sea foam in the 
!  visible and near-infrared; In situ measurements and remote sensing 
!  implications.  JGR 101: 14361-14371.)

      real a0,a1,a2,a3
      real b0,b1,b2,b3

      data a0,a1,a2,a3 /0.9976, 0.21940,5.554E-2, 6.700E-3/
      data b0,b1,b2,b3 /5.0260,-0.01138,9.552E-6,-2.698E-9/

      real, parameter :: d2r = 3.14159265358979323846/180. 
      real, parameter :: r2d = 1.0/d2r

      real cos40

!  To be set by initialization routine when negative.

      real :: wfac(nlt)=-1.

      contains

      subroutine setlte

!  Sets parameters for ocean irradiance.  33 spectral bands across the 
!  entire solar spectrum are optimized for ocean biological and heat 
!  research, based on Gregg and Casey, 2008.  Skill Assessment of a 
!  Spectral Ocean-Atmosphere Radiative Model. Journal of Marine Systems 
!  http://dx.doi.org/10.1016/j.jmarsys.2008.05.007

      integer nl
      real    tlog, fac

      cos40 = cos(40.*d2r)

!  Spectral reflectance factor 

      do nl = 1,nlt
         if (lam(nl) .lt. 900)then
            tlog = alog( 1.0E-36 + exp(-(aw(nl)+0.5*bw(nl))) )
            fac = a0 + tlog*(a1 + tlog*(a2 + a3*tlog))
         else
            tlog = float(lam(nl))
            fac = b0 + tlog*(b1 + tlog*(b2 + b3*tlog))
         endif
         wfac(nl) = min(fac,1.0)
         wfac(nl) = max(fac,0.0)
      enddo

      return
      end subroutine setlte


      subroutine sfcrfl(cosz,ws,rod,ros,OASIM)

!  Computes surface reflectance for direct (rod) and diffuse (ros)
!  components separately, as a function of theta, wind speed or
!  stress.  From Gregg and Carder, 1989 A simple spectral solar 
!  irradiance model for cloudless maritime atmospheres. Limnology 
!  and Oceanography 35: 1657-1675

!  Includes spectral dependence of foam reflectance derived from Frouin
!  et al., 1996 (JGR)

!  cosz = cosine of solar zenith angle
!  ws = wind speed in m s-1
!  rod, ros = direct and diffuse spectral reflectance for OASIM bands

      real,    intent(OUT) :: rod(nlt),ros(nlt) !direct and diffuse reflectance
      real,    intent(IN)  :: ws, cosz
      logical, intent(IN)  :: OASIM

      real, parameter :: rn=1.0/1.341, roair=1.2E3
      real, parameter :: CERESFAC   = 0.068/0.077
      real, parameter :: OCNALBDF   = .08*CERESFAC
      real, parameter :: A0         = 0.40670980*CERESFAC
      real, parameter :: A1         =-1.23236340*CERESFAC
      real, parameter :: A2         = 1.42240510*CERESFAC
      real, parameter :: A3         =-0.55573341*CERESFAC

      real cn, rof,  b, rtheta, rthetar
      real rospd, rosps
      real sintr, rmin, rpls, sinp, tanp
      real sinrmin, sinrpls, tanrmin, tanrpls
      integer nl

      if(wfac(1)<0.0) call setlte

!  Foam (rof) and diffuse (rosps) reflectance


      if (ws .gt. 4.0)then
         if (ws .le. 7.0)then
            cn  = 6.2E-4 + 1.56E-3/ws
            rof = roair*cn*2.2E-5*ws*ws - 4.0E-4
         else
            cn  = 0.49E-3 + 0.065E-3*ws
            rof = (roair*cn*4.5E-5 - 4.0E-5)*ws*ws
         endif
         rosps = 0.057
      else
         rof   = 0.0
         rosps = 0.066
      endif

      if(.not.OASIM) rosps = ocnalbdf

      if(OASIM) then


!  Direct
!   Fresnel reflectance for theta < 40, ws < 2 m/s

         if (cosz .ge. cos40 .or. ws .lt. 2.0)then
            if (cosz .eq. 1.0)then
               rospd = 0.0211
            else
               rtheta  = acos(cosz)
               sintr   = sqrt(1.0-cosz**2)/rn
               rthetar = asin(sintr)
               print *, rtheta,rthetar,cosz,sintr
               rmin    = rtheta - rthetar
               rpls    = rtheta + rthetar
               sinrmin = sin(rmin)
               sinrpls = sin(rpls)
               tanrmin = tan(rmin)
               tanrpls = tan(rpls)

               sinp = (sinrmin*sinrmin)/(sinrpls*sinrpls)
               tanp = (tanrmin*tanrmin)/(tanrpls*tanrpls)
               rospd = 0.5*(sinp + tanp)

            endif

         else

!   Empirical fit otherwise

            rtheta  = acos(cosz)*r2d
            b = -7.14E-4*ws + 0.0618
            rospd = 0.0253*exp(b*(rtheta-40.0))

         endif

      else

         rosps = ocnalbdf
         rospd = A0+(A1+(A2+A3*cosz)*cosz)*cosz
      endif

!  Reflectance totals

      do nl = 1,nlt
         rod(nl) = rospd + rof*wfac(nl)
         ros(nl) = rosps + rof*wfac(nl)
      enddo

      return
      end subroutine sfcrfl


      subroutine ocalbedo(cosz,ws,albd,albs,OASIM)

!  Computes ocean surface albedo from solar zenith angle (cosz) 
!  and wind speed (ws, m/s).  Computation is done fo rthe 33 bands
!  of OASIM (Gregg and Casey 2008), but albedo is converted to the 8
!  bands of the GEOS-5 model.
!  Albedo is provided as direct (albd) and diffuse (albs).     

!  cosz = cosine of solar zenith angle
!  ws   = wind speed in m s-1
!  albd, albs = direct and diffuse spectral reflectance for GEOS-5 bands

      logical, intent(IN)  :: OASIM
      real,    intent(IN)  :: cosz, ws
      real,    intent(OUT) :: albd(nltgmao),albs(nltgmao)

      real rod(nlt),ros(nlt)
      integer n, nl

!  Derive surface reflectance as a function of cosz and ws
!  for OASIM Bands

      call sfcrfl(cosz,ws,rod,ros,OASIM)

 !  Albedo at GMAO Bands

      albd(1) = rod(1)
      albs(1) = ros(1)

      albd(2) = rod(1)
      albs(2) = ros(1)

      albd(3) = rod(2)
      albs(3) = ros(2)

      albd(4) = (rod(2)*0.52 + rod(3) + rod(4))/2.52
      albs(4) = (ros(2)*0.52 + ros(3) + ros(4))/2.52

      albd(5:) = 0.0
      albs(5:) = 0.0

      n = 0
      do nl = 5,17
         albd(5) = albd(5) + rod(nl)
         albs(5) = albs(5) + ros(nl)
         n = n+1
      enddo
      albd(5) = albd(5)/float(n)
      albs(5) = albs(5)/float(n)

      n = 0
      do nl = 18,23
         albd(6) = albd(6) + rod(nl)
         albs(6) = albs(6) + ros(nl)
         n = n+1
      enddo
      albd(6) = albd(6)/float(n)
      albs(6) = albs(6)/float(n)

      n = 0
      do nl = 24,31
         albd(7) = albd(7) + rod(nl)
         albs(7) = albs(7) + ros(nl)
         n = n+1
      enddo
      albd(7) = albd(7)/float(n)
      albs(7) = albs(7)/float(n)

      n = 0
      do nl = 32,nlt
         albd(8) = albd(8) + rod(nl)
         albs(8) = albs(8) + ros(nl)
         n = n+1
      enddo
      albd(8) = albd(8)/float(n)
      albs(8) = albs(8)/float(n)

      return
      end subroutine ocalbedo


      end module OCalbedoMod
