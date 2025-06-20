      subroutine radmod(zd,Edtop,Estop,rmud,a,bt,bb,Dmax,Edz,Esz,Euz)
 
!  Model of irradiance in the water column.  Accounts for three 
!  irradiance streams:
 
!  Edz = direct downwelling irradiance
!  Esz = diffuse downwelling irradiance
!  Euz = diffuse upwelling irradiance
 
!  Uses Ackelson's (1994, JGR) mod's to the Aas (1987, AO) 
!  two-stream model.
 
!  Propagation is done in energy units, tests are done in quanta,
!  final is quanta for phytoplankton growth.
 
!  Commented out terms produce a max error of 
!  0.8% in Esz for a > 0.004 and bb > 0.0001 and
!  3.9% in Euz for a > 0.004 and bb > 0.00063
 
#define SPINUP 0
#if SPINUP
#define EUZ 0
#else
#define EUZ 1
#endif
      data rbot /0.0/ !bottom reflectance
      data rd /1.5/   !these are taken from Ackleson, et al. 1994 (JGR)
      data ru /3.0/
 
!  Constants
      rmus = 1.0/0.83            !avg cosine diffuse down
      rmuu = 1.0/0.4             !avg cosine diffuse up
 
!  Downwelling irradiance: Edz, Esz
!  Compute irradiance components at depth
      cd = (a+bt)*rmud
      Edz = Edtop*exp(-cd*zd)
      au = a*rmuu
      Bu = ru*bb*rmuu
      Cu = au+Bu
      as = a*rmus
      Bs = rd*bb*rmus
      Cs = as+Bs
      Bd = bb*rmud
      Fd = (bt-bb)*rmud
      bquad = Cs - Cu
      cquad = Bs*Bu - Cs*Cu
      sqarg = bquad*bquad - 4.0*cquad
      a1 = 0.5*(-bquad + sqrt(sqarg))
      a2 = 0.5*(-bquad - sqrt(sqarg))
      S = -(Bu*Bd + Cu*Fd)
      SEdz = S*Edz
      a2ma1 = a2 - a1
      rM = SEdz/(a1*a2ma1)
      rN = SEdz/(a2*a2ma1)
!      ea2Dmax = exp(a2ma1*Dmax)
!      c1 = (rN-rM)*exp(-a1*Dmax) - (Estop-rM+rN)*ea2Dmax
!     *                             /(1.0-ea2Dmax)
!      c2 = Estop - rM + rN - c1
      c2 = Estop - rM + rN 
!      a1arg = a1*zd
!      a1arg = min(a1arg,82.8)
!      Ta1z = exp(a1arg)
      Ta2z = exp(a2*zd)
!      Esz = c1*Ta1z + c2*Ta2z + rM - rN
      Esz = c2*Ta2z + rM - rN
      Esz = max(Esz,0.0)
!      Eutmp = ((a1+Cs)*c1)*Ta1z + ((a2+Cs)*c2)*Ta2z + Cs*rM
!     *             - Cs*rN - Fd*Edz
#if EUZ
      Eutmp = ((a2+Cs)*c2)*Ta2z + Cs*rM - Cs*rN - Fd*Edz
      Euz = Eutmp/Bu
      Euz = max(Euz,0.0)
#else
      Euz = 0.0
#endif
 
      return
      end
 
 
!! ******************************************************************
!      subroutine upeu 
  
!!  Computes surface upwelling irradiance.  The approximation used for
!!  upwelling irradiance a depth, i.e., that rmud = rmus, is not
!!  valid at the surface, and a full treatment of diffuse and direct
!!  path lengths is required.  
  
!      save
!#include "gloparam.F.h"
!#include "radpth.h"
!#include "comlte.h"
!#include "ostate.F.h"
!      real Edtop(nlt),Estop(nlt)
!      real sfceu(nlt)
!      real bbc(10)
!      data ifst /0/
!      data bbc / 0.002, 0.00071, 0.0032, 0.00071, 0.0029,
!     *           0.0,   0.0,     0.0,    0.0,     0.0/
  
!!  Constants
!      if (ifst .eq. 0)then
!       rmus = 1.0/0.83            !avg cosine diffuse down
!       rmuu = 1.0/0.4             !avg cosine diffuse up
!       bbw = 0.5    !backscattering to forward scattering ratio
!       rbot = 0.0                 !bottom reflectance
!       rd = 1.5   !these are taken from Ackleson, et al. 1994 (JGR)
!       ru = 3.0
!       Dmax = 500.0    !depth at which Ed = 0
!       do nl = 1,nlt
!        sfceu(nl) = 0.0
!       enddo
!       ifst = 1
!      endif
  
!      m = indext2
!      k = 1
!      do i = inwst,inwnd
!       do nl = 1,nlt
!        Edtop(nl) = Ed(i,nl)
!        Estop(nl) = Es(i,nl)
!       enddo
!       do nl = 1,nlt
!        a = aw(nl) + acdom(nl)
!        bt = bw(nl)
!        bb = bbw*bw(nl)
!        do n = nnut+1,ntyp-nzoo
!         a = a + P(i,k,m,n)*ac(n-nnut,nl)
!         bt = bt + P(i,k,m,n)*bc(n-nnut,nl)
!         bb = bb + P(i,k,m,n)*bbc(n-nnut)*bc(n-nnut,nl)
!        enddo
!!  Compute Eu
!!   Eu from Ed
!        ad = a*rmud(i)
!        bd = rd*bb*rmud(i)
!        au = a*rmuu
!        bu = ru*bb*rmuu
!        cd = ad+bd
!        cu = au+bu
!        bquad = cd - cu
!        cquad = bd*bu - cd*cu
!        sqarg = bquad*bquad - 4.0*cquad
!        a2 = 0.5*(-bquad - sqrt(sqarg))
!        sfceu1 = (a2+cd)/bu
!        EuEd = Etmp*Ta2z
  
!!   Eu from Es
!        as = a*rmus
!        bs = rd*bb*rmus
!        cs = as+bs
!        bquad = cs - cu
!        cquad = bs*bu - cs*cu
!        sqarg = bquad*bquad - 4.0*cquad
!        a2 = 0.5*(-bquad - sqrt(sqarg))
!        sfceu2 = (a2+cs)/bu
!        sfceu(nl) = Edtop(nl)*sfceu1 + Estop(nl)*sfceu2
!       enddo
!!       if (i .eq. 334)then
!!        write(6,*)ihr,sfceu(10)
!!       endif
!      enddo
  
!      return
!      end
