module cfmip_data_mod

  use MAPL
  use GEOS_UtilsMod, only : GEOS_Qsat, GEOS_DQsat

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public CFMIP_IC


logical, parameter :: debug=.FALSE.

contains

SUBROUTINE CFMIP_IC                                          &
              ( FILENAME, NT, NLEVEL, NLAYR,                 &
                PREF_MODEL_E,                                &
                TIME,                                        &
		YY,                                          &
		MO,                                          &
                DD,                                          &
                HH,                                          &
                MM,                                          &
                PCP_OBS,                                     &
		TS_AIR,                                      &
		TG_SOIL,                                     &
                TSKIN,                                       &
                QSFCAIR,                                     &
                QSKIN,                                       &
                PSFC,                                        &
		LHF,                                         & 
		SHF,                                         & 
                ptop,                                        &
                tt,                                          &
                qq,                                          &
                uu,                                          &
                vv,                                          &
                oo,                                          &
                T_H_adv,                                     &
                T_V_adv,                                     &
                Q_H_adv,                                     &
                Q_V_adv,                                     &           
                P_MODEL_E                                    )

      

        IMPLICIT NONE

!---------------------------------------------------------------
! ARGUMENTS:
!---------------------------------------------------------------
! INPUT :
!-------
! FILENAME:     file name of data
! NT:           time slice number of data
! NLEVEL:       number of vertical levels of data
! NLAYR:        number of vertical layers of Single column model
! PREF_MODEL_E: Model edge-pressures (hPa) 
!	
! OUTPUT
!-------
! TIME
! YY
! MO
! DD
! HH
! MM
! PCP_OBS
! TS_AIR
! TG_SOIL
! TSKIN
! QSFCAIR
! QSKIN
! PSFC
! LHF
! SHF
! ptop
! tt
! qq
! uu
! vv
! oo
! T_H_adv
! T_V_adv
! Q_H_adv
! Q_V_adv
!---------------------------------------------------------------
! Declarations for input args
	INTEGER,  INTENT(IN) :: Nt, NLEVEL, NLAYR
	CHARACTER(LEN=*), INTENT(IN) :: FILENAME
	REAL,  DIMENSION(0:NLAYR),     INTENT(IN   ) :: PREF_MODEL_E

! Declarations for output args - scalars and single level
        real ptop
        real, dimension(nt), intent(out) :: time, yy, mo, dd, hh, mm
        real, dimension(nt), intent(out) :: pcp_obs, tskin, qsfcair
        real, dimension(nt), intent(out) :: qskin, psfc, shf, lhf, ts_air, tg_soil

! Declarations for output args - multiple-level
        real, dimension(nt,nlayr),intent(out):: tt, qq, uu, vv   
        real, dimension(nt,nlayr),intent(out):: T_H_adv, T_V_adv
        real, dimension(nt,nlayr),intent(out):: Q_H_adv, Q_V_adv
        real, dimension(nt,0:nlayr), intent(out):: P_MODEL_E, OO

! local vars

	integer :: IRD, UNIT
        real, dimension(nlayr) :: p_model

! declaration for stuff in data set
        integer bdate
        real phis, lon, lat
        real, dimension(64) :: lev_new
        real, dimension(nlevel) :: lev
        integer, dimension(nt) :: tsec
        real, dimension(nt) :: lhflx, shflx, ps, ptend, tsair, tg, solin
        real, dimension(nt) :: zenith, dayfrac, delta, ecc, u_srf, v_srf, rh_srf, alb_srf
        real, dimension(nt,nlevel) :: tdata, qdata, udata, vdata, omegadata
        real, dimension(nt,nlevel) :: divT, divq, vertdivT, vertdivq
        real, dimension(nt,nlevel) :: divT3d, divq3d, div, o3mmr

       !! real ,external :: InterpInPres
!
        ird = 1
        pcp_obs=0.
        tskin=0.
        qsfcair=0.
        qskin=0.
        psfc=0. 
        shf=0.
        lhf=0.
        ts_air=0.
        tg_soil=0.

        tt=0. 
        qq=0. 
        uu=0. 
        vv=0.  
        oo=0.  
        T_H_adv=0.
        T_V_adv=0.
        Q_H_adv=0.
        Q_V_adv=0.
        P_MODEL_E=0.
        
        time=0.
 
        UNIT = GETFILE(FILENAME, form = "unformatted")

        read(UNIT) bdate
        read(UNIT) phis
        read(UNIT) lon
        read(UNIT) lat
        read(UNIT) lev
        read(UNIT) tsec
        read(UNIT) lhflx 
        read(UNIT) shflx 
        read(UNIT) ps 
        read(UNIT) ptend 
        read(UNIT) tsair 
        read(UNIT) tg 
        read(UNIT) solin
        read(UNIT) zenith
        read(UNIT) dayfrac
        read(UNIT) delta
        read(UNIT) ecc
        read(UNIT) u_srf
        read(UNIT) v_srf
        read(UNIT) rh_srf
        read(UNIT) alb_srf

        read(UNIT) tdata 
        read(UNIT) qdata 
        read(UNIT) udata 
        read(UNIT) vdata 
        read(UNIT) omegadata
        read(UNIT) divT 
        read(UNIT) divq 
        read(UNIT) vertdivT 
        read(UNIT) vertdivq
        read(UNIT) divT3d 
        read(UNIT) divq3d 
        read(UNIT) div 
        read(UNIT) o3mmr

        yy(ird) = 1.0*1998
        mo(ird) = 1.0*7
        dd(ird) = 1.0*15
        hh(ird) = 1.0*0
        mm(ird) = 1.0*0

        pcp_obs(ird) = -999.9
        ts_air(ird) = tsair(ird)
        tg_soil(ird) = tg(ird)
        tskin(ird) = tg(ird)
        psfc(ird) = ps(ird)
        qskin(ird) = 1000.*geos_qsat(tskin(ird),psfc(ird)/100.)*rh_srf(ird)
        qsfcair(ird) = qdata(ird,nlevel)

        ptop = minval(lev)

        P_MODEL_E(ird,0:) =  pref_model_e(0:)*psfc(ird)/pref_model_e(nlayr)
        p_model(:) = 0.5* ( p_model_e(ird,0:nlayr-1) + p_model_e(ird,1:nlayr) )

               if (debug) then
                  write(*,*) "THIS IS PREF in CFMIP" 
                  write(*,*) pref_model_e
                  write(*,*) "THIS IS P_MODEL_E in CFMIP" 
                  write(*,*) p_model_e
                  write(*,*) "THIS IS P_MODEL in CFMIP" 
                  write(*,*) p_model
               endif

!        tdata(ird,:) = (/239.096, 228.1121, 221.6722, 217.8232, 215.1689, 212.6234, 210.2998, 208.3557, 206.8936, 205.8506, 205.318, 205.2797, 206.054, 207.3618, 209.0052, 210.8095, 212.7212, 214.7583, 216.9371, 219.2846, 221.7891, 224.4648, 227.2775, 230.2074, 233.2128, 236.3087, 239.45, 242.6231, 245.8229, 248.9953,  252.1724, 255.4494, 258.7444, 261.9404, 265.0162, 267.8944, 270.5506, 272.9888, 275.3893, 277.7449, 279.9416, 281.9912, 284.0184, 285.9781, 287.74, 289.2596, 290.5135, 291.64, 292.6471, 293.5434, 294.3362, 295.032, 295.1397, 283.9792, 285.074, 286.2524, 287.2851, 288.1753, 288.929, 289.5704, 290.1137, 290.551, 290.8827, 291.1213/)

!        qdata(ird,:) = (/2.738439e-06, 2.386476e-06, 2.167775e-06, 2.022134e-06, 1.943735e-06, 1.91291e-06, 1.922756e-06, 1.994845e-06, 2.177456e-06, 2.538252e-06, 3.02783e-06, 3.38639e-06, 3.708099e-06, 4.088135e-06, 4.774295e-06, 5.848523e-06, 7.479741e-06, 9.693273e-06, 1.269689e-05, 1.668489e-05, 2.230726e-05, 3.001641e-05, 4.038344e-05, 5.391814e-05, 7.088484e-05, 9.229892e-05, 0.0001177767, 0.0001470888, 0.000184272, 0.0002343657, 0.000310628, 0.0004146358, 0.0005423331, 0.0006852034, 0.0008001997, 0.0009065709, 0.001097643, 0.001395561, 0.00175284, 0.002025659, 0.002289154, 0.002710893, 0.003151285, 0.003462769, 0.003690755, 0.003811365, 0.003788855, 0.003731875, 0.003731875, 0.003731875, 0.003731875, 0.003731875, 0.003731875, 0.01021351, 0.01021351, 0.01021351, 0.01021351, 0.01021351, 0.01021351, 0.01021351, 0.01021351, 0.01021351, 0.01021351, 0.01021351/)

!        lev_new = (/494.4198, 1483.253, 2472.099, 3460.937, 4449.765, 5438.624, 6427.467, 7416.284, 8405.158, 9394.112, 10383.19, 11372.46, 12362.02, 13351.92, 14342.19, 15341.24, 16380.03, 17482.86, 18653.28, 19894.93, 21211.71, 22607.61, 24086.83, 25653.76, 27312.95, 29069.2, 30927.47, 32892.95, 34971.04, 37167.38, 39487.86, 41937.78, 44520.95, 47227.08, 50031.8, 52908.93, 55830.45, 58767.62, 61697.88, 64606.13, 67479.39, 70305.65, 73070.53, 75759.24, 78359.75, 80861.59, 83252.66, 85520.44, 87655.99, 89652.63, 91502.81, 93198.96, 93468.15, 93468.34, 94736.95, 96116.12, 97336.52, 98397.48, 99302.08, 100076.6, 100736.0, 101269.1, 101674.7, 101967.2 /)

!        tt(ird,:)  = InterpInPres ( tdata(ird,:) , lev_new , p_model )
!        qq(ird,:)  = InterpInPres ( qdata(ird,:) , lev_new , p_model )*1000.
        tt(ird,:)  = InterpInPres ( tdata(ird,:) , lev , p_model )
        qq(ird,:)  = InterpInPres ( qdata(ird,:) , lev , p_model )*1000.
        uu(ird,:)  = InterpInPres ( udata(ird,:) , lev , p_model )
        vv(ird,:)  = InterpInPres ( vdata(ird,:) , lev , p_model )
        oo(ird,0:) = InterpInPres ( omegadata(ird,:) , lev , p_model_e(ird,0: ) )

        T_H_adv(ird,:) = -1.0 * InterpInPres ( divT(ird,:) , lev , p_model )
        T_V_adv(ird,:) = -1.0 * InterpInPres ( vertdivT(ird,:) , lev , p_model )
        Q_H_adv(ird,:) = -1.0 * InterpInPres ( divQ(ird,:) , lev , p_model )*1000.
        Q_V_adv(ird,:) = -1.0 * InterpInPres ( vertdivQ(ird,:) , lev , p_model )*1000.

               if (debug) then
                  write(*,*) nt, nlayr
                  write(*,*) p_model
                  write(*,*) tt
                  write(*,*) qq
                  write(*,*) oo
                  write(*,*) T_H_adv
                  write(*,*) Q_H_adv
                  write(*,*) T_V_adv
                  write(*,*) Q_V_adv
               endif

end SUBROUTINE CFMIP_IC

function  InterpInPres ( v , p, px ) result(vx)

! PX - target (GEOS5) midlevel pressures
! P  - dataset pressures

!!integer,  intent(in) :: lo,lx
real,    dimension(:), intent(in)  :: v,p
real,    dimension(:), intent(in)  :: px
!!real,    dimension(lx)       :: vx
real,    dimension(size(px))       :: vx

integer  :: l1,l2,l,lo,lx

logical  :: pflip

lo = size( p )
lx = size( px )

if (p(lo) > p(lo-1)) pflip = .false. ! top-down like GEOS-5
if (p(lo) < p(lo-1)) pflip = .true.  ! bottom up, reverse of GEOS-5


if (pflip) then

   where( px >= maxval(p))
      vx = v(1)
   endwhere

   do l=1,lo-1
      where( ( px >= p(l+1))  .and.( px < p(l) ) )
         vx = v(l+1)+(px-p(l+1))*(v(l)-v(l+1))/(p(l)-p(l+1))
      endwhere
   end do

   where( px < minval(p))
      vx = v(lo)
   endwhere

else

   where( px >= maxval(p))
      vx = v(lo)
   endwhere

   do l=1,lo-1
      where( ( px <= p(l+1))  .and.( px > p(l) ) )
         vx = v(l+1)+(px-p(l+1))*(v(l)-v(l+1))/(p(l)-p(l+1))
      endwhere
   end do

   where( px < minval(p))
      vx = v(1)
   endwhere

endif

end function InterpInPres

end module cfmip_data_mod


