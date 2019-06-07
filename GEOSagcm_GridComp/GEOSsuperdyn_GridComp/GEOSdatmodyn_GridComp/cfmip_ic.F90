module cfmip_data_mod

  use MAPL_Mod
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


