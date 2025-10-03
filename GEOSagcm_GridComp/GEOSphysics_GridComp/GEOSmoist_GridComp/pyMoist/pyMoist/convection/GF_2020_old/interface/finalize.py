from ndsl.dsl.gt4py import computation, PARALLEL, interval

def finalize(
    
):
    from __externals__ import FEED_3DMIDEL

    with computation(PARALLEL), interval(...):
       if FEED_3DMIDEL:
            # update GEOS-5 model state with the feedback from cumulus convection
            # to include the tendencies from the convection,  update the vars th1,q1,v1 and u1
          
       





IF(FEED_3DMODEL)THEN
      !-- update GEOS-5 model state with the feedback from cumulus convection
      !- to include the tendencies from the convection,  update the vars th1,q1,v1 and u1
      DO j=1,myp
        DO i=1,mxp
          IF(do_this_column(i,j) == 0) CYCLE
          !- conv precip rate: mm/s = kg m-2 s-1
          CNPCPRATE (i,j) =  CONPRR (i,j)
          IF(ITEST==0) CNPCPRATE(i,j) =  0.

          !--- sublimation/evaporation tendencies (kg/kg/s)
          DO k=1,mzp
          !--- sublimation/evaporation tendencies (kg/kg/s)
               REVSU (i,j,k) = REVSU_GF(flip(k),i,j)
          !--- preciptation fluxes (kg/kg/s)
               PRFIL (i,j,k) = PRFIL_GF(flip(k),i,j)
          ENDDO

        ENDDO
      ENDDO
!-----
      IF(USE_MOMENTUM_TRANSP > 0) THEN
        DO j=1,myp
          DO i=1,mxp
            IF(do_this_column(i,j) == 0) CYCLE
          ENDDO
        ENDDO
      ENDIF


      IF(APPLY_SUB_MP == 1) THEN
        DO j=1,myp
          DO i=1,mxp
            IF(do_this_column(i,j) == 0) CYCLE
          ENDDO
        ENDDO
      ENDIF

      IF(USE_TRACER_TRANSP==1) THEN

        DO j=1,myp
          DO i=1,mxp
            IF(do_this_column(i,j) == 0) CYCLE
            DO k=1,mzp
              !- update tracer mass mixing ratios
              DO ispc=1,mtp

                 CNV_Tracers(ispc)%Q(i,j,k) = CNV_Tracers(ispc)%Q(i,j,k) + DT_moist * SRC_CHEM(ispc,flip(k),i,j)

                 !-- final check for negative tracer mass mixing ratio
                 CNV_Tracers(ispc)%Q(i,j,k) = max(CNV_Tracers(ispc)%Q(i,j,k), mintracer)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

      ENDIF

      entr_dp = MAPL_UNDEF
      entr_md = MAPL_UNDEF
      entr_sh = MAPL_UNDEF

      DO IENS=1, maxiens
        if(icumulus_gf(IENS) == ON) then
           DO j=1,myp
            DO i=1,mxp
             if(ierr4d(i,j,IENS) .ne. 0) cycle
              DO k=mzp,flip(ktop4d(i,j,IENS))-1,-1

                !- Export entrainment rates used by GF
                if (IENS==DEEP) entr_dp(i,j,k) = entr5d(i,flip(k),j,IENS)
                if (IENS==MID ) entr_md(i,j,k) = entr5d(i,flip(k),j,IENS)
                if (IENS==SHAL) entr_sh(i,j,k) = entr5d(i,flip(k),j,IENS)

                !- special treatment for CNV_DQCDT: 'convective_condensate_source',  UNITS     = 'kg m-2 s-1',
                !- SRC_CI contains contributions from deep, shallow,... . So, do not need to be accumulated over  CNV_DQCDT
                !- note the SRC_CI has different array structure (k,i,j) _not_ (i,j,k)
                CNV_DQCDT(i,j,k)=  SRC_CI(flip(k),i,j) * DZ(i,j,k) * AIR_DEN(i,j,k) !units: kg[w]/(kg[air] s) * m * kg[air]/m^3 = kg[w]/(m^2 s)

                !
                !-'detraining_mass_flux',UNITS     = 'kg m-2 s-1',
                CNV_MFD (i,j,k) = CNV_MFD (i,j,k) + up_massdetr5d(i,flip(k),j,IENS)


                !-'cloud_base_mass_flux',    units = 'kg m-2 s-1',
                CNV_MF0 (i,j,k) = CNV_MF0 (i,j,k) + ( zup5d(i,flip(k),j,IENS) )

                !-convective mass flux [kg m^{-2} s^{-1}] - with downdrafts
           !!!  CNV_MFC (i,j,k) = CNV_MFC (i,j,k) + ( zup5d(i,flip(k),j,IENS) + &
           !!!                                        zdn5d(i,flip(k),j,IENS)*edt4d(i,j,IENS) )

                !---only updraft
                CNV_MFC (i,j,k) = CNV_MFC (i,j,k) + ( zup5d(i,flip(k),j,IENS) )

                ! deep convective total water flux. assumes .033 fractional area.
                qsatup = MAPL_EQsat(tup5d(i,flip(k),j,IENS),press(flip(k),i,j),dtqw) + clwup5d(i,flip(k),j,IENS)/0.033
                WQT_DC (i,j,k) = WQT_DC (i,j,k) + zup5d(i,flip(k),j,IENS)*(qsatup - rvap(flip(k),i,j))

                if(zup5d(i,flip(k),j,IENS) > 1.0e-6) then
                   !-'entrainment parameter',  UNITS     ='m-1',
                   ENTLAM  (i,j,k) =  ENTLAM (i,j,k) + (up_massentr5d(i,flip(k),j,IENS)/(DZ(i,j,k)*zup5d(i,flip(k),j,IENS)))

                   !-'updraft_vertical_velocity',            UNITS     = 'hPa s-1',
                   CNV_CVW (i,j,k) = -0.2 ! hPa/s =>  4 m/s
                endif

                !-'grid_mean_convective_condensate', UNITS     ='kg kg-1'
                CNV_QC  (i,j,k) = CNV_QC (i,j,k) + clwup5d(i,flip(k),j,IENS)
                !
                !
                !~ !----------------------------------------------------------------------------------------------------
                !- not using progno-cloud to calculate the precip from the convective column
                !- if CNV_PRC3 will be send to progno-cloud, set CNPCPRATE = zero
                !-'convective_precipitation_from_GF',UNITS     = 'kg m-2 s-1',
                !- JAN/17/2017 : the units above are wrong. The correct are kg[precip water]/kg[air]
                CNV_PRC3(i,j,k) = CNV_PRC3 (i,j,k) +                 (prup5d(i,flip(k),j,IENS) + &
                                                     edt4d(i,j,IENS)* prdn5d(i,flip(k),j,IENS) ) &
                                                     * DT_moist/(DZ(i,j,k)*AIR_DEN(i,j,k))

                !-'updraft_areal_fraction',
                if(zup5d(i,flip(k),j,IENS) > 1.0e-6) CNV_UPDF(i,j,k) = 0.033

             ENDDO
            ENDDO
           ENDDO
        endif
      ENDDO
  ENDIF


  !
  !--- cold pool/"convection tracer"
  IF(CONVECTION_TRACER==1) THEN
           DO j=1,myp
            DO i=1,mxp

              tau_cp=FRLAND(i,j)*tau_land_cp + (1.0-FRLAND(i,j))*tau_ocea_cp

              DO k=1,mzp

                !- sink term (exp decay 1h)
                snk_cnvtr =  dt_moist * abs(CNV_TR (i,j,k))/tau_cp

                !- downdraft convective mass flux [kg m^{-2} s^{-1}]
                ! iens =?
                !src_cnvtr =  edt4d(i,j,iens)*zdn5d(i,flip(k),j,iens)

                !- downdraft detrainment mass flux [kg m^{-2} s^{-1}]
                ! iens =?
                !src_cnvtr =  edt4d(i,j,iens)*dd_massdetr5d(i,flip(k),j,iens)


                !- source term
                !- downdraft detrainment of buoyancy [ J/kg s^{-1}]
                !- negative sign => source for updraft lifting
                src_cnvtr = - dt_moist * min(0.,SRC_BUOY(flip(k),i,j))

                !- 'continuity' equation = ADV + SRC - SINK
                CNV_TR (i,j,k) = CNV_TR (i,j,k)  + src_cnvtr  - snk_cnvtr

              ENDDO
            ENDDO
           ENDDO
          !print*,"buoy_exc2=",maxval(SRC_BUOY),minval(SRC_BUOY)
   ENDIF

!
  CNV_TOPP_DP = MAPL_UNDEF
  CNV_TOPP_MD = MAPL_UNDEF
  CNV_TOPP_SH = MAPL_UNDEF
  IF(maxval(icumulus_gf)>0) then
      DO IENS=1, maxiens
        if(icumulus_gf(IENS) == ON .and. IENS== DEEP) then
         DO j=1,myp
           DO i=1,mxp
            if(ierr4d(i,j,DEEP) /= 0) cycle
            CNV_TOPP_DP(i,j)     = press(ktop4d(i,j,DEEP),i,j)
            MFDP      (i,j)      =   xmb4d(i,j,DEEP)
            SIGMA_DEEP(i,j)      = sigma4d(i,j,DEEP)
            MUPDP     (i,j,1:mzp)=zup5d(i,flip(1):flip(mzp):-1,j,DEEP)
            MDNDP     (i,j,1:mzp)=zdn5d(i,flip(1):flip(mzp):-1,j,DEEP) * edt4d(i,j,IENS)
           ENDDO
         ENDDO
        elseif(icumulus_gf(IENS) == ON .and. IENS== SHAL) then
         DO j=1,myp
           DO i=1,mxp
            if(ierr4d(i,j,SHAL) /= 0) cycle
            CNV_TOPP_SH(i,j)=press(ktop4d(i,j,SHAL),i,j)
            MFSH (i,j)      =xmb4d(i,j,SHAL)
            MUPSH(i,j,1:mzp)=zup5d(i,flip(1):flip(mzp):-1,j,SHAL)
           ENDDO
         ENDDO
        elseif(icumulus_gf(IENS) == ON .and. IENS== MID) then
         DO j=1,myp
           DO i=1,mxp
             if(ierr4d(i,j,MID) /= 0) cycle
             CNV_TOPP_MD(i,j)     = press(ktop4d(i,j,MID),i,j)
             MFMD      (i,j)      = xmb4d(i,j,MID)
             SIGMA_MID (i,j)      = sigma4d(i,j,MID )
             MUPMD     (i,j,1:mzp)=zup5d(i,flip(1):flip(mzp):-1,j,MID)
           ENDDO
         ENDDO
        endif
      ENDDO
      !- for output purposes, ierr=0 (convection is ON) will be changed to 1
      where(ierr4d==0)ierr4d=1
      where(ierr4d>1)ierr4d=0
      DO j=1,myp
       DO i=1,mxp
        !- Tendencies
         DTDT_GF(i,j,1:mzp)=SRC_T(flip(1):flip(mzp):-1,i,j)
         DQDT_GF(i,j,1:mzp)=SRC_Q(flip(1):flip(mzp):-1,i,j)
         DUDT_GF(i,j,1:mzp)=SRC_U(flip(1):flip(mzp):-1,i,j)
         DVDT_GF(i,j,1:mzp)=SRC_V(flip(1):flip(mzp):-1,i,j)
        !- Error codes
         ERRDP(i,j)=float(ierr4d(i,j,DEEP))
         ERRSH(i,j)=float(ierr4d(i,j,SHAL))
         ERRMD(i,j)=float(ierr4d(i,j,MID ))
       ENDDO
      ENDDO
  ENDIF
  !
  if( allocated(src_chem))  deallocate(src_chem,stat=alloc_stat) !tendency   from convection

  !- for debugging purposes only
  if(wrtgrads) call alloc_grads_arr(1,mzp,2,jl)