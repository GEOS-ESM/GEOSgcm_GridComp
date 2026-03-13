  IF(itest==5 .and. draft == "mid_up" ) then

      !--- part 1 GAMMA format
      csum =0.
      zubeg=0.
      lev_start=min(.9,.1+csum*.013)
      kb_adj=max(kb,2)
      kb_adj=min(kb_adj,kt-1)
      if(kb_adj==kt) stop "kb_adj==kt"

      tunning = 0.30
      alpha2  = (tunning*(beta_deep -2.)+1.)/(1.-tunning)

      do k=27,3,-1
            if(x_alpha(k) >= alpha2)exit
      enddo
      k1=k+1
      if(x_alpha(k1) .ne. x_alpha(k1-1))then
        a=x_alpha(k1)-x_alpha(k1-1)
        b=x_alpha(k1-1)*(k1) -(k1-1)*x_alpha(k1)
        x1= (alpha2-b)/a
        y1=a*x1+b
        g_a=g_alpha(k1)-g_alpha(k1-1)
        g_b=g_alpha(k1-1)*k1 - (k1-1)*g_alpha(k1)
        g_alpha2=g_a*x1+g_b
      else
        g_alpha2=g_alpha(k1)
      endif

      fzu = gamma(alpha2 + beta_deep)/(g_alpha2*g_beta_deep)
      fzu=0.01*fzu
      do k=kb_adj,min(kte,kt)
         kratio= (po_cup(k)-po_cup(kb_adj))/(po_cup(kt)-po_cup(kb_adj))
         zu(k) = zubeg+FZU*kratio**(alpha2-1.0) * (1.0-kratio)**(beta_deep-1.0)

      enddo
      !- normalize ZU
      zu(kts:min(kte,kt+1))= zu(kts:min(kte,kt+1))/ (1.e-12+maxval(zu(kts:min(kte,kt+1))))

      !--- part 2: BETA format
      pmaxzu=psur-px*(psur-po_cup(kt))
      kb_adj=minloc(abs(po_cup(kts:kt)-pmaxzu),1)
      kb_adj=max(kb,kb_adj)
      kb_adj=min(kb_adj,kt)
      !beta=4.  !=> must be larger than 1
                !=> higher makes the profile sharper
                !=> around the maximum zu
      !- 2nd approach for beta and alpha parameters
      !- the tunning parameter must be between 0.5 (low  level max zu)
      !-                                   and 1.5 (high level max zu)
      !tunning= 1.0
      tunning = 0.6 !-- clo-X and tune0.6 experiment.
      !
      beta    = 2.0/tunning
      alpha   = tunning*beta
      !
      !- this alpha constrains the location of the maximun ZU to be at
      !- "kb_adj" vertical level
      alpha=1. + (beta-1.0)*(float(kb_adj)/float(kt+1))/(1.0-(float(kb_adj)/float(kt+1)))
      !
      ! imposing zu(ktop) = 0
      do k=klcl-1,min(kte,kt)
          kratio= float(k)/float(kt+1)
          zuh(k) = kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
      enddo
      !- normalize ZU
      zuh(kts:min(kte,kt+1))= zuh(kts:min(kte,kt+1))/ (1.e-12+maxval(zuh(kts:min(kte,kt+1))))

      !--- part 3: BETA format from sfc to max zuh, then GAMMA format
      Do k=kts,max(kts,maxloc(zuh(:),1)-2)
          zu(k)=zuh(k)
      ENDDO
      Do k=max(kts,maxloc(zuh(:),1)-1),min(maxloc(zuh(:),1)+1,kt)
          zu(k)=0.5*(zu(k)+zuh(k))
      ENDDO

      !-- special treatment below k22/klcl
      DO k=klcl,kts+1,-1
        zu(k)=zu(k+1)*0.5
      enddo
     !-- smooth section
     IF(do_smooth) then
          !--from surface
          zul(kts+1)=zu(kts+1)*0.25
          do k=kts+2,maxloc(zu,1)
             zul(k)=(zu(k-1)+zu(k))*0.5
          enddo
          do k=kts+1,maxloc(zu,1)
             zu(k)=(zul(k)+zu(k))*0.5
          enddo

          !--from ktop
          zul(kt-1)=zu(kt-1)*0.1
	  !print*,"ZUMD=",kt,zu(kt),zul(kt)
          do k=kt-2,kt-min(maxloc(zu,1),5),-1
             zul(k)=(zul(k+1)+zu(k))*0.5
          enddo
	  wgty=0.
          do k=kt,kt-min(maxloc(zu,1),5),-1
	     wgty=wgty+1./(float(min(maxloc(zu,1),5))+1)
	     zu(k)=zul(k)*(1.-wgty)+ zu(k)*wgty
             !print*,"zuMD=",k,zu(k),zul(k),(zul(k)+zu(k))*0.5,min(maxloc(zu,1),5),wgty
          enddo
      ENDIF
      zu(kts)=0.
  !---------------------------------------------------------

