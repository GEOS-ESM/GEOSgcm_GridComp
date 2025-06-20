;It sounds like ATAU and BTAU don’t differ?  If not, great!  Those (like tsa1, tsa2, tsb1, tsb2) would be hard to check…
 
;Based on the Catchment code, the following fitted functions have to
;be within tolerance of each other (say, within 1%),
;even if the different parameter values themselves differ (i.e., because different compilers find different combinations of parameters that work):
 
;(1)        AR1W= AMIN1(1.,AMAX1(0.,(1.+ars1(n)*CATDEFW)                       &
;                 /(1.+ars2(n)*CATDEFW+ars3(n)*CATDEFW*CATDEFW)))     ! for 0<catdef<cdcr1
; 
;(2)        WMIN=AMIN1(1.,AMAX1(0.,arw4(n)+(1.-arw4(n))*(1.+arw1(n)*CATDEFX)       &
;                 /(1.+arw2(n)*CATDEFX+arw3(n)*CATDEFX*CATDEFX)))       ! for 0<catdef<cdcr1
; 
;(3)       if (ara1(n) .ne. ara3(n)) then
;            cdi=(ara4(n)-ara2(n))/(ara1(n)-ara3(n))
;          else
;            cdi=0.
;          endif
; 
;(4)         ZBAR=SQRT(1.e-20+catdef(n)/bf1(n))-bf2(n)              ! for 0<catdef<cdcr1
; 
;(5)      bf3
; 
;(6)       rzx=rzexc(n)/vgwmax(n)
;           if(rzx .gt. .01) then
;            ax=tsa1(n)
;            bx=tsb1(n)
;         elseif(rzx .lt. -.01) then
;            ax=tsa2(n)
;            bx=tsb2(n)
;         else
;            ax=tsa2(n)+(rzx+.01)*(tsa1(n)-tsa2(n))/.02
;            bx=tsb2(n)+(rzx+.01)*(tsb1(n)-tsb2(n))/.02
;         endif
;         tsc2=exp(ax+bx*catdef(n))   ! for 0<catdef<cdcr2  and -20<rzexc<20
 
; compute cdcr1, cdcr2

FUNCTION asat, ars1_1,ars2_1,ars3_1,ars1_2,ars2_2,ars3_2, cdcr1

  max_err = 0.
  for i = 0, fix(cdcr1),2 do begin
     X = 1. *i
     AR1W_1= MIN([1.,MAX([0.,(1.+ ars1_1*X)/(1.+ ars2_1*X + ars3_1*X*X)])])
     AR1W_2= MIN([1.,MAX([0.,(1.+ ars1_2*X)/(1.+ ars2_2*X + ars3_2*X*X)])])
     max_err = max ([max_err, abs(AR1W_1 - AR1W_2)])
  endfor
  return, max_err
  
END

FUNCTION wmin, arw1_1, arw2_1, arw3_1, arw4_1,arw1_2, arw2_2, arw3_2, arw4_2, cdcr1

  max_err = 0.
  for i = 0, fix(cdcr1),2 do begin
     X = 1. *i
     wmin_1= MIN([1.,MAX([0.,arw4_1 + (1.+arw1_1*X)/(1.+ arw2_1*X + arw3_1*X*X)])])
     wmin_2= MIN([1.,MAX([0.,arw4_2 + (1.+arw1_2*X)/(1.+ arw2_2*X + arw3_2*X*X)])])
     max_err = max ([max_err, abs(wmin_1 - wmin_2)])
  endfor
  return, max_err
  
END

FUNCTION  zbar,bf1_1, bf2_1, bf1_2, bf2_2, cdcr2
  
  max_err = 0.
  for i = 0, fix(cdcr2),2 do begin
     X = 1. *i
     zbar1 = (1.e-20 + X/bf1_1)^(0.5) - bf2_1
     zbar2 = (1.e-20 + X/bf1_2)^(0.5) - bf2_2
     max_err = max ([max_err, abs(zbar1 - zbar2)])
  endfor
  return, max_err
  
END


FUNCTION TSC2, TSA1_1, TSB1_1, TSA2_1, TSB2_1, TSA1_2, TSB1_2, TSA2_2, TSB2_2, VGWMAX, CDCR2

  max_err = 0.
  for j = -20, 20, 5 do begin
     R = 1.*j
     rzx = R/VGWMAX
     if (rzx gt .01) then begin
        ax_1 = TSA1_1
        bx_1 = TSB1_1
        ax_2 = TSA1_2
        bx_2 = TSB1_2        
     endif else if (rzx lt -0.01) then begin
        ax_1 = TSA2_1
        bx_1 = TSB2_1
        ax_2 = TSA2_2
        bx_2 = TSB2_2                
     endif else if ((rzx ge -0.01) and (rzx le 0.01)) then begin
        ax_1 = tsa2_1 +(rzx+.01)*(tsa1_1-tsa2_1)/.02
        bx_1 = tsb2_1 +(rzx+.01)*(tsb1_1-tsb2_1)/.02        
        ax_2 = tsa2_2 +(rzx+.01)*(tsa1_2-tsa2_2)/.02
        bx_2 = tsb2_2 +(rzx+.01)*(tsb1_2-tsb2_2)/.02                
     endif
     
     for i = 0,fix(cdcr2),10 do begin
        X = 1.*i
        tsc2_1=exp(ax_1 + bx_1 * X)
        tsc2_2=exp(ax_2 + bx_2 * X)
        max_err = max ([max_err, abs(tsc2_1 - tsc2_2)])
     endfor
  endfor
  return, max_err
END

FUNCTION cdi, ara1_1, ara2_1, ara3_1, ara4_1, ara1_2, ara2_2, ara3_2, ara4_2

  if (ara1_1 ne ara3_1) then begin
     cdi_1 = (ara4_1 - ara2_1)/(ara1_1 - ara3_1)
  endif else begin
     cdi_1 = 0.
  endelse
  
  if (ara1_2 ne ara3_2) then begin
     cdi_2 = (ara4_2 - ara2_2)/(ara1_2 - ara3_2)
  endif else begin
     cdi_2 = 0.
  endelse

  return, abs (cdi_1 - cdi_2)
    
END
;

pro prob_analy
  
  tol = 1000000l
  path = '/discover/nobackup/smahanam/bcs/comp_zerodiff/'
  file1 = 'NLv4_SLES11/SMAP_EASEv2_M09/clsm/catch_params.nc4'
  file2 = 'NLv4_SLES12/SMAP_EASEv2_M09/clsm/catch_params.nc4'
  file3 = 'NLv4_SLES12_i18_dble/SMAP_EASEv2_M09/clsm/catch_params.nc4'
  file4 = 'NLv4_SLES12_i19_dble/SMAP_EASEv2_M09/clsm/catch_params.nc4'

  temp = read_ascii(path + 'NLv4_SLES11/SMAP_EASEv2_M09/clsm/bad_sat_param.tiles')
  btiles = long (temp.field1(0,*))
  btiles= reform(btiles, [149],/overwrite)
  temp = read_ascii(path + 'NLv4_SLES12/SMAP_EASEv2_M09/clsm/bad_sat_param.tiles')
  btiles2 = long (temp.field1(0,*))
  btiles2= reform(btiles2, [68],/overwrite)
  btiles = [btiles,  btiles2]
 
  a = sort(btiles) 
  btiles = btiles(a)
  btiles = btiles - 1
  ncid1 = NCDF_OPEN(path + file1,/NOWRITE)
  ncid2 = NCDF_OPEN(path + file2,/NOWRITE)
  ncid3 = NCDF_OPEN(path + file3,/NOWRITE)
  ncid4 = NCDF_OPEN(path + file4,/NOWRITE)

  ; compute CDCR1, CDCR2, VGWMAX
  NCDF_VARGET, ncid1, 'POROS', POROS
  NCDF_VARGET, ncid1, 'DP2BR', DP2BR
  NCDF_VARGET, ncid1, 'BEE'  , BEE
  NCDF_VARGET, ncid1, 'PSIS' , PSIS
  NCDF_VARGET, ncid1, 'WPWET', WPWET

  VGWMAX = POROS * 1000.
  NTILES = N_ELEMENTS (POROS)
  CDCR1 = fltarr (NTILES)
  CDCR2 = fltarr (NTILES)
  
  for n = 0l, NTILES -1l do begin
     
     CDCR2(n)  = (1.-WPWET(n)) * POROS(n) * DP2BR(n)
     term1     = -1.+((PSIS(n)-DP2BR(n)/1000.)/PSIS(n))^((BEE(n)-1.)/BEE(n))
     term2     = PSIS(n)*BEE(n)/(BEE(n)-1)
     CDCR1(n)  = 1000.*POROS(n)*(DP2BR(n)/1000.-(-term2*term1))   
        
  endfor
  goto, plot_err
  
  ; ar.new
  ; ======

  ; Asat
  NCDF_VARGET, ncid1, 'ARS1', ARS1_1
  NCDF_VARGET, ncid1, 'ARS2', ARS2_1
  NCDF_VARGET, ncid1, 'ARS3', ARS3_1
  NCDF_VARGET, ncid2, 'ARS1', ARS1_2
  NCDF_VARGET, ncid2, 'ARS2', ARS2_2
  NCDF_VARGET, ncid2, 'ARS3', ARS3_2
  NCDF_VARGET, ncid3, 'ARS1', ARS1_3
  NCDF_VARGET, ncid3, 'ARS2', ARS2_3
  NCDF_VARGET, ncid3, 'ARS3', ARS3_3
  NCDF_VARGET, ncid4, 'ARS1', ARS1_4
  NCDF_VARGET, ncid4, 'ARS2', ARS2_4
  NCDF_VARGET, ncid4, 'ARS3', ARS3_4

  ; shave off < 1.e-6
  ARS1_1 = round (ARS1_1 * tol) / double(tol)
  ARS2_1 = round (ARS2_1 * tol) / double(tol)
  ARS3_1 = round (ARS3_1 * tol) / double(tol)
  ARS1_2 = round (ARS1_2 * tol) / double(tol)
  ARS2_2 = round (ARS2_2 * tol) / double(tol)
  ARS3_2 = round (ARS3_2 * tol) / double(tol)
  ARS1_3 = round (ARS1_3 * tol) / double(tol)
  ARS2_3 = round (ARS2_3 * tol) / double(tol)
  ARS3_3 = round (ARS3_3 * tol) / double(tol)
  ARS1_4 = round (ARS1_4 * tol) / double(tol)
  ARS2_4 = round (ARS2_4 * tol) / double(tol)
  ARS3_4 = round (ARS3_4 * tol) / double(tol)

  ARS1_1 (btiles) = 0.
  ARS2_1 (btiles) = 0.
  ARS3_1 (btiles) = 0.
  ARS1_2 (btiles) = 0.
  ARS2_2 (btiles) = 0.
  ARS3_2 (btiles) = 0.
  
  ;
  asat_err1 = 0.
  asat_err2 = 0.
  asat_loc  = 1
  for n = 0l, NTILES -1l do begin
     diff = (ARS1_1(n) + ARS2_1(n) + ARS3_1(n)) - (ARS1_2(n) + ARS2_2(n) + ARS3_2(n))
     if (diff ne 0.) then begin
        asat_err1 = [asat_err1, $
                     asat (ARS1_1(n), ARS2_1(n), ARS3_1(n),ARS1_2(n), ARS2_2(n), ARS3_2(n), cdcr1(n))]
        asat_loc = [asat_loc, n + 1]
     endif
        diff = (ARS1_3(n) + ARS2_3(n) + ARS3_3(n)) - (ARS1_4(n) + ARS2_4(n) + ARS3_4(n))
        
     if (diff ne 0.) then $
        asat_err2 = [asat_err2, $
                     asat (ARS1_3(n), ARS2_3(n), ARS3_3(n),ARS1_4(n), ARS2_4(n), ARS3_4(n), cdcr1(n))]
     
  endfor
  
  asat_err1 =  asat_err1(1:*)
  asat_err2 =  asat_err2(1:*)
  asat_loc  =  asat_loc (1:*)

  print , "DONE ASAT", n_elements (asat_err1), n_elements (asat_err2)
  
  ARS1_1 = 0.
  ARS2_1 = 0.
  ARS3_1 = 0.
  ARS1_2 = 0.
  ARS2_2 = 0.
  ARS3_2 = 0.
  ARS1_3 = 0.
  ARS2_3 = 0.
  ARS3_3 = 0.
  ARS1_4 = 0.
  ARS2_4 = 0.
  ARS3_4 = 0.

  ; WMIN
  NCDF_VARGET, ncid1, 'ARW1',ARW1_1 
  NCDF_VARGET, ncid1, 'ARW2',ARW2_1 
  NCDF_VARGET, ncid1, 'ARW3',ARW3_1
  NCDF_VARGET, ncid1, 'ARW4',ARW4_1
  NCDF_VARGET, ncid2, 'ARW1',ARW1_2 
  NCDF_VARGET, ncid2, 'ARW2',ARW2_2 
  NCDF_VARGET, ncid2, 'ARW3',ARW3_2
  NCDF_VARGET, ncid2, 'ARW4',ARW4_2
  NCDF_VARGET, ncid3, 'ARW1',ARW1_3 
  NCDF_VARGET, ncid3, 'ARW2',ARW2_3 
  NCDF_VARGET, ncid3, 'ARW3',ARW3_3
  NCDF_VARGET, ncid3, 'ARW4',ARW4_3
  NCDF_VARGET, ncid4, 'ARW1',ARW1_4 
  NCDF_VARGET, ncid4, 'ARW2',ARW2_4 
  NCDF_VARGET, ncid4, 'ARW3',ARW3_4
  NCDF_VARGET, ncid4, 'ARW4',ARW4_4

  ARW1_1 = round (ARW1_1 * tol) / double(tol) 
  ARW2_1 = round (ARW2_1 * tol) / double(tol)
  ARW3_1 = round (ARW3_1 * tol) / double(tol)
  ARW4_1 = round (ARW4_1 * tol) / double(tol)
  ARW1_2 = round (ARW1_2 * tol) / double(tol)
  ARW2_2 = round (ARW2_2 * tol) / double(tol)
  ARW3_2 = round (ARW3_2 * tol) / double(tol)
  ARW4_2 = round (ARW4_2 * tol) / double(tol)
  ARW1_3 = round (ARW1_3 * tol) / double(tol)
  ARW2_3 = round (ARW2_3 * tol) / double(tol)
  ARW3_3 = round (ARW3_3 * tol) / double(tol)
  ARW4_3 = round (ARW4_3 * tol) / double(tol)
  ARW1_4 = round (ARW1_4 * tol) / double(tol)
  ARW2_4 = round (ARW2_4 * tol) / double(tol)
  ARW3_4 = round (ARW3_4 * tol) / double(tol)
  ARW4_4 = round (ARW4_4 * tol) / double(tol)

  ARW1_1(btiles) = 0.
  ARW2_1(btiles) = 0.
  ARW3_1(btiles) = 0.
  ARW4_1(btiles) = 0.
  ARW1_2(btiles) = 0.
  ARW2_2(btiles) = 0.
  ARW3_2(btiles) = 0.
  ARW4_2(btiles) = 0.

  wmin_err1 = 0.
  wmin_err2 = 0.
  wmin_loc  = 0
  
  for n = 0l, NTILES -1l do begin
     diff = (ARW1_1(n) + ARW2_1(n) + ARW3_1(n) + ARW4_1(n)) - (ARW1_2(n) + ARW2_2(n) + ARW3_2(n) + ARW4_2(n))
     if (diff ne 0.) then begin
        wmin_err1 = [wmin_err1, $
                     wmin(ARW1_1(n), ARW2_1(n), ARW3_1(n), ARW4_1(n), ARW1_2(n), ARW2_2(n), ARW3_2(n), ARW4_2(n), cdcr1(n))]
        wmin_loc = [wmin_loc, n + 1]
     endif
     diff = (ARW1_3(n) + ARW2_3(n) + ARW3_3(n) + ARW4_3(n)) - (ARW1_4(n) + ARW2_4(n) + ARW3_4(n) + ARW4_4(n))
     if (diff ne 0.) then $
        wmin_err2 = [wmin_err2, $
                     wmin(ARW1_3(n), ARW2_3(n), ARW3_3(n), ARW4_3(n), ARW1_4(n), ARW2_4(n), ARW3_4(n), ARW4_4(n), cdcr1(n))]     
  endfor

  wmin_err1 = wmin_err1(1:*)
  wmin_err2 = wmin_err2(1:*)
  wmin_loc  = wmin_loc (1:*)
  
    PRINT ,'DONE WMIN', n_elements (wmin_err1), n_elements (wmin_err1)
  
  ARW1_1 = 0. 
  ARW2_1 = 0.
  ARW3_1 = 0.
  ARW4_1 = 0.
  ARW1_2 = 0.
  ARW2_2 = 0.
  ARW3_2 = 0.
  ARW4_2 = 0.
  ARW1_3 = 0.
  ARW2_3 = 0.
  ARW3_3 = 0.
  ARW4_3 = 0.
  ARW1_4 = 0.
  ARW2_4 = 0.
  ARW3_4 = 0.
  ARW4_4 = 0.

  ; shape
  NCDF_VARGET, ncid1, 'ARA1',ARA1_1 
  NCDF_VARGET, ncid1, 'ARA2',ARA2_1 
  NCDF_VARGET, ncid1, 'ARA3',ARA3_1
  NCDF_VARGET, ncid1, 'ARA4',ARA4_1
  NCDF_VARGET, ncid2, 'ARA1',ARA1_2 
  NCDF_VARGET, ncid2, 'ARA2',ARA2_2 
  NCDF_VARGET, ncid2, 'ARA3',ARA3_2
  NCDF_VARGET, ncid2, 'ARA4',ARA4_2
  NCDF_VARGET, ncid3, 'ARA1',ARA1_3 
  NCDF_VARGET, ncid3, 'ARA2',ARA2_3 
  NCDF_VARGET, ncid3, 'ARA3',ARA3_3
  NCDF_VARGET, ncid3, 'ARA4',ARA4_3
  NCDF_VARGET, ncid4, 'ARA1',ARA1_4 
  NCDF_VARGET, ncid4, 'ARA2',ARA2_4 
  NCDF_VARGET, ncid4, 'ARA3',ARA3_4
  NCDF_VARGET, ncid4, 'ARA4',ARA4_4

  ARA1_1 = round (ARA1_1 * tol) / double(tol) 
  ARA2_1 = round (ARA2_1 * tol) / double(tol)
  ARA3_1 = round (ARA3_1 * tol) / double(tol)
  ARA4_1 = round (ARA4_1 * tol) / double(tol)
  ARA1_2 = round (ARA1_2 * tol) / double(tol)
  ARA2_2 = round (ARA2_2 * tol) / double(tol)
  ARA3_2 = round (ARA3_2 * tol) / double(tol)
  ARA4_2 = round (ARA4_2 * tol) / double(tol)
  ARA1_3 = round (ARA1_3 * tol) / double(tol)
  ARA2_3 = round (ARA2_3 * tol) / double(tol)
  ARA3_3 = round (ARA3_3 * tol) / double(tol)
  ARA4_3 = round (ARA4_3 * tol) / double(tol)
  ARA1_4 = round (ARA1_4 * tol) / double(tol)
  ARA2_4 = round (ARA2_4 * tol) / double(tol)
  ARA3_4 = round (ARA3_4 * tol) / double(tol)
  ARA4_4 = round (ARA4_4 * tol) / double(tol)


   ARA1_1 (btiles) = 0.
   ARA2_1 (btiles) = 0.
   ARA3_1 (btiles) = 0.
   ARA4_1 (btiles) = 0.
   ARA1_2 (btiles) = 0.
   ARA2_2 (btiles) = 0.
   ARA3_2 (btiles) = 0.
   ARA4_2 (btiles) = 0.
  
  ara_err1 = 0.
  ara_err2 = 0.
  ara_loc  = 1
  
  for n = 0l, NTILES -1l do begin
     diff = (ARA1_1(n) + ARA2_1(n) + ARA3_1(n) + ARA4_1(n)) - (ARA1_2(n) + ARA2_2(n) + ARA3_2(n) + ARA4_2(n))
     if (diff ne 0.) then begin
        ara_err1 = [ara_err1, $
                    cdi(ARA1_1(n), ARA2_1(n), ARA3_1(n), ARA4_1(n), ARA1_2(n), ARA2_2(n), ARA3_2(n), ARA4_2(n))]
        ara_loc = [ara_loc, n + 1]
     endif
     diff = (ARA1_3(n) + ARA2_3(n) + ARA3_3(n) + ARA4_3(n)) - (ARA1_4(n) + ARA2_4(n) + ARA3_4(n) + ARA4_4(n))
     if (diff ne 0.) then $
        ara_err2 = [ara_err2, $
                     cdi(ARA1_3(n), ARA2_3(n), ARA3_3(n), ARA4_3(n), ARA1_4(n), ARA2_4(n), ARA3_4(n), ARA4_4(n))]     
  endfor
  
  ara_err1 = ara_err1(1:*)
  ara_err2 = ara_err2(1:*)
  ara_loc  = ara_loc (1:*)
  
  print , 'DONE ARA', n_elements (ara_err1), n_elements (ara_err2)
  
  ARA1_1 = 0. 
  ARA2_1 = 0.
  ARA3_1 = 0.
  ARA4_1 = 0.
  ARA1_2 = 0.
  ARA2_2 = 0.
  ARA3_2 = 0.
  ARA4_2 = 0.
  ARA1_3 = 0.
  ARA2_3 = 0.
  ARA3_3 = 0.
  ARA4_3 = 0.
  ARA1_4 = 0.
  ARA2_4 = 0.
  ARA3_4 = 0.
  ARA4_4 = 0.  
    
  ; bf.dat
  ; ======

  NCDF_VARGET, ncid1, 'BF1', BF1_1
  NCDF_VARGET, ncid1, 'BF2', BF2_1
  NCDF_VARGET, ncid1, 'BF3', BF3_1
  NCDF_VARGET, ncid2, 'BF1', BF1_2
  NCDF_VARGET, ncid2, 'BF2', BF2_2
  NCDF_VARGET, ncid2, 'BF3', BF3_2
  NCDF_VARGET, ncid3, 'BF1', BF1_3
  NCDF_VARGET, ncid3, 'BF2', BF2_3
  NCDF_VARGET, ncid3, 'BF3', BF3_3
  NCDF_VARGET, ncid4, 'BF1', BF1_4
  NCDF_VARGET, ncid4, 'BF2', BF2_4
  NCDF_VARGET, ncid4, 'BF3', BF3_4
  
 ; shave off < 1.e-6
  BF1_1 = round (BF1_1 * tol) / double(tol)
  BF2_1 = round (BF2_1 * tol) / double(tol)
  BF3_1 = round (BF3_1 * tol) / double(tol)
  BF1_2 = round (BF1_2 * tol) / double(tol)
  BF2_2 = round (BF2_2 * tol) / double(tol)
  BF3_2 = round (BF3_2 * tol) / double(tol)
  BF1_3 = round (BF1_3 * tol) / double(tol)
  BF2_3 = round (BF2_3 * tol) / double(tol)
  BF3_3 = round (BF3_3 * tol) / double(tol)
  BF1_4 = round (BF1_4 * tol) / double(tol)
  BF2_4 = round (BF2_4 * tol) / double(tol)
  BF3_4 = round (BF3_4 * tol) / double(tol)


  BF1_1(btiles) = 0.
  BF2_1(btiles) = 0.
  BF3_1(btiles) = 0.
  BF1_2(btiles) = 0.
  BF2_2(btiles) = 0.
  BF3_2(btiles) = 0.
  
  ; zbar
  zbar_err1 = 0.
  zbar_err2 = 0.
  topo_err1 = 0.
  topo_err2 = 0.
  zbar_loc  = 1
  topo_loc  = 1
  
  for n = 0l, NTILES -1l do begin

     diff = (BF1_1(n) + BF2_1(n)) - (BF1_2(n) + BF2_2(n))
     if (diff ne 0.) then begin
        zbar_err1 = [zbar_err1, $
                     zbar(BF1_1(n),BF2_1(n), BF1_2(n), BF2_2(n), CDCR2(n))]
        zbar_loc = [zbar_loc, n + 1]

     endif
     diff = (BF1_3(n) + BF2_3(n)) - (BF1_4(n) + BF2_4(n))
     if (diff ne 0.) then $
        zbar_err2 = [zbar_err2, $
                     zbar(BF1_3(n),BF2_3(n), BF1_4(n), BF2_4(n), CDCR2(n))]
     if (BF3_1(n) - BF3_2(n)) then begin
        topo_err1 = [topo_err1, abs (BF3_1(n) - BF3_2(n))]
        topo_loc = [topo_loc, n + 1]
     endif
     if (BF3_3(n) - BF3_4(n)) then $
        topo_err2 = [topo_err2, abs (BF3_3(n) - BF3_4(n))]
     
  endfor

  zbar_err1 = zbar_err1(1:*)
  zbar_err2 = zbar_err2(1:*)
  topo_err1 = topo_err1(1:*)
  topo_err2 = topo_err2(1:*)
  zbar_loc  = zbar_loc (1:*)
  topo_loc  = topo_loc (1:*)
  
  print , 'DONE ZBAR', n_elements (zbar_err1), n_elements (zbar_err2)
  BF1_1  = 0.
  BF2_1  = 0.
  BF3_1  = 0.
  BF1_2  = 0.
  BF2_2  = 0.
  BF3_2  = 0.
  BF1_3  = 0.
  BF2_3  = 0.
  BF3_3  = 0.
  BF1_4  = 0.
  BF2_4  = 0.
  BF3_4  = 0.
  
  ; ts.dat
  ; ======

  ;goto, jump
  
  NCDF_VARGET, ncid1, 'TSA1',TSA1_1 
  NCDF_VARGET, ncid1, 'TSB1',TSB1_1 
  NCDF_VARGET, ncid1, 'TSA2',TSA2_1
  NCDF_VARGET, ncid1, 'TSB2',TSB2_1
  NCDF_VARGET, ncid2, 'TSA1',TSA1_2 
  NCDF_VARGET, ncid2, 'TSB1',TSB1_2 
  NCDF_VARGET, ncid2, 'TSA2',TSA2_2
  NCDF_VARGET, ncid2, 'TSB2',TSB2_2
  NCDF_VARGET, ncid3, 'TSA1',TSA1_3 
  NCDF_VARGET, ncid3, 'TSB1',TSB1_3 
  NCDF_VARGET, ncid3, 'TSA2',TSA2_3
  NCDF_VARGET, ncid3, 'TSB2',TSB2_3
  NCDF_VARGET, ncid4, 'TSA1',TSA1_4 
  NCDF_VARGET, ncid4, 'TSB1',TSB1_4 
  NCDF_VARGET, ncid4, 'TSA2',TSA2_4
  NCDF_VARGET, ncid4, 'TSB2',TSB2_4
  
  TSA1_1 = round (TSA1_1 * tol) / double(tol) 
  TSB1_1 = round (TSB1_1 * tol) / double(tol)
  TSA2_1 = round (TSA2_1 * tol) / double(tol)
  TSB2_1 = round (TSB2_1 * tol) / double(tol)
  TSA1_2 = round (TSA1_2 * tol) / double(tol)
  TSB1_2 = round (TSB1_2 * tol) / double(tol)
  TSA2_2 = round (TSA2_2 * tol) / double(tol)
  TSB2_2 = round (TSB2_2 * tol) / double(tol)
  TSA1_3 = round (TSA1_3 * tol) / double(tol)
  TSB1_3 = round (TSB1_3 * tol) / double(tol)
  TSA2_3 = round (TSA2_3 * tol) / double(tol)
  TSB2_3 = round (TSB2_3 * tol) / double(tol)
  TSA1_4 = round (TSA1_4 * tol) / double(tol)
  TSB1_4 = round (TSB1_4 * tol) / double(tol)
  TSA2_4 = round (TSA2_4 * tol) / double(tol)
  TSB2_4 = round (TSB2_4 * tol) / double(tol)  

  TSA1_1(btiles) = 0.
  TSB1_1(btiles) = 0.
  TSA2_1(btiles) = 0.
  TSB2_1(btiles) = 0.
  TSA1_2(btiles) = 0.
  TSB1_2(btiles) = 0.
  TSA2_2(btiles) = 0.
  TSB2_2(btiles) = 0.

  tsc_err1 = 0.
  tsc_err2 = 0.
  tsc_loc  = 1
  
  for n = 0l, NTILES -1l do begin

     diff = (TSA1_1(n) + TSB1_1(n) + TSA2_1(n) + TSB2_1(n)) - (TSA1_2(n) + TSB1_2(n) + TSA2_2(n) + TSB2_2(n))
     if (diff ne 0.) then begin
        tsc_err1 = [tsc_err1,TSC2(TSA1_1(n), TSB1_1(n), TSA2_1(n), TSB2_1(n), TSA1_2(n), TSB1_2(n), TSA2_2(n), TSB2_2(n), VGWMAX(n), CDCR2(n))]
        tsc_loc = [tsc_loc, n + 1]
     endif
        diff = (TSA1_3(n) + TSB1_3(n) + TSA2_3(n) + TSB2_3(n)) - (TSA1_4(n) + TSB1_4(n) + TSA2_4(n) + TSB2_4(n))
        
     if (diff ne 0.) then $
        tsc_err2 = [tsc_err2,TSC2(TSA1_3(n), TSB1_3(n), TSA2_3(n), TSB2_3(n), TSA1_4(n), TSB1_4(n), TSA2_4(n), TSB2_4(n), VGWMAX(n), CDCR2(n))]     
     
  endfor
  
  tsc_err1 = tsc_err1(1:*) 
  tsc_err2 = tsc_err2(1:*)
  tsc_loc  = tsc_loc (1:*)
  jump:
  
  save, asat_err1, asat_err2, wmin_err1, wmin_err2, ara_err1, ara_err2,  zbar_err1,zbar_err2, topo_err1, topo_err2, tsc_err1, tsc_err2,asat_loc, wmin_loc, ara_loc, zbar_loc, topo_loc, tsc_loc, file = 'error_file_excl_bad'
  stop
  plot_err:

  restore, 'error_file_excl_bad'
  
  load_colors
  sp,1,1,/f,/h,/color,/land
  !P.Multi = [0, 1, 2, 0, 0]
  !p.background = 255
  !P.position   = 0
  !p.charsize   = 1.2
  Erase,255

  ; ASAT = 0.1; WMIN 1; ARA = 200.,
  ; zbar = 0.25, topo = 0.1, tsc2 = 0.08
  div = 50.
  xmax = 0.1                              ; ASAT
  histo1 = histogram(asat_err1,min = 0., max = xmax, nbins=100, location =loc)
  histo2 = histogram(asat_err1,min = 0., max = xmax/div, nbins=100, location =loc2)
  ymax = max ([histo1, histo2])
  binz = xmax / (99)

  plot,indgen(100),replicate(0.,100),xrange = [0.,xmax],yrange =[0., ymax], color =0, charsize = 1.,XSTYLE=1,YSTYLE=1,title = 'Distribution of Max (abs(difference)) of A!Dsat!N (NLv4 vs NLv4p) ' + strtrim(string (n_elements(asat_err1)),2) + ' tiles affected', xtitle = 'Max (Abs (Diff)) of A!Dsat!N for CATDEF [0,CDCR1]', ytitle = 'Frequency'

  for i = 0, 98 do begin
     xbox = [loc(i), loc(i), loc(i+1), loc(i+1), loc(i)]
     ybox = [0, histo1(i), histo1(i), 0, 0]
     polyfill, xbox,ybox, color = 0
  endfor
  loc = loc2
  
  plot,indgen(100),replicate(0.,100),xrange = [0.,xmax/div],yrange =[0., ymax], color =0, charsize = 1.,XSTYLE=1,YSTYLE=1, title = 'Distribution of Max (abs(difference)) of A!Dsat!N (NLv4 vs NLv4p) ' + strtrim(string (n_elements(asat_err1)),2) + ' tiles affected', xtitle = 'Max (Abs (Diff)) of A!Dsat!N for CATDEF [0,CDCR1]', ytitle = 'Frequency'

  for i = 0, 98 do begin
     xbox = [loc(i), loc(i), loc(i+1), loc(i+1), loc(i)]
     ybox = [0, histo2(i), histo2(i), 0, 0]
     polyfill, xbox,ybox, color = 0
  endfor

  ; ############################################################################
  div = 100.
  xmax = 1.                              ; WMIN
  histo1 = histogram(wmin_err1,min = 0., max = xmax, nbins=100, location =loc)
  histo2 = histogram(wmin_err1,min = 0., max = xmax/div, nbins=100, location =loc2)
  ymax = max ([histo1, histo2])
  binz = xmax / (99)

  plot,indgen(100),replicate(0.,100),xrange = [0.,xmax],yrange =[0., ymax], color =0, charsize = 1.,XSTYLE=1,YSTYLE=1,title = 'Distribution of Max (abs(difference)) of !9q!X!D0!N (NLv4 vs NLv4p) ' + strtrim(string (n_elements(wmin_err1)),2) + ' tiles affected', xtitle = 'Max (Abs (Diff)) of !9q!X for CATDEF [0,CDCR1]', ytitle = 'Frequency'

  for i = 0, 98 do begin
     xbox = [loc(i), loc(i), loc(i+1), loc(i+1), loc(i)]
     ybox = [0, histo1(i), histo1(i), 0, 0]
     polyfill, xbox,ybox, color = 0
  endfor

  loc = loc2
  plot,indgen(100),replicate(0.,100),xrange = [0.,xmax/div],yrange =[0., ymax], color =0, charsize = 1.,XSTYLE=1,YSTYLE=1, title = 'Distribution of Max (abs(difference)) of !9q!X!D0!N (NLv4 vs NLv4p) ' + strtrim(string (n_elements(wmin_err1)),2) + ' tiles affected', xtitle = 'Max (Abs (Diff)) of !9q!X for CATDEF [0,CDCR1]', ytitle = 'Frequency'

  for i = 0, 98 do begin
     xbox = [loc(i), loc(i), loc(i+1), loc(i+1), loc(i)]
     ybox = [0, histo2(i), histo2(i), 0, 0]
     polyfill, xbox,ybox, color = 0
  endfor
  
  ; ############################################################################
  div = 1000.
  xmax = 200.                              ; ARA
  histo1 = histogram(ara_err1,min = 0., max = xmax, nbins=100, location =loc)
  histo2 = histogram(ara_err1,min = 0., max = xmax/div, nbins=100, location =loc2)
  ymax = max ([histo1, histo2])
  binz = xmax / (99)

  plot,indgen(100),replicate(0.,100),xrange = [0.,xmax],yrange =[0., ymax], color =0, charsize = 1.,XSTYLE=1,YSTYLE=1,title = 'Distribution of Max (abs(difference)) of Shape Param (NLv4 vs NLv4p) ' + strtrim(string (n_elements(ara_err1)),2) + ' tiles affected', xtitle = 'Max (Abs (Diff)) of ARA for CATDEF [0,CDCR1]', ytitle = 'Frequency'

  for i = 0, 98 do begin
     xbox = [loc(i), loc(i), loc(i+1), loc(i+1), loc(i)]
     ybox = [0, histo1(i), histo1(i), 0, 0]
     polyfill, xbox,ybox, color = 0
  endfor
  loc = loc2
  
  plot,indgen(100),replicate(0.,100),xrange = [0.,xmax/div],yrange =[0., ymax], color =0, charsize = 1.,XSTYLE=1,YSTYLE=1, title = 'Distribution of Max (abs(difference)) of Shape Param (NLv4 vs NLv4p) ' + strtrim(string (n_elements(ara_err1)),2) + ' tiles affected', xtitle = 'Max (Abs (Diff)) of ARA for CATDEF [0,CDCR1]', ytitle = 'Frequency'

  for i = 0, 98 do begin
     xbox = [loc(i), loc(i), loc(i+1), loc(i+1), loc(i)]
     ybox = [0, histo2(i), histo2(i), 0, 0]
     polyfill, xbox,ybox, color = 0
  endfor

  ; ############################################################################
  
  xmax = 0.25                              ; ZBAR
  histo1 = histogram(zbar_err1,min = 0., max = xmax, nbins=100, location =loc)
  histo2 = histogram(zbar_err1,min = 0., max = xmax/div, nbins=100, location =loc2)
  ymax = max ([histo1, histo2])
  binz = xmax / (99)

  plot,indgen(100),replicate(0.,100),xrange = [0.,xmax],yrange =[0., ymax], color =0, charsize = 1.,XSTYLE=1,YSTYLE=1,title = 'Distribution of Max (abs(difference)) of ZBAR (NLv4 vs NLv4p) ' + strtrim(string (n_elements(zbar_err1)),2) + ' tiles affected', xtitle = 'Max (Abs (Diff)) of ZBAR for CATDEF [0,CDCR2]', ytitle = 'Frequency'

  for i = 0, 98 do begin
     xbox = [loc(i), loc(i), loc(i+1), loc(i+1), loc(i)]
     ybox = [0, histo1(i), histo1(i), 0, 0]
     polyfill, xbox,ybox, color = 0
  endfor

  loc = loc2
  plot,indgen(100),replicate(0.,100),xrange = [0.,xmax/div],yrange =[0., ymax], color =0, charsize = 1.,XSTYLE=1,YSTYLE=1, title = 'Distribution of Max (abs(difference)) of ZBAR (NLv4 vs NLv4p) ' + strtrim(string (n_elements(zbar_err1)),2) + ' tiles affected', xtitle = 'Max (Abs (Diff)) of ZBAR for CATDEF [0,CDCR2]', ytitle = 'Frequency'

  for i = 0, 98 do begin
     xbox = [loc(i), loc(i), loc(i+1), loc(i+1), loc(i)]
     ybox = [0, histo2(i), histo2(i), 0, 0]
     polyfill, xbox,ybox, color = 0
  endfor

  ; ############################################################################
  
  xmax = .1                             ; topo mean
  histo1 = histogram(topo_err1,min = 0., max = xmax, nbins=100, location =loc)
  histo2 = histogram(topo_err1,min = 0., max = xmax/div, nbins=100, location =loc2)
  ymax = max ([histo1, histo2])
  binz = xmax / (99)

  plot,indgen(100),replicate(0.,100),xrange = [0.,xmax],yrange =[0., ymax], color =0, charsize = 1.,XSTYLE=1,YSTYLE=1,title = 'Distribution of Max (abs(difference)) of mean CTI (NLv4 vs NLv4p) ' + strtrim(string (n_elements(topo_err1)),2) + ' tiles affected', xtitle = 'Max (Abs (Diff)) of mean CTI', ytitle = 'Frequency'

  for i = 0, 98 do begin
     xbox = [loc(i), loc(i), loc(i+1), loc(i+1), loc(i)]
     ybox = [0, histo1(i), histo1(i), 0, 0]
     polyfill, xbox,ybox, color = 0
  endfor
  loc = loc2
  
  plot,indgen(100),replicate(0.,100),xrange = [0.,xmax/div],yrange =[0., ymax], color =0, charsize = 1.,XSTYLE=1,YSTYLE=1, title = 'Distribution of Max (abs(difference)) of mean CTI (NLv4 vs NLv4p) ' + strtrim(string (n_elements(topo_err1)),2) + ' tiles affected', xtitle = 'Max (Abs (Diff)) of mean CTI', ytitle = 'Frequency'

  for i = 0, 98 do begin
     xbox = [loc(i), loc(i), loc(i+1), loc(i+1), loc(i)]
     ybox = [0, histo2(i), histo2(i), 0, 0]
     polyfill, xbox,ybox, color = 0
  endfor

  ; ############################################################################
  
  xmax = 0.08                              ; TSC2
  histo1 = histogram(tsc_err1,min = 0., max = xmax, nbins=100, location =loc)
  histo2 = histogram(tsc_err1,min = 0., max = xmax/div, nbins=100, location =loc2)
  ymax = max ([histo1, histo2])
  binz = xmax / (99)

  plot,indgen(100),replicate(0.,100),xrange = [0.,xmax],yrange =[0., ymax], color =0, charsize = 1.,XSTYLE=1,YSTYLE=1,title = 'Distribution of Max (abs(difference)) of TSC2 (NLv4 vs NLv4p) ' + strtrim(string (n_elements(tsc_err1)),2) + ' tiles affected', xtitle = 'Max (Abs (Diff)) of TSC for RZEXC [-20,20], CATDEF [0., CDCR2]', ytitle = 'Frequency'

  for i = 0, 98 do begin
     xbox = [loc(i), loc(i), loc(i+1), loc(i+1), loc(i)]
     ybox = [0, histo1(i), histo1(i), 0, 0]
     polyfill, xbox,ybox, color = 0
  endfor

  loc = loc2
  plot,indgen(100),replicate(0.,100),xrange = [0.,xmax/div],yrange =[0., ymax], color =0, charsize = 1.,XSTYLE=1,YSTYLE=1, title = 'Distribution of Max (abs(difference)) of TSC2 (NLv4 vs NLv4p) ' + strtrim(string (n_elements(tsc_err1)),2) + ' tiles affected', xtitle = 'Max (Abs (Diff)) of TSC for RZEXC [-20,20], CATDEF [0., CDCR2]', ytitle = 'Frequency'

  for i = 0, 98 do begin
     xbox = [loc(i), loc(i), loc(i+1), loc(i+1), loc(i)]
     ybox = [0, histo2(i), histo2(i), 0, 0]
     polyfill, xbox,ybox, color = 0
  endfor

  !P.Multi = [0, 1, 2, 0, 0]
  !p.background = 255
  !P.position   = 0
  !p.charsize   = 1.2
  Erase,255

  ; -------
  ; max_loc
  ; -------
    
  ;; Asat
  NCDF_VARGET, ncid1, 'ARS1', ARS1_1
  NCDF_VARGET, ncid1, 'ARS2', ARS2_1
  NCDF_VARGET, ncid1, 'ARS3', ARS3_1
  NCDF_VARGET, ncid2, 'ARS1', ARS1_2
  NCDF_VARGET, ncid2, 'ARS2', ARS2_2
  NCDF_VARGET, ncid2, 'ARS3', ARS3_2
  
  ymin = min(asat_err1, subscript_max = kloc)
  n = asat_loc(kloc) -1
  y1 = fltarr (fix(cdcr1(n)))
  y2 = fltarr (fix(cdcr1(n)))
  for i = 0, fix(cdcr1(n)) -1 do begin
     X = 1. *i
     y1(i)= MIN([1.,MAX([0.,(1.+ ars1_1(n)*X)/(1.+ ars2_1(n)*X + ars3_1(n)*X*X)])])
     y2(i)= MIN([1.,MAX([0.,(1.+ ars1_2(n)*X)/(1.+ ars2_2(n)*X + ars3_2(n)*X*X)])])
  endfor
  
  plot,indgen(fix(cdcr1(n))),y1,xrange = [0.,cdcr1(n)],yrange =[0., max([y1,y2])], color =0, charsize = 1.,XSTYLE=1,$
              YSTYLE=1,title = 'A!Dsat!N in tile '+ strtrim(string (asat_loc(kloc)),2) , xtitle = 'CATDEF', ytitle = 'A!Dsat!N'
  
  oplot,indgen(fix(cdcr1(n))),y1, color = 0, thick =2
  oplot,indgen(fix(cdcr1(n))),y2, color = 156, thick =2
  ARS1_1 = 0.
  ARS2_1 = 0.
  ARS3_1 = 0.
  ARS1_2 = 0.
  ARS2_2 = 0.
  ARS3_2 = 0.
    
   ; WMIN
  NCDF_VARGET, ncid1, 'ARW1',ARW1_1 
  NCDF_VARGET, ncid1, 'ARW2',ARW2_1 
  NCDF_VARGET, ncid1, 'ARW3',ARW3_1
  NCDF_VARGET, ncid1, 'ARW4',ARW4_1
  NCDF_VARGET, ncid2, 'ARW1',ARW1_2 
  NCDF_VARGET, ncid2, 'ARW2',ARW2_2 
  NCDF_VARGET, ncid2, 'ARW3',ARW3_2
  NCDF_VARGET, ncid2, 'ARW4',ARW4_2

  ymin = min(wmin_err1, subscript_max = kloc)
  n = wmin_loc(kloc) -1
  y1 = fltarr (fix(cdcr1(n)))
  y2 = fltarr (fix(cdcr1(n)))
  for i = 0, fix(cdcr1(n)) -1 do begin
     X = 1. *i
     y1(i) = MIN([1.,MAX([0.,arw4_1(n) + (1.+arw1_1(n)*X)/(1.+ arw2_1(n)*X + arw3_1(n)*X*X)])])
     y2(i) = MIN([1.,MAX([0.,arw4_2(n) + (1.+arw1_2(n)*X)/(1.+ arw2_2(n)*X + arw3_2(n)*X*X)])])
  endfor  
  plot,indgen(fix(cdcr1(n))),y1,xrange = [0.,cdcr1(n)],yrange =[0., max([y1,y2])], color =0, charsize = 1.,XSTYLE=1,YSTYLE=1,title = '!9q!X!D0!N in tile ' + strtrim(string (wmin_loc(kloc)),2), xtitle = 'CATDEF', ytitle = '!9q!X!D0!N'
  oplot,indgen(fix(cdcr1(n))),y1, color = 0, thick =2
  oplot,indgen(fix(cdcr1(n))),y2, color = 156, thick =2
  ARW1_1 = 0. 
  ARW2_1 = 0.
  ARW3_1 = 0.
  ARW4_1 = 0.
  ARW1_2 = 0.
  ARW2_2 = 0.
  ARW3_2 = 0.
  ARW4_2 = 0.
  
  ; shape
  NCDF_VARGET, ncid1, 'ARA1',ARA1_1 
  NCDF_VARGET, ncid1, 'ARA2',ARA2_1 
  NCDF_VARGET, ncid1, 'ARA3',ARA3_1
  NCDF_VARGET, ncid1, 'ARA4',ARA4_1
  NCDF_VARGET, ncid2, 'ARA1',ARA1_2 
  NCDF_VARGET, ncid2, 'ARA2',ARA2_2 
  NCDF_VARGET, ncid2, 'ARA3',ARA3_2
  NCDF_VARGET, ncid2, 'ARA4',ARA4_2

  ymin = min(ARA_err1, subscript_max = kloc)
  n = ARA_loc(kloc) -1
  y1 = fltarr (fix(cdcr1(n)))
  y2 = fltarr (fix(cdcr1(n)))
  
  for i = 0, fix(cdcr1(n)) -1 do begin
     if (ara1_1(N) ne ara3_1(N)) then begin
        y1(i) = (ara4_1(n) - ara2_1(n))/(ara1_1(n) - ara3_1(n))
     endif else begin
        y1(i) = 0.
     endelse
  
     if (ara1_2(n) ne ara3_2(n)) then begin
        y2(i) = (ara4_2(n) - ara2_2(n))/(ara1_2(n) - ara3_2(n))
     endif else begin
        y2(i) = 0.
     endelse
  endfor
  plot,indgen(fix(cdcr1(n))),y1,xrange = [0.,cdcr1(n)],yrange =[0., max([y1,y2])], color =0, charsize = 1.,XSTYLE=1,YSTYLE=1,title = 'Shape param in tile ' + strtrim(string (ara_loc(kloc)),2), xtitle = 'CATDEF', ytitle = 'ARA'
  oplot,indgen(fix(cdcr1(n))),y1, color = 0, thick =2
  oplot,indgen(fix(cdcr1(n))),y2, color = 156, thick =2  
  ARA1_1 = 0. 
  ARA2_1 = 0.
  ARA3_1 = 0.
  ARA4_1 = 0.
  ARA1_2 = 0.
  ARA2_2 = 0.
  ARA3_2 = 0.
  ARA4_2 = 0.
  
  ; bf.dat
  ; ======

  NCDF_VARGET, ncid1, 'BF1', BF1_1
  NCDF_VARGET, ncid1, 'BF2', BF2_1
  NCDF_VARGET, ncid1, 'BF3', BF3_1
  NCDF_VARGET, ncid2, 'BF1', BF1_2
  NCDF_VARGET, ncid2, 'BF2', BF2_2
  NCDF_VARGET, ncid2, 'BF3', BF3_2

  ymin = min(zbar_err1, subscript_max = kloc)
  n = zbar_loc(kloc) -1
  y1 = fltarr (fix(cdcr2(n)))
  y2 = fltarr (fix(cdcr2(n)))

  for i = 0, fix(cdcr2(n)) -1 do begin
     X = 1. *i
     y1(i) = (1.e-20 + X/bf1_1(n))^(0.5) - bf2_1(n)
     y2(i) = (1.e-20 + X/bf1_2(n))^(0.5) - bf2_2(n)     
  endfor
  ymin = min(topo_err1, subscript_max = kloc2)
  print, topo_err1(kloc2), BF3_1(topo_loc(kloc2)-1),BF3_2(topo_loc(kloc2)-1), topo_loc(kloc)
  plot,indgen(fix(cdcr2(n))),y1,xrange = [0.,cdcr2(n)],yrange =[0., max([y1,y2])], color =0, charsize = 1.,XSTYLE=1,YSTYLE=1,title = 'ZBAR in tile ' + strtrim(string (zbar_loc(kloc)),2)+ ' and Topo diff of ' + strtrim(string (BF3_1(topo_loc(kloc2)-1)),2) + ' & ' + strtrim(string (BF3_2(topo_loc(kloc2)-1)),2) +' at ' +  strtrim(string (topo_loc(kloc2)),2), xtitle = 'CATDEF', ytitle = 'ZBAR'
  
  oplot,indgen(fix(cdcr2(n))),y1, color = 0, thick =2
  oplot,indgen(fix(cdcr2(n))),y2, color = 156, thick =2
  
  BF1_1  = 0.
  BF2_1  = 0.
  BF3_1  = 0.
  BF1_2  = 0.
  BF2_2  = 0.
  BF3_2  = 0.

  ; ts.dat
  ; ======
  
  NCDF_VARGET, ncid1, 'TSA1',TSA1_1 
  NCDF_VARGET, ncid1, 'TSB1',TSB1_1 
  NCDF_VARGET, ncid1, 'TSA2',TSA2_1
  NCDF_VARGET, ncid1, 'TSB2',TSB2_1
  NCDF_VARGET, ncid2, 'TSA1',TSA1_2 
  NCDF_VARGET, ncid2, 'TSB1',TSB1_2 
  NCDF_VARGET, ncid2, 'TSA2',TSA2_2
  NCDF_VARGET, ncid2, 'TSB2',TSB2_2


  
  closeps,'idl.ps'
  spawn, 'ps2pdf idl.ps '+ 'histo3.pdf'
  stop  
end

;_________________________________________________________________________________________

pro load_colors

R = intarr (256)
G = intarr (256)
B = intarr (256)

R (*) = 255
G (*) = 255
B (*) = 255

r_drought = [0,   0,   0,   0,  47, 200, 255, 255, 255, 255, 249, 197]
g_drought = [0, 115, 159, 210, 255, 255, 255, 255, 219, 157,   0,   0]
b_drought = [0,   0,   0,   0,  67, 130, 255,   0,   0,   0,   0,   0]

colors = indgen (11) + 1
R (0:11) = r_drought
G (0:11) = g_drought
B (0:11) = b_drought

r_green = [200, 150,  47,  60,   0,   0,   0,   0]
g_green = [255, 255, 255, 230, 219, 187, 159, 131]
b_green = [200, 150,  67,  15,   0,   0,   0,   0]

r_blue  = [ 55,   0,   0,   0,   0,   0,   0,   0,   0,   0]
g_blue  = [255, 255, 227, 195, 167, 115,  83,   0,   0,   0]
b_blue  = [199, 255, 255, 255, 255, 255, 255, 255, 200, 130]

r_red   = [255, 240, 255, 255, 255, 255, 255, 233, 197]
g_red   = [255, 255, 219, 187, 159, 131,  51,  23,   0]
b_red   = [153,  15,   0,   0,   0,   0,   0,   0,   0]

r_grey  = [245, 225, 205, 185, 165, 145, 125, 105,  85]
g_grey  = [245, 225, 205, 185, 165, 145, 125, 105,  85]
b_grey  = [245, 225, 205, 185, 165, 145, 125, 105,  85]

r_type  = [255,106,202,251,  0, 29, 77,109,142,233,255,255,255,127,164,164,217,217,204,104,  0]
g_type  = [245, 91,178,154, 85,115,145,165,185, 23,131,131,191, 39, 53, 53, 72, 72,204,104, 70]
b_type  = [215,154,214,153,  0,  0,  0,  0, 13,  0,  0,200,  0,  4,  3,200,  1,200,204,200,200]

r_lct2  = [  0,   0,   0,   0,   0,   0,   0,   0,   0,  55, 120, 190, 240, 255, 255, 255, 255, 255, 233, 197, 158]
g_lct2  = [  0,   0,   0,  83, 115, 167, 195, 227, 255, 255, 255, 255, 255, 219, 187, 159, 131,  51,  23,   0,   0]
b_lct2  = [130, 200, 255, 255, 255, 255, 255, 255, 255, 199, 135,  67,  15,   0,   0,   0,   0 ,  0,   0,   0,   0]

r_grads_rb = [160, 110, 30,   0,   0,   0,   0, 160, 230, 230, 240, 250, 240]
g_grads_rb = [  0,  0,  60, 150, 200, 210, 220, 230, 220, 175, 130,  60,   0] 
b_grads_rb = [200, 220, 255, 255, 200, 140,  0,  50,  50,  45,  40,  60, 130] 

r_veg  = [233,255,255,255,210,  0,  0,  0,204,170,255,220,205,  0,  0,170,  0, 40,120,140,190,150,255,255,  0,  0,  0,195,255,  0]
g_veg  = [ 23,131,191,255,255,255,155,  0,204,240,255,240,205,100,160,200, 60,100,130,160,150,100,180,235,120,150,220, 20,245, 70]
b_veg  = [  0,  0,  0,178,255,255,255,200,204,240,100,100,102,  0,  0,  0,  0,  0,  0,  0,  0,  0, 50,175, 90,120,130,  0,215,200]

R (20:27) = r_green
G (20:27) = g_green
B (20:27) = b_green

R (30:39) = r_blue
G (30:39) = g_blue
B (30:39) = b_blue

R (40:48) = r_red
G (40:48) = g_red
B (40:48) = b_red

R (50:58) = r_grey
G (50:58) = g_grey
B (50:58) = b_grey

R (60:80) = r_type
G (60:80) = g_type
B (60:80) = b_type

R (90:119) = r_veg
G (90:119) = g_veg
B (90:119) = b_veg

R (120:132) = r_grads_rb
G (120:132) = g_grads_rb
B (120:132) = b_grads_rb

R (140:160) = r_lct2
G (140:160) = g_lct2
B (140:160) = b_lct2
TVLCT,R ,G ,B

end
