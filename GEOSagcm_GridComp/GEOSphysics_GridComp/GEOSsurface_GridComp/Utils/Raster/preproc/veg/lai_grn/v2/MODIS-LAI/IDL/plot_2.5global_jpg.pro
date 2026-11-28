limit_global = [-60,-180,90,180]
limit_conus  = [29.,-125.,54.,-82.]
limit_east   = [29,-84.,36.,-77]
limit_scan   = [55.,0.,75.,55.]
limits = limit_global

n_levels = 21
loadct2,0
colors = indgen(n_levels)+ 80
upval = 7.
lwval = 0.
levels = [lwval,lwval+(upval-lwval)/(n_levels -1) +indgen(n_levels -1)*(upval-lwval)/(n_levels -1)]
;n_levels = 256
;mg_loadct, 0, file='/home/sarith/IDL/lib/brewer/brewer.tbl',rgb_table=ct
;      RED   = CT(0:N_levels-1,0)
;      GREEN = CT(0:N_levels-1,1)
;      BLUE  = CT(0:N_levels-1,2)
;
;colors = indgen(n_levels)
;upval = 7.
;lwval = 0.
;levels = [lwval,lwval+(upval-lwval)/(n_levels -1) +indgen(n_levels -1)*(upval-lwval)/(n_levels -1)]

i_sib=8640
j_sib=4320

x25=indgen(i_sib)*2.5/60. -180.+ 1.25/60.
y25=indgen(j_sib)*2.5/60. -90. + 1.25/60.

;for dofm =0,365,8 do begin  
;for dofm =0,365,8 do begin  
for dofm =152,152,8 do begin
    mytime=timegen(start=julday(1,1,2006),final=julday(1,1,2006)+dofm,unit='d')
    caldat,mytime,month,day,yyyy
    ii=n_elements(mytime)-1
    date_stamp=string(yyyy[ii],month[ii],day[ii],format='(i4.4,".",i2.2,".",i2.2)')
    MMDD= string(month[ii],day[ii],format='(i2.2,".",i2.2)')   

    fout = '/l_data/model_parameters/LAI/MODIS-LAI/MODIS_8-DayClim/MOD15A2.YYYY.' + MMDD + '.global_2.5min.data'
    print,fout
    openr, lun,fout,/get_lun,/F77_UNFORMATTED

    if file_test(fout) then begin
            lai = fltarr(i_sib,j_sib)
            strip = fltarr (8640)
            for j = 0,4320 - 1 do begin
            readu,lun,strip
            lai (*,j) = strip(*)
        endfor
        close,lun

        lai (where (lai eq 9999.)) =  !VALUES.F_NAN
      !P.position=0
      !P.Multi = [0, 1, 1, 0, 1]
      Erase, Color= 255
      MAP_SET,/ROBINSON,/hires,color= 0,/NoErase,limit=limits,title =mmdd
      contour,lai,x25,y25,levels = levels,c_colors=colors,/overplot

   alpha=fltarr(n_levels,2)
   alpha(*,0)=levels
   alpha(*,1)=levels
   h=[0,1]
   !P.position=[0.30, 0.0+0.005, 0.70, 0.015+0.005]
   clev = levels
   clev (*) = 1
  contour,alpha,levels,h,levels=levels,c_colors=colors,/fill,/xstyle,/ystyle, $
           /noerase,yticks=1,ytickname=[' ',' '] ,xrange=[min(levels),max(levels)], $
           xtitle=' ', color=0,xtickv=levels, $
           C_charsize=1.0, charsize=0.5 ,xticks=10
   contour,alpha,levels,h,levels=levels,color=0,/overplot,c_label=clev
   for k = 0,n_elements(levels) -1 do xyouts,levels(k),1.1, string(levels(k),'(f4.2)'),orientation=90,color=0
;thisDevice = !D.NAME
;SET_PLOT, 'Z'
;
;TVLCT, red, green, blue
;SET_PLOT, thisDevice
;
;image = TVRD()
;thisImage = BYTSCL(image)
;s = SIZE(thisImage)
;image3d = BYTARR(s(1), s(2),3)
;image3d( *, *,0) = red(thisImage)
;image3d( *, *,1) = green(thisImage)
;image3d( *, *,2) = blue(thisImage)
;write_jpeg,'MOD15A2_YYYY' + MMDD + '.jpg',image3d,quality=75,true=3
;wdelete
      write_to_jpeg,file_name= 'MOD15A2_YYYY' + MMDD + '.jpg'
    endif
endfor
end
