;
; This IDL code is good to plot global maps of MODIS scale factors
; Contact - Sarith (sarith.p.mahanama@nasa.gov)
;
; Usage -Set path to rst and til files & output path
gfile = '../DC144_91/DC0144xPC0091_DE0360xPE0180-Pfaf.notiny'
pathout = '../DC144_91/'
poles  = 'PC'
dateline = 'DC'
; End user defined variables

tilefile=gfile
ip =0
nc =0
nr =0
chard ='   '
openr,10,strtrim(tilefile,2)+'.til'
readf,10,ip
readf,10,ip
readf,10,chard
readf,10,nc
readf,10,nr
close,10

resol=strtrim((string(nc)),2) +'x'+strtrim((string(nr)),2)

dx = 360./fix(nc)
dy = 180./fix(nr-1)
if (poles eq 'PE') then dy = 180./fix(nr)
xylim=[-90.,-179.999,90.,180.]

y = [-90.0+ 0.25*dy,(indgen(nr-2))*dy -90.+dy,90.-0.25*dy]
if (poles eq 'PE') then y = [indgen(nr)*dy -90.+0.5*dy]
loadct2,0

if (dateline eq 'DE') then begin
x = indgen(nc)*dx -180. + dx/2.
endif else begin
x = indgen(nc)*dx -180. 
endelse

tilefile=strtrim(gfile,2)+'.til'
get_frac, til=tilefile, ntiles=ntiles ,ii=ii, jj=jj,fr=fr
m_scale_factor = fltarr(ntiles)

sibname =['albvr','albnr',  'albvf','albnf']
sibname = reform(sibname,2,2,/overwrite)
idum =0l
month = ['JAN',  'FEB',  'MAR',  'APR',  'MAY',   $
'JUN',  'JUL',  'AUG',  'SEP',  'OCT',  'NOV',  'DEC']


for  ialbt = 0,1 do begin
    for ialbs = 0,1 do begin
        sp,1,4,/f,/h,/color,landscape=1
        ofile = strtrim(pathout,2) +   'modis_scale_factor.' + sibname(ialbs,ialbt) + '.clim'

        psfile =pathout+'modis_scale_factor.' + sibname(ialbs,ialbt) + '.clim.ps'

        openr,unit,ofile,/GET_LUN,/xdr

        print,ofile
        
        for mon = 1,12 do begin
            readu,unit,idum,m_scale_factor,idum
            ctitle = 'Scale Factor : '+ sibname(ialbs,ialbt) + ' - ' + month(mon-1)
            
            plot_geos5_grid, nc = nc, nr=nr,ntiles=ntiles ,ii=ii, jj=jj,fr=fr,   $
                     x=x,y=y,m_scale_factor,ctitle=ctitle                      
            
        endfor
        close,unit
        FREE_LUN,unit
        wait,30
        closeps,psfile
    endfor
endfor

end
