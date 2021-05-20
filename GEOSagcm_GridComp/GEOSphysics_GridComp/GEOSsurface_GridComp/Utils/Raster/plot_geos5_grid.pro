pro plot_geos5_grid, nc = nc, nr=nr,ntiles=ntiles ,ii=ii, jj=jj,fr=fr,   $
                     x=x,y=y,m_scale_factor,ctitle=ctitle

m_scale_grid = replicate(0.,nc,nr)
land_frac = replicate(0.,nc,nr)
colors=80+indgen(21)
upval = floor(max(m_scale_factor)) +1.
lwval=0.
levm = [lwval,lwval+(upval-lwval)/21+indgen(20)*(upval-lwval)/21.,upval]

for n = 0l,ntiles-1l do begin
    land_frac(ii(n)-1,jj(n)-1) = land_frac(ii(n)-1,jj(n)-1) + fr(n)
    m_scale_grid(ii(n)-1,jj(n)-1) =  m_scale_grid(ii(n)-1,jj(n)-1) + fr(n)*m_scale_factor(n)
endfor
for j = 0, nr -1 do begin
    for i = 0, nc-1 do begin
        if(land_frac(i,j) gt 0.)  then m_scale_grid(i,j) =  m_scale_grid(i,j) /land_frac(i,j)
    endfor
endfor

gplot, m_scale_grid,x,y, /shade, lev=levm, c_color=colors, title=ctitle,oc_color=255,oceanmask=1

end   
