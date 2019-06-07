pro readtilefile,f=f,til=til

close,1
openr,1,f
ntil=1L & ngrd=1L 
i1=0L & j1=0L
n1$=''
i2=0L & j2=0L
n2$=''
r=fltarr(12)



readf,1,ntil
readf,1,ngrd

readf,1,n1$
readf,1,i1
readf,1,j1

readf,1,n2$
readf,1,i2
readf,1,j2


til=fltarr(12,ntil)


for i=0L,ntil-1 do begin

    readf,1,r
    til(*,i)=r

endfor


close,1

return
end
