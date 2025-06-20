pro readlaigrn,f=f,nt=nt,lai=lai,grn=grn

a=fltarr(nt)

lai=fltarr(nt,12)
grn=fltarr(nt,12)

close,1
openr,1,/f77_u,f

for i=0L,11Ldo begin

print,' month=',i+1
   readu,1,a
    lai(*,i)=a
   readu,1,a
    grn(*,i)=a

endfor

close,1

return
end
