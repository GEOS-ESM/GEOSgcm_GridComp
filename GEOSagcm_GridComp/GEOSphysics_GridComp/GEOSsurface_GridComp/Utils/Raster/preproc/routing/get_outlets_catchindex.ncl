begin
nc=291284
ns=22612
ng=525

msk=asciiread("inputs/Pfaf_msk.txt",nc,"integer")
outid=new(ns,integer)
k=0
do i=0,nc-1
  if(msk(i).eq.2)then
    outid(k)=i+1
    k=k+1
  end if
end do
do i=k,ns-1
  outid(i)=nc+i-k+1
end do
asciiwrite("outputs/outlet_catchindex.txt",outid)

mskall=new(nc+ng,integer)
mskall(0:nc-1)=msk
mskall(nc:)=2
asciiwrite("outputs/Pfaf_msk_all.txt",mskall)

final=asciiread("inputs/Pfaf_finalID.txt",nc,"integer")
finalall=new(nc+ng,integer)
finalall(0:nc-1)=final
do i=nc,nc+ng-1
  finalall(i)=i+1
end do
asciiwrite("outputs/Pfaf_finalID_all.txt",finalall)

end
