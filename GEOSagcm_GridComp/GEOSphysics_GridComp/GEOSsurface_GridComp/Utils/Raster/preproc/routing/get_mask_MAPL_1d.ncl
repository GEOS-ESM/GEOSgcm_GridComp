begin

nt=2592

f=addfile("inputs/MAPL_Tripolar.nc","r")
msk_MAPL=f->mask

t2lati=asciiread("inputs/TM0072xTM0036_tile_to_MAPL_lati.txt",nt,"integer")
t2loni=asciiread("inputs/TM0072xTM0036_tile_to_MAPL_loni.txt",nt,"integer")

msk_tile=new(nt,integer)
msk_tile=msk_tile
do i=0,nt-1
  msk_tile(i)=toint(msk_MAPL(t2lati(i)-1,t2loni(i)-1))
end do

asciiwrite("outputs/mask_MAPL_1d_TM0072xTM0036.txt",msk_tile)


end
