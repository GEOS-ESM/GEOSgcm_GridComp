pro cleanunder,a

bad=where( abs( a ) lt 1.0e-7 )

if min(bad) gt -1 then a(bad)=0.00


return
end
