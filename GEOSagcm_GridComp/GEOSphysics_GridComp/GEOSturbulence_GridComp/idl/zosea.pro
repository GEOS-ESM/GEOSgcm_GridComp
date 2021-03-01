

ustar=(indgen(200)+1)*.002

Z0 = (0.05*1.533e-5)/USTAR + (0.018/9.81)*USTAR^2


  cn=1000.*(0.4/alog(10./z0))^2
  u=ustar/sqrt(cn*1.e-3)

      RE  = Z0*UStar/1.533E-5
      ZT = Z0 * (0.4/RE)
zq = 1.5*zt


ct = 1000.*(0.4/alog(10./z0))*(0.4/alog(10./zT))
cq = 1000.*(0.4/alog(10./z0))*(0.4/alog(10./zq))


  plot,u,cn
  oplot,u,ct
;  oplot,u,cq

ZT = Z0 * exp(2.0 - 2.48*RE^.25)
zq = 1.5*zt

ct = 1000.*(0.4/alog(10./z0))*(0.4/alog(10./zT))
cq = 1000.*(0.4/alog(10./z0))*(0.4/alog(10./zq))


  oplot,u,ct
;  oplot,u,cq
 
ZT = Z0 * exp(0.4*(1.935 - 4.0*sqrt(RE-0.1)))
zq = 1.5*zt

ct = 1000.*(0.4/alog(10./z0))*(0.4/alog(10./zT))
cq = 1000.*(0.4/alog(10./z0))*(0.4/alog(10./zq))


  oplot,u,ct
  oplot,u,re

end
