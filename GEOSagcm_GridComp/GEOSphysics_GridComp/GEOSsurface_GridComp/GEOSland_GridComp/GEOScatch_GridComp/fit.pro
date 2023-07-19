

; type 2

 vgdd=[13.66, 13.66, 14.62, 15.70, 16.33, 16.62, 16.66, 16.60, 16.41, 15.73, 14.62, 13.66]
 vgrd=[211.32, 211.32, 218.78, 243.40, 294.87, 345.90, 355.18, 341.84, 307.22, 244.84, 218.78, 211.32]
 lai =[0.520, 0.520, 0.867, 2.107, 4.507, 6.773, 7.173, 6.507, 5.040, 2.173, 0.867, 0.520]

 ll=.1*indgen(110)

 ygdd=  16.6 - 4*exp(-.75*ll)
 ygrd=  200. + 21.*ll

 plot,ll,ygdd
 oplot,lai,vgdd,psym=6

; plot,ll,ygrd
; oplot,lai,vgrd,psym=6


; type 3

 vgdd=[13.76, 13.80, 13.86, 13.88, 13.90, 13.93, 13.91, 13.89, 13.88, 13.86, 13.80, 13.76   ]
 vgrd=[565.41, 587.05, 623.46, 638.13, 652.86, 675.04,    660.24, 645.49, 638.13, 623.46, 587.05, 565.41]
 lai =[8.760,9.160, 9.827,10.093,10.360,10.760,10.493,10.227,10.093, 9.827, 9.160, 8.760]

 ll=8+.03*indgen(100)

 ygdd=  14.14 - 5.3*exp(-.3*ll)
 ygrd=  135. + 50.*ll

;plot,lai,vgdd,psym=6,yrange=[13.5,14]
;oplot,ll,ygdd

 ;plot,ll,ygrd
 ;oplot,lai,vgrd,psym=6

; type 4

 vgdd=[    0.218,   0.227,   0.233,   0.239,   0.260,   0.299,    0.325,   0.313,   0.265,   0.244,   0.233,   0.227]
 vgrd=[    24.43,  24.63,  24.80,  24.96,  25.72,  27.74,    30.06,  28.86,  25.90,  25.11,  24.80,  24.63]
 lai =[0.782,0.893, 1.004, 1.116, 1.782, 3.671, 4.782, 4.227, 2.004, 1.227, 1.004, 0.893]


 ll=.05*indgen(110)

 ygdd=  .355 - .175*exp(-.35*ll)
 ygrd=  23.5 + 1.30*ll

 ;plot,lai,vgdd,psym=6
 ;oplot,ll,ygdd

; plot,ll,ygrd
; oplot,lai,vgrd,psym=6

; type 5

 vgdd=[    2.813,   2.813,   2.662,   2.391,   2.391,   2.391,     2.391,   2.975,   3.138,   3.062,   2.907,   2.813]
 vgrd=[    103.60, 103.60, 102.35, 100.72, 100.72, 100.72,    100.72, 105.30, 107.94, 106.59, 104.49, 103.60]
 lai =[3.760,3.760, 2.760, 1.760, 1.760, 1.760, 1.760, 5.760,10.760, 7.760, 4.760, 3.760]

ll=.1*indgen(110)

 ygdd=  3.15 - 1.7*exp(-0.45*ll)
 ygrd=  100. + .95*ll

;plot,lai,vgdd,psym=6
;oplot,ll,ygdd

; plot,lai,vgrd,psym=6
; oplot,ll,ygrd

; type 6

 vgdd=[    0.10629, 0.10629, 0.10629, 0.10629, 0.10629, 0.12299,    0.21521, 0.22897, 0.19961, 0.10629, 0.10629, 0.10629]
 vgrd=[    22.86,  22.86,  22.86,  22.86,  22.86,  23.01,    24.36,  24.69,  24.04,  22.86,  22.86,  22.86]
 lai =[0.739,0.739, 0.739, 0.739, 0.739, 1.072, 5.072, 5.739, 4.405, 0.739, 0.739, 0.739]

 ll=.06*indgen(110)

 ygdd=  .27 - .20*exp(-.25*ll)
 ygrd=  22.5 + .4*ll

; plot,lai,vgdd,psym=6
; oplot,ll,ygdd

; plot,ll,ygrd
; oplot,lai,vgrd,psym=6


end
