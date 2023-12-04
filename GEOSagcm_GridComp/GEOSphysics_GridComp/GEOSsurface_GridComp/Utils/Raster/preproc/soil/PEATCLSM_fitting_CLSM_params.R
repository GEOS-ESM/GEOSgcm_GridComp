
######################################################################
# Fitting of CLSM catchment functions (Ducharne et al. 2000)
# to theoretical relationships of PEATCLSM 
# June 2019 (Michel Bechtold, KU Leuven, michel.bechtold@kuleuven.be)
# Reference: Bechtold et al. 2019 (https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018MS001574)
######################################################################
library(stats)
library(weights)

# PEATCLSM Campbell parameters 
phydr = data.frame(Ks=NA, theta_s=NA, b=NA, psi_s=NA)
phydr$Ks        = 2.8*10^(-5) *100*84000       # [m3.m-3]
phydr$theta_s   = 0.930                        # [m3.m-3]
phydr$b         = 3.5                          # [-]       # shape parameter
phydr$psi_s     = -0.03                        # [m]       # air entry pressure
# stdev of microtography elevation 
# --> taken into account for catdef (hummocks are still unsaturated at zbar = 0)
microtopo_std   = 0.12                         # [m]

# Campbell Functions
f_ret_Campbell <- function(theta_s, b, psi_s, psi){
  theta = theta_s * (psi/psi_s)^(-1/b)   # psi = psi_s * (theta/theta_s)^(-b)
  theta[theta>theta_s] = theta_s
  return(theta)
}
f_K_Campbell <- function(ks, b, psi_s, psi){
  k = ks* (psi_s/psi)^(2+2/b)
  k[k>ks] = ks
  return(k)
}
Campbell_1d <- function(z_){
  theta = theta_s * ((z_*100)/(psi_s*100))^(-1/b)
  theta[(z_*100)>=(psi_s*100)]=theta_s
  return(theta)
}
#### Azl, Azu, and Sy_soil 
# see Dettmann and Bechtold, 2015, Hydrol Processes
Campbell_1d_Azl <- function(z_,zl,theta_s,psi_s,b,microtopo_std){
  Fs = pnorm(z_, mean = 0, sd = microtopo_std, log = FALSE)
  #Fs = 0
  theta = theta_s * (((zl-z_)*100)/(psi_s*100))^(-1/b)
  theta[((zl-z_)*100)>=(psi_s*100)]=theta_s
  theta_Fs = (1-Fs)*theta
  return(theta_Fs)
}
Campbell_1d_Azu <- function(z_,zu,theta_s,psi_s,b,microtopo_std){
  Fs = pnorm(z_, mean = 0, sd = microtopo_std, log = FALSE)
  #Fs = 0
  theta = theta_s * (((zu-z_)*100)/(psi_s*100))^(-1/b)
  theta[((zu-z_)*100)>=(psi_s*100)]=theta_s
  theta_Fs = (1-Fs)*theta
  return(theta_Fs)
}
get_Sy_soil <- function(zl_,zu_,phydr,microtopo_std){
theta_s = phydr$theta_s
b = phydr$b
psi_s = phydr$psi_s
Sy = rep(NA,length(zl_))
for (i in 1:length(zl_)){
  zl=zl_[i]
  zu=zu_[i]
  A = 0
  for (j in 1:200){
    zm=0.5*(zl_[j]+zu_[j])
    Azl = Campbell_1d_Azl(zm,zl,theta_s,psi_s,b,microtopo_std)
    Azu = Campbell_1d_Azu(zm,zu,theta_s,psi_s,b,microtopo_std)
    A = A + (zu_[j]-zl_[j]) *  (Azu-Azl)
  }
  Sy[i] = 1/(1*(zu-zl)) * A
}
return(Sy)
}

# calculate from -1 to 1 m
zl_ = seq(-1,1,0.01)
zu_ = seq(-0.99,1.01,0.01)
########################################################### 
# Fit bf1 bf2
########################################################### 
### Theoretical Sy - zbar curve
Sy_soil = get_Sy_soil(zl_,zu_,phydr,microtopo_std)
x=(zl_+zu_)/2
Sy_surface = pnorm(0.5*(zu_+zl_), mean = 0, sd = microtopo_std, log = FALSE)
Sy = Sy_soil + Sy_surface

### derive catdef - zbar curve from Sy  
catdef = rep(NA,length(Sy))
for (i in 1:length(zl_)){
  zl=zl_[i]
  zu=zu_[i]
  catdef[i] = sum(Sy_soil[i:length(zl_)])*0.01
}
catdef = catdef*1000

######## FIT bf parameters to theoretical curve
## fit limited for zbar (negative = below ground) range from -0.6 m  to 0.1 m
controllist = nls.control(maxiter = 100,minFactor=0.0000000000000001)
result=1
while (result==1){
  result <- tryCatch({
    bf1_limits = c(10,1000)
    bf2_limits = c(0.01,10.0)
    bf1_ini = runif(1,99,101)
    bf2_ini = runif(1,0.29,0.31)
    #t = nls(y ~ ((mean(y) + c)+ a*(x - alf_r) + 0.5 * b * (x - alf_r)^2), data=df3, start = list(a = a_ini, b = b_ini, c = c_ini), algorithm="port", lower = c(-0.2,0.0,-4), upper= c(0.0,0.01,4),
    #        control=controllist )
    df = data.frame(catdef=catdef[x>(-0.6) & x<(0.1)])
    df$x = -x[x>(-0.6) & x<(0.1)]
    t = nls(x ~ ((1.0E-20 + (catdef)/a)^0.5 - b), data=df, 
            start = list(a = bf1_ini, b = bf2_ini), algorithm="port", 
            lower = c(bf1_limits[1],bf2_limits[1]), upper= c(bf1_limits[2],bf2_limits[2]),control=controllist )
    ll=2
  }, error = function(e) {
    ll=1
  })
}
bf1 = coef(summary(t))[1]
bf2 = coef(summary(t))[2]
plot(x,catdef,xlim=c(-0.6,0.1),ylim=c(0,300))
zbar = ((1.0E-20 + (catdef)/bf1)^0.5 - bf2)
points(-zbar,catdef,pch=19,col='red')

########################################################### 
# Area partioniong parameters
# FIT ars* parameters to theoretical curve
########################################################### 
# AR1: saturated area --> in PEATCLSM AR1: area with surface water
ar1_theoretical = pnorm(-zbar, mean = 0, sd = microtopo_std, log = FALSE)
fitted_range = -zbar>(-0.6) & -zbar<(0.0)
controllist = nls.control(maxiter = 100,minFactor=0.0000000000000001)
result=1
while (result==1){
  result <- tryCatch({
    ars1_limits = c(-0.01,-0.001)
    ars2_limits = c(0.001,0.1)
    ars3_limits = c(0.0001,0.1)
    ars1_ini = runif(1,-0.006,-0.004)
    ars2_ini = runif(1,0.0019,0.0021)
    ars3_ini = runif(1,0.0019,0.0021)
    #t = nls(y ~ ((mean(y) + c)+ a*(x - alf_r) + 0.5 * b * (x - alf_r)^2), data=df3, start = list(a = a_ini, b = b_ini, c = c_ini), algorithm="port", lower = c(-0.2,0.0,-4), upper= c(0.0,0.01,4),
    #        control=controllist )
    df = data.frame(catdef=catdef[fitted_range])
    df$ar1 = ar1_theoretical[fitted_range]
    t = nls(ar1 ~ (1+ars1*(catdef))/(1+ars2*(catdef)+ars3*(catdef)^2), data=df, 
            start = list(ars1 = ars1_ini, ars2 = ars2_ini, ars3 = ars3_ini), algorithm="port", 
            lower = c(ars1_limits[1],ars2_limits[1],ars3_limits[1]), upper= c(ars1_limits[2],ars2_limits[2],ars3_limits[2]),control=controllist )
    ll=2
  }, error = function(e) {
    ll=1
  })
}
ars1 = coef(summary(t))[1]
ars2 = coef(summary(t))[2]
ars3 = coef(summary(t))[3]
cond = fitted_range
ar1 = (1+ars1*(catdef))/(1+ars2*(catdef)+ars3*(catdef)^2)
# ar1 vs catdef
plot(catdef[cond],ar1_theoretical[cond],ylim=c(0,0.6))
points(catdef[cond],ar1[cond],pch=19,col='red')
# ar1 vs zbar
plot(zbar[cond],ar1_theoretical[cond],ylim=c(0,0.6))
points(zbar[cond],ar1[cond],pch=19,col='red')


###########################################3
# mean root zone and surface soil moisture from water table distribution
zbar=seq(-4,1,0.01)
catdef = ((-zbar + bf2)^2 - 1.0E-20) * bf1 
catdef[(-zbar + bf2)<=0] = 0
#plot(zbar,catdef,xlim=c(-4,1))
rzmean = rep(NA,length(catdef))
srfmean = rep(NA,length(catdef))
#zbar_plot = c(-1.0)
### root zone moisture in dependence of water table 
theta = rep(NA,length(catdef))
theta_srf = rep(NA,length(catdef))
# input for integrate(Campbell_1d...)
theta_s = phydr$theta_s
psi_s = phydr$psi_s 
b = phydr$b
for (i in 1:length(zbar)){
  theta[i] = integrate(Campbell_1d,zbar[i],zbar[i]+1.0)$value
  theta_srf[i] = 20*integrate(Campbell_1d,zbar[i],zbar[i]+0.05)$value
}
###
for (k in 1:length(zbar)){
  pd = dnorm(-zbar[k]+zbar, mean = 0, sd = microtopo_std, log = FALSE)
  # weighted mean
  rzmean[k]=wtd.mean(theta/theta_s,weights=pd)
  # calculate mean surface moisture only for ar2 and ar4
  # theta < theta_s
  cond_srf = theta_srf < theta_s-0.0001
  srfmean[k]=wtd.mean(theta_srf[cond_srf]/theta_s,weights=pd[cond_srf])
}

# compare with and without water table distribution
#plot(zbar,rzmean,xlim=c(-0.5,0.1))
#lines(zbar,theta/theta_s,col="red")
#plot(catdef,rzmean,xlim=c(0,1000),ylim=c(0,1))
#lines(catdef,theta/theta_s,col="red")
#plot(zbar,srfmean,xlim=c(-0.9,0.0))
#lines(zbar,theta_srf/theta_s,col="red")

### Coefficients of polynomial function in PEATCLSM
# for integration of equilibrium surface soil moisture in ar2 and ar4
cond_fit=zbar>(-1.2) & zbar<0.0
m_srf <- lm(srfmean[cond_fit] ~ poly(-zbar[cond_fit], 4, raw=TRUE))
lines(zbar[cond_fit],predict(m_srf),col='blue',lwd=4)
m_srf$coefficients

############## MEAN Rootzone moisture as function of catdef
f_RZEQXI <- function(catdef,ara1,ara2,arw1,arw2,arw3,arw4) {
   # for PEATCLSM ara3=ara1 and ara4=ara2
   ara3=ara1
   ara4=ara2
   if (ara1 != ara3) {
   cdi=(ara4-ara2)/(ara1-ara3)
   } else {
     cdi=0.
   }
   ara3=ara1
   ara4=ara2
   
   RZEQXI = catdef
   RZEQXI[] = NA
   AR1=RZEQXI
   for (k in 1:length(catdef)){
     AR1[k]=min(1.,max(0.,(1.+ars1*catdef[k])/(1.+ars2*catdef[k]+ars3*catdef[k]*catdef[k]))) 
     #AR1 = (1+ars1*(catdef))/(1+ars2*(catdef)+ars3*(catdef)^2)
     
     if (catdef[k] >= cdi) {
       ax=ara3*catdef[k]+ara4
     } else {
       ax=ara1*catdef[k]+ara2
     }
     
     WMIN_TEST = 1+arw2*catdef[k]+arw3*catdef[k]*catdef[k]
     if (abs(WMIN_TEST)>1.e-6){
       WMIN = arw4 + (1-arw4) * (1+arw1*catdef[k])/(1+arw2*catdef[k]+arw3*(catdef[k])^2)
     } else {
       WMIN = min(1,max(0,arw4+(1-arw4)*(1+arw1*catdef[k])/1.e-6))
     }
     
     ARG1=max(-40., min(40., -ax*(1.-WMIN)))
     EXPARG1=exp(ARG1)
     RZEQX=(WMIN-1.-(2./ax))*EXPARG1 + WMIN + (2./ax)
     RZEQXI[k]=ax*EXPARG1*( -1. -(2./ax) - (2./(ax*ax)) + WMIN + (WMIN/ax) )  + WMIN + 2./ax
     AR20=1.+(-ax-1.+ax*WMIN)*EXPARG1
     RZEQXI[k]=RZEQXI[k]/(AR20+1.E-20)
   }
   return(RZEQXI)
}
 
#### FIT RZEQXI
fitted_range = zbar>(-0.7) & zbar<0.0
controllist = nls.control(maxiter = 100,minFactor=0.0000000000000001)
result=1
while (result==1){
  result <- tryCatch({
    ara1_limits = c(0.1,100.0)
    ara2_limits = c(0.001,3000)
    arw1_limits = c(0.001,0.1)
    arw2_limits = c(0.001,0.1)
    arw3_limits = c(1e-5,1e-3)
    arw4_limits = c(0.0001,0.8)
    ara1_ini=runif(1,ara1_limits[1],ara1_limits[2])   
    ara2_ini=runif(1,ara2_limits[1],ara2_limits[2])
    arw1_ini=runif(1,arw1_limits[1],arw1_limits[2])
    arw2_ini=runif(1,arw2_limits[1],arw2_limits[2])
    arw3_ini=runif(1,arw3_limits[1],arw3_limits[2])
    arw4_ini=runif(1,arw4_limits[1],arw4_limits[2])
    #t = nls(y ~ ((mean(y) + c)+ a*(x - alf_r) + 0.5 * b * (x - alf_r)^2), data=df3, start = list(a = a_ini, b = b_ini, c = c_ini), algorithm="port", lower = c(-0.2,0.0,-4), upper= c(0.0,0.01,4),
    #        control=controllist )
    df = data.frame(catdef=catdef[fitted_range])
    df$rzmean = rzmean[fitted_range]
    t = nls(rzmean ~ f_RZEQXI(catdef,ara1,ara2,arw1,arw2,arw3,arw4), data=df, 
            start = list(ara1 = ara1_ini, ara2 = ara2_ini, arw1 = arw1_ini, arw2 = arw2_ini, arw3 = arw3_ini, arw4 = arw4_ini), algorithm="port", 
            lower = c(ara1_limits[1],ara2_limits[1],arw1_limits[1],arw2_limits[1],arw3_limits[1],arw4_limits[1]), 
            upper= c(ara1_limits[2],ara2_limits[2],arw1_limits[2],arw2_limits[2],arw3_limits[2],arw4_limits[2]),control=controllist )
    ll=2
  }, error = function(e) {
    ll=1
  })
}

ara1 = coef(summary(t))[1]
ara2 = coef(summary(t))[2]
arw1 = coef(summary(t))[3]
arw2 = coef(summary(t))[4]
arw3 = coef(summary(t))[5]
arw4 = coef(summary(t))[6]

RZEQXI = f_RZEQXI(catdef,ara1,ara2,arw1,arw2,arw3,arw4)
plot(catdef[fitted_range],RZEQXI[fitted_range],ylim=c(0,1))
points(catdef[fitted_range],rzmean[fitted_range],col='red',pch=19)
plot(zbar,RZEQXI,xlim=c(-0.7,0),ylim=c(0,1))
points(zbar,rzmean,col='red',pch=19)

# for ar* file
sprintf("%1.7e   %1.7e   %1.7e   %1.7e   %1.7e   %1.7e   %1.7e   %1.7e   %1.7e   %1.7e   %1.7e", ars1, ars2, ars3, ara1, ara2, ara1 ,ara2, arw1, arw2, arw3, arw4)
# for bf* file
sprintf("%1.7e  %1.7e", bf1, bf2)
# permanent wilting point here defined as surface soil moisture at 100 cm wtd
PWPW = f_ret_Campbell(phydr$theta_s, phydr$b, phydr$psi_s, -1)
PWPW = f_ret_Campbell(phydr$theta_s, phydr$b, phydr$psi_s, -1)/0.93
sprintf("%1.3f", PWPW)
