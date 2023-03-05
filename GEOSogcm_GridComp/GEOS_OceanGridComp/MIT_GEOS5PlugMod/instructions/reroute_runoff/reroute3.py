#!/usr/bin/env python

#f2py -c -m reroute3 reroute3.f90
#module load other/comp/gcc-5.3-sp3
#module load other/SSSO_Ana-PyD/SApd_4.2.0_py2.7_gcc-5.3-sp3
import os
import scipy as sp
import netCDF4 as nc
import math as math
from reroute3 import reroute3

#dir=os.environ['NOBACKUP']+'/workdir/move_runoff'

# Read ocean land mask
print 'Reading ocean land mask'
gfile='./TEST.mit_mask.20000415_0015z.nc4'
ff=nc.Dataset(gfile)
wet=ff.variables['MASKO'][0]
wet_orig=wet.copy()

# Mask Uruguay, Amazon, Neva for 1/4 degree grid
#wet[340:355, 880:895]=0
#wet[494:505, 915:926]=0
#wet[796:805, 1232:1241]=0

sh=wet.shape
nps=360
ssize=sh[0]

# list of blank squeres
# cat data.exch2 | tail -n 1333 | head -n 1332 | tr -d "\n" | tr -d " "
#f = open('blanklist.txt', 'r+')
#BL = sp.array(f.read().split(','))
#f.close()
BlankList=sp.array([1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,21,22,23,24,\
  65,71,75,76,90,95,96,101,102,109,110,111,112,113,114,115,116,117,118,119,\
  120,121,122,123,124,125,126,127,128,129,130,131,132,\
  189,190,193,194,195,196,199,\
  200,201,202,203,205,206,207,208,209,211,212,213,214,215,216,247,253,\
  267,268,269,270,287,288,305,306,323,324,341,342,359,360,362,377,378,\
  380,381,382,395,396,400,412,413,414,430,432,450,468])-1

nbl=BlankList.shape[0]

# dimentions of the five gcmfaces:
f=sp.zeros((5,2))
f[0,:]=[math.sqrt((nps+nbl)/13),3*math.sqrt((nps+nbl)/13)]
f[1,:]=[math.sqrt((nps+nbl)/13),3*math.sqrt((nps+nbl)/13)]
f[2,:]=[math.sqrt((nps+nbl)/13),math.sqrt((nps+nbl)/13)]
f[3,:]=[3*math.sqrt((nps+nbl)/13),math.sqrt((nps+nbl)/13)]
f[4,:]=[3*math.sqrt((nps+nbl)/13),math.sqrt((nps+nbl)/13)]

f=sp.int_(f)
f1=sp.int_(f*ssize)

# reshape wet to nps squers
x=sp.reshape(wet,(ssize,sh[1]/ssize,ssize))

# convesion to llc faces
print 'convesion to llc faces'
nface=0;
nface2=0;
y=[sp.zeros((f1[0,1],f1[0,0])),\
   sp.zeros((f1[1,1],f1[1,0])),\
   sp.zeros((f1[2,1],f1[2,0])),\
   sp.zeros((f1[3,1],f1[3,0])),\
   sp.zeros((f1[4,1],f1[4,0])),]

for face in range(5):
  for j in range(f[face,1]):
    for i in range(f[face,0]):
      if nface not in BlankList:
        y[face][j*ssize:(j+1)*ssize,i*ssize:(i+1)*ssize]=x[:,nface2,:]
        nface2=nface2+1;
      else:
        y[face][j*ssize:(j+1)*ssize,i*ssize:(i+1)*ssize]=0;
      nface=nface+1;


# Building new composite array
print 'Building new composite array'
h=f1[0,0]

y2=sp.zeros((h*5+2,h*5+2))
y2[1:h*3+1,1:h+1]=y[0]
y2[1:h*3+1,h+1:h*2+1]=y[1]
y2[h*3+1:h*4+1,h+1:h*2+1]=y[2]
y2[h*3+1:h*4+1,h*2+1:h*5+1]=y[3]
y2[h*4+1:h*5+1,h*2+1:h*5+1]=y[4]

y2[1:h*3+1,0]=y[4][h-1,::-1]
y2[1:h*3+1,h*2+1]=y[3][0,::-1]
y2[h*3,h*2+1:h*5+1]=y[1][::-1,h-1]
y2[h*3+1,1:h+1]=y[2][::-1,0]
y2[h*3+1:h*4+1,h+1]=y[0][h*3-1,::-1]
y2[h*4+1,h+1:2*h+1]=y[4][::-1,0]
y2[h*4+1:h*5+1,2*h+1]=y[2][h-1,::-1]
y2[h*5+1,2*h+1:5*h+1]=y[0][::-1,0]

# calculation number of wet neigbors
print 'calculation number of wet neigbors'
tmp=sp.zeros((h*5+2,h*5+2))
ny2=sp.zeros((h*5+2,h*5+2))
tmp=y2
ny2[1:-1,1:-1]=tmp[1:-1,1:-1]+tmp[:-2,:-2]+tmp[1:-1,:-2]\
    +tmp[2:,:-2]+tmp[2:,1:-1]+tmp[2:,2:]\
    +tmp[1:-1,2:]+tmp[:-2,2:]+tmp[:-2,1:-1]

# converting back to llc faces from composite grid
print 'converting back to llc faces from composite grid'
ny=[sp.zeros((f1[0,1],f1[0,0])),\
   sp.zeros((f1[1,1],f1[1,0])),\
   sp.zeros((f1[2,1],f1[2,0])),\
   sp.zeros((f1[3,1],f1[3,0])),\
   sp.zeros((f1[4,1],f1[4,0])),]

ny[0]=ny2[1:h*3+1,1:h+1]
ny[1]=ny2[1:h*3+1,h+1:h*2+1]
ny[2]=ny2[h*3+1:h*4+1,h+1:h*2+1]
ny[3]=ny2[h*3+1:h*4+1,h*2+1:h*5+1]
ny[4]=ny2[h*4+1:h*5+1,h*2+1:h*5+1]

# converting back to the stacked format
print 'converting back to the stacked format'
nface=0
nface2=0;
nx=sp.zeros((sh[0],sh[1]/sh[0],sh[0]))
for face in range(5):
  for j in range(f[face,1]):
    for i in range(f[face,0]):
      if nface not in BlankList:
        nx[:,nface2,:]=ny[face][j*ssize:(j+1)*ssize,i*ssize:(i+1)*ssize]
        nface2=nface2+1;
      nface=nface+1;

nwet=sp.reshape(nx,(sh[0],sh[1]))

# Make slmask
slmask=wet.copy()
slmask[nwet<6.0]=0.0

#import matplotlib.pyplot as plt

#plt.figure(1)
#plt.pcolormesh(nwet[:,0:30])
#ax = plt.gca()
#ax.set_aspect('equal')
#plt.colorbar()
#plt.show(block=False)

#plt.figure(2)
#plt.pcolormesh(wet[:,0:30])
#ax = plt.gca()
#ax.set_aspect('equal')
#plt.colorbar()
#plt.show(block=False)

#plt.figure(3)
#plt.pcolormesh(slmask[:,0:30])
#ax = plt.gca()
#ax.set_aspect('equal')
#plt.colorbar()
#plt.show(block=False)

# Read routing table
print 'Reading routing data'
runoff_table='./CF0090x6C_LL5400xLL0015-Pfafstetter.TRN'

Nt=sp.fromfile(runoff_table,dtype='i4',count=3)[1]

vlist=[('','i4'),('Nt','i4'),('','i4'),
       ('','i4'),('sind','i4',Nt),('','i4'),
       ('','i4'),('dind','i4',Nt),('','i4'),
       ('','i4'),('weight','f4',Nt),('','i4')]

route=sp.fromfile(runoff_table,dtype=vlist);

# Read ocean i,j corresponding to tmp from tile file
print 'Reading tile data '
tfile='./CF0090x6C_LL5400xLL0015-Pfafstetter.til'
cols=(0,1,2,3,8,9)
rec=[('type','i4'),('area','f4'),('lon','f4'),('lat','f4'),
     ('ii','i4'),('jj','i4')]
tdata=sp.loadtxt(tfile, skiprows=8,usecols=cols,dtype=rec); 
#sp.save(dir+'/tdata.npy',tdata)
#tdata=sp.load(dir+'/tdata.npy')

# Append tile number
rec.append(('tnum','i4'))
tdata=sp.array(tdata,dtype=rec); tdata['tnum']=sp.arange(tdata.size)+1

# Subset ocean only tiles 
otdata=tdata[sp.where(tdata['type']==0)[0]]

# Remove tiles which are over ocean land (slmask[i,j]==0)
tslmask=sp.zeros(otdata.size)
for tt, ind in enumerate(otdata[['jj','ii']]):
    tslmask[tt]=slmask[ind[0]-1,ind[1]-1]
otdata=otdata[sp.where(tslmask > 0)[0]]

# Rerouting discharge
print 'Rerouting discharge...'
ind=(route['sind']!=route['dind'])
Nt=route['dind'][ind].size
vlist=[('','i4'),('Nt','i4'),('','i4'),
       ('','i4'),('sind','i4',Nt),('','i4'),
       ('','i4'),('dind','i4',Nt),('','i4'),
       ('','i4'),('weight','f4',Nt),('','i4')]

route_new=sp.array(sp.zeros(1,),dtype=vlist)
route_new['f0']=4
route_new['Nt']=Nt
route_new['f2']=4
route_new['f3']=Nt*4
route_new['sind']=route['sind'][ind]
route_new['f5']=Nt*4
route_new['f6']=Nt*4
route_new['dind'],route_new['weight']=reroute3(sp.ascontiguousarray(route['dind'][ind]),
                                              sp.ascontiguousarray(route['weight'][ind]),                    
                                              sp.ascontiguousarray(slmask),
                                              sp.ascontiguousarray(tdata['ii']), 
                                              sp.ascontiguousarray(tdata['jj']),
                                              sp.ascontiguousarray(tdata['lon']), 
                                              sp.ascontiguousarray(tdata['lat']),
                                              sp.ascontiguousarray(tdata['area']),
                                              sp.ascontiguousarray(otdata['tnum']))
route_new['f8']=Nt*4
route_new['f9']=Nt*4
route_new['f11']=Nt*4

route_new.tofile('./runoff_new.bin')
print '...done'

dstind_old=route['dind'][ind]
gind=sp.zeros((Nt,4),dtype='i4')
gind[:,0]=tdata[dstind_old-1]['ii']-1
gind[:,1]=tdata[dstind_old-1]['jj']-1
gind[:,2]=tdata[route_new['dind'][0]-1]['ii']-1
gind[:,3]=tdata[route_new['dind'][0]-1]['jj']-1
