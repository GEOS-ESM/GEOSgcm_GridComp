import devpath
import os
ISSM_DIR = os.getenv('ISSM_DIR') # for binaries

import numpy as np
from triangle import triangle
from model import *
from netCDF4 import Dataset
from InterpFromGridToMesh import InterpFromGridToMesh
from bamg import bamg
from xy2ll import xy2ll
from ll2xy import ll2xy
from clusters.discover_geos import export_discover

if not os.path.exists('./netcdfs'):
    os.mkdir('./netcdfs')

# Step 1: Mesh generation
#Generate initial uniform mesh (resolution = 20000 m)
# project mesh onto new coordinate system
md = triangle(model(), '/discover/nobackup/agstubbl/ISSM/data/GRIS/GreenlandOutline.exp', 20000)
[md.mesh.lat, md.mesh.long] = xy2ll(md.mesh.x, md.mesh.y, + 1, 39, 71)
[xi, yi] = ll2xy(md.mesh.lat, md.mesh.long, + 1, 45, 70)
md.mesh.x = xi
md.mesh.y = yi

ncdata_vx = Dataset('/discover/nobackup/agstubbl/ISSM/data/GRIS/greenland_vel_mosaic200_2017_2018_vx_v02.1.nc', mode='r')
ncdata_vy = Dataset('/discover/nobackup/agstubbl/ISSM/data/GRIS/greenland_vel_mosaic200_2017_2018_vy_v02.1.nc', mode='r')

# Get velocities (Note: You can use ncprint('file') to see an ncdump)
x1 = np.squeeze(ncdata_vx.variables['x'][:].data)
y1 = np.squeeze(ncdata_vx.variables['y'][:].data)
velx = np.squeeze(ncdata_vx.variables['Band1'][:].data)
vely = np.squeeze(ncdata_vy.variables['Band1'][:].data)
ncdata_vx.close()
ncdata_vy.close()

velx[np.abs(velx)>1e9] = 0
vely[np.abs(vely)>1e9] = 0

vx = InterpFromGridToMesh(x1, y1, velx, md.mesh.x, md.mesh.y, 0)
vy = InterpFromGridToMesh(x1, y1, vely, md.mesh.x, md.mesh.y, 0)
speed = np.sqrt(vx**2 + vy**2)

# Mesh Greenland (refine according to flow speed)
md = bamg(md, 'hmax', 75000, 'hmin', 5000, 'gradation', 1.4, 'field', speed, 'err', 8)

# export mesh
export_discover(md, './netcdfs/GRIS_mesh.nc')