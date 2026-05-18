import devpath
import os
ISSM_DIR = os.getenv('ISSM_DIR') # for binaries

import numpy as np
from triangle import triangle
from model import *
from netCDF4 import Dataset
from InterpFromGridToMesh import InterpFromGridToMesh
from bamg import bamg
from clusters.discover_geos import export_discover

if not os.path.exists('./netcdfs'):
    os.mkdir('./netcdfs')

# Step 1: Mesh generation
# Generate initial uniform mesh (resolution = 60000 m)
# project mesh onto new coordinate system
md = triangle(model(), '/discover/nobackup/agstubbl/ISSM/data/AIS/AntarcticaOutline.exp', 60000)

print('   Loading velocities data from NetCDF')
nsidc_vel = Dataset('/discover/nobackup/agstubbl/ISSM/data/AIS/Antarctica_ice_velocity.nc')
xmin = nsidc_vel.xmin
xmin = float(xmin.lstrip()[0:10])
ymax = nsidc_vel.ymax
ymax = float(ymax.lstrip()[0:10])
spacing = nsidc_vel.spacing
spacing = float((spacing.lstrip())[0:4])
nx = nsidc_vel.nx
ny = nsidc_vel.ny
velx = nsidc_vel['vx'][:].data
vely = nsidc_vel['vy'][:].data
# Build coordinates
x2 = xmin + np.arange(nx + 1) * spacing
y2 = (ymax - ny * spacing) + np.arange(ny + 1) * spacing

# print('   Set observed velocities')
md.initialization.vx = InterpFromGridToMesh(x2, y2, np.flipud(velx), md.mesh.x, md.mesh.y, 0)
md.initialization.vy = InterpFromGridToMesh(x2, y2, np.flipud(vely), md.mesh.x, md.mesh.y, 0)
md.initialization.vz = np.zeros(md.mesh.numberofvertices)
md.initialization.vel = np.sqrt(md.initialization.vx**2 + md.initialization.vy**2)
del velx, vely

md = bamg(md, 'hmax', 75000, 'hmin', 5000, 'gradation', 1.4, 'field', md.initialization.vel, 'err', 8)

# export mesh
export_discover(md, './netcdfs/AIS_mesh.nc')