import devpath
import sys,os
import numpy as np
from triangle import triangle
from model import *
from netCDF4 import Dataset
from InterpFromGridToMesh import InterpFromGridToMesh
from bamg import bamg
from clusters.discover_geos import export_discover

h_max = sys.argv[1] if len(sys.argv) > 1 else 24000
h_min = sys.argv[2] if len(sys.argv) > 2 else 2000
h_max = float(h_max)
h_min = float(h_min)

print(f'h_max: {h_max}')
print(f'h_min: {h_min}')

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
vx = nsidc_vel['vx'][:].data
vy = nsidc_vel['vy'][:].data
# Build coordinates
x2 = xmin + np.arange(nx + 1) * spacing
y2 = (ymax - ny * spacing) + np.arange(ny + 1) * spacing

vx_ = InterpFromGridToMesh(x2, y2, np.flipud(vx), md.mesh.x, md.mesh.y, 0)
vy_ = InterpFromGridToMesh(x2, y2, np.flipud(vy), md.mesh.x, md.mesh.y, 0)
speed = np.sqrt(vx_**2 + vy_**2)
del vx, vy, vx_, vy_

md = bamg(md, 'hmax', h_max, 'hmin', h_min, 'gradation', 1.4, 'field', speed, 'err', 8)

# export mesh
export_discover(md, './netcdfs/AIS_mesh.nc')
