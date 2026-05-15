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
from loadmodel import loadmodel
from setmask import setmask
from parameterize import parameterize
from setflowequation import setflowequation
from ll2xy import ll2xy
from clusters.discover_geos import export_discover
from marshall import marshall

if not os.path.exists('./netcdfs'):
    os.mkdir('./netcdfs')

# Step 1: Mesh generation
#Generate initial uniform mesh (resolution = 20000 m)
# project mesh onto new coordinate system
md = triangle(model(), '/discover/nobackup/agstubbl/ISSM/data/GreenlandOutline.exp', 20000)
[md.mesh.lat, md.mesh.long] = xy2ll(md.mesh.x, md.mesh.y, + 1, 39, 71)
[xi, yi] = ll2xy(md.mesh.lat, md.mesh.long, + 1, 45, 70)
md.mesh.x = xi
md.mesh.y = yi

ncdata_vx = Dataset('/discover/nobackup/agstubbl/ISSM/data/greenland_vel_mosaic200_2017_2018_vx_v02.1.nc', mode='r')
ncdata_vy = Dataset('/discover/nobackup/agstubbl/ISSM/data/greenland_vel_mosaic200_2017_2018_vy_v02.1.nc', mode='r')

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
export_discover(md, './netcdfs/Greenland_mesh.nc')

# Step 2: parameterize model
md = loadmodel('./netcdfs/Greenland_mesh.nc')
md = setmask(md, '', '')
md = parameterize(md, './Greenland_parameterize.py')
md = setflowequation(md, 'SSA', 'all')

# export parameterization
export_discover(md, "./netcdfs/Greenland_parameterization.nc",delete_rundir=True)

# Step 3: basal friction inversion 
#Control general
md.inversion.iscontrol = 1
md.inversion.nsteps = 200
md.inversion.step_threshold = 0.99 * np.ones((md.inversion.nsteps))
md.inversion.maxiter_per_step = 10 * np.ones((md.inversion.nsteps))

#Cost functions
md.inversion.cost_functions = [101, 103, 501]
md.inversion.cost_functions_coefficients = np.ones((md.mesh.numberofvertices, 3))
md.inversion.cost_functions_coefficients[:, 0] = 350
md.inversion.cost_functions_coefficients[:, 1] = 0.2
md.inversion.cost_functions_coefficients[:, 2] = 2e-6

#Controls
md.inversion.control_parameters = ['FrictionCoefficient']
md.inversion.gradient_scaling = 50 * np.ones((md.inversion.nsteps, 1))
md.inversion.min_parameters = 1 * np.ones((md.mesh.numberofvertices, 1))
md.inversion.max_parameters = 200 * np.ones((md.mesh.numberofvertices, 1))

#Additional parameters
md.stressbalance.restol = 0.01
md.stressbalance.reltol = 0.1

md.private.solution = 'Stressbalance'
md.settings.waitonlock = 0
md.toolkits.ToolkitsFile(md.miscellaneous.name + '.toolkits')
marshall(md,md.miscellaneous.name+'.bin') # create .bin file

# export configuration for loading solution in next step
export_discover(md,'./netcdfs/Greenland_inversion.nc',delete_rundir=True)
