import devpath
import os
ISSM_DIR = os.getenv('ISSM_DIR') # for binaries

import numpy as np
from triangle import triangle
from model import *
from netCDF4 import Dataset
from InterpFromGridToMesh import InterpFromGridToMesh
from bamg import bamg
from loadmodel import loadmodel
from setmask import setmask
from parameterize import parameterize
from clusters.discover_geos import export_discover
from marshall import marshall

if not os.path.exists('./netcdfs'):
    os.mkdir('./netcdfs')

# Step 1: Mesh generation
# Generate initial uniform mesh (resolution = 60000 m)
# project mesh onto new coordinate system
md = triangle(model(), '/discover/nobackup/agstubbl/ISSM/data/AntarcticaOutline.exp', 60000)

print('   Loading velocities data from NetCDF')
nsidc_vel = Dataset('/discover/nobackup/agstubbl/ISSM/data/Antarctica_ice_velocity.nc')
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
export_discover(md, './netcdfs/Antarctica_mesh.nc')

# Step 2: parameterize model
md = loadmodel('./netcdfs/Antarctica_mesh.nc')

md = setmask(md, '', '')
md = parameterize(md, './Antarctica_parameterize.py')


# export parameterization
export_discover(md, "./netcdfs/Antarctica_parameterization.nc",delete_rundir=True)

# Step 3: basal friction inversion 
	# Control general
md.inversion.nsteps = 400
md.inversion.iscontrol=1
md.inversion.maxsteps=400
md.inversion.maxiter=400
md.inversion.dxmin=0.01
md.inversion.gttol=1.0e-8

md.inversion.step_threshold = 0.99 * np.ones((md.inversion.nsteps))
md.inversion.maxiter_per_step = 40 * np.ones((md.inversion.nsteps))

md.inversion.gradient_scaling = 50 * np.ones((md.inversion.nsteps, 1))
md.inversion.min_parameters = 1 * np.ones((md.mesh.numberofvertices, 1))
md.inversion.max_parameters = 200 * np.ones((md.mesh.numberofvertices, 1))

#Cost functions
md.inversion.cost_functions = [101, 103, 501]
md.inversion.cost_functions_coefficients = np.ones((md.mesh.numberofvertices, 3))
md.inversion.cost_functions_coefficients[:, 0] = 1
md.inversion.cost_functions_coefficients[:, 1] = 1
md.inversion.cost_functions_coefficients[:, 2] = 2e-10

# Controls
md.inversion.control_parameters = ['FrictionCoefficient']
md.inversion.min_parameters=1*np.ones(np.shape(md.mesh.x))
md.inversion.max_parameters=200*np.ones(np.shape(md.mesh.x))

# Additional parameters
md.stressbalance.restol=0.0000001
md.stressbalance.reltol=0.0000001
md.stressbalance.abstol=np.nan

# Solve
md.private.solution = 'Stressbalance'
md.settings.waitonlock = 0
md.toolkits.ToolkitsFile(md.miscellaneous.name + '.toolkits')
marshall(md,md.miscellaneous.name+'.bin')

# export configuration for loading solution in next step
export_discover(md,'./netcdfs/Antarctica_inversion.nc',delete_rundir=True)
