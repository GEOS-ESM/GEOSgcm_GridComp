import devpath
import os
ISSM_DIR = os.getenv('ISSM_DIR') # for binaries

import numpy as np
from model import *
from loadmodel import loadmodel
from setmask import setmask
from parameterize import parameterize
from clusters.discover_geos import export_discover
from marshall import marshall

# Step 2: parameterize model
md = loadmodel('./netcdfs/AIS_mesh.nc')

md = setmask(md, '', '')
md = parameterize(md, './AIS_parameterize.py')


# export parameterization
export_discover(md, "./netcdfs/AIS_parameterization.nc",delete_rundir=True)

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
md.toolkits.ToolkitsFile('ISSM_'+md.miscellaneous.name + '.toolkits')
marshall(md,'ISSM_'+md.miscellaneous.name+'.bin')

# export configuration for loading solution in next step
export_discover(md,'./netcdfs/AIS_inversion.nc',delete_rundir=True)
