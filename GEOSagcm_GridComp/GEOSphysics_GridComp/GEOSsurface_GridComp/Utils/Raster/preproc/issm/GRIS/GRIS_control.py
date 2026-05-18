import devpath
import os
ISSM_DIR = os.getenv('ISSM_DIR') # for binaries

import numpy as np
from model import *
from loadmodel import loadmodel
from setmask import setmask
from parameterize import parameterize
from setflowequation import setflowequation
from clusters.discover_geos import export_discover
from marshall import marshall


# Step 2: parameterize model
md = loadmodel('./netcdfs/GRIS_mesh.nc')
md = setmask(md, '', '')
md = parameterize(md, './GRIS_parameterize.py')
md = setflowequation(md, 'SSA', 'all')

# export parameterization
export_discover(md, "./netcdfs/GRIS_parameterization.nc",delete_rundir=True)

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
md.toolkits.ToolkitsFile('ISSM_'+md.miscellaneous.name + '.toolkits')
marshall(md,'ISSM_'+md.miscellaneous.name+'.bin') # create .bin file

# export configuration for loading solution in next step
export_discover(md,'./netcdfs/GRIS_inversion.nc',delete_rundir=True)
