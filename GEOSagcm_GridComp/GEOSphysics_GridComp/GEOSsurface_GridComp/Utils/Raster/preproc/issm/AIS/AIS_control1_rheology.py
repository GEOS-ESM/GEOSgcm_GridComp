import devpath
import os
ISSM_DIR = os.getenv('ISSM_DIR') # for binaries

import numpy as np
from model import *
from marshall import marshall
from m1qn3inversion import m1qn3inversion
from loadmodel import loadmodel
from reinitializelevelset import reinitializelevelset
from verbose import verbose
from cuffey import cuffey
from clusters.discover_geos import export_discover

# Load model parameterization
md = loadmodel("./netcdfs/AIS_parameterization.nc")

# We need to make some adjustments to the mask to make sure it is compatible with the observations
pos = np.where(md.inversion.vel_obs==0)[0]
levelset_initial = md.mask.ice_levelset.copy()
md.mask.ice_levelset[pos] = +1 # set these areas to be ice-free
md.mask.ice_levelset=reinitializelevelset(md,md.mask.ice_levelset)

#Control general
md.inversion=m1qn3inversion(md.inversion)
md.inversion.iscontrol=1
md.verbose=verbose('solution',True,'control',True)

#Cost functions
md.inversion.cost_functions= [101, 103, 502]
md.inversion.cost_functions_coefficients=np.zeros((md.mesh.numberofvertices,len(md.inversion.cost_functions)))
md.inversion.cost_functions_coefficients[:,0]=5000
md.inversion.cost_functions_coefficients[:,1]=100
md.inversion.cost_functions_coefficients[:,2]=1e-18

# where there's no ice or now flow, don't invert:
pos=np.where( (md.mask.ice_levelset>0) | (md.inversion.vel_obs==0))[0]
md.inversion.cost_functions_coefficients[:,0][pos]=0
md.inversion.cost_functions_coefficients[:,1][pos]=0

#Controls
md.inversion.control_parameters=['MaterialsRheologyBbar']
md.inversion.maxsteps=200
md.inversion.maxiter =200
pos=np.where( (md.mask.ocean_levelset<0) )[0]
md.inversion.min_parameters=cuffey(273)*np.ones(md.mesh.numberofvertices)  
md.inversion.max_parameters=cuffey(200)*np.ones(md.mesh.numberofvertices) 
md.inversion.control_scaling_factors=1e8
md.inversion.dfmin_frac = .4
md.inversion.dxmin = 1e-5

#Additional parameters
md.stressbalance.abstol=np.nan

#Extract model over ice shelves
mds= md.extract((md.mask.ocean_levelset<0))

# marshall binary for solve via command line 
mds.private.solution = 'Stressbalance'
mds.settings.waitonlock = 0
mds.toolkits.ToolkitsFile('ISSM_'+mds.miscellaneous.name + '.toolkits')
marshall(mds,'ISSM_'+mds.miscellaneous.name+'.bin')

export_discover(mds,'./netcdfs/AIS_inversion_rheology.nc',delete_rundir=True)