import devpath
import os
ISSM_DIR = os.getenv('ISSM_DIR') # for binaries

from model import *
from toolkits import toolkits
from loadresultsfromdisk import loadresultsfromdisk
from marshall import marshall
from verbose import verbose 
from clusters.discover_geos import export_discover
from loadmodel import loadmodel
import numpy as np

md = loadmodel('./netcdfs/AIS_inversion.nc')
md = loadresultsfromdisk(md,'ISSM_AIS.outbin')

# update friction and velocity initialization based on inversion
md.friction.coefficient=md.results.StressbalanceSolution.FrictionCoefficient
md.initialization.vx=md.results.StressbalanceSolution.Vx
md.initialization.vy=md.results.StressbalanceSolution.Vy
md.initialization.vel=np.sqrt(md.initialization.vx**2 + md.initialization.vy**2)

# zero frontal melt rate and calving rate, placeholder zero smb
md.frontalforcings.meltingrate = np.zeros(md.mesh.numberofvertices)
md.calving.calvingrate = np.zeros(md.mesh.numberofvertices)
md.smb.mass_balance = np.zeros(md.mesh.numberofvertices) 

# Write the binary input file
# Additional options
md.inversion.iscontrol = 0
md.stressbalance.restol=1.0e-12
md.stressbalance.reltol=1.0e-12
md.transient.requested_outputs = ['default']
md.transient.isthermal=0
md.settings.waitonlock = 0
md.private.solution = 'Transient'
md.verbose = verbose('000000000')
md.toolkits = toolkits()
marshall(md,'ISSM_'+md.miscellaneous.name+'.bin') # create .bin file
md.toolkits.ToolkitsFile('ISSM_'+md.miscellaneous.name + '.toolkits')
export_discover(md,'./netcdfs/AIS_initialization.nc',delete_rundir=True)