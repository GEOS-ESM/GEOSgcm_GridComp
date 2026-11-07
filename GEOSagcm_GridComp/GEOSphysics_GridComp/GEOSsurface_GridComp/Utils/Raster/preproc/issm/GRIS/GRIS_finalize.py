import devpath
import os
ISSM_DIR = os.getenv('ISSM_DIR') # for binaries

from model import *
from loadresultsfromdisk import loadresultsfromdisk
from marshall import marshall
from verbose import verbose 
from cluster.discover_geos import export_discover
from loadmodel import loadmodel
import numpy as np
from toolkits import toolkits

md = loadmodel('./netcdfs/GRIS_inversion.nc')
md = loadresultsfromdisk(md,'ISSM_GRIS.outbin')
md.friction.coefficient=md.results.StressbalanceSolution.FrictionCoefficient
md.initialization.vx=md.results.StressbalanceSolution.Vx
md.initialization.vy=md.results.StressbalanceSolution.Vy
md.initialization.vel=np.sqrt(md.initialization.vx**2 + md.initialization.vy**2)

# Write the binary input file
# Additional options
md.inversion.iscontrol = 0
md.stressbalance.restol=1.0e-12
md.stressbalance.reltol=1.0e-12

# Transient options
md.transient.ismovingfront = 0
md.transient.isthermal = 0
md.transient.isstressbalance = 1
md.transient.ismasstransport = 1
md.transient.isgroundingline = 1

nv = md.mesh.numberofvertices
onboundary = np.asarray(md.mesh.vertexonboundary).flatten().astype(bool)

# BC for mass transport
spcthickness = np.full((nv, 1), np.nan)
spcthickness[onboundary, 0] = np.asarray(md.geometry.thickness).flatten()[onboundary]
md.masstransport.spcthickness = spcthickness
md.masstransport.stabilization = 5  # SUPG

# Levelset options (https://gmd.copernicus.org/articles/17/6227/2024/)
spclevelset = np.full((nv, 1), np.nan)
spclevelset[onboundary, 0] = np.asarray(md.mask.ice_levelset).flatten()[onboundary]
md.levelset.spclevelset = spclevelset
md.levelset.reinit_frequency = 20  # Should be more than 10
md.levelset.stabilization = 5      # SUPG

# No frontal melt by default
md.frontalforcings.meltingrate = np.zeros((nv, 1))

md.private.solution = 'Transient'
md.verbose = verbose('000000000')
md.toolkits = toolkits()
marshall(md,'ISSM_'+md.miscellaneous.name+'.bin') # create .bin file
md.toolkits.ToolkitsFile('ISSM_'+md.miscellaneous.name + '.toolkits')
export_discover(md,'./netcdfs/GRIS_initialization.nc',delete_rundir=True)

